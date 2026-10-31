#include "fx/conv.h"

#include "base/fft.h"
#include "base/simd.h"

#include <string.h>

size_t udsp_conv_parts(size_t frames, size_t block) {
  size_t p = (frames + block - 1) / block;
  return p ? p : 1;
}

void udsp_conv_ir(const float *ir, size_t frames, size_t stride, size_t block,
                  float *spectra, const float *table, float *scratch) {
  const size_t parts = udsp_conv_parts(frames, block),
               s = UDSP_CONV_SPECTRUM(block);
  const float scale = 1.0f / (float)(2 * block);
  for (size_t p = 0; p < parts; p++) {
    /* a partition in the frame's first half, zeros after it: overlap-save
     * keeps the circular product's second half, which then holds no wrap */
    memset(scratch, 0, 2 * block * sizeof(float));
    for (size_t i = 0; i < block && p * block + i < frames; i++) {
      scratch[i] = ir[(p * block + i) * stride] * scale;
    }
    float *re = spectra + p * s;
    udsp_fft_real(scratch, re, re + block + 1, 2 * block, table);
  }
}

void udsp_conv_in_init(UDSP_CONV_IN *in, size_t block, size_t parts,
                       float *line, float *frame) {
  in->block = block;
  in->parts = parts;
  in->pos = 0;
  in->line = line;
  in->frame = frame;
  memset(line, 0, parts * UDSP_CONV_SPECTRUM(block) * sizeof(float));
  memset(frame, 0, 2 * block * sizeof(float));
}

void udsp_conv_in_push(UDSP_CONV_IN *in, const float *x, size_t stride,
                       const float *table) {
  const size_t b = in->block;
  memmove(in->frame, in->frame + b, b * sizeof(float));
  for (size_t i = 0; i < b; i++) {
    in->frame[b + i] = x[i * stride];
  }
  in->pos = in->pos + 1 == in->parts ? 0 : in->pos + 1;
  float *re = in->line + in->pos * UDSP_CONV_SPECTRUM(b);
  udsp_fft_real(in->frame, re, re + b + 1, 2 * b, table);
}

void udsp_conv_mac(const UDSP_CONV_IN *in, const float *spectra, size_t parts,
                   float *acc) {
  udsp_conv_mac_range(in, spectra, 0, parts, 0, acc);
}

void udsp_conv_mac_range(const UDSP_CONV_IN *in, const float *spectra,
                         size_t from, size_t to, size_t lead, float *acc) {
  const size_t bins = in->block + 1, s = UDSP_CONV_SPECTRUM(in->block);
  float *restrict ar = acc, *restrict ai = acc + bins;
  size_t slot = (in->pos + in->parts - (from - lead) % in->parts) % in->parts;
  for (size_t p = from; p < to; p++) {
    const float *restrict xr = in->line + slot * s, *restrict xi = xr + bins;
    const float *restrict hr = spectra + p * s, *restrict hi = hr + bins;
    size_t k = 0;
    for (; k + UDSP_W <= bins; k += UDSP_W) {
      UDSP_VF a = udsp_vload(xr + k), b = udsp_vload(xi + k);
      UDSP_VF c = udsp_vload(hr + k), d = udsp_vload(hi + k);
      udsp_vstore(ar + k, udsp_vload(ar + k) + (a * c - b * d));
      udsp_vstore(ai + k, udsp_vload(ai + k) + (a * d + b * c));
    }
    for (; k < bins; k++) {
      ar[k] += xr[k] * hr[k] - xi[k] * hi[k];
      ai[k] += xr[k] * hi[k] + xi[k] * hr[k];
    }
    slot = slot ? slot - 1 : in->parts - 1;
  }
}

void udsp_conv_out(float *acc, size_t block, float *out, size_t stride,
                   const float *table, float *scratch) {
  udsp_fft_real_inverse(acc, acc + block + 1, scratch, 2 * block, table);
  for (size_t i = 0; i < block; i++) {
    out[i * stride] = scratch[block + i];
  }
}

void udsp_conv_head(const float *ir, size_t frames, size_t stride, size_t taps,
                    float *hr) {
  for (size_t j = 0; j < taps; j++) {
    hr[taps - 1 - j] = j < frames ? ir[j * stride] : 0.0f;
  }
}

void udsp_conv_direct(const float *hr, size_t taps, const float *x, float *out,
                      size_t n) {
  const float *x0 = x - (taps - 1);
  size_t i = 0;
  /* four outputs' vectors side by side: each one's sum is a chain of
   * adds over the taps (in their order, as one at a time), and four in
   * flight hide its latency */
  for (; i + 4 * UDSP_W <= n; i += 4 * UDSP_W) {
    UDSP_VF a0 = udsp_vload(out + i), a1 = udsp_vload(out + i + UDSP_W);
    UDSP_VF a2 = udsp_vload(out + i + 2 * UDSP_W);
    UDSP_VF a3 = udsp_vload(out + i + 3 * UDSP_W);
    const float *p = x0 + i;
    for (size_t j = 0; j < taps; j++, p++) {
      const UDSP_VF h = udsp_vsplat(hr[j]);
      a0 = a0 + h * udsp_vload(p);
      a1 = a1 + h * udsp_vload(p + UDSP_W);
      a2 = a2 + h * udsp_vload(p + 2 * UDSP_W);
      a3 = a3 + h * udsp_vload(p + 3 * UDSP_W);
    }
    udsp_vstore(out + i, a0);
    udsp_vstore(out + i + UDSP_W, a1);
    udsp_vstore(out + i + 2 * UDSP_W, a2);
    udsp_vstore(out + i + 3 * UDSP_W, a3);
  }
  for (; i + UDSP_W <= n; i += UDSP_W) {
    UDSP_VF acc = udsp_vload(out + i);
    for (size_t j = 0; j < taps; j++) {
      acc = acc + udsp_vsplat(hr[j]) * udsp_vload(x0 + i + j);
    }
    udsp_vstore(out + i, acc);
  }
  for (; i < n; i++) {
    float acc = out[i];
    for (size_t j = 0; j < taps; j++) {
      acc = acc + hr[j] * x0[i + j];
    }
    out[i] = acc;
  }
}
