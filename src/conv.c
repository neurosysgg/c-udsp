#include "conv.h"

#include "fft.h"
#include "simd.h"

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
  const size_t bins = in->block + 1, s = UDSP_CONV_SPECTRUM(in->block);
  float *restrict ar = acc, *restrict ai = acc + bins;
  size_t slot = in->pos;
  for (size_t p = 0; p < parts; p++) {
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
