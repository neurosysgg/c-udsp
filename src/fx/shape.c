#include "fx/shape.h"

#include "base/simd.h"
#include "synth/osc.h"

void udsp_softclip(float *buf, size_t n, float drive) {
  const float inv = 1.0f / (1.0f + drive * 0.25f);
  size_t i = 0;
#if UDSP_W > 1
  const UDSP_VF vd = udsp_vsplat(drive), one = udsp_vsplat(1.0f);
  for (; i + UDSP_W <= n; i += UDSP_W) {
    const UDSP_VF x = udsp_vload(buf + i) * vd;
    const UDSP_VF y = x / (one + udsp_vabs(x));
    udsp_vstore(buf + i, y * inv * 2.0f);
  }
#endif
  for (; i < n; i++) {
    const float x = buf[i] * drive;
    const float y = x / (1.0f + (x < 0.0f ? -x : x));
    buf[i] = y * inv * 2.0f;
  }
}

float udsp_shape1(UDSP_SHAPE kind, float x) {
  switch (kind) {
  case UDSP_SHAPE_HARD:
    return x > 1.0f ? 1.0f : (x < -1.0f ? -1.0f : x);
  case UDSP_SHAPE_FOLD: {
    /* sin(x), x's turns wrapped to 0..1 */
    float t = 0.15915494f * x;
    t -= (float)(int)t;
    t += t < 0.0f ? 1.0f : 0.0f;
    return udsp_osc_sine(t);
  }
  case UDSP_SHAPE_TUBE:
    return x >= 0.0f ? x / (1.0f + x) : x / (1.0f - 2.0f * x);
  default:
    if (x >= 3.0f || x <= -3.0f) {
      return x > 0.0f ? 1.0f : -1.0f;
    }
    return x * (27.0f + x * x) / (27.0f + 9.0f * x * x);
  }
}

void udsp_waveshape(float *buf, size_t n, UDSP_SHAPE kind, float drive) {
  for (size_t i = 0; i < n; i++) {
    buf[i] = udsp_shape1(kind, buf[i] * drive);
  }
}

#if UDSP_W > 1
/* udsp_shape1() per lane: each branch computed, the lanes picked (the
 * same operations, so the same bits) */
static UDSP_VF vshape(UDSP_SHAPE kind, UDSP_VF x) {
  const UDSP_VF one = udsp_vsplat(1.0f), neg = udsp_vsplat(-1.0f);
  switch (kind) {
  case UDSP_SHAPE_HARD:
    return udsp_vsel(x > one, one, udsp_vsel(x < neg, neg, x));
  case UDSP_SHAPE_TUBE:
    return udsp_vsel(x >= 0.0f, x / (1.0f + x), x / (1.0f - 2.0f * x));
  default: {
    /* one expression, as udsp_shape1()'s: a compiler that contracts
     * (clang) fuses the same multiplies and adds in both */
    const UDSP_VF p = x * (27.0f + x * x) / (27.0f + 9.0f * x * x);
    const UDSP_VI big = (x >= 3.0f) | (x <= -3.0f);
    return udsp_vsel(big, udsp_vsel(x > 0.0f, one, neg), p);
  }
  }
}
#endif

#define MIX_RUN 256

void udsp_shape_mix(float *io, size_t n, UDSP_SHAPE kind, float drive,
                    float mix, float out) {
  const float dry = 1.0f - mix;
  size_t i = 0;
#if UDSP_W > 1
  float w[MIX_RUN];
  while (n - i >= UDSP_W) {
    size_t k = n - i < MIX_RUN ? n - i : MIX_RUN;
    k -= k % UDSP_W;
    float *x = io + i;
    if (kind == UDSP_SHAPE_FOLD) {
      /* sin(x), x's turns wrapped to 0..1: the turns, then the sines */
      for (size_t j = 0; j < k; j += UDSP_W) {
        UDSP_VF t = 0.15915494f * (udsp_vload(x + j) * drive);
        t -= __builtin_convertvector(__builtin_convertvector(t, UDSP_VI),
                                     UDSP_VF);
        t += udsp_vsel(t < 0.0f, udsp_vsplat(1.0f), udsp_vsplat(0.0f));
        udsp_vstore(w + j, t);
      }
      udsp_osc_sines(w, w, k);
    } else {
      for (size_t j = 0; j < k; j += UDSP_W) {
        udsp_vstore(w + j, vshape(kind, udsp_vload(x + j) * drive));
      }
    }
    for (size_t j = 0; j < k; j += UDSP_W) {
      const UDSP_VF v = udsp_vload(x + j);
      udsp_vstore(x + j, (v * dry + udsp_vload(w + j) * mix) * out);
    }
    i += k;
  }
#endif
  for (; i < n; i++) {
    const float x = io[i];
    io[i] = (x * dry + udsp_shape1(kind, x * drive) * mix) * out;
  }
}
