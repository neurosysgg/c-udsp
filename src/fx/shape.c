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
