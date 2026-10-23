#include "fx/gain.h"

#include <math.h>

#include "base/simd.h"

void udsp_gain_apply(float *buf, size_t n, float gain) {
  for (size_t i = 0; i < n; i++) {
    buf[i] *= gain;
  }
}

void udsp_gain_stereo(float *buf, size_t frames, float gl, float gr) {
  size_t i = 0, n = frames * 2;
#if UDSP_W > 1
  UDSP_VF g;
  for (int k = 0; k < UDSP_W; k++) {
    g[k] = k & 1 ? gr : gl;
  }
  for (; i + UDSP_W <= n; i += UDSP_W) {
    udsp_vstore(buf + i, udsp_vload(buf + i) * g);
  }
#endif
  for (; i < n; i++) {
    buf[i] *= i & 1 ? gr : gl;
  }
}

void udsp_mix_add(float *dst, const float *src, size_t n, float gain) {
  size_t i = 0;
#if UDSP_W > 1
  const UDSP_VF g = udsp_vsplat(gain);
  for (; i + UDSP_W <= n; i += UDSP_W) {
    udsp_vstore(dst + i, udsp_vload(dst + i) + udsp_vload(src + i) * g);
  }
#endif
  for (; i < n; i++) {
    dst[i] += src[i] * gain;
  }
}

void udsp_pan_gains(float pan, UDSP_PAN_LAW law, float *gl, float *gr) {
  pan = pan < -1.0f ? -1.0f : (pan > 1.0f ? 1.0f : pan);
  const float quarter = 1.5707963267948966f;
  if (law == UDSP_PAN_EQUAL_POWER) {
    const float th = (pan + 1.0f) * 0.5f * quarter;
    *gl = cosf(th);
    *gr = sinf(th);
    return;
  }
  /* the side panned away from: cos of how far, so exactly 1 at centre */
  *gl = pan > 0.0f ? cosf(pan * quarter) : 1.0f;
  *gr = pan < 0.0f ? cosf(-pan * quarter) : 1.0f;
}
