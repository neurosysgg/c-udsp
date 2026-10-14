#include "shape.h"

#include "simd.h"

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
