#include "synth/ramp.h"

#include "base/simd.h"

void udsp_ramp_init(UDSP_RAMP *r, unsigned frames) {
  r->value = 0.0f;
  r->target = 0.0f;
  r->step = 1.0f / (float)(frames ? frames : 1);
}

void udsp_ramp_stereo(UDSP_RAMP *r, float *io, size_t frames) {
  size_t f = 0;
  float e = r->value;
  /* the moving part one frame at a time: each step is one float add, and
   * k steps are not one k * step multiply (a reference repeating the
   * adds has to land on the same values) */
  for (; f < frames && e != r->target; f++) {
    if (e < r->target) {
      e += r->step;
      if (e > r->target) {
        e = r->target;
      }
    } else {
      e -= r->step;
      if (e < r->target) {
        e = r->target;
      }
    }
    io[f * 2] *= e;
    io[f * 2 + 1] *= e;
  }
  r->value = e;
  if (f < frames && e != 1.0f) {
    float *p = io + f * 2;
    size_t n = (frames - f) * 2, i = 0;
#if UDSP_W > 1
    const UDSP_VF ve = udsp_vsplat(e);
    for (; i + UDSP_W <= n; i += UDSP_W) {
      udsp_vstore(p + i, udsp_vload(p + i) * ve);
    }
#endif
    for (; i < n; i++) {
      p[i] *= e;
    }
  }
}
