#include "svf.h"

static float sin_poly(float x) {
  const float x2 = x * x;
  return x * (1.0f - x2 * (1.0f / 6.0f) * (1.0f - x2 * (1.0f / 20.0f)));
}

static float clampf(float v, float lo, float hi) {
  return v < lo ? lo : (v > hi ? hi : v);
}

void udsp_svf_reset(UDSP_SVF *s) {
  s->lp[0] = s->lp[1] = s->bp[0] = s->bp[1] = 0.0f;
}

void udsp_svf_set(UDSP_SVF *s, float cutoff_hz, float res, UDSP_SVF_MODE mode,
                  unsigned rate) {
  float x = 3.14159265f * cutoff_hz / (float)(rate ? rate : 48000);
  s->f = 2.0f * sin_poly(clampf(x, 0.0f, 0.7f));
  s->q = 1.0f - res;
  s->mode = mode;
}

/* one sample through one channel's filter: state in, state out */
static inline float svf_step(float x, float f, float q, UDSP_SVF_MODE mode,
                             float *lp, float *bp) {
  const float hp = x - *lp - q * *bp;
  *bp += f * hp;
  *lp += f * *bp;
  const float out =
      mode == UDSP_SVF_HIGH ? hp : (mode == UDSP_SVF_BAND ? *bp : *lp);
  *lp = clampf(*lp, -4.0f, 4.0f);
  *bp = clampf(*bp, -4.0f, 4.0f);
  return clampf(out, -4.0f, 4.0f);
}

/* both channels in one loop: each channel is one serial chain (~20
 * cycles a sample), two of them overlap -- 1.5-1.8x over a channel at a
 * time, the same operations per channel, so the same bits */
void udsp_svf_stereo(UDSP_SVF *s, float *io, size_t frames) {
  const float f = s->f, q = s->q;
  const UDSP_SVF_MODE mode = s->mode;
  float lp0 = s->lp[0], bp0 = s->bp[0], lp1 = s->lp[1], bp1 = s->bp[1];
  for (size_t i = 0; i < frames; i++) {
    io[i * 2] = svf_step(io[i * 2], f, q, mode, &lp0, &bp0);
    io[i * 2 + 1] = svf_step(io[i * 2 + 1], f, q, mode, &lp1, &bp1);
  }
  s->lp[0] = lp0;
  s->bp[0] = bp0;
  s->lp[1] = lp1;
  s->bp[1] = bp1;
}
