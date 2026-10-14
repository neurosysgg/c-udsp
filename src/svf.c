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

void udsp_svf_stereo(UDSP_SVF *s, float *io, size_t frames) {
  const float f = s->f, q = s->q;
  for (int c = 0; c < 2; c++) {
    float lp = s->lp[c], bp = s->bp[c];
    for (size_t i = 0; i < frames; i++) {
      const float hp = io[i * 2 + c] - lp - q * bp;
      bp += f * hp;
      lp += f * bp;
      const float out =
          s->mode == UDSP_SVF_HIGH ? hp : (s->mode == UDSP_SVF_BAND ? bp : lp);
      lp = clampf(lp, -4.0f, 4.0f);
      bp = clampf(bp, -4.0f, 4.0f);
      io[i * 2 + c] = clampf(out, -4.0f, 4.0f);
    }
    s->lp[c] = lp;
    s->bp[c] = bp;
  }
}
