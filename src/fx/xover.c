#include "fx/xover.h"

#include <string.h>

#include "base/simd.h"
#include "synth/osc.h"

#define BUTTER_K 1.41421356f /* 1 / Q, Q = 1 / sqrt(2) */

void udsp_xover_reset(UDSP_XOVER *x) {
  memset(x->sec, 0, sizeof x->sec);
  memset(x->ap, 0, sizeof x->ap);
}

void udsp_xover_set(UDSP_XOVER *x, int bands, const float *hz, unsigned rate) {
  bands = bands < 1 ? 1 : (bands > UDSP_XOVER_BANDS ? UDSP_XOVER_BANDS : bands);
  x->bands = bands;
  const float lo = 10.0f / (float)rate;
  for (int s = 0; s < bands - 1; s++) {
    float t = hz[s] / (float)rate;
    t = t < lo ? lo : (t > 0.49f ? 0.49f : t);
    /* tan(pi t): eq.c's */
    const float g = udsp_osc_sine(0.5f * t) / udsp_osc_sine(0.5f * t + 0.25f);
    UDSP_XOVER_SPLIT *c = &x->c[s];
    c->k = BUTTER_K;
    c->a1 = 1.0f / (1.0f + g * (g + BUTTER_K));
    c->a2 = g * c->a1;
    c->a3 = g * c->a2;
  }
}

/* one trapezoidal SVF tick: v1 the band, v2 the low output */
static inline void tick(const UDSP_XOVER_SPLIT *c, float x, float *ic,
                        float *v1, float *v2) {
  const float v3 = x - ic[1];
  *v1 = c->a1 * ic[0] + c->a2 * v3;
  *v2 = ic[1] + c->a2 * ic[0] + c->a3 * v3;
  ic[0] = 2.0f * *v1 - ic[0];
  ic[1] = 2.0f * *v2 - ic[1];
}

/* rest's frames split at c: its low side into low, its high side kept;
 * both channels in one loop (two chains overlap: svf.c's finding) */
UDSP_NO_SLP static void split(const UDSP_XOVER_SPLIT *c, float (*sec)[2][2],
                              float *rest, float *low, size_t frames) {
  float s[3][2][2];
  memcpy(s, sec, sizeof s);
  for (size_t i = 0; i < frames; i++) {
    for (int ch = 0; ch < 2; ch++) {
      const float in = rest[2 * i + ch];
      float v1, v2, w1, w2;
      tick(c, in, s[0][ch], &v1, &v2);
      const float hp = in - c->k * v1 - v2;
      tick(c, v2, s[1][ch], &w1, &w2);
      low[2 * i + ch] = w2;
      tick(c, hp, s[2][ch], &w1, &w2);
      rest[2 * i + ch] = hp - c->k * w1 - w2;
    }
  }
  memcpy(sec, s, sizeof s);
}

/* io's frames through c's allpass: x - 2k band */
UDSP_NO_SLP static void allpass(const UDSP_XOVER_SPLIT *c, float (*st)[2],
                                float *io, size_t frames) {
  float s[2][2];
  memcpy(s, st, sizeof s);
  for (size_t i = 0; i < frames; i++) {
    for (int ch = 0; ch < 2; ch++) {
      const float in = io[2 * i + ch];
      float v1, v2;
      tick(c, in, s[ch], &v1, &v2);
      io[2 * i + ch] = in - 2.0f * c->k * v1;
    }
  }
  memcpy(st, s, sizeof s);
}

void udsp_xover_split(UDSP_XOVER *x, const float *in, float *const *out,
                      size_t frames) {
  const int n = x->bands;
  float *rest = out[n - 1];
  if (rest != in) {
    memcpy(rest, in, frames * 2 * sizeof(float));
  }
  for (int s = 0; s < n - 1; s++) {
    split(&x->c[s], x->sec[s], rest, out[s], frames);
  }
  for (int b = 0; b < n - 2; b++) {
    for (int s = b + 1; s < n - 1; s++) {
      allpass(&x->c[s], x->ap[b][s], out[b], frames);
    }
  }
}
