#include "fx/xover.h"

#include <string.h>

#include "base/simd.h"
#include "synth/osc.h"

#define BUTTER_K 1.41421356f /* 1 / Q, Q = 1 / sqrt(2) */
#define LR2_K 2.0f           /* Q = 1 / 2: two first-order poles */

void udsp_xover_reset(UDSP_XOVER *x) {
  memset(x->sec, 0, sizeof x->sec);
  memset(x->ap, 0, sizeof x->ap);
}

static void set(UDSP_XOVER *x, int bands, const float *hz, unsigned rate,
                int lr2) {
  const float k = lr2 ? LR2_K : BUTTER_K;
  bands = bands < 1 ? 1 : (bands > UDSP_XOVER_BANDS ? UDSP_XOVER_BANDS : bands);
  x->bands = bands;
  x->lr2 = lr2;
  const float lo = 10.0f / (float)rate;
  for (int s = 0; s < bands - 1; s++) {
    float t = hz[s] / (float)rate;
    t = t < lo ? lo : (t > 0.49f ? 0.49f : t);
    /* tan(pi t): eq.c's */
    const float g = udsp_osc_sine(0.5f * t) / udsp_osc_sine(0.5f * t + 0.25f);
    UDSP_XOVER_SPLIT *c = &x->c[s];
    c->k = k;
    c->a1 = 1.0f / (1.0f + g * (g + k));
    c->a2 = g * c->a1;
    c->a3 = g * c->a2;
  }
}

void udsp_xover_set(UDSP_XOVER *x, int bands, const float *hz, unsigned rate) {
  set(x, bands, hz, rate, 0);
}

void udsp_xover_set_lr2(UDSP_XOVER *x, int bands, const float *hz,
                        unsigned rate) {
  set(x, bands, hz, rate, 1);
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

/* split() at 12 dB: one section, its low pass low, its high pass
 * inverted kept */
UDSP_NO_SLP static void split_lr2(const UDSP_XOVER_SPLIT *c, float (*sec)[2],
                                  float *rest, float *low, size_t frames) {
  float s[2][2];
  memcpy(s, sec, sizeof s);
  for (size_t i = 0; i < frames; i++) {
    for (int ch = 0; ch < 2; ch++) {
      const float in = rest[2 * i + ch];
      float v1, v2;
      tick(c, in, s[ch], &v1, &v2);
      low[2 * i + ch] = v2;
      rest[2 * i + ch] = c->k * v1 + v2 - in;
    }
  }
  memcpy(sec, s, sizeof s);
}

/* allpass() at 12 dB: low minus high, 2 band + 2 low - x */
UDSP_NO_SLP static void allpass_lr2(const UDSP_XOVER_SPLIT *c, float (*st)[2],
                                    float *io, size_t frames) {
  float s[2][2];
  memcpy(s, st, sizeof s);
  for (size_t i = 0; i < frames; i++) {
    for (int ch = 0; ch < 2; ch++) {
      const float in = io[2 * i + ch];
      float v1, v2;
      tick(c, in, s[ch], &v1, &v2);
      io[2 * i + ch] = c->k * v1 + 2.0f * v2 - in;
    }
  }
  memcpy(st, s, sizeof s);
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

#if defined(UDSP_HAVE_V4)
/* tick() in four lanes */
typedef struct {
  UDSP_V4 a1, a2, a3, ic0, ic1;
} LANES;

static inline void vtick(LANES *s, UDSP_V4 x, UDSP_V4 *v1, UDSP_V4 *v2,
                         UDSP_V4 *ic0, UDSP_V4 *ic1) {
  const UDSP_V4 v3 = x - s->ic1;
  *v1 = s->a1 * s->ic0 + s->a2 * v3;
  *v2 = s->ic1 + s->a2 * s->ic0 + s->a3 * v3;
  *ic0 = 2.0f * *v1 - s->ic0;
  *ic1 = 2.0f * *v2 - s->ic1;
}

static inline void advance(LANES *s, UDSP_V4 ic0, UDSP_V4 ic1, UDSP_V4I on,
                           int edge) {
  s->ic0 = edge ? udsp_v4sel(on, ic0, s->ic0) : ic0;
  s->ic1 = edge ? udsp_v4sel(on, ic1, s->ic1) : ic1;
}

/* Three bands as one wavefront (eq.c's chain, as a graph): a channel's
 * seven sections in four-lane vectors, each a frame behind the one
 * feeding it --
 *   a  split 0's first section, split 1's (L R L R)  frames t, t - 2
 *   b  their low sides' second                     t - 1, t - 3
 *   c  their high sides' second                    t - 1, t - 3
 *   g  band 0's allpass at split 1 (L R, two idle)  t - 2
 * b and c take a's own lanes, g b's split 0 lanes, a's split 1 lanes c's
 * split 0 lanes: one shuffle a step. The same operations per section as
 * split() and allpass(), so the same bits. */
typedef struct {
  LANES a, b, c, g;
  UDSP_V4 ka, kg; /* k for a, b, c (split 0's, split 1's); g's 2 k */
  UDSP_V4 av2, ahp, bv2, cout; /* the last step's outputs */
} WAVE3;

static inline void step3(WAVE3 *w, const float *in, float *const *out, int t,
                         int frames, int edge) {
  UDSP_V4 x = {0.0f, 0.0f, 0.0f, 0.0f};
  if (!edge || t < frames) {
    x[0] = in[2 * t];
    x[1] = in[2 * t + 1];
  }
  x = UDSP_SHUF4(x, w->cout, 0, 1, 4, 5);
  UDSP_V4 v1, v2, i0, i1, w1, w2, j0, j1;
  UDSP_V4I on = {0};
  if (edge) {
    const UDSP_V4I tt = (UDSP_V4I){0} + t, lag01 = {0, 0, 2, 2};
    on = (tt >= lag01) & (tt - lag01 < frames);
  }
  vtick(&w->a, x, &v1, &v2, &i0, &i1);
  const UDSP_V4 hp = x - w->ka * v1 - v2;
  advance(&w->a, i0, i1, on, edge);

  if (edge) {
    const UDSP_V4I tt = (UDSP_V4I){0} + t, lag13 = {1, 1, 3, 3};
    on = (tt >= lag13) & (tt - lag13 < frames);
  }
  vtick(&w->b, w->av2, &w1, &w2, &j0, &j1);
  advance(&w->b, j0, j1, on, edge);
  const UDSP_V4 low = w2;
  UDSP_V4 c1, c2;
  vtick(&w->c, w->ahp, &c1, &c2, &j0, &j1);
  const UDSP_V4 high = w->ahp - w->ka * c1 - c2;
  advance(&w->c, j0, j1, on, edge);

  if (edge) {
    const UDSP_V4I tt = (UDSP_V4I){0} + t, lag2 = {2, 2, 2, 2};
    on = (tt >= lag2) & (tt - lag2 < frames);
  }
  UDSP_V4 g1, g2;
  vtick(&w->g, w->bv2, &g1, &g2, &j0, &j1);
  const UDSP_V4 ap = w->bv2 - w->kg * g1;
  advance(&w->g, j0, j1, on, edge);

  if (!edge || (t >= 2 && t - 2 < frames)) {
    out[0][2 * (t - 2)] = ap[0];
    out[0][2 * (t - 2) + 1] = ap[1];
  }
  if (!edge || t >= 3) {
    out[1][2 * (t - 3)] = low[2];
    out[1][2 * (t - 3) + 1] = low[3];
    out[2][2 * (t - 3)] = high[2];
    out[2][2 * (t - 3) + 1] = high[3];
  }
  w->av2 = v2;
  w->ahp = hp;
  w->bv2 = low;
  w->cout = high;
}

static void lanes_in(LANES *l, const UDSP_XOVER_SPLIT *lo,
                     const UDSP_XOVER_SPLIT *hi, const float (*s0)[2],
                     const float (*s1)[2]) {
  l->a1 = (UDSP_V4){lo->a1, lo->a1, hi->a1, hi->a1};
  l->a2 = (UDSP_V4){lo->a2, lo->a2, hi->a2, hi->a2};
  l->a3 = (UDSP_V4){lo->a3, lo->a3, hi->a3, hi->a3};
  l->ic0 = (UDSP_V4){s0[0][0], s0[1][0], s1[0][0], s1[1][0]};
  l->ic1 = (UDSP_V4){s0[0][1], s0[1][1], s1[0][1], s1[1][1]};
}

static void lanes_out(const LANES *l, float (*s0)[2], float (*s1)[2]) {
  for (int ch = 0; ch < 2; ch++) {
    s0[ch][0] = l->ic0[ch];
    s0[ch][1] = l->ic1[ch];
    if (s1) {
      s1[ch][0] = l->ic0[2 + ch];
      s1[ch][1] = l->ic1[2 + ch];
    }
  }
}

static void split3(UDSP_XOVER *x, const float *in, float *const *out,
                   size_t frames) {
  const UDSP_XOVER_SPLIT *c0 = &x->c[0], *c1 = &x->c[1];
  static const float idle[2][2] = {{0.0f, 0.0f}, {0.0f, 0.0f}};
  WAVE3 w;
  lanes_in(&w.a, c0, c1, x->sec[0][0], x->sec[1][0]);
  lanes_in(&w.b, c0, c1, x->sec[0][1], x->sec[1][1]);
  lanes_in(&w.c, c0, c1, x->sec[0][2], x->sec[1][2]);
  lanes_in(&w.g, c1, c1, x->ap[0][1], idle);
  w.ka = (UDSP_V4){c0->k, c0->k, c1->k, c1->k};
  w.kg = (UDSP_V4){0.0f, 0.0f, 0.0f, 0.0f} + 2.0f * c1->k;
  w.av2 = w.ahp = w.bv2 = w.cout = (UDSP_V4){0.0f, 0.0f, 0.0f, 0.0f};
  const int n = (int)frames, steps = n + 3;
  const int mid = n > 3 ? 3 : n;
  int t = 0;
  for (; t < mid; t++) {
    step3(&w, in, out, t, n, 1);
  }
  for (; t < n; t++) {
    step3(&w, in, out, t, n, 0);
  }
  for (; t < steps; t++) {
    step3(&w, in, out, t, n, 1);
  }
  lanes_out(&w.a, x->sec[0][0], x->sec[1][0]);
  lanes_out(&w.b, x->sec[0][1], x->sec[1][1]);
  lanes_out(&w.c, x->sec[0][2], x->sec[1][2]);
  lanes_out(&w.g, x->ap[0][1], NULL);
}
#endif

void udsp_xover_split(UDSP_XOVER *x, const float *in, float *const *out,
                      size_t frames) {
  const int n = x->bands;
#if defined(UDSP_HAVE_V4)
  if (n == 3 && !x->lr2) {
    if (frames) {
      split3(x, in, out, frames);
    }
    return;
  }
#endif
  float *rest = out[n - 1];
  if (rest != in) {
    memcpy(rest, in, frames * 2 * sizeof(float));
  }
  for (int s = 0; s < n - 1; s++) {
    if (x->lr2) {
      split_lr2(&x->c[s], x->sec[s][0], rest, out[s], frames);
    } else {
      split(&x->c[s], x->sec[s], rest, out[s], frames);
    }
  }
  for (int b = 0; b < n - 2; b++) {
    for (int s = b + 1; s < n - 1; s++) {
      if (x->lr2) {
        allpass_lr2(&x->c[s], x->ap[b][s], out[b], frames);
      } else {
        allpass(&x->c[s], x->ap[b][s], out[b], frames);
      }
    }
  }
}
