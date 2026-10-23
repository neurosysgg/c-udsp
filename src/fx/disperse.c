#include "fx/disperse.h"

#include <string.h>

#include "base/simd.h"
#include "synth/osc.h"
#include "synth/pitch.h"

void udsp_disperse_reset(UDSP_DISPERSE *d) {
  memset(d->ic, 0, sizeof d->ic);
  d->used = 0;
}

/* stage i's place in the fan, -0.5 .. 0.5: the base-2 van der Corput
 * sequence from its second term (0, -1/4, 1/4, -3/8, 1/8, ...) */
static float fan(unsigned i) {
  float v = 0.0f, b = 0.5f;
  for (unsigned n = i + 1; n; n >>= 1, b *= 0.5f) {
    v += (n & 1) ? b : 0.0f;
  }
  return v - 0.5f;
}

void udsp_disperse_set(UDSP_DISPERSE *d, float amount, float hz, float pinch,
                       float spread, unsigned rate) {
  amount = amount<0.0f ? 0.0f : amount>(float) UDSP_DISPERSE_MAX
               ? (float)UDSP_DISPERSE_MAX
               : amount;
  pinch = pinch < 0.1f ? 0.1f : (pinch > 40.0f ? 40.0f : pinch);
  spread = spread < 0.0f ? 0.0f : (spread > 8.0f ? 8.0f : spread);
  d->whole = (int)amount;
  d->frac = amount - (float)d->whole;
  const int n = d->whole + (d->frac > 0.0f && d->whole < UDSP_DISPERSE_MAX);
  const float k = 1.0f / pinch, lo = 10.0f / (float)rate;
  d->k = k;
  for (int i = 0; i < n; i++) {
    float t = hz * (float)udsp_exp2(spread * fan((unsigned)i)) / (float)rate;
    t = t < lo ? lo : (t > 0.49f ? 0.49f : t);
    const float g = udsp_osc_sine(0.5f * t) / udsp_osc_sine(0.5f * t + 0.25f);
    UDSP_DISPERSE_STAGE *s = &d->st[i];
    s->a1 = 1.0f / (1.0f + g * (g + k));
    s->a2 = g * s->a1;
    s->a3 = g * s->a2;
  }
}

/* one frame of a stage's channel: x in, x - w * band out */
static inline float tick(const UDSP_DISPERSE_STAGE *s, float w, float *ic,
                         float x) {
  const float v3 = x - ic[1];
  const float v1 = s->a1 * ic[0] + s->a2 * v3;
  const float v2 = ic[1] + s->a2 * ic[0] + s->a3 * v3;
  ic[0] = 2.0f * v1 - ic[0];
  ic[1] = 2.0f * v2 - ic[1];
  return x - w * v1;
}

/* m stages over the block, a frame through all of them at a time: each
 * stage's two channels are a chain of ~15 cycles a frame, so 2m chains
 * overlap instead of two (svf.c's finding, session 221) -- the same
 * operations each, the same bits as a stage at a time */
#define GROUP 4
static inline __attribute__((always_inline)) void
stages(const UDSP_DISPERSE_STAGE *st, const float *w, float (*ic)[2][2],
       float *io, size_t frames, const int m) {
  float s[GROUP][2][2];
  memcpy(s, ic, (size_t)m * sizeof s[0]);
  for (size_t i = 0; i < frames; i++) {
    float l = io[2 * i], r = io[2 * i + 1];
    UDSP_UNROLL
    for (int j = 0; j < m; j++) {
      l = tick(&st[j], w[j], s[j][0], l);
      r = tick(&st[j], w[j], s[j][1], r);
    }
    io[2 * i] = l;
    io[2 * i + 1] = r;
  }
  memcpy(ic, s, (size_t)m * sizeof s[0]);
}

#if UDSP_W > 1
/* WAVE stages in vector lanes, a wavefront: lane k is stage k at the
 * frame k steps behind lane 0's, its input the output lane k - 1 made a
 * step before. Each stage still runs its frames in order through tick()'s
 * operations, so the bits are a stage at a time's; but the chains of a
 * step are all independent, WAVE of them a channel. At the block's ends
 * (the first and last WAVE - 1 steps) the lanes without a frame of
 * their own keep their state. */
#define WAVE_G 2 /* vectors of stages: two keep enough chains in flight */
#define WAVE (WAVE_G * UDSP_W)
static void wavefront(const UDSP_DISPERSE_STAGE *st, const float *w,
                      float (*ic)[2][2], float *io, size_t frames) {
  UDSP_VF a1[WAVE_G], a2[WAVE_G], a3[WAVE_G], wv[WAVE_G];
  UDSP_VF s1[2][WAVE_G], s2[2][WAVE_G], out[2][WAVE_G];
  UDSP_VI lane[WAVE_G];
  for (int g = 0; g < WAVE_G; g++) {
    for (int j = 0; j < UDSP_W; j++) {
      const int k = g * UDSP_W + j;
      a1[g][j] = st[k].a1;
      a2[g][j] = st[k].a2;
      a3[g][j] = st[k].a3;
      wv[g][j] = w[k];
      lane[g][j] = k;
      for (int c = 0; c < 2; c++) {
        s1[c][g][j] = ic[k][c][0];
        s2[c][g][j] = ic[k][c][1];
        out[c][g][j] = 0.0f;
      }
    }
  }
  const size_t steps = frames + WAVE - 1;
  for (size_t n = 0; n < steps; n++) {
    const int edge = n < WAVE - 1 || n >= frames;
    for (int c = 0; c < 2; c++) {
      UDSP_VF x[WAVE_G];
      x[0] = udsp_vshift_in(out[c][0], n < frames ? io[2 * n + c] : 0.0f);
      UDSP_UNROLL
      for (int g = 1; g < WAVE_G; g++) {
        x[g] = udsp_vshift_in(out[c][g], out[c][g - 1][UDSP_W - 1]);
      }
      UDSP_UNROLL
      for (int g = 0; g < WAVE_G; g++) {
        const UDSP_VF v3 = x[g] - s2[c][g];
        const UDSP_VF v1 = a1[g] * s1[c][g] + a2[g] * v3;
        const UDSP_VF v2 = s2[c][g] + a2[g] * s1[c][g] + a3[g] * v3;
        const UDSP_VF n1 = 2.0f * v1 - s1[c][g], n2 = 2.0f * v2 - s2[c][g];
        out[c][g] = x[g] - wv[g] * v1;
        if (edge) {
          const UDSP_VI f = (int)n - lane[g];
          const UDSP_VI ok = (f >= 0) & (f < (int)frames);
          s1[c][g] = udsp_vsel(ok, n1, s1[c][g]);
          s2[c][g] = udsp_vsel(ok, n2, s2[c][g]);
        } else {
          s1[c][g] = n1;
          s2[c][g] = n2;
        }
      }
      if (n >= WAVE - 1) {
        io[2 * (n - (WAVE - 1)) + c] = out[c][WAVE_G - 1][UDSP_W - 1];
      }
    }
  }
  for (int g = 0; g < WAVE_G; g++) {
    for (int j = 0; j < UDSP_W; j++) {
      for (int c = 0; c < 2; c++) {
        ic[g * UDSP_W + j][c][0] = s1[c][g][j];
        ic[g * UDSP_W + j][c][1] = s2[c][g][j];
      }
    }
  }
}
#endif

void udsp_disperse_stereo(UDSP_DISPERSE *d, float *io, size_t frames) {
  const int n = d->whole + (d->frac > 0.0f && d->whole < UDSP_DISPERSE_MAX);
  float w[UDSP_DISPERSE_MAX];
  for (int i = 0; i < n; i++) {
    /* the last one partly in: its output faded in over its input */
    w[i] = (i < d->whole ? 1.0f : d->frac) * 2.0f * d->k;
  }
  int i = 0;
#if UDSP_W > 1
  for (; i + WAVE <= n; i += WAVE) {
    wavefront(d->st + i, w + i, d->ic + i, io, frames);
  }
#endif
  for (; i + GROUP <= n; i += GROUP) {
    stages(d->st + i, w + i, d->ic + i, io, frames, GROUP);
  }
  for (; i < n; i++) {
    stages(d->st + i, w + i, d->ic + i, io, frames, 1);
  }
  /* a stage out keeps no state: in again, it starts from rest */
  if (d->used > n) {
    memset(d->ic[n], 0, (size_t)(d->used - n) * sizeof d->ic[0]);
  }
  d->used = n;
}
