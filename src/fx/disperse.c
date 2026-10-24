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
/* G vectors of stages in lanes, a wavefront: lane k is stage k at the
 * frame k steps behind lane 0's, its input the output lane k - 1 made a
 * step before. Each stage still runs its frames in order through tick()'s
 * operations, so the bits are a stage at a time's; but the chains of a
 * step are all independent, G * UDSP_W of them a channel. At the block's
 * ends (the first and last G * UDSP_W - 1 steps) the lanes without a
 * frame of their own keep their state. Lanes past the m stages pass
 * their input on as it is (zero coefficients and state: v1 = +0, out =
 * x - +0 = x, -0 too) and keep no state, so a remainder short of a
 * wavefront still runs as one instead of a stage at a time. */
#define WAVE_G 2 /* vectors of stages: two keep enough chains in flight */
#define WAVE (WAVE_G * UDSP_W)
static inline __attribute__((always_inline)) void
wavefront(const UDSP_DISPERSE_STAGE *st, const float *w, float (*ic)[2][2],
          float *io, size_t frames, int m, const int G) {
  const int lanes = G * UDSP_W;
  UDSP_VF a1[WAVE_G], a2[WAVE_G], a3[WAVE_G], wv[WAVE_G];
  UDSP_VF s1[2][WAVE_G], s2[2][WAVE_G], out[2][WAVE_G];
  UDSP_VI lane[WAVE_G];
  /* the lanes as floats first, then a vector load each (eq.c's: lanes
   * stored one at a time into a vector in memory don't forward to its
   * load) */
  float f[8][WAVE];
  int ln[WAVE];
  for (int k = 0; k < lanes; k++) {
    const UDSP_DISPERSE_STAGE z = {0}, *sk = k < m ? &st[k] : &z;
    f[0][k] = sk->a1;
    f[1][k] = sk->a2;
    f[2][k] = sk->a3;
    f[3][k] = k < m ? w[k] : 0.0f;
    for (int c = 0; c < 2; c++) {
      f[4 + 2 * c][k] = k < m ? ic[k][c][0] : 0.0f;
      f[5 + 2 * c][k] = k < m ? ic[k][c][1] : 0.0f;
    }
    ln[k] = k;
  }
  for (int g = 0; g < G; g++) {
    const int at = g * UDSP_W;
    a1[g] = udsp_vload(f[0] + at);
    a2[g] = udsp_vload(f[1] + at);
    a3[g] = udsp_vload(f[2] + at);
    wv[g] = udsp_vload(f[3] + at);
    for (int c = 0; c < 2; c++) {
      s1[c][g] = udsp_vload(f[4 + 2 * c] + at);
      s2[c][g] = udsp_vload(f[5 + 2 * c] + at);
      out[c][g] = (UDSP_VF){0};
    }
    memcpy(&lane[g], ln + at, sizeof lane[g]);
  }
  const size_t steps = frames + (size_t)lanes - 1;
  for (size_t n = 0; n < steps; n++) {
    const int edge = n < (size_t)lanes - 1 || n >= frames;
    for (int c = 0; c < 2; c++) {
      UDSP_VF x[WAVE_G];
      x[0] = udsp_vshift_in(out[c][0], n < frames ? io[2 * n + c] : 0.0f);
      UDSP_UNROLL
      for (int g = 1; g < G; g++) {
        x[g] = udsp_vshift_in(out[c][g], out[c][g - 1][UDSP_W - 1]);
      }
      UDSP_UNROLL
      for (int g = 0; g < G; g++) {
        const UDSP_VF v3 = x[g] - s2[c][g];
        const UDSP_VF v1 = a1[g] * s1[c][g] + a2[g] * v3;
        const UDSP_VF v2 = s2[c][g] + a2[g] * s1[c][g] + a3[g] * v3;
        const UDSP_VF n1 = 2.0f * v1 - s1[c][g], n2 = 2.0f * v2 - s2[c][g];
        out[c][g] = x[g] - wv[g] * v1;
        if (edge) {
          const UDSP_VI fr = (int)n - lane[g];
          const UDSP_VI ok = (fr >= 0) & (fr < (int)frames);
          s1[c][g] = udsp_vsel(ok, n1, s1[c][g]);
          s2[c][g] = udsp_vsel(ok, n2, s2[c][g]);
        } else {
          s1[c][g] = n1;
          s2[c][g] = n2;
        }
      }
      if (n >= (size_t)lanes - 1) {
        io[2 * (n - ((size_t)lanes - 1)) + c] = out[c][G - 1][UDSP_W - 1];
      }
    }
  }
  for (int g = 0; g < G; g++) {
    for (int c = 0; c < 2; c++) {
      udsp_vstore(f[4 + 2 * c] + g * UDSP_W, s1[c][g]);
      udsp_vstore(f[5 + 2 * c] + g * UDSP_W, s2[c][g]);
    }
  }
  for (int k = 0; k < m; k++) {
    for (int c = 0; c < 2; c++) {
      ic[k][c][0] = f[4 + 2 * c][k];
      ic[k][c][1] = f[5 + 2 * c][k];
    }
  }
}

/* a remainder of up to this many stages a stage at a time instead: one
 * stage's two chains beat a vector of pads under gcc (17 stages at -v3
 * 17.6 -> 13.0 us per 1024 frames, SSE2 23.6 -> 20.6), not under clang
 * (SSE2's 9 stages 12.2 -> 14.1); two or more always lose.
 * -DWAVE_SCALAR=N forces it (bench/README.md "session 263") */
#ifndef WAVE_SCALAR
#ifdef __clang__
#define WAVE_SCALAR 0
#else
#define WAVE_SCALAR 1
#endif
#endif

static void wave_full(const UDSP_DISPERSE_STAGE *st, const float *w,
                      float (*ic)[2][2], float *io, size_t frames, int m) {
  wavefront(st, w, ic, io, frames, m, WAVE_G);
}

static void wave_one(const UDSP_DISPERSE_STAGE *st, const float *w,
                     float (*ic)[2][2], float *io, size_t frames, int m) {
  wavefront(st, w, ic, io, frames, m, 1);
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
  for (; i < n && n - i > WAVE_SCALAR; i += WAVE) {
    const int m = n - i < WAVE ? n - i : WAVE;
    (m > UDSP_W ? wave_full : wave_one)(d->st + i, w + i, d->ic + i, io, frames,
                                        m);
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
