#include "fx/eq.h"

#include <math.h>

#include "base/simd.h"
#include "synth/osc.h"
#include "synth/pitch.h"

/* an 8th-order Butterworth's sections: 1 / (2 cos((2i + 1) pi / 16)) */
static const float BUTTER8_Q[UDSP_EQ_STAGES] = {0.50979558f, 0.60134489f,
                                                0.89997622f, 2.56291545f};

void udsp_eq_reset(UDSP_EQ_BAND *b) {
  b->ramp = 0;
  for (int s = 0; s < UDSP_EQ_STAGES; s++) {
    b->ic1[s][0] = b->ic1[s][1] = b->ic2[s][0] = b->ic2[s][1] = 0.0f;
  }
}

/* tan(pi t), t a fraction of the rate: a half turn's angle as turns */
static float tan_rate(float t) {
  return udsp_osc_sine(0.5f * t) / udsp_osc_sine(0.5f * t + 0.25f);
}

static void stage(UDSP_EQ_STAGE *s, float g, float k, float m0, float m1,
                  float m2) {
  s->g = g;
  s->k = k;
  s->a1 = 1.0f / (1.0f + g * (g + k));
  s->a2 = g * s->a1;
  s->a3 = g * s->a2;
  s->m0 = m0;
  s->m1 = m1;
  s->m2 = m2;
}

void udsp_eq_set(UDSP_EQ_BAND *b, UDSP_EQ_SHAPE shape, float hz, float q,
                 float gain_db, unsigned rate) {
  float t = hz / (float)rate;
  const float lo = 10.0f / (float)rate;
  t = t < lo ? lo : (t > 0.49f ? 0.49f : t);
  q = q < 0.1f ? 0.1f : (q > 40.0f ? 40.0f : q);
  b->ramp = 0;
  gain_db = gain_db < -30.0f ? -30.0f : (gain_db > 30.0f ? 30.0f : gain_db);
  const float g = tan_rate(t);
  const float k = 1.0f / q;
  /* the amplitude's square root: 10^(dB / 40) */
  const float a = (float)udsp_exp2(gain_db * (3.321928094887362 / 40.0));
  b->stages = 1;
  switch (shape) {
  case UDSP_EQ_LOW_CUT_48:
  case UDSP_EQ_HIGH_CUT_48: {
    const int high = shape == UDSP_EQ_HIGH_CUT_48;
    /* the band's q over a Butterworth's 0.71 sharpens the last section */
    const float res = q * 1.41421356f;
    for (int i = 0; i < UDSP_EQ_STAGES; i++) {
      float sq = BUTTER8_Q[i] * (i == UDSP_EQ_STAGES - 1 ? res : 1.0f);
      float sk = 1.0f / sq;
      if (high) {
        stage(&b->st[i], g, sk, 0.0f, 0.0f, 1.0f);
      } else {
        stage(&b->st[i], g, sk, 1.0f, -sk, -1.0f);
      }
    }
    b->stages = UDSP_EQ_STAGES;
    break;
  }
  case UDSP_EQ_LOW_CUT_12:
    stage(&b->st[0], g, k, 1.0f, -k, -1.0f);
    break;
  case UDSP_EQ_HIGH_CUT_12:
    stage(&b->st[0], g, k, 0.0f, 0.0f, 1.0f);
    break;
  case UDSP_EQ_LOW_SHELF:
    stage(&b->st[0], g / sqrtf(a), k, 1.0f, k * (a - 1.0f), a * a - 1.0f);
    break;
  case UDSP_EQ_HIGH_SHELF:
    stage(&b->st[0], g * sqrtf(a), k, a * a, k * (1.0f - a) * a, 1.0f - a * a);
    break;
  case UDSP_EQ_NOTCH:
    stage(&b->st[0], g, k, 1.0f, -k, 0.0f);
    break;
  default: {
    /* the bell: its bandwidth kept as the gain changes (Q / A) */
    const float kb = 1.0f / (q * a);
    stage(&b->st[0], g, kb, 1.0f, kb * (a * a - 1.0f), 0.0f);
    break;
  }
  }
}

void udsp_eq_ramp(UDSP_EQ_BAND *b, const UDSP_EQ_BAND *to, unsigned frames) {
  if (to->stages != b->stages || frames == 0) {
    for (int k = 0; k < UDSP_EQ_STAGES; k++) {
      b->st[k] = to->st[k];
    }
    b->stages = to->stages;
    b->ramp = 0;
    return;
  }
  const float n = (float)frames;
  for (int k = 0; k < b->stages; k++) {
    const UDSP_EQ_STAGE *x = &b->st[k], *y = &to->st[k];
    b->step[k][0] = (y->g - x->g) / n;
    b->step[k][1] = (y->k - x->k) / n;
    b->step[k][2] = (y->m0 - x->m0) / n;
    b->step[k][3] = (y->m1 - x->m1) / n;
    b->step[k][4] = (y->m2 - x->m2) / n;
    b->end[k] = *y;
  }
  b->ramp = frames;
}

/* one sample through a section's channel: state in, state out */
static inline float step(const UDSP_EQ_STAGE *s, float v0, float *ic1,
                         float *ic2) {
  const float v3 = v0 - *ic2;
  const float v1 = s->a1 * *ic1 + s->a2 * v3;
  const float v2 = *ic2 + s->a2 * *ic1 + s->a3 * v3;
  *ic1 = 2.0f * v1 - *ic1;
  *ic2 = 2.0f * v2 - *ic2;
  return s->m0 * v0 + s->m1 * v1 + s->m2 * v2;
}

/* n frames of a ramp through section k: its coefficients a frame on,
 * then the frame */
UDSP_NO_SLP static void ramp_stage(UDSP_EQ_BAND *b, int k, float *io,
                                   size_t n) {
  UDSP_EQ_STAGE s = b->st[k];
  const float *d = b->step[k];
  float l1 = b->ic1[k][0], l2 = b->ic2[k][0];
  float r1 = b->ic1[k][1], r2 = b->ic2[k][1];
  for (size_t i = 0; i < n; i++) {
    s.g += d[0];
    s.k += d[1];
    s.m0 += d[2];
    s.m1 += d[3];
    s.m2 += d[4];
    s.a1 = 1.0f / (1.0f + s.g * (s.g + s.k));
    s.a2 = s.g * s.a1;
    s.a3 = s.g * s.a2;
    io[i * 2] = step(&s, io[i * 2], &l1, &l2);
    io[i * 2 + 1] = step(&s, io[i * 2 + 1], &r1, &r2);
  }
  b->st[k] = s;
  b->ic1[k][0] = l1;
  b->ic2[k][0] = l2;
  b->ic1[k][1] = r1;
  b->ic2[k][1] = r2;
}

/* a section at a time over the block, both channels in one loop (svf.c's
 * reasoning: two serial chains overlap) */
UDSP_NO_SLP void udsp_eq_stereo(UDSP_EQ_BAND *b, float *io, size_t frames) {
  if (b->ramp) {
    const size_t n = frames < b->ramp ? frames : b->ramp;
    for (int k = 0; k < b->stages; k++) {
      ramp_stage(b, k, io, n);
    }
    b->ramp -= (unsigned)n;
    if (b->ramp == 0) {
      /* the last frame's on the target exactly */
      for (int k = 0; k < b->stages; k++) {
        b->st[k] = b->end[k];
      }
    }
    io += 2 * n;
    frames -= n;
  }
  for (int k = 0; k < b->stages; k++) {
    const UDSP_EQ_STAGE s = b->st[k];
    float l1 = b->ic1[k][0], l2 = b->ic2[k][0];
    float r1 = b->ic1[k][1], r2 = b->ic2[k][1];
    for (size_t i = 0; i < frames; i++) {
      io[i * 2] = step(&s, io[i * 2], &l1, &l2);
      io[i * 2 + 1] = step(&s, io[i * 2 + 1], &r1, &r2);
    }
    b->ic1[k][0] = l1;
    b->ic2[k][0] = l2;
    b->ic1[k][1] = r1;
    b->ic2[k][1] = r2;
  }
}

#if UDSP_W > 1
/* A wavefront over sections (disperse.c's): lane k is section k at the
 * frame k steps behind lane 0's, its input the output lane k - 1 made a
 * step before. Each section runs its frames in order through step()'s
 * operations, so the bits are a section at a time's, but a step's chains
 * are all independent: one recurrence's latency a frame for the whole
 * chain instead of one a section. Lanes past the last section have zero
 * coefficients and state and stay zero; at the block's ends (the first
 * and last m - 1 steps) the sections without a frame of their own keep
 * their state. */
#define CHAIN_LANES 16 /* a wavefront's sections at most */
#define CHAIN_MAX_G (CHAIN_LANES / UDSP_W)

typedef struct {
  const UDSP_EQ_STAGE *st;
  float *ic1, *ic2; /* [2]: left, right */
} SECTION;

/* step n: the sections a frame on each; edge: some lanes have no frame */
static inline __attribute__((always_inline)) void
wave_step(const UDSP_VF *a1, const UDSP_VF *a2, const UDSP_VF *a3,
          const UDSP_VF *m0, const UDSP_VF *m1, const UDSP_VF *m2,
          UDSP_VF (*s1)[CHAIN_MAX_G], UDSP_VF (*s2)[CHAIN_MAX_G],
          UDSP_VF (*out)[CHAIN_MAX_G], const UDSP_VI *lane, float *io,
          size_t frames, size_t n, int m, int last_g, int last_j,
          const int edge, const int G) {
  UDSP_UNROLL
  for (int c = 0; c < 2; c++) {
    UDSP_VF x[CHAIN_MAX_G];
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
      out[c][g] = m0[g] * x[g] + m1[g] * v1 + m2[g] * v2;
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
    if (n >= (size_t)m - 1) {
      io[2 * (n - ((size_t)m - 1)) + c] = out[c][last_g][last_j];
    }
  }
}

static inline __attribute__((always_inline)) void
wavefront(const SECTION *sec, int m, float *io, size_t frames, const int G) {
  UDSP_VF a1[CHAIN_MAX_G], a2[CHAIN_MAX_G], a3[CHAIN_MAX_G];
  UDSP_VF m0[CHAIN_MAX_G], m1[CHAIN_MAX_G], m2[CHAIN_MAX_G];
  UDSP_VF s1[2][CHAIN_MAX_G], s2[2][CHAIN_MAX_G], out[2][CHAIN_MAX_G];
  UDSP_VI lane[CHAIN_MAX_G];
  for (int g = 0; g < G; g++) {
    for (int j = 0; j < UDSP_W; j++) {
      const int k = g * UDSP_W + j;
      const UDSP_EQ_STAGE z = {0}, *s = k < m ? sec[k].st : &z;
      a1[g][j] = s->a1;
      a2[g][j] = s->a2;
      a3[g][j] = s->a3;
      m0[g][j] = s->m0;
      m1[g][j] = s->m1;
      m2[g][j] = s->m2;
      lane[g][j] = k;
      for (int c = 0; c < 2; c++) {
        s1[c][g][j] = k < m ? sec[k].ic1[c] : 0.0f;
        s2[c][g][j] = k < m ? sec[k].ic2[c] : 0.0f;
        out[c][g][j] = 0.0f;
      }
    }
  }
  const int last_g = (m - 1) / UDSP_W, last_j = (m - 1) % UDSP_W;
  const size_t steps = frames + (size_t)m - 1;
  /* the edges' steps, then the middle's (every lane its own frame) */
  size_t n = 0;
  const size_t mid = frames > (size_t)m - 1 ? (size_t)m - 1 : frames;
  for (; n < mid; n++) {
    wave_step(a1, a2, a3, m0, m1, m2, s1, s2, out, lane, io, frames, n, m,
              last_g, last_j, 1, G);
  }
  for (; n < frames; n++) {
    wave_step(a1, a2, a3, m0, m1, m2, s1, s2, out, lane, io, frames, n, m,
              last_g, last_j, 0, G);
  }
  for (; n < steps; n++) {
    wave_step(a1, a2, a3, m0, m1, m2, s1, s2, out, lane, io, frames, n, m,
              last_g, last_j, 1, G);
  }
  for (int k = 0; k < m; k++) {
    for (int c = 0; c < 2; c++) {
      sec[k].ic1[c] = s1[c][k / UDSP_W][k % UDSP_W];
      sec[k].ic2[c] = s2[c][k / UDSP_W][k % UDSP_W];
    }
  }
}

/* the sections gathered so far, one wavefront (or a few: CHAIN_LANES at
 * a time) */
static void run_sections(const SECTION *sec, int m, float *io, size_t frames) {
  while (m > 0) {
    const int k = m < CHAIN_LANES ? m : CHAIN_LANES;
    const int g = (k + UDSP_W - 1) / UDSP_W;
    if (g == 1) {
      wavefront(sec, k, io, frames, 1);
    } else if (g == 2) {
      wavefront(sec, k, io, frames, 2);
#if CHAIN_MAX_G > 2
    } else if (g == 3) {
      wavefront(sec, k, io, frames, 3);
#endif
    } else {
      wavefront(sec, k, io, frames, CHAIN_MAX_G);
    }
    sec += k;
    m -= k;
  }
}
#endif

void udsp_eq_chain(UDSP_EQ_BAND *const *bands, int n, float *io,
                   size_t frames) {
#if UDSP_W > 1
  if (frames == 0) {
    return;
  }
  SECTION sec[CHAIN_LANES];
  int m = 0;
  for (int i = 0; i < n; i++) {
    UDSP_EQ_BAND *b = bands[i];
    if (b->ramp) {
      run_sections(sec, m, io, frames);
      m = 0;
      udsp_eq_stereo(b, io, frames);
      continue;
    }
    for (int k = 0; k < b->stages; k++) {
      if (m == CHAIN_LANES) {
        run_sections(sec, m, io, frames);
        m = 0;
      }
      sec[m].st = &b->st[k];
      sec[m].ic1 = b->ic1[k];
      sec[m].ic2 = b->ic2[k];
      m++;
    }
  }
  run_sections(sec, m, io, frames);
#else
  for (int i = 0; i < n; i++) {
    udsp_eq_stereo(bands[i], io, frames);
  }
#endif
}

/* The trapezoidal SVF is its analog prototype at a warped frequency:
 * H = m0 + (m1 s + m2) / (s^2 + k s + 1), s = j tan(pi f / fs) / g */
float udsp_eq_response(const UDSP_EQ_BAND *b, float hz, unsigned rate) {
  float t = hz / (float)rate;
  t = t < 0.0f ? 0.0f : (t > 0.4999f ? 0.4999f : t);
  const double w = tan_rate(t);
  double mag = 1.0;
  for (int k = 0; k < b->stages; k++) {
    const UDSP_EQ_STAGE *s = &b->st[k];
    const double x = w / s->g; /* s = j x */
    /* denominator 1 - x^2 + j k x, numerator m2 + j m1 x */
    const double dr = 1.0 - x * x, di = s->k * x;
    const double nr = s->m2, ni = s->m1 * x;
    const double den = dr * dr + di * di;
    /* m0 + n / d */
    const double hr = s->m0 + (nr * dr + ni * di) / den;
    const double hi = (ni * dr - nr * di) / den;
    mag *= sqrt(hr * hr + hi * hi);
  }
  return (float)mag;
}
