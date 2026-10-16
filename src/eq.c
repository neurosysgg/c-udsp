#include "eq.h"

#include <math.h>

#include "osc.h"
#include "pitch.h"

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
    stage(&b->st[0], g * sqrtf(a), k, a * a, k * (1.0f - a) * a,
          1.0f - a * a);
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
static void ramp_stage(UDSP_EQ_BAND *b, int k, float *io, size_t n) {
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
void udsp_eq_stereo(UDSP_EQ_BAND *b, float *io, size_t frames) {
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
