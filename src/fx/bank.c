#include "fx/bank.h"

#include <string.h>

#include "base/simd.h"
#include "synth/osc.h"
#include "synth/pitch.h"

void udsp_bank_reset(UDSP_BANK *b) {
  memset(b->ic1, 0, sizeof b->ic1);
  memset(b->ic2, 0, sizeof b->ic2);
}

/* a stage's Q for a band of Q q out of order stages alike: a stage's
 * -3 dB width w makes order of them -3 dB at the width where each is
 * down 3 / order dB, sqrt(2^(1 / order) - 1) of w */
static float stage_q(float q, uint32_t order) {
  static const float WIDEN[UDSP_BANK_ORDER] = {1.0f, 0.643594f, 0.509825f,
                                               0.435017f};
  return q * WIDEN[order - 1];
}

/* the integral of (1 + x^2)^-order over x: sqrt(pi) G(order - 1/2) /
 * G(order) */
float udsp_bank_noise_width(float hz, float q, uint32_t order) {
  static const float AREA[UDSP_BANK_ORDER] = {3.14159265f, 1.57079633f,
                                              1.17809725f, 0.98174770f};
  order = order < 1 ? 1 : order > UDSP_BANK_ORDER ? UDSP_BANK_ORDER : order;
  return hz / (2.0f * stage_q(q > 0.05f ? q : 0.05f, order)) * AREA[order - 1];
}

void udsp_bank_set(UDSP_BANK *b, uint32_t n, const float *hz, const float *q,
                   const float *gain, uint32_t order, unsigned rate) {
  n = n < 1 ? 1 : n > UDSP_BANK_MAX ? UDSP_BANK_MAX : n;
  order = order < 1 ? 1 : order > UDSP_BANK_ORDER ? UDSP_BANK_ORDER : order;
  const uint32_t lanes = (n + UDSP_W - 1) / UDSP_W * UDSP_W;
  if (n != b->n || order != b->order) {
    udsp_bank_reset(b);
  }
  b->n = n;
  b->lanes = lanes;
  b->order = order;
  const float lo = 10.0f / (float)rate;
  for (uint32_t i = 0; i < UDSP_BANK_MAX; i++) {
    if (i >= n) {
      b->a1[i] = b->a2[i] = b->a3[i] = b->k[i] = b->gain[i] = 0.0f;
      continue;
    }
    float t = hz[i] / (float)rate;
    t = t < lo ? lo : (t > 0.49f ? 0.49f : t);
    const float g = udsp_osc_sine(0.5f * t) / udsp_osc_sine(0.5f * t + 0.25f);
    const float sq = stage_q(q[i] > 0.05f ? q[i] : 0.05f, order);
    const float k = 1.0f / sq;
    b->k[i] = k;
    b->a1[i] = 1.0f / (1.0f + g * (g + k));
    b->a2[i] = g * b->a1[i];
    b->a3[i] = g * b->a2[i];
    b->gain[i] = gain ? gain[i] : 1.0f;
  }
}

/* the prototype's k j x / (1 - x^2 + j k x), x = tan(pi hz / rate) / g
 * (g = a2 / a1), to the power of the stages */
void udsp_bank_response(const UDSP_BANK *b, uint32_t i, float hz, unsigned rate,
                        float *re, float *im) {
  float t = hz / (float)rate;
  t = t < 0.0f ? 0.0f : (t > 0.4999f ? 0.4999f : t);
  const double w =
      (double)udsp_osc_sine(0.5f * t) / (double)udsp_osc_sine(0.5f * t + 0.25f);
  if (i >= b->n || b->a1[i] == 0.0f) {
    *re = *im = 0.0f;
    return;
  }
  const double x = w * b->a1[i] / b->a2[i], k = b->k[i];
  const double dr = 1.0 - x * x, di = k * x, den = dr * dr + di * di;
  /* k j x / (dr + j di) */
  const double sr = k * x * di / den, si = k * x * dr / den;
  double hr = b->gain[i], hi = 0.0;
  for (uint32_t s = 0; s < b->order; s++) {
    const double r = hr * sr - hi * si;
    hi = hr * si + hi * sr;
    hr = r;
  }
  *re = (float)hr;
  *im = (float)hi;
}

/* one stage over a vector of bands: v1 (the band-pass) times k */
static inline UDSP_VF stage(UDSP_VF x, UDSP_VF a1, UDSP_VF a2, UDSP_VF a3,
                            UDSP_VF k, float *ic1p, float *ic2p) {
  const UDSP_VF ic1 = udsp_vload(ic1p), ic2 = udsp_vload(ic2p);
  const UDSP_VF v3 = x - ic2;
  const UDSP_VF v1 = a1 * ic1 + a2 * v3;
  const UDSP_VF v2 = ic2 + a2 * ic1 + a3 * v3;
  udsp_vstore(ic1p, v1 + v1 - ic1);
  udsp_vstore(ic2p, v2 + v2 - ic2);
  return k * v1;
}

/* a vector of bands' coefficients stay in registers over the frames;
 * each frame's stages in order */
void udsp_bank_run(UDSP_BANK *b, const float *x, size_t stride, float *out,
                   size_t frames) {
  const uint32_t lanes = b->lanes, order = b->order;
  for (uint32_t j = 0; j < lanes; j += UDSP_W) {
    const UDSP_VF a1 = udsp_vload(b->a1 + j), a2 = udsp_vload(b->a2 + j),
                  a3 = udsp_vload(b->a3 + j), k = udsp_vload(b->k + j),
                  gain = udsp_vload(b->gain + j);
    for (size_t f = 0; f < frames; f++) {
      UDSP_VF v = udsp_vsplat(x[f * stride]);
      for (uint32_t s = 0; s < order; s++) {
        v = stage(v, a1, a2, a3, k, b->ic1[s] + j, b->ic2[s] + j);
      }
      udsp_vstore(out + f * lanes + j, v * gain);
    }
  }
}

void udsp_follow_set(UDSP_FOLLOW *f, uint32_t lanes, float attack,
                     float release, float hold) {
  /* 1 - e^(-1 / t) */
  const double LOG2E = 1.4426950408889634;
  f->lanes = lanes;
  f->attack = attack < 1.0f ? 1.0f : (float)(1.0 - udsp_exp2(-LOG2E / attack));
  f->release =
      release < 1.0f ? 1.0f : (float)(1.0 - udsp_exp2(-LOG2E / release));
  f->hold = hold < 0.0f ? 0.0f : hold;
}

/* rising: held restarts; falling: held counts up, the release only once
 * it's past hold. Each lane's operations the scalar build's */
#if UDSP_W > 1
void udsp_follow_run(UDSP_FOLLOW *f, float *io, size_t frames) {
  const UDSP_VF att = udsp_vsplat(f->attack), rel = udsp_vsplat(f->release),
                hold = udsp_vsplat(f->hold), one = udsp_vsplat(1.0f),
                zero = udsp_vsplat(0.0f);
  const uint32_t lanes = f->lanes;
  for (uint32_t j = 0; j < lanes; j += UDSP_W) {
    UDSP_VF env = udsp_vload(f->env + j), held = udsp_vload(f->held + j);
    for (size_t i = 0; i < frames; i++) {
      float *p = io + i * lanes + j;
      const UDSP_VF r = udsp_vabs(udsp_vload(p));
      const UDSP_VI up = r > env;
      held = udsp_vsel(up, zero, held + one);
      const UDSP_VF step =
          udsp_vsel(up, att, udsp_vsel(held > hold, rel, zero));
      env += step * (r - env);
      udsp_vstore(p, env);
    }
    udsp_vstore(f->env + j, env);
    udsp_vstore(f->held + j, held);
  }
}
#else
void udsp_follow_run(UDSP_FOLLOW *f, float *io, size_t frames) {
  const uint32_t lanes = f->lanes;
  for (uint32_t j = 0; j < lanes; j++) {
    float env = f->env[j], held = f->held[j];
    for (size_t i = 0; i < frames; i++) {
      float *p = io + i * lanes + j;
      const float r = fabsf(*p);
      const int up = r > env;
      held = up ? 0.0f : held + 1.0f;
      const float step = up ? f->attack : held > f->hold ? f->release : 0.0f;
      env += step * (r - env);
      *p = env;
    }
    f->env[j] = env;
    f->held[j] = held;
  }
}
#endif
