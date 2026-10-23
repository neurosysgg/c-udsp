#include "fx/disperse.h"

#include <string.h>

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

/* both channels in one loop: two independent chains overlap (svf.c's
 * finding, session 221) -- the same operations each, the same bits */
static void stage(const UDSP_DISPERSE_STAGE *s, float k, float (*ic)[2],
                  float *io, size_t frames, float mix) {
  const float a1 = s->a1, a2 = s->a2, a3 = s->a3, w = mix * 2.0f * k;
  float l1 = ic[0][0], l2 = ic[0][1], r1 = ic[1][0], r2 = ic[1][1];
  for (size_t i = 0; i < frames; i++) {
    const float xl = io[2 * i], xr = io[2 * i + 1];
    const float lv3 = xl - l2, rv3 = xr - r2;
    const float lv1 = a1 * l1 + a2 * lv3, rv1 = a1 * r1 + a2 * rv3;
    const float lv2 = l2 + a2 * l1 + a3 * lv3, rv2 = r2 + a2 * r1 + a3 * rv3;
    l1 = 2.0f * lv1 - l1;
    r1 = 2.0f * rv1 - r1;
    l2 = 2.0f * lv2 - l2;
    r2 = 2.0f * rv2 - r2;
    io[2 * i] = xl - w * lv1;
    io[2 * i + 1] = xr - w * rv1;
  }
  ic[0][0] = l1;
  ic[0][1] = l2;
  ic[1][0] = r1;
  ic[1][1] = r2;
}

void udsp_disperse_stereo(UDSP_DISPERSE *d, float *io, size_t frames) {
  for (int i = 0; i < d->whole; i++) {
    stage(&d->st[i], d->k, d->ic[i], io, frames, 1.0f);
  }
  if (d->frac > 0.0f && d->whole < UDSP_DISPERSE_MAX) {
    /* the next stage's output faded in over its input */
    stage(&d->st[d->whole], d->k, d->ic[d->whole], io, frames, d->frac);
  }
  /* a stage out keeps no state: in again, it starts from rest */
  const int n = d->whole + (d->frac > 0.0f && d->whole < UDSP_DISPERSE_MAX);
  if (d->used > n) {
    memset(d->ic[n], 0, (size_t)(d->used - n) * sizeof d->ic[0]);
  }
  d->used = n;
}
