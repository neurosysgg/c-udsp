#include "onset.h"

#include "synth/pitch.h"

/* a one-pole's coefficient for a time constant of ms: exp(-1 / frames),
 * without libm (the same bits everywhere) */
static float pole_of(float rate, float ms) {
  const double frames = (double)rate * ms / 1000.0;
  return frames <= 0.0 ? 0.0f : (float)udsp_exp2(-1.4426950408889634 / frames);
}

/* dB -> linear, without libm */
static float lin(float db) {
  return (float)udsp_exp2(db * (3.321928094887362 / 20.0));
}

void udsp_onset_set(UDSP_ONSET *o, float rate, float threshold_db,
                    float rise_db, float guard_ms) {
  o->down = pole_of(rate, 20.0f);
  o->pole = pole_of(rate, 20.0f);
  o->thr = lin(threshold_db);
  o->rise = lin(rise_db);
  const double g = (double)rate * guard_ms / 1000.0;
  o->guard = g < 1.0 ? 1u : g > 1e9 ? 1000000000u : (uint32_t)g;
}

void udsp_onset_init(UDSP_ONSET *o, float rate) {
  o->fast = o->slow = 0.0f;
  o->since = UINT32_MAX;
  udsp_onset_set(o, rate, -30.0f, 6.0f, 40.0f);
}

size_t udsp_onset_run(UDSP_ONSET *o, const float *key, size_t frames,
                      uint32_t *out, size_t cap) {
  size_t n = 0;
  float fast = o->fast, slow = o->slow;
  uint32_t since = o->since;
  for (size_t i = 0; i < frames; i++) {
    const float l = key[2 * i] < 0.0f ? -key[2 * i] : key[2 * i];
    const float r = key[2 * i + 1] < 0.0f ? -key[2 * i + 1] : key[2 * i + 1];
    const float x = l > r ? l : r;
    const float before = slow;
    const int up = x > fast;
    fast = up ? x : fast * o->down;
    slow = fast + (slow - fast) * o->pole;
    since = since == UINT32_MAX ? since : since + 1;
    if (up && fast > o->thr && fast > o->rise * before && since >= o->guard) {
      if (n < cap) {
        out[n] = (uint32_t)i;
      }
      n++;
      since = 0;
    }
  }
  o->fast = fast;
  o->slow = slow;
  o->since = since;
  return n < cap ? n : cap;
}
