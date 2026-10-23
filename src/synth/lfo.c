#include "synth/lfo.h"

#include "synth/osc.h"

#define TURNS (1.0f / 4294967296.0f)

static float next_random(UDSP_LFO *l) {
  uint32_t x = l->rng;
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  l->rng = x;
  return (float)(int32_t)x * (1.0f / 2147483648.0f);
}

void udsp_lfo_init(UDSP_LFO *l, UDSP_LFO_SHAPE shape, uint32_t seed) {
  l->shape = shape;
  l->phase = 0;
  l->inc = 0;
  l->rng = seed ? seed : 0x9e3779b9u;
  l->held = next_random(l);
}

float udsp_lfo_value(const UDSP_LFO *l) {
  const float t = (float)l->phase * TURNS;
  switch (l->shape) {
  case UDSP_LFO_TRIANGLE:
    return t < 0.25f ? 4.0f * t
                     : (t < 0.75f ? 2.0f - 4.0f * t : 4.0f * t - 4.0f);
  case UDSP_LFO_SAW_UP:
    return 2.0f * t - 1.0f;
  case UDSP_LFO_SAW_DOWN:
    return 1.0f - 2.0f * t;
  case UDSP_LFO_SQUARE:
    return l->phase < 0x80000000u ? 1.0f : -1.0f;
  case UDSP_LFO_RANDOM:
    return l->held;
  default:
    return udsp_osc_sine(t);
  }
}

void udsp_lfo_advance(UDSP_LFO *l, uint32_t frames) {
  uint64_t to = (uint64_t)l->phase + (uint64_t)frames * l->inc;
  if (to >> 32) {
    l->held = next_random(l);
  }
  l->phase = (uint32_t)to;
}

void udsp_lfo_reset(UDSP_LFO *l, uint32_t phase) {
  l->phase = phase;
  l->held = next_random(l);
}
