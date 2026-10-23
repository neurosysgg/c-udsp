#include "offline/dither.h"

void udsp_dither_init(UDSP_DITHER *d, uint32_t seed) {
  d->state = seed ? seed : 1u;
}

/* uniform in [-0.5, 0.5) from the top 24 bits of the next state */
static float uniform(UDSP_DITHER *d) {
  uint32_t x = d->state;
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  d->state = x;
  return (float)(x >> 8) * (1.0f / 16777216.0f) - 0.5f;
}

void udsp_dither_tpdf(UDSP_DITHER *d, float *buf, size_t n, int bits) {
  /* a step of a quantizer over [-1, 1): 2^-(bits - 1), exact in a float */
  float step = 1.0f;
  for (int b = 1; b < bits; b++) {
    step *= 0.5f;
  }
  for (size_t i = 0; i < n; i++) {
    const float a = uniform(d), b = uniform(d);
    buf[i] += (a + b) * step;
  }
}
