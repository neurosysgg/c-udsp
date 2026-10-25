#include "fx/velvet.h"

static uint32_t velvet_next(uint32_t *s) {
  uint32_t x = *s;
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  *s = x;
  return x;
}

void udsp_velvet_make(UDSP_VELVET *v, const UDSP_VELVET_SEG *seg, size_t n,
                      uint32_t seed) {
  uint32_t s = seed ? seed : 1;
  v->n = 0;
  for (size_t k = 0; k < n; k++) {
    const uint32_t span = seg[k].to - seg[k].from;
    for (uint32_t c = 0; c < seg[k].taps && v->n < UDSP_VELVET_MAX; c++) {
      /* cell c is [c span / taps, (c + 1) span / taps), in 64 bits: no
       * overflow at any length */
      const uint64_t lo = (uint64_t)c * span / seg[k].taps;
      const uint64_t hi = (uint64_t)(c + 1) * span / seg[k].taps;
      const uint32_t r = velvet_next(&s);
      const uint64_t w = hi > lo ? hi - lo : 1;
      v->at[v->n] = seg[k].from + (uint32_t)(lo + (uint64_t)r % w);
      v->g[v->n] = velvet_next(&s) & 1 ? -seg[k].gain : seg[k].gain;
      v->n++;
    }
  }
}

void udsp_velvet_read(const UDSP_VELVET *v, const UDSP_LINE *l, float *out,
                      size_t n) {
  for (size_t i = 0; i < n; i++) {
    const uint32_t back = (uint32_t)(n - i);
    float y = 0.0f;
    for (uint32_t t = 0; t < v->n; t++) {
      y += v->g[t] * udsp_line_tap(l, back + v->at[t]);
    }
    out[i] = y;
  }
}
