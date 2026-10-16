#ifndef UDSP_LINE_H
#define UDSP_LINE_H

#include <stddef.h>
#include <stdint.h>

/* A delay line, mono, over a caller-owned ring a power of two long, read
 * at any delay -- whole frames, or fractional through a cubic (Hermite)
 * between the four around it: what a delay whose time moves, a chorus
 * and a reverb's modulated lines are made of. Read, then push: a delay
 * of 1 is the frame pushed last. */

typedef struct {
  float *buf;
  uint32_t mask; /* the ring's length - 1 */
  uint32_t w;    /* where the next frame goes */
} UDSP_LINE;

/* the ring's length (floats) for delays up to max frames */
uint32_t udsp_line_size(uint32_t max);

/* takes buf (udsp_line_size()'s floats) and clears it */
void udsp_line_init(UDSP_LINE *l, float *buf, uint32_t size);
void udsp_line_clear(UDSP_LINE *l);

static inline void udsp_line_push(UDSP_LINE *l, float x) {
  l->buf[l->w] = x;
  l->w = (l->w + 1) & l->mask;
}

/* the frame d frames back, d >= 1 */
static inline float udsp_line_tap(const UDSP_LINE *l, uint32_t d) {
  return l->buf[(l->w - d) & l->mask];
}

/* Catmull-Rom between y0 and y1 at x (0..1), ym1 and y2 around them */
static inline float udsp_hermite(float ym1, float y0, float y1, float y2,
                                 float x) {
  const float c1 = 0.5f * (y1 - ym1);
  const float c2 = ym1 - 2.5f * y0 + 2.0f * y1 - 0.5f * y2;
  const float c3 = 0.5f * (y2 - ym1) + 1.5f * (y0 - y1);
  return ((c3 * x + c2) * x + c1) * x + y0;
}

/* the line d frames back, d >= 2 and fractional */
static inline float udsp_line_read(const UDSP_LINE *l, float d) {
  const uint32_t i = (uint32_t)d;
  const float f = d - (float)i;
  return udsp_hermite(udsp_line_tap(l, i - 1), udsp_line_tap(l, i),
                      udsp_line_tap(l, i + 1), udsp_line_tap(l, i + 2), f);
}

#endif
