#ifndef UDSP_SMOOTH_H
#define UDSP_SMOOTH_H

#include <stddef.h>
#include <stdint.h>

/* A parameter's value gliding to a new target in a fixed time, so a
 * knob's move or an automation step doesn't zipper (a gain jumping
 * clicks, a cutoff jumping thumps). Linear, and the same length whatever
 * the distance: a new target takes `frames` frames from wherever the
 * value is, and the last one lands on the target exactly. The state is
 * the caller's, as everywhere here. */

typedef struct {
  float value;
  float target;
  float step;
  uint32_t left;   /* frames until the target */
  uint32_t frames; /* a glide's length */
} UDSP_SMOOTH;

/* at value, settled; glides of `frames` (0: jumps) */
void udsp_smooth_init(UDSP_SMOOTH *s, float value, uint32_t frames);

/* glide to target from where it is */
void udsp_smooth_to(UDSP_SMOOTH *s, float target);

/* the value now, then a frame on */
static inline float udsp_smooth_next(UDSP_SMOOTH *s) {
  float v = s->value;
  if (s->left) {
    s->value = --s->left ? s->value + s->step : s->target;
  }
  return v;
}

/* n values (udsp_smooth_next's, in order) into out */
void udsp_smooth_render(UDSP_SMOOTH *s, float *out, size_t n);

/* a frame on without the value: a caller using it once a block skips
 * the rest */
void udsp_smooth_skip(UDSP_SMOOTH *s, size_t n);

static inline int udsp_smooth_settled(const UDSP_SMOOTH *s) {
  return s->left == 0;
}

#endif
