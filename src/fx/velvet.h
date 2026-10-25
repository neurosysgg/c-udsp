#ifndef UDSP_VELVET_H
#define UDSP_VELVET_H

#include <stddef.h>
#include <stdint.h>

#include "fx/line.h"

/* A velvet-noise filter (Alary, Politis and Välimäki's decorrelator, DAFx
 * 2017): a sparse FIR of +-g taps, each segment of its length cut into as
 * many equal cells as it has taps, one tap at a random place in each, of a
 * random sign. A mono signal through it, added to the side, widens it
 * without the comb a single delay leaves; the segments' gains stepping
 * down are the response's envelope. The draw is the seed's (xorshift32):
 * the same seed and segments, the same taps. */

#define UDSP_VELVET_MAX 64

typedef struct {
  uint32_t from, to; /* frames back, from < to */
  uint32_t taps;
  float gain;
} UDSP_VELVET_SEG;

typedef struct {
  uint32_t n;
  uint32_t at[UDSP_VELVET_MAX]; /* frames back, ascending */
  float g[UDSP_VELVET_MAX];
} UDSP_VELVET;

/* the taps of n segments in turn (past UDSP_VELVET_MAX in all, the rest
 * left out); a seed of 0 is taken as 1 */
void udsp_velvet_make(UDSP_VELVET *v, const UDSP_VELVET_SEG *seg, size_t n,
                      uint32_t seed);

/* the filter over the last n frames written to l (udsp_line_write): out[i]
 * the sum of g x[i - at], x[i] the i-th of them; l holds the last tap
 * plus n frames */
void udsp_velvet_read(const UDSP_VELVET *v, const UDSP_LINE *l, float *out,
                      size_t n);

#endif
