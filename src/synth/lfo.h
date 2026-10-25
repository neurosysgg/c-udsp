#ifndef UDSP_LFO_H
#define UDSP_LFO_H

#include <stdint.h>

/* A low-frequency oscillator for modulation, read at control rate: its
 * value now, then a number of frames on. A 32-bit phase as osc.h's (a
 * rate exact and drift-free); bipolar, -1..1. The random shape is a
 * sample and hold: a new xorshift value at every turn; the smooth one
 * goes from the last value to the new one over the turn, a smoothstep
 * (value noise: still at every turn's start, never past either value --
 * Kilohearts' Ensemble's voices wander so). */

typedef enum {
  UDSP_LFO_SINE,
  UDSP_LFO_TRIANGLE, /* 0 rising at phase 0, as the sine */
  UDSP_LFO_SAW_UP,   /* -1 .. 1 */
  UDSP_LFO_SAW_DOWN, /* 1 .. -1 */
  UDSP_LFO_SQUARE,   /* 1 the first half */
  UDSP_LFO_RANDOM,   /* held for a turn */
  UDSP_LFO_SMOOTH,   /* random values a turn apart, smoothstep between */
  UDSP_LFO_SHAPES,
} UDSP_LFO_SHAPE;

typedef struct {
  UDSP_LFO_SHAPE shape;
  uint32_t phase; /* a turn is 2^32 */
  uint32_t inc;   /* per frame: osc.h's udsp_osc_inc_hz */
  uint32_t rng;   /* never 0 */
  float held;     /* the random shape's value, the smooth one's next */
  float from;     /* the smooth one's last */
} UDSP_LFO;

void udsp_lfo_init(UDSP_LFO *l, UDSP_LFO_SHAPE shape, uint32_t seed);

/* the value at the phase */
float udsp_lfo_value(const UDSP_LFO *l);

/* frames on (a new random value for every turn passed) */
void udsp_lfo_advance(UDSP_LFO *l, uint32_t frames);

/* to a phase (a note's retrigger, a bar's start): a new random value
 * (the smooth shape's from where it is) */
void udsp_lfo_reset(UDSP_LFO *l, uint32_t phase);

#endif
