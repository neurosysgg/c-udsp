#ifndef UDSP_SHAPE_H
#define UDSP_SHAPE_H

#include <stddef.h>

/* Waveshapers: memoryless curves over a flat sample stream. */

/* x / (1 + |x|) of the input times drive (lifted from µOS's utracker,
 * fx.c): a rational soft clip, monotonic, +-1 in the limit, no tanh
 * (all but identical here at a quarter of the cost); the 2 / (1 + drive
 * / 4) makeup keeps it from being a volume knob too. drive >= 1. */
void udsp_softclip(float *buf, size_t n, float drive);

/* A saturator's curves, each passing a small signal as it is (slope 1
 * at 0) and none a volume knob of its own (no makeup: the caller's
 * output gain is):
 * - soft: tanh's shape, a rational (x (27 + x^2) / (27 + 9 x^2), Padé's)
 *   reaching +-1 at |x| = 3 and staying there;
 * - hard: clipped at +-1;
 * - fold: a sine, sin(x): past pi / 2 it folds back down (a wave
 *   folder's bright overtones);
 * - tube: asymmetric, x / (1 + x) above 0 and x / (1 - 2 x) under: the
 *   halves flatten unevenly (even harmonics, and a DC offset the caller
 *   filters out). */
typedef enum {
  UDSP_SHAPE_SOFT,
  UDSP_SHAPE_HARD,
  UDSP_SHAPE_FOLD,
  UDSP_SHAPE_TUBE,
  UDSP_SHAPES,
} UDSP_SHAPE;

float udsp_shape1(UDSP_SHAPE kind, float x);

/* buf[i] = shape(buf[i] * drive) */
void udsp_waveshape(float *buf, size_t n, UDSP_SHAPE kind, float drive);

/* a saturator's block with its knobs still: io[i] = (x * (1 - mix) +
 * shape(x * drive) * mix) * out, x = io[i] -- the same bits as that
 * expression over udsp_shape1(), a sample at a time */
void udsp_shape_mix(float *io, size_t n, UDSP_SHAPE kind, float drive,
                    float mix, float out);

#endif
