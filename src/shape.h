#ifndef UDSP_SHAPE_H
#define UDSP_SHAPE_H

#include <stddef.h>

/* Waveshapers: memoryless curves over a flat sample stream. */

/* x / (1 + |x|) of the input times drive (lifted from µOS's utracker,
 * fx.c): a rational soft clip, monotonic, +-1 in the limit, no tanh
 * (all but identical here at a quarter of the cost); the 2 / (1 + drive
 * / 4) makeup keeps it from being a volume knob too. drive >= 1. */
void udsp_softclip(float *buf, size_t n, float drive);

#endif
