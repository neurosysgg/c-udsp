#ifndef UDSP_DITHER_H
#define UDSP_DITHER_H

#include <stddef.h>
#include <stdint.h>

/* TPDF dither for a float signal about to be quantized to `bits` bits
 * (an export at 16 or 24): the sum of two independent uniform values of
 * half a step each way -- a triangle over one step either side, which
 * leaves the quantization error independent of the signal (no
 * distortion on a quiet fade, only a constant hiss). The noise comes
 * from a xorshift32 whose state is the caller's: the same seed, the same
 * bits on every build. */

typedef struct {
  uint32_t state;
} UDSP_DITHER;

/* seed 0 is taken as 1 (a xorshift's state can't be 0) */
void udsp_dither_init(UDSP_DITHER *d, uint32_t seed);

/* n samples (any channels, interleaved: each sample its own noise) in
 * place, for a quantizer of bits bits over [-1, 1) (a step 2^-(bits-1)) */
void udsp_dither_tpdf(UDSP_DITHER *d, float *buf, size_t n, int bits);

#endif
