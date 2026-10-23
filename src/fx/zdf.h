#ifndef UDSP_ZDF_H
#define UDSP_ZDF_H

#include <stddef.h>

/* A zero-delay-feedback state-variable filter (Simper's trapezoidal
 * SVF), mono: a synth voice's filter. Unlike Chamberlin's (svf.h) it is
 * stable up to just under half the rate and under any modulation -- a
 * cutoff swept every few frames, resonance up to self-oscillation's
 * edge -- at a few more operations a sample. Its tan() is the
 * oscillator's polynomial sine over its cosine: no libm, the same bits
 * everywhere. Low-pass has unity gain at DC; at the cutoff every mode
 * but the notch has gain 1 / k (k = 2 - 2 res: -6 dB with no
 * resonance). */

typedef enum {
  UDSP_ZDF_LOW,
  UDSP_ZDF_BAND,
  UDSP_ZDF_HIGH,
  UDSP_ZDF_NOTCH,
} UDSP_ZDF_MODE;

typedef struct {
  float k, a1, a2, a3; /* coefficients: udsp_zdf_set() */
  UDSP_ZDF_MODE mode;
  float ic1, ic2; /* state */
} UDSP_ZDF;

/* zeroes the state (the coefficients are left alone) */
void udsp_zdf_reset(UDSP_ZDF *f);

/* cutoff in Hz (clamped to 10 Hz .. 0.49 of the rate), resonance 0..1
 * (1: Q 50, a ring short of oscillating) */
void udsp_zdf_set(UDSP_ZDF *f, float cutoff_hz, float res, UDSP_ZDF_MODE mode,
                  unsigned rate);

/* filters n samples in place */
void udsp_zdf_mono(UDSP_ZDF *f, float *io, size_t n);

#endif
