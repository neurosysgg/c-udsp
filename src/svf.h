#ifndef UDSP_SVF_H
#define UDSP_SVF_H

#include <stddef.h>

/* Chamberlin's state-variable filter, stereo (lifted from µOS's
 * utracker, fx.c): one structure for low-, high- and band-pass, cheap
 * to retune per block. Its frequency coefficient 2 sin(pi fc / fs) is
 * stable only well below fs / 4, so the cutoff is clamped there (past
 * it the topology oscillates instead of rolling off), and the states are
 * clamped to +-4 so a step into high resonance can't run away. The sine
 * is a 3-term polynomial, 0.01 % off at 10 kHz / 48 kHz: no libm. */

typedef enum {
  UDSP_SVF_LOW,
  UDSP_SVF_HIGH,
  UDSP_SVF_BAND,
} UDSP_SVF_MODE;

typedef struct {
  float f, q; /* coefficients: udsp_svf_set() */
  UDSP_SVF_MODE mode;
  float lp[2], bp[2]; /* state, per channel */
} UDSP_SVF;

/* zeroes the state (the coefficients are left alone) */
void udsp_svf_reset(UDSP_SVF *s);

/* cutoff in Hz, resonance 0 (none) .. 0.99 */
void udsp_svf_set(UDSP_SVF *s, float cutoff_hz, float res, UDSP_SVF_MODE mode,
                  unsigned rate);

/* filters frames interleaved stereo frames in place */
void udsp_svf_stereo(UDSP_SVF *s, float *io, size_t frames);

#endif
