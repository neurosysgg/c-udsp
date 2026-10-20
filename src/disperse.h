#ifndef UDSP_DISPERSE_H
#define UDSP_DISPERSE_H

#include <stddef.h>

/* A disperser, stereo (Kilohearts' Disperser): up to UDSP_DISPERSE_MAX
 * second-order allpasses in series -- each Simper's trapezoidal SVF as
 * x - 2k band (eq.h's), so it retunes per block under modulation --
 * smearing the phase around their frequency: a transient spread into a
 * downward chirp, the level untouched. The amount is the stages, a
 * fraction fading the last one in (so it glides); the pinch their Q;
 * the spread fans their frequencies over that many octaves around the
 * centre, stage i always at the same place in the fan (a van der Corput
 * sequence), so more stages never move the ones there. No libm. */

#define UDSP_DISPERSE_MAX 128

typedef struct {
  float a1, a2, a3;
} UDSP_DISPERSE_STAGE;

typedef struct {
  UDSP_DISPERSE_STAGE st[UDSP_DISPERSE_MAX];
  float k;
  int whole;  /* stages wholly in */
  float frac; /* of the next one */
  float ic[UDSP_DISPERSE_MAX][2][2]; /* ic1, ic2 a channel */
  int used; /* stages whose state may not be at rest */
} UDSP_DISPERSE;

/* zeroes the state (the coefficients are left alone) */
void udsp_disperse_reset(UDSP_DISPERSE *d);

/* amount 0 .. UDSP_DISPERSE_MAX stages, the centre in Hz, pinch the
 * stages' Q (0.1 .. 40), spread in octaves (0 .. 8); every frequency
 * clamped to 10 Hz .. 0.49 of the rate */
void udsp_disperse_set(UDSP_DISPERSE *d, float amount, float hz, float pinch,
                       float spread, unsigned rate);

/* frames interleaved stereo frames in place */
void udsp_disperse_stereo(UDSP_DISPERSE *d, float *io, size_t frames);

#endif
