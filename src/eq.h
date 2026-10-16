#ifndef UDSP_EQ_H
#define UDSP_EQ_H

#include <stddef.h>

/* A parametric EQ's band, stereo: Simper's trapezoidal SVF (zdf.h's)
 * with its outputs mixed into each shape (Simper's "SvfLinearTrap"
 * paper: out = m0 x + m1 band + m2 low), so a band retunes per block
 * under any modulation without the thumps a biquad's coefficient jump
 * makes. The cuts come at 12 dB an octave (one section, its Q the
 * band's) and 48 (four sections: an 8th-order Butterworth whose last,
 * sharpest section takes the band's Q as a resonance), Ableton's EQ
 * Eight's eight shapes. tan() is the oscillator's sine over its cosine
 * and the dB the pitch module's udsp_exp2: no libm (sqrt aside, which
 * IEEE rounds the same everywhere), the same bits on every build. */

typedef enum {
  UDSP_EQ_LOW_CUT_48,
  UDSP_EQ_LOW_CUT_12,
  UDSP_EQ_LOW_SHELF,
  UDSP_EQ_BELL,
  UDSP_EQ_NOTCH,
  UDSP_EQ_HIGH_SHELF,
  UDSP_EQ_HIGH_CUT_12,
  UDSP_EQ_HIGH_CUT_48,
  UDSP_EQ_SHAPES,
} UDSP_EQ_SHAPE;

#define UDSP_EQ_STAGES 4

typedef struct {
  float g, k;          /* the prototype's: for udsp_eq_response() */
  float a1, a2, a3;    /* the section's recurrence */
  float m0, m1, m2;    /* its mix */
} UDSP_EQ_STAGE;

typedef struct {
  UDSP_EQ_STAGE st[UDSP_EQ_STAGES];
  int stages;
  float ic1[UDSP_EQ_STAGES][2], ic2[UDSP_EQ_STAGES][2]; /* state */
  /* a ramp (udsp_eq_ramp): per frame steps of g, k, m0, m1, m2, the
   * frames left, and where it ends */
  float step[UDSP_EQ_STAGES][5];
  unsigned ramp;
  UDSP_EQ_STAGE end[UDSP_EQ_STAGES];
} UDSP_EQ_BAND;

/* zeroes the state and ends a ramp (the coefficients are left alone) */
void udsp_eq_reset(UDSP_EQ_BAND *b);

/* hz clamped to 10 Hz .. 0.49 of the rate, q to 0.1 .. 40 (0.71 a
 * Butterworth: no resonance on the cuts), gain in dB (the bell and the
 * shelves; +-30) */
void udsp_eq_set(UDSP_EQ_BAND *b, UDSP_EQ_SHAPE shape, float hz, float q,
                 float gain_db, unsigned rate);

/* glide to `to`'s coefficients (a udsp_eq_set() of its own) over the
 * next `frames` frames filtered: g, k and the mix move linearly a frame
 * at a time and the recurrence is worked out from them each frame, so a
 * retuned band neither steps (a coefficient jump's zipper) nor leaves
 * the filter's stable region on the way. A change of shape that changes
 * the sections' count jumps instead. Any ramp under way is replaced. */
void udsp_eq_ramp(UDSP_EQ_BAND *b, const UDSP_EQ_BAND *to, unsigned frames);

/* filters frames interleaved stereo frames in place */
void udsp_eq_stereo(UDSP_EQ_BAND *b, float *io, size_t frames);

/* the band's gain (linear, not dB) at hz: an EQ's curve on a screen */
float udsp_eq_response(const UDSP_EQ_BAND *b, float hz, unsigned rate);

#endif
