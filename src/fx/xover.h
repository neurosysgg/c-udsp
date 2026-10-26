#ifndef UDSP_XOVER_H
#define UDSP_XOVER_H

#include <stddef.h>

/* A Linkwitz-Riley crossover, stereo: a signal split into up to
 * UDSP_XOVER_BANDS bands at 24 dB an octave -- each split two
 * Butterworth sections of Simper's trapezoidal SVF (eq.h's) each way, so
 * the low and high halves meet at -6 dB in phase and sum to an allpass.
 * The splits go up from the lowest, each on what's above the one before;
 * a lower band then goes through the allpass of every split above it,
 * so all the bands sum to the input through those allpasses: flat in
 * level, the phase turned (Ableton's Multiband Dynamics', OTT's). Retuned
 * per block without thumps (the trapezoidal SVF's state survives a
 * coefficient change). No libm.
 *
 * Or at 12 dB an octave (udsp_xover_set_lr2: Linkwitz-Riley's second
 * order, LFOTool's split): one section each way at Q 1/2 -- the low half
 * its low pass, the high half its high pass inverted, so the two meet at
 * -6 dB in phase and sum to a first-order allpass; the lower bands go
 * through that allpass of each split above (low minus high). */

#define UDSP_XOVER_BANDS 5

typedef struct {
  float k, a1, a2, a3;
} UDSP_XOVER_SPLIT;

typedef struct {
  int bands, lr2;
  UDSP_XOVER_SPLIT c[UDSP_XOVER_BANDS - 1];
  /* each split's three sections (the first on its input, then the low
   * and the high side's second) and each band's allpasses (band b's at
   * split s > b): ic1, ic2 a channel */
  float sec[UDSP_XOVER_BANDS - 1][3][2][2];
  float ap[UDSP_XOVER_BANDS - 1][UDSP_XOVER_BANDS - 1][2][2];
} UDSP_XOVER;

/* zeroes the state (the coefficients are left alone) */
void udsp_xover_reset(UDSP_XOVER *x);

/* bands 1 .. UDSP_XOVER_BANDS, split at hz[0 .. bands - 2] (ascending;
 * each clamped to 10 Hz .. 0.49 of the rate) */
void udsp_xover_set(UDSP_XOVER *x, int bands, const float *hz, unsigned rate);

/* the same at 12 dB an octave (Linkwitz-Riley's second order) */
void udsp_xover_set_lr2(UDSP_XOVER *x, int bands, const float *hz,
                        unsigned rate);

/* in (frames interleaved stereo frames) into out[0 .. bands - 1], the
 * lowest band first; in may be out[bands - 1] */
void udsp_xover_split(UDSP_XOVER *x, const float *in, float *const *out,
                      size_t frames);

/* The Linkwitz-Riley low pass at linear phase, an FIR: the bilinear
 * Butterworth's (second order, warped to hz) magnitude squared -- 1 / (1 +
 * (tan(pi f / rate) / tan(pi hz / rate))^4), a filter run forward and
 * backward -- at zero phase, sampled on an n-point grid (n a power of
 * two, table udsp_fft_table()'s for it, scratch 2 n + 2 floats) and the
 * taps (odd) about its middle cut out, centred on (taps - 1) / 2: iZotope
 * Ozone 11's crossovers (its Imager's, its Dynamics' "type 2"; ctui-wm's
 * notes/vst-replace.md, measured within 0.001-0.1 dB). Bands made of such
 * low passes' differences (the first one, then each minus the one below,
 * the centre's impulse minus the last) sum to that delay exactly. n wants
 * to outlast the response: 32768 holds 20 Hz at 48 kHz to -100 dB. */
void udsp_xover_fir(float *h, size_t taps, double hz, unsigned rate, size_t n,
                    const float *table, float *scratch);

#endif
