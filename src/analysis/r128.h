#ifndef UDSP_R128_H
#define UDSP_R128_H

#include <stddef.h>
#include <stdint.h>

/* A program's integrated loudness (EBU R128 / ITU-R BS.1770-4): what a
 * ReplayGain 2.0 scanner measures a track or an album by. Unlike
 * loudness.h's momentary meter, this one is the full spec: the
 * K-weighting derived for the caller's sample rate, 400 ms blocks every
 * 100 ms, the channel weights, and both gates (absolute -70 LUFS,
 * relative -10 LU). Plus the sample peak, which ReplayGain's clipping
 * prevention wants alongside.
 *
 * Stateful, which DESIGN.md asks a reason for: a program's loudness is
 * over all of it, and a track decoded whole into memory would be the
 * allocation this library never makes. So the caller feeds it chunk by
 * chunk into a UDSP_R128 it owns, and reads the result at the end.
 *
 * The gated blocks go into a fixed histogram, not a list: a block's
 * energy is summed exactly into its 0.1 LU bin, so the only rounding is
 * which side of the relative gate the one bin straddling it falls on.
 * Two histograms add up, which is what an album's loudness is (the gates
 * over all its tracks' blocks together, not an average of the tracks'). */

#define UDSP_R128_MAX_CHANNELS 8
/* -70 to +10 LUFS in 0.1 LU bins; a louder block (only a test signal gets
 * there) counts in the top one, its energy still exact */
#define UDSP_R128_BINS 800

typedef struct {
  uint64_t count[UDSP_R128_BINS];
  double energy[UDSP_R128_BINS];
  float peak; /* the largest |sample|, any channel */
} UDSP_R128_HIST;

typedef struct {
  /* set by udsp_r128_init(), read-only after (weight may be changed
   * before the first udsp_r128_add(), for a layout the default misreads) */
  double k[2][5]; /* see udsp_r128_kweight() */
  float weight[UDSP_R128_MAX_CHANNELS];
  int channels;
  size_t hop; /* frames per 100 ms */

  double z[UDSP_R128_MAX_CHANNELS][4]; /* both stages' filter state */
  double sub[4];                       /* the last four hops' weighted sums */
  double acc;                          /* the current hop's, so far */
  size_t acc_n;
  uint64_t hops;
  UDSP_R128_HIST hist;
} UDSP_R128;

/* channels is 1..UDSP_R128_MAX_CHANNELS, in WAVE / FLAC / Vorbis order;
 * the weights follow BS.1770 for that count: 1.0 for the front ones,
 * 1.41 for the surrounds, 0 for an LFE (4: L R Ls Rs; 5: L R C Ls Rs;
 * 6: L R C LFE Ls Rs; 7: L R C LFE Cs Ls Rs; 8: L R C LFE Lb Rb Ls Rs). */
void udsp_r128_init(UDSP_R128 *m, int sample_rate, int channels);

/* frames of interleaved audio, in m's channel count */
void udsp_r128_add(UDSP_R128 *m, const float *interleaved, size_t frames);

/* the integrated loudness in LUFS, or -HUGE_VAL when no block passed the
 * absolute gate (silence, or under 400 ms of audio) */
double udsp_r128_integrated(const UDSP_R128_HIST *h);

/* dst += src: an album out of its tracks */
void udsp_r128_hist_add(UDSP_R128_HIST *dst, const UDSP_R128_HIST *src);

/* BS.1770's two K-weighting stages for sample_rate, as {b0, b1, b2, a1,
 * a2} each (a0 = 1): the spec's analog prototypes put through the
 * bilinear transform, which at 48 kHz gives its published table */
void udsp_r128_kweight(int sample_rate, double k[2][5]);

#endif
