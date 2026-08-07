#ifndef UDSP_SPECTROGRAM_H
#define UDSP_SPECTROGRAM_H

/* fixed grid size for the scrolling waterfall -- a caller maps this onto
 * whatever cell/pixel size it actually renders at (nearest-neighbor
 * lookup, same as any other fixed-resolution-source/variable-resolution-
 * display mismatch), so this only needs to be "enough columns of history
 * to look like a waterfall, enough rows to look like a spectrum," not
 * tied to any specific display size. */
#define UDSP_SPECTROGRAM_COLS 96
#define UDSP_SPECTROGRAM_BINS 48

/* one column per udsp_spectrogram_push() call, oldest at index 0, newest
 * at COLS-1 -- intensity is 0-255, not a float dBFS value: quantizing
 * here rather than at render time keeps this struct an order of
 * magnitude smaller (COLS*BINS bytes vs. floats), useful for a caller
 * that copies/broadcasts it often (e.g. through a lock-free snapshot). If
 * a consumer ever needs the pre-quantized float value, add it as a
 * second field rather than changing this one back to float. */
typedef struct UDSP_SPECTROGRAM {
  unsigned char intensity[UDSP_SPECTROGRAM_COLS][UDSP_SPECTROGRAM_BINS];
} UDSP_SPECTROGRAM;

void udsp_spectrogram_init(UDSP_SPECTROGRAM *sg);

/* mono_tail/sample_rate/scratch_re/scratch_im: same contract as
 * udsp_spectrum_bars() (spectrum.h), which this calls with n_bars =
 * UDSP_SPECTROGRAM_BINS to get one new column, then scrolls the history
 * left by one column and appends it -- this is the one stateful module
 * in the library (a spectrogram is inherently "spectrum bars plus memory
 * of the recent past"), which is exactly why it's its own module built
 * on top of the stateless udsp_spectrum_bars() rather than a flag/mode on
 * it. */
void udsp_spectrogram_push(UDSP_SPECTROGRAM *sg, const float *mono_tail,
                           int sample_rate, float *scratch_re, float *scratch_im);

/* same as udsp_spectrogram_push(), but takes an already-computed
 * UDSP_FFT_N/2-length magnitude spectrum (udsp_spectrum_magnitude()'s
 * output) instead of a raw mono tail -- for a caller that's already run
 * the FFT for something else (e.g. a spectrum-bars widget) off the same
 * window and doesn't want to pay for a second windowing+FFT pass over
 * identical input just to get this module's own binning. */
void udsp_spectrogram_push_mag(UDSP_SPECTROGRAM *sg, const float *mag, int sample_rate);

#endif
