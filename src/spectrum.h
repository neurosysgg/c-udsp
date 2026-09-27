#ifndef UDSP_SPECTRUM_H
#define UDSP_SPECTRUM_H

#include <stddef.h>

/* FFT size for every spectrum-bar computation in this library -- fixed
 * (not derived from sample rate/caller timing) and must stay a power of
 * two (udsp_fft_forward()'s precondition). 8192 samples is ~186ms at
 * 44.1kHz -- still short enough to feel responsive in a live visualizer,
 * and gives ~5.4Hz bin resolution, plenty of detail even in the sub-bass
 * end of a log-spaced bar layout. A caller with different
 * responsiveness/resolution needs can fork this constant locally; it
 * isn't a runtime argument because every buffer this module touches
 * (scratch_re/scratch_im, mono_tail) is sized from it at compile time. */
#define UDSP_FFT_N 8192

/* lowest frequency any bar/bin edge is ever placed at -- below this, bar
 * width per the log scale collapses to near-nothing anyway (see
 * udsp_spectrum_bars()'s doc comment), so this is also the log scale's
 * effective bottom. */
#define UDSP_SPECTRUM_MIN_HZ 20.0f

/* default bar count for a "spectrum bars" style visualizer -- just a
 * reasonable default; udsp_spectrum_bars() itself takes n_bars as a
 * runtime argument, so nothing stops a caller from asking for a
 * different count (e.g. many more, narrower bins for a spectrogram-style
 * waterfall). */
#define UDSP_SPECTRUM_BARS 128

/* mono_tail must point at exactly UDSP_FFT_N contiguous, most-recent
 * mono samples. scratch_re/scratch_im must each be UDSP_FFT_N floats:
 * caller-owned working buffers, clobbered on every call. mag_out gets
 * UDSP_FFT_N/2 magnitude values (bin 0 is DC) -- it may alias scratch_re,
 * same as udsp_fft_magnitude()'s own in-place contract, which is exactly
 * how udsp_spectrum_bars() below uses it.
 *
 * This is the windowing+FFT half of what udsp_spectrum_bars() does,
 * split out so a caller that needs the same FFT tail binned more than
 * one way (e.g. a spectrum-bars widget and a spectrogram column reading
 * the same window) runs it once, not once per binning -- see
 * udsp_spectrum_bin() for the other half. */
void udsp_spectrum_magnitude(const float *mono_tail, float *scratch_re,
                             float *scratch_im, float *mag_out);

/* mag must be a UDSP_FFT_N/2-length magnitude spectrum, as produced by
 * udsp_spectrum_magnitude(). Writes n_bars dBFS values into bars_out, one
 * per logarithmically-spaced frequency band from UDSP_SPECTRUM_MIN_HZ to
 * sample_rate/2 (Nyquist) -- log spacing because that's how both music's
 * harmonic content and human pitch perception are actually distributed,
 * unlike a linear FFT bin index (a linear layout would spend most of its
 * bars on treble ultrasonics nobody needs and almost none on the octaves
 * below ~500Hz where most of a mix's energy lives). Each bar's value is
 * the loudest (not averaged -- a quiet wide band shouldn't get diluted
 * by a few silent bins) bin magnitude within its band, converted via
 * udsp_analyze_dbfs(). Stateless -- one call, one answer -- so the same
 * function backs both a "few wide bars" widget and a "many narrow bins"
 * spectrogram column with nothing but a different n_bars. */
void udsp_spectrum_bin(const float *mag, int sample_rate, float *bars_out,
                       int n_bars);

/* udsp_spectrum_magnitude() + udsp_spectrum_bin() in one call -- the
 * convenience shape for a caller that only needs one binning of one FFT
 * (most callers; see udsp_spectrum_magnitude()'s own doc comment for
 * when the split versions are worth calling directly instead). */
void udsp_spectrum_bars(const float *mono_tail, int sample_rate,
                        float *scratch_re, float *scratch_im, float *bars_out,
                        int n_bars);

/* --- the planned path: udsp_spectrum_bars() with nothing recomputed ---
 *
 * Everything that only depends on the FFT size, the rate and the bar count
 * (the window, the FFT's twiddles, each bar's bin range) is worked out once
 * into a caller-owned plan; a spectrum is then a real FFT that windows its
 * input as it loads it, squared magnitudes, and a max per band: ~20x
 * faster than udsp_spectrum_bars() at the same size, same bars
 * (bench/README.md). The FFT size is the plan's, not UDSP_FFT_N: a live
 * meter wants short windows (2048: ~43 ms at 48 kHz), a detailed display
 * long ones (8192). The plan is read-only after udsp_spectrum_plan():
 * share it freely. */

/* the largest FFT a plan takes, and the most bars it holds (more than a
 * few thousand is pointless on any display) */
#define UDSP_SPECTRUM_MAX_N 8192
#define UDSP_SPECTRUM_MAX_BARS 1024

/* floats of scratch the planned functions take for an n-point plan: bins
 * 0..n/2, real and imaginary halves, the second on a 64-byte boundary
 * past the first (a vector load never straddles two cache lines there
 * when the scratch itself is 64-byte aligned). UDSP_SPECTRUM_SCRATCH fits
 * any plan */
#define UDSP_SPECTRUM_SCRATCH_N(n) (2 * ((n) / 2 + 16))
#define UDSP_SPECTRUM_SCRATCH UDSP_SPECTRUM_SCRATCH_N(UDSP_SPECTRUM_MAX_N)

/* the tables first: in a plan that starts on a 64-byte boundary (the
 * caller's to arrange: aligned_alloc(), _Alignas) they do too. Sized for
 * UDSP_SPECTRUM_MAX_N; an n-point plan uses the first n window floats and
 * 2n table floats, so a smaller plan touches (and keeps in cache) only
 * that much */
typedef struct UDSP_SPECTRUM_PLAN {
  float window[UDSP_SPECTRUM_MAX_N]; /* Hann, times the bars' 2/n */
  float fft[2 * UDSP_SPECTRUM_MAX_N]; /* UDSP_FFT_TABLE_FLOATS(n) used */
  int bin0[UDSP_SPECTRUM_MAX_BARS], bin1[UDSP_SPECTRUM_MAX_BARS];
  int n;
  int sample_rate;
  int n_bars;
} UDSP_SPECTRUM_PLAN;

/* fills plan for n-point FFTs (a power of two, 16..UDSP_SPECTRUM_MAX_N;
 * anything else: UDSP_FFT_N) of sample_rate audio into n_bars (clamped to
 * 1..MAX_BARS) bars: the log-spaced bands udsp_spectrum_bin() uses at that
 * size */
void udsp_spectrum_plan(UDSP_SPECTRUM_PLAN *plan, int n, int sample_rate,
                        int n_bars);

/* the power spectrum of the plan->n samples at mono_tail, Hann-windowed:
 * afterwards scratch[k] (k = 0..n/2) holds bin k's squared amplitude,
 * scaled like the bars (a full-scale sine's bin near 1, i.e. 0 dBFS).
 * scratch: UDSP_SPECTRUM_SCRATCH_N(plan->n) floats. Returns scratch. For
 * bins rather than bars (band energies, a spectrogram of one's own), or
 * one spectrum binned several ways (udsp_spectrum_bin_planned()). */
float *udsp_spectrum_power(const UDSP_SPECTRUM_PLAN *plan,
                           const float *mono_tail, float *scratch);

/* plan's n_bars dBFS bars out of a power spectrum from
 * udsp_spectrum_power() (of any plan with the same n and rate): the
 * planned udsp_spectrum_bin() */
void udsp_spectrum_bin_planned(const UDSP_SPECTRUM_PLAN *plan,
                               const float *power, float *bars_out);

/* udsp_spectrum_power() + udsp_spectrum_bin_planned(): udsp_spectrum_bars()
 * through plan */
void udsp_spectrum_bars_planned(const UDSP_SPECTRUM_PLAN *plan,
                                const float *mono_tail, float *scratch,
                                float *bars_out);

#endif
