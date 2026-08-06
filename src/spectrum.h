#ifndef UDSP_SPECTRUM_H
#define UDSP_SPECTRUM_H

#include <stddef.h>

/* FFT size for every spectrum-bar computation in this library -- fixed
 * (not derived from sample rate/caller timing) and must stay a power of
 * two (udsp_fft_forward()'s precondition). 2048 samples is ~46ms at
 * 44.1kHz -- short enough to feel responsive in a live visualizer, long
 * enough for ~21Hz bin resolution, i.e. usable detail even in the
 * sub-bass end of a log-spaced bar layout. A caller with different
 * responsiveness/resolution needs can fork this constant locally; it
 * isn't a runtime argument because every buffer this module touches
 * (scratch_re/scratch_im, mono_tail) is sized from it at compile time. */
#define UDSP_FFT_N 2048

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
#define UDSP_SPECTRUM_BARS 28

/* mono_tail must point at exactly UDSP_FFT_N contiguous, most-recent
 * mono samples. scratch_re/scratch_im must each be UDSP_FFT_N floats:
 * caller-owned working buffers, clobbered on every call, so this
 * function itself never allocates.
 *
 * Writes n_bars dBFS values into bars_out, one per logarithmically-spaced
 * frequency band from UDSP_SPECTRUM_MIN_HZ to sample_rate/2 (Nyquist) --
 * log spacing because that's how both music's harmonic content and human
 * pitch perception are actually distributed, unlike a linear FFT bin
 * index (a linear layout would spend most of its bars on treble
 * ultrasonics nobody needs and almost none on the octaves below ~500Hz
 * where most of a mix's energy lives). Each bar's value is the loudest
 * (not averaged -- a quiet wide band shouldn't get diluted by a few
 * silent bins) bin magnitude within its band, converted via
 * udsp_analyze_dbfs(). Stateless -- one call, one answer -- so the same
 * function backs both a "few wide bars" widget and a "many narrow bins"
 * spectrogram column with nothing but a different n_bars. */
void udsp_spectrum_bars(const float *mono_tail, int sample_rate,
                        float *scratch_re, float *scratch_im,
                        float *bars_out, int n_bars);

#endif
