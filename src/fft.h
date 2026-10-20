#ifndef UDSP_FFT_H
#define UDSP_FFT_H

#include <stddef.h>

/* generic power-of-two FFT machinery -- no audio-domain concepts (no Hz,
 * no dB, no channels) live here at all, deliberately: this is the one
 * numerical primitive spectrum.c builds on, kept isolated so it stays a
 * pure, independently-verifiable transform (see DESIGN.md's pipeline
 * layering). Every function below trusts its n argument is already a
 * power of two -- callers only ever pass a compile-time constant
 * (UDSP_FFT_N in spectrum.c), so there is no runtime check, per this
 * project's "trust internal callers" convention. */

/* in-place radix-2 Cooley-Tukey, decimation-in-time. re/im are n-length
 * arrays holding the signal's real/imaginary parts on entry and its
 * (unnormalized) DFT on return. For a real-valued input, caller sets
 * im[i] = 0 for all i and only the first n/2+1 output bins are distinct
 * (the rest mirror them, standard real-FFT symmetry) -- spectrum.c is the
 * one place that relies on this. */
void udsp_fft_forward(float *re, float *im, size_t n);

/* multiplies in[i] by a Hann window and writes the result to out[i], for
 * i in [0, n) -- reduces spectral leakage from analyzing a finite, non-
 * periodic slice of a continuous signal. in and out may be the same
 * buffer. */
void udsp_window_hann(const float *in, float *out, size_t n);

/* out[i] = sqrt(re[i]^2 + im[i]^2) for i in [0, n) -- magnitude of each
 * complex bin, discarding phase (nothing downstream of spectrum.c ever
 * needs phase). */
void udsp_fft_magnitude(const float *re, const float *im, size_t n, float *out);

/* --- the table-driven transforms (the fast path) ---
 *
 * Everything udsp_fft_forward()/udsp_window_hann() recompute per call
 * (a twiddle recurrence per stage, a cos() per sample) comes from a
 * caller-owned table instead, filled once. The tables are read-only after
 * that: nothing here is state, one table serves any number of callers and
 * threads. */

/* floats a twiddle table for an n-point transform holds */
#define UDSP_FFT_TABLE_FLOATS(n) (2 * (n))

/* fills table (UDSP_FFT_TABLE_FLOATS(n) floats) for transforms of size n.
 * Each stage's twiddles are contiguous (stage h's e^(-i*pi*k/h) at [h + k]
 * and [n + h + k]), so a butterfly loads them as plain vectors, and each
 * is computed directly (no recurrence: exact to float rounding). */
void udsp_fft_table(float *table, size_t n);

/* udsp_fft_forward() with twiddles from table (filled for n) */
void udsp_fft_forward_table(float *re, float *im, size_t n, const float *table);

/* the DFT of n real samples (n >= 4, a power of two), through an n/2-point
 * complex transform: about half the work of udsp_fft_forward() on the same
 * input with im = 0. Writes bins 0..n/2 (n/2 + 1 of them: the rest mirror
 * them) to re/im, which each hold n/2 + 1 floats and mustn't overlap in.
 * table is filled for n. Same values as udsp_fft_forward() gives those
 * bins, to float rounding. */
void udsp_fft_real(const float *in, float *re, float *im, size_t n,
                   const float *table);

/* udsp_fft_real() of in[i] * w[i] (a window table): the window applied as
 * the transform loads its input, so the windowed samples are never stored
 * -- one pass and n floats of memory less than udsp_window_apply() first */
void udsp_fft_real_windowed(const float *in, const float *w, float *re,
                            float *im, size_t n, const float *table);

/* udsp_fft_real()'s inverse: bins 0..n/2 (re/im, n/2 + 1 each; im[0] and
 * im[n/2] read as 0) back to n real samples in out, unnormalized -- n
 * times the signal, so a forward then an inverse gives n * in (a caller
 * folds the 1/n into what it multiplies the spectrum by). re/im are used
 * as scratch and left changed. table is filled for n; n >= 4. */
void udsp_fft_real_inverse(float *re, float *im, float *out, size_t n,
                           const float *table);

/* out[i] = re[i]^2 + im[i]^2: the squared magnitude, which is all a peak
 * search, a band energy or a dB value needs (10*log10 of it), without a
 * sqrt per bin. out may be re or im. */
void udsp_fft_power(const float *re, const float *im, size_t n, float *out);

/* fills w (n floats) with udsp_window_hann()'s window */
void udsp_window_hann_table(float *w, size_t n);

/* out[i] = in[i] * w[i]: any window table applied. in and out may be the
 * same buffer */
void udsp_window_apply(const float *in, const float *w, float *out, size_t n);

#endif
