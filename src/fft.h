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

#endif
