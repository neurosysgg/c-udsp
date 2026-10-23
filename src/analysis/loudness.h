#ifndef UDSP_LOUDNESS_H
#define UDSP_LOUDNESS_H

#include <stddef.h>

/* floor for udsp_loudness_lufs_momentary() below -- LUFS has no natural
 * "meter bottom" the way dBFS's 0dBFS ceiling gives it one, so this picks
 * a floor low enough that real silence always reads as "the bottom of
 * the bar" rather than a wildly negative number. */
#define UDSP_LUFS_FLOOR -70.0f

/* Simplified, momentary-only loudness estimate loosely following ITU-R
 * BS.1770's K-weighting + mean-square approach -- see DESIGN.md for
 * exactly which parts of the real spec this does and doesn't implement
 * (short version: the K-weighting filter cascade and the -0.691 + 10*log10
 * formula are real; the integrated, gated loudness BS.1770 actually
 * defines is not attempted, and the filter coefficients are the spec's
 * fixed 48kHz design applied at whatever the caller's real sample rate
 * is, not re-derived per rate). Good enough for a live "how loud does
 * this sound" meter bar; not a certified loudness measurement.
 *
 * mono is n samples of mono audio (see udsp_mono_downmix()). scratch is
 * no longer used (the filter runs in one pass now; it was an n-float
 * working buffer) and may be NULL; the parameter stays for the callers.
 * Returns LUFS, clamped to UDSP_LUFS_FLOOR. */
float udsp_loudness_lufs_momentary(const float *mono, size_t n, float *scratch);

#endif
