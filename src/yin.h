#ifndef UDSP_YIN_H
#define UDSP_YIN_H

#include <stddef.h>
#include <stdint.h>

/* A note's fundamental frequency from its sound: YIN (de Cheveigne and
 * Kawahara, 2002) -- the squared difference between the signal and
 * itself delayed, normalized by its running mean, its first dip under a
 * threshold refined by a parabola. For a sampler mapping a folder of
 * multi-samples whose roots aren't written anywhere (a file's index for
 * its name): offline, on a stretch of the sustain, not per block.
 * Stateless, the caller's scratch, no allocation. */

/* the fundamental of x[0..n) in Hz, searched between fmin and fmax, or 0
 * when no dip goes under threshold (0.1-0.2: lower is stricter; noise, a
 * drum). n must hold the longest period twice: n >= 2 * rate / fmin.
 * scratch: rate / fmin + 2 floats. *aperiodicity (unless NULL): the
 * dip's depth, 0 perfectly periodic, toward 1 noise. */
float udsp_yin(const float *x, size_t n, uint32_t rate, float fmin, float fmax,
               float threshold, float *scratch, float *aperiodicity);

/* the MIDI note (69 = 440 Hz) nearest hz, and how far off it is in cents
 * (-50..50) unless NULL; -1 for hz <= 0 */
int udsp_yin_note(float hz, float *cents);

#endif
