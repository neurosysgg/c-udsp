#ifndef UDSP_STRETCH_H
#define UDSP_STRETCH_H

#include <stddef.h>
#include <stdint.h>

/* A phase vocoder's time stretch, offline: a signal made longer or
 * shorter at its own pitch (read faster afterwards, its pitch moved at
 * its length: ctui-daw's convolver's response). Hann frames of win
 * samples, the synthesis hop win / 4, each analysis frame where the
 * stretch puts it; identity phase locking (Laroche and Dolson): each
 * spectral peak's phase advanced by its own frequency, the bins around
 * it rotated with it, so a partial keeps its shape. A transient smears
 * over about a window (a response's first reflection a little); noise --
 * a reverb's tail -- stretches cleanly. Its angles are its own
 * arithmetic (no libm's atan2 / sin / cos), so it's the same bits on
 * every build and libc. */

/* the window for a rate: the power of two nearest 21 ms (1024 at 44.1
 * and 48 kHz) */
size_t udsp_stretch_window(uint32_t rate);

/* the frames n stretched by factor give: n * factor rounded, at least 1 */
size_t udsp_stretch_frames(size_t n, double factor);

/* the scratch floats udsp_stretch() needs for win */
size_t udsp_stretch_scratch(size_t win);

/* x's n frames (at x[i * stride]) stretched by factor (1/8 .. 8) into
 * udsp_stretch_frames(n, factor) frames at y[i * ystride] */
void udsp_stretch(const float *x, size_t n, size_t stride, double factor,
                  size_t win, float *y, size_t ystride, float *scratch);

#endif
