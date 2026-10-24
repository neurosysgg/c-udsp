#ifndef UDSP_IRFB_H
#define UDSP_IRFB_H

#include <stddef.h>

/* An impulse response fed back through itself, offline (Kilohearts'
 * Convolver's feedback; ctui-daw's convolver bakes it into its
 * response): what comes out of y = h * (x + g * y delayed by d frames)
 * for an impulse x -- H / (1 - g e^(-iwd) H) -- computed once in the
 * frequency domain over an n-point transform. The loop's gain is
 * g times |H| at each frequency, so it settles only where |g H| < 1:
 * udsp_irfb_peak() gives the largest |H| for a caller to scale g by.
 * The response's tail past the transform folds back onto its start; a
 * transform of twice the frames kept puts what folds at the level the
 * tail has decayed to by then (udsp_irfb_frames' -90 dB, or what the
 * cap leaves). The rotation e^(-iwd) is the transform's own twiddle
 * table, so it's the same bits on every build and libc. */

/* the frames the response of h (len frames) takes to fall 90 dB with a
 * loop gain of g (|g| < 1; at the peak frequency) through d frames, at
 * most cap */
size_t udsp_irfb_frames(size_t len, size_t d, float g, size_t cap);

/* the transform for a response of frames: the power of two at least
 * twice it (4 at least) */
size_t udsp_irfb_size(size_t frames);

/* the scratch floats udsp_irfb() and udsp_irfb_peak() need for an
 * n-point transform */
#define UDSP_IRFB_SCRATCH(n) (2 * ((n) + 2))

/* the largest |H| of h's len frames (at h[i * stride], len <= n) over an
 * n-point transform (table filled for n by udsp_fft_table) */
float udsp_irfb_peak(const float *h, size_t len, size_t stride, size_t n,
                     const float *table, float *scratch);

/* h (len frames at h[i * stride], len <= n) fed back with gain g through
 * d frames: the response's first frames frames (<= n) into y[i *
 * ystride]; table filled for n */
void udsp_irfb(const float *h, size_t len, size_t stride, float g, size_t d,
               float *y, size_t frames, size_t ystride, size_t n,
               const float *table, float *scratch);

#endif
