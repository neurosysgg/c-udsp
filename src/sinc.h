#ifndef UDSP_SINC_H
#define UDSP_SINC_H

#include <stddef.h>
#include <stdint.h>

/* An offline resampler: a Kaiser-windowed sinc (120 dB down past the
 * lower rate's Nyquist), polyphase over the rates' ratio reduced to
 * up / down -- 44.1 -> 48 kHz is 160 phases, each its own row of taps
 * in a caller-owned table, normalised so a constant passes exactly.
 * A ratio with more than UDSP_SINC_PHASES_MAX phases (an odd rate:
 * 47999 Hz) interpolates linearly between that many rows instead.
 *
 * For audio imported at another rate than a project's (a DAW's clips:
 * converted once, not as they play). Its zero crossings each side set
 * the transition band: UDSP_SINC_ZEROS puts the passband's edge at
 * 0.4695 of the lower rate (20.7 kHz at 44.1). The table is built
 * without libm (its own sine and I0), and the dot product is four
 * lanes in a fixed order whatever the target, so every build gives the
 * same bits. The input comes in windows -- one channel, any part of it
 * -- and the output is the same however it's cut. */

#define UDSP_SINC_ZEROS 128
#define UDSP_SINC_PHASES_MAX 4096
/* rates above this are refused (keeps an output's index times down
 * within 64 bits for any length a file can have) */
#define UDSP_SINC_RATE_MAX 1000000u

typedef struct {
  uint32_t up, down; /* dst / gcd, src / gcd */
  uint32_t taps;     /* per row, a multiple of 4 */
  uint32_t rows;     /* up, or UDSP_SINC_PHASES_MAX + 1 interpolated */
  int interp;
  const float *table; /* rows * taps */
} UDSP_SINC;

/* the floats udsp_sinc_init's table needs; 0 for a rate of 0 or past
 * UDSP_SINC_RATE_MAX, or zeros 0 */
size_t udsp_sinc_floats(uint32_t src_rate, uint32_t dst_rate, uint32_t zeros);

/* fills table (udsp_sinc_floats of it) and s; 0 for what floats
 * refuses. The same rates make a single row passing the input through. */
int udsp_sinc_init(UDSP_SINC *s, uint32_t src_rate, uint32_t dst_rate,
                   uint32_t zeros, float *table);

/* the output's frames for in_frames of input: ceil(in * up / down) */
uint64_t udsp_sinc_length(const UDSP_SINC *s, uint64_t in_frames);

/* the input frames [*lo, *hi) that outputs [first, first + n) read (lo
 * may be negative, hi past the input's end: those read as silence) */
void udsp_sinc_span(const UDSP_SINC *s, uint64_t first, size_t n, int64_t *lo,
                    int64_t *hi);

/* outputs [first, first + n) of one channel into out, from the input's
 * frames [in_first, in_first + in_n) at in; frames outside that window
 * read as silence -- so a window holding udsp_sinc_span's frames (or as
 * many of them as the input has) gives every output exactly */
void udsp_sinc_run(const UDSP_SINC *s, const float *in, int64_t in_first,
                   size_t in_n, uint64_t first, float *out, size_t n);

#endif
