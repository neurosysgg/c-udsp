#ifndef UDSP_DELAY_H
#define UDSP_DELAY_H

#include <stddef.h>
#include <stdint.h>

/* A stereo feedback delay over a caller-owned ring (lifted from µOS's
 * utracker, fx.c, which kept a static pool of 500 ms rings: here the
 * caller sizes and owns it, e.g. UDSP_DELAY_RING_FLOATS(max_frames)). */

#define UDSP_DELAY_RING_FLOATS(max_frames) ((size_t)(max_frames) * 2)

typedef struct {
  float *ring;  /* 2 * cap floats */
  uint32_t cap; /* frames; the longest delay */
  uint32_t write;
} UDSP_DELAY;

/* takes ring (2 * cap floats) and clears it */
void udsp_delay_init(UDSP_DELAY *d, float *ring, uint32_t cap);
void udsp_delay_reset(UDSP_DELAY *d);

/* in place: out = dry + mix * the signal len frames ago (1..cap), which
 * was fed dry + feedback * itself */
void udsp_delay_stereo(UDSP_DELAY *d, float *io, size_t frames, uint32_t len,
                       float feedback, float mix);

#endif
