#ifndef UDSP_TRUEPEAK_H
#define UDSP_TRUEPEAK_H

#include <stddef.h>

/* A signal's true peak (ITU-R BS.1770-4 annex 2): the largest value of the
 * waveform between its samples, read at 4x through a 48-tap polyphase
 * low-pass -- what a sample peak misses by up to 3 dB and a converter or an
 * MP3 encoder doesn't. Stateful as r128.c is: the filter's history runs
 * across blocks, and the same frames give the same peak at any block size. */

#define UDSP_TRUEPEAK_MAX_CHANNELS 2
#define UDSP_TRUEPEAK_PHASES 4
#define UDSP_TRUEPEAK_TAPS 12 /* a phase's */

typedef struct {
  float h[UDSP_TRUEPEAK_PHASES][UDSP_TRUEPEAK_TAPS];
  /* twice the taps: a window of them always contiguous */
  float z[UDSP_TRUEPEAK_MAX_CHANNELS][2 * UDSP_TRUEPEAK_TAPS];
  int at;
  int channels;
  float peak; /* the largest |value|, any channel, since the last reset */
} UDSP_TRUEPEAK;

/* channels 1..UDSP_TRUEPEAK_MAX_CHANNELS; the tables are worked out here,
 * so not on a realtime path */
void udsp_truepeak_init(UDSP_TRUEPEAK *t, int channels);
void udsp_truepeak_add(UDSP_TRUEPEAK *t, const float *interleaved,
                       size_t frames);

#endif
