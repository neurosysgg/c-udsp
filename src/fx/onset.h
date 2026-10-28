#ifndef UDSP_ONSET_H
#define UDSP_ONSET_H

#include <stddef.h>
#include <stdint.h>

/* Onsets in a stereo signal, to the frame: a kick's start, for a ducker
 * keyed by it (ctui-daw's daw.duck, ctui-wm's notes/vst-replace.md
 * "The Duck"). Two followers: a fast one of the louder channel's |x| (up
 * at once, down over 20 ms: a 50 Hz tail's half cycle) and a slow one of
 * the fast one (a one-pole over 20 ms, both ways). An onset is a frame
 * where the fast one, rising, is over the threshold and over `rise`
 * times the slow one -- a jump from what the signal had been doing, not
 * its level (a kick's body, its tail and a held tone stay under it) --,
 * at least `guard` frames after the last. No allocation; the same bits at
 * any block size. */

typedef struct {
  float fast, slow; /* the followers */
  float down;       /* the fast one's release a frame */
  float pole;       /* the slow one's */
  float thr, rise;  /* linear */
  uint32_t guard;   /* frames */
  uint32_t since;   /* frames since the last onset (capped) */
} UDSP_ONSET;

/* at rest at rate: a threshold of -30 dBFS, a rise of 2 (6 dB), 40 ms
 * between onsets */
void udsp_onset_init(UDSP_ONSET *o, float rate);
/* the threshold (dBFS), the rise (dB) and the guard (ms) */
void udsp_onset_set(UDSP_ONSET *o, float rate, float threshold_db,
                    float rise_db, float guard_ms);

/* frames stereo frames of key (interleaved): the onsets' offsets into
 * out, at most cap of them (the rest still counted against the guard);
 * how many */
size_t udsp_onset_run(UDSP_ONSET *o, const float *key, size_t frames,
                      uint32_t *out, size_t cap);

#endif
