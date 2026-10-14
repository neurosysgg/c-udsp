#ifndef UDSP_SAMPLER_H
#define UDSP_SAMPLER_H

#include <stddef.h>
#include <stdint.h>

/* A sample played back at a pitch (lifted from µOS's utracker,
 * sampler.c): read at a variable rate with linear interpolation, where
 * pitch and speed are one control (the classic sampler). The sample
 * stays at its own rate: the rate conversion is the same multiply the
 * pitch is, so doing it at load too would round twice. The voice's
 * position is Q32 fixed point, so a loop played for minutes doesn't
 * drift against anything else. Frames are the caller's (decoded however
 * it likes); the voice is a caller-owned struct. */

typedef struct {
  const float *frames; /* interleaved, channels per frame */
  uint32_t frame_count;
  uint32_t sample_rate; /* the sample's own */
  uint32_t channels;    /* 1 or 2 */
  uint32_t loop_start;  /* loop_end <= loop_start: one-shot */
  uint32_t loop_end;
  uint8_t root_note; /* the MIDI note that plays it unpitched */
} UDSP_SAMPLE;

typedef struct {
  const UDSP_SAMPLE *smp;
  uint64_t pos_q32;
  uint64_t step_q32;
  int playing;
} UDSP_SAMPLER_VOICE;

/* starts v on s at note, out_rate being the output's rate; a step past
 * 64 source frames per output frame is clamped (a malformed header's
 * rate, not a pitch) */
void udsp_sampler_note_on(UDSP_SAMPLER_VOICE *v, const UDSP_SAMPLE *s, int note,
                          uint32_t out_rate);

/* stops at once (a release is the caller's envelope) */
void udsp_sampler_stop(UDSP_SAMPLER_VOICE *v);

/* overwrites frames interleaved STEREO frames (a mono sample on both);
 * silence once a one-shot runs off its end, which clears playing */
void udsp_sampler_render(UDSP_SAMPLER_VOICE *v, float *out, size_t frames);

#endif
