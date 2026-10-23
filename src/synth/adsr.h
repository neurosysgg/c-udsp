#ifndef UDSP_ADSR_H
#define UDSP_ADSR_H

#include <stddef.h>

/* An ADSR envelope, per frame: a linear attack to 1 (a synth's punch:
 * the time given from 0, less from where a retrigger finds it), then
 * exponential decay to the sustain level and release to 0 -- the analog
 * shape, each segment covering its time from full scale down to
 * UDSP_ADSR_FLOOR (-80 dB) of its distance, then landing on its target.
 * A gate on attacks from wherever the value is, a gate off releases from
 * wherever it is: no jumps, so no clicks. The coefficients come from
 * udsp_exp2 (no libm): the same bits on every build. */

#define UDSP_ADSR_FLOOR 1e-4f

typedef enum {
  UDSP_ADSR_IDLE, /* 0, done */
  UDSP_ADSR_ATTACK,
  UDSP_ADSR_DECAY,
  UDSP_ADSR_SUSTAIN,
  UDSP_ADSR_RELEASE,
} UDSP_ADSR_STAGE;

typedef struct {
  float value;
  UDSP_ADSR_STAGE stage;
  float attack;  /* the step per frame */
  float decay;   /* the distance's factor per frame */
  float release; /* the same */
  float sustain;
} UDSP_ADSR;

/* the times in ms (0: a frame), sustain 0..1; the stage and value kept,
 * so a change while it runs takes over from there */
void udsp_adsr_set(UDSP_ADSR *e, float attack_ms, float decay_ms, float sustain,
                   float release_ms, unsigned rate);

/* note on (attack) or off (release) */
void udsp_adsr_gate(UDSP_ADSR *e, int on);

/* n values, a frame on each */
void udsp_adsr_render(UDSP_ADSR *e, float *out, size_t n);

/* n frames on, as udsp_adsr_render() would leave it, nothing written:
 * an envelope only read now and then (a held one costs nothing) */
void udsp_adsr_skip(UDSP_ADSR *e, size_t n);

static inline int udsp_adsr_idle(const UDSP_ADSR *e) {
  return e->stage == UDSP_ADSR_IDLE;
}

#endif
