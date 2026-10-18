#ifndef UDSP_OPENV_H
#define UDSP_OPENV_H

#include <stddef.h>
#include <stdint.h>

/* Operator's envelope (ctui-daw's FM synth, d8): from an initial level
 * to a peak (attack), to a sustain level (decay), held, then to an end
 * level (release) -- each segment its time long and bent by its slope
 * (0 a line; towards 1 it moves late, towards -1 early: an exponential
 * of up to 2^SLOPE_RANGE between its ends), landing on its target
 * exactly. Retriggered it attacks from wherever it is (no click), idle
 * from its initial level. The loop modes:
 *
 *   hold     held at sustain until the gate's off
 *   loop     at sustain, back to the initial level over the loop time,
 *            then the attack again, until the gate's off
 *   trigger  the gate's end ignored: attack, decay, then the release
 *
 * Operator's beat and sync modes are hold plus udsp_openv_restart() at
 * the beats the caller counts (it knows the song's). Per frame, in
 * double where it accumulates and from udsp_exp2: the same bits on every
 * build. */

#define UDSP_OPENV_SLOPE_RANGE 6.0

typedef enum {
  UDSP_OPENV_HOLD,
  UDSP_OPENV_LOOP,
  UDSP_OPENV_TRIGGER,
} UDSP_OPENV_MODE;

typedef enum {
  UDSP_OPENV_IDLE, /* at its end level, done */
  UDSP_OPENV_ATTACK,
  UDSP_OPENV_DECAY,
  UDSP_OPENV_SUSTAIN,
  UDSP_OPENV_RETURN, /* loop: from sustain back to the initial level */
  UDSP_OPENV_RELEASE,
} UDSP_OPENV_STAGE;

typedef struct {
  float initial, peak, sustain, end;              /* levels */
  float attack_ms, decay_ms, release_ms;          /* 0: a frame */
  float loop_ms;                                  /* loop mode's return */
  float attack_slope, decay_slope, release_slope; /* -1 .. 1 */
  UDSP_OPENV_MODE mode;
} UDSP_OPENV_SHAPE;

typedef struct {
  UDSP_OPENV_SHAPE shape;
  unsigned rate;
  float value;
  UDSP_OPENV_STAGE stage;
  /* the segment running: pos of len frames from `from` to `to`, its
   * curve u = r^pos over den = r^len - 1 (den 0: a line) */
  uint32_t pos, len;
  float from, to;
  double u, r, den;
} UDSP_OPENV;

/* idle at the shape's end level; set it before gating */
void udsp_openv_init(UDSP_OPENV *e, const UDSP_OPENV_SHAPE *s, unsigned rate);

/* a new shape: the segment running keeps its course, the next ones take
 * it (sustain and the end level at once when they're where it is) */
void udsp_openv_set(UDSP_OPENV *e, const UDSP_OPENV_SHAPE *s);

/* note on (the attack) or off (the release; ignored in trigger mode) */
void udsp_openv_gate(UDSP_OPENV *e, int on);

/* beat / sync: the attack again from where it is, while it's held */
void udsp_openv_restart(UDSP_OPENV *e);

/* n values, a frame on each */
void udsp_openv_render(UDSP_OPENV *e, float *out, size_t n);

static inline int udsp_openv_idle(const UDSP_OPENV *e) {
  return e->stage == UDSP_OPENV_IDLE;
}

#endif
