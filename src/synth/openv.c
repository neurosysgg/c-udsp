#include "synth/openv.h"

#include "synth/pitch.h"

static uint32_t frames_of(float ms, unsigned rate) {
  const double f = (double)ms * rate / 1000.0 + 0.5;
  return f < 1.0 ? 1u : f > 4294967295.0 ? 0xffffffffu : (uint32_t)f;
}

static float clamp_slope(float s) {
  return s < -1.0f ? -1.0f : s > 1.0f ? 1.0f : s;
}

/* a segment of ms from the value to `to`, bent by slope */
static void segment(UDSP_OPENV *e, UDSP_OPENV_STAGE st, float to, float ms,
                    float slope) {
  e->stage = st;
  e->pos = 0;
  e->len = frames_of(ms, e->rate);
  e->from = e->value;
  e->to = to;
  const double k = (double)clamp_slope(slope) * UDSP_OPENV_SLOPE_RANGE;
  e->u = 1.0;
  if (k == 0.0) {
    e->r = 1.0;
    e->den = 0.0;
  } else {
    e->r = udsp_exp2(k / e->len);
    e->den = udsp_exp2(k) - 1.0;
  }
}

static void attack(UDSP_OPENV *e) {
  const UDSP_OPENV_SHAPE *s = &e->shape;
  segment(e, UDSP_OPENV_ATTACK, s->peak, s->attack_ms, s->attack_slope);
}

static void release(UDSP_OPENV *e) {
  const UDSP_OPENV_SHAPE *s = &e->shape;
  segment(e, UDSP_OPENV_RELEASE, s->end, s->release_ms, s->release_slope);
}

/* the segment after the one that just landed */
static void next(UDSP_OPENV *e) {
  const UDSP_OPENV_SHAPE *s = &e->shape;
  switch (e->stage) {
  case UDSP_OPENV_ATTACK:
    segment(e, UDSP_OPENV_DECAY, s->sustain, s->decay_ms, s->decay_slope);
    break;
  case UDSP_OPENV_DECAY:
    if (s->mode == UDSP_OPENV_LOOP) {
      segment(e, UDSP_OPENV_RETURN, s->initial, s->loop_ms, 0.0f);
    } else if (s->mode == UDSP_OPENV_TRIGGER) {
      release(e);
    } else {
      e->stage = UDSP_OPENV_SUSTAIN;
    }
    break;
  case UDSP_OPENV_RETURN:
    attack(e);
    break;
  default:
    e->stage = UDSP_OPENV_IDLE;
    break;
  }
}

void udsp_openv_init(UDSP_OPENV *e, const UDSP_OPENV_SHAPE *s, unsigned rate) {
  e->shape = *s;
  e->rate = rate;
  e->value = s->end;
  e->stage = UDSP_OPENV_IDLE;
  e->pos = e->len = 0;
  e->from = e->to = s->end;
  e->u = e->r = 1.0;
  e->den = 0.0;
}

void udsp_openv_set(UDSP_OPENV *e, const UDSP_OPENV_SHAPE *s) {
  const UDSP_OPENV_SHAPE old = e->shape;
  e->shape = *s;
  if (e->stage == UDSP_OPENV_SUSTAIN) {
    e->value = s->sustain;
  } else if (e->stage == UDSP_OPENV_IDLE && e->value == old.end) {
    e->value = s->end;
  }
}

void udsp_openv_gate(UDSP_OPENV *e, int on) {
  if (on) {
    if (e->stage == UDSP_OPENV_IDLE) {
      e->value = e->shape.initial;
    }
    attack(e);
  } else if (e->stage != UDSP_OPENV_IDLE && e->stage != UDSP_OPENV_RELEASE &&
             e->shape.mode != UDSP_OPENV_TRIGGER) {
    release(e);
  }
}

void udsp_openv_restart(UDSP_OPENV *e) {
  if (e->stage != UDSP_OPENV_IDLE && e->stage != UDSP_OPENV_RELEASE) {
    attack(e);
  }
}

/* a segment's frames short of its landing, the stage tested once */
static void run(UDSP_OPENV *e, float *out, size_t n) {
  const double from = e->from, d = e->to - e->from;
  uint32_t pos = e->pos;
  float value = e->value;
  if (e->den == 0.0) {
    const double len = e->len;
    for (size_t i = 0; i < n; i++) {
      out[i] = value;
      pos++;
      value = (float)(from + d * ((double)pos / len));
    }
  } else {
    const double r = e->r, den = e->den;
    double u = e->u;
    for (size_t i = 0; i < n; i++) {
      out[i] = value;
      pos++;
      u *= r;
      value = (float)(from + d * ((u - 1.0) / den));
    }
    e->u = u;
  }
  e->pos = pos;
  e->value = value;
}

void udsp_openv_render(UDSP_OPENV *e, float *out, size_t n) {
  size_t i = 0;
  while (i < n) {
    if (e->stage == UDSP_OPENV_IDLE || e->stage == UDSP_OPENV_SUSTAIN) {
      for (; i < n; i++) {
        out[i] = e->value;
      }
      return;
    }
    const size_t left = e->pos < e->len ? e->len - e->pos - 1 : 0;
    const size_t k = n - i < left ? n - i : left;
    run(e, out + i, k);
    i += k;
    if (i < n) { /* the landing */
      out[i++] = e->value;
      e->pos++;
      e->value = e->to;
      next(e);
    }
  }
}
