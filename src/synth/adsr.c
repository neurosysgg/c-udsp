#include "synth/adsr.h"

#include "synth/pitch.h"

/* log2(UDSP_ADSR_FLOOR) */
#define FLOOR_LOG2 -13.287712379549449

static float frames_of(float ms, unsigned rate) {
  float f = ms * (float)rate / 1000.0f;
  return f < 1.0f ? 1.0f : f;
}

/* the factor per frame that takes a distance of 1 to the floor in the
 * segment's frames */
static float factor(float ms, unsigned rate) {
  return (float)udsp_exp2(FLOOR_LOG2 / frames_of(ms, rate));
}

void udsp_adsr_set(UDSP_ADSR *e, float attack_ms, float decay_ms, float sustain,
                   float release_ms, unsigned rate) {
  e->attack = 1.0f / frames_of(attack_ms, rate);
  e->decay = factor(decay_ms, rate);
  e->release = factor(release_ms, rate);
  e->sustain = sustain < 0.0f ? 0.0f : (sustain > 1.0f ? 1.0f : sustain);
}

void udsp_adsr_gate(UDSP_ADSR *e, int on) {
  if (on) {
    e->stage = UDSP_ADSR_ATTACK;
  } else if (e->stage != UDSP_ADSR_IDLE) {
    e->stage = UDSP_ADSR_RELEASE;
  }
}

/* a loop per stage, n frames on: per frame the operations of one switch
 * over the stage, in its order (the same bits), without its dispatch --
 * and a held sustain or idle is a fill. out NULL: only the state moves */
static inline __attribute__((always_inline)) void walk(UDSP_ADSR *e, float *out,
                                                       size_t n) {
  float v = e->value;
  UDSP_ADSR_STAGE st = e->stage;
  size_t i = 0;
  while (i < n) {
    switch (st) {
    case UDSP_ADSR_ATTACK: {
      /* within half a step: the sum's rounding doesn't cost a frame */
      const float a = e->attack, top = 1.0f - 0.5f * e->attack;
      for (; i < n; i++) {
        if (out) {
          out[i] = v;
        }
        v += a;
        if (v >= top) {
          v = 1.0f;
          st = UDSP_ADSR_DECAY;
          i++;
          break;
        }
      }
      break;
    }
    case UDSP_ADSR_DECAY: {
      const float s = e->sustain, k = e->decay;
      for (; i < n; i++) {
        if (out) {
          out[i] = v;
        }
        const float d = (v - s) * k;
        if (d < UDSP_ADSR_FLOOR && d > -UDSP_ADSR_FLOOR) {
          v = s;
          st = UDSP_ADSR_SUSTAIN;
          i++;
          break;
        }
        v = s + d;
      }
      break;
    }
    case UDSP_ADSR_RELEASE: {
      const float k = e->release;
      for (; i < n; i++) {
        if (out) {
          out[i] = v;
        }
        v *= k;
        if (v < UDSP_ADSR_FLOOR) {
          v = 0.0f;
          st = UDSP_ADSR_IDLE;
          i++;
          break;
        }
      }
      break;
    }
    case UDSP_ADSR_SUSTAIN:
    case UDSP_ADSR_IDLE:
      /* the first frame what it was (a sustain moved while held), then
       * the level */
      if (out) {
        out[i++] = v;
      }
      v = st == UDSP_ADSR_SUSTAIN ? e->sustain : 0.0f;
      if (out) {
        for (; i < n; i++) {
          out[i] = v;
        }
      }
      i = n;
      break;
    }
  }
  e->value = v;
  e->stage = st;
}

void udsp_adsr_render(UDSP_ADSR *e, float *out, size_t n) { walk(e, out, n); }

void udsp_adsr_skip(UDSP_ADSR *e, size_t n) { walk(e, NULL, n); }
