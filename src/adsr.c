#include "adsr.h"

#include "pitch.h"

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

void udsp_adsr_render(UDSP_ADSR *e, float *out, size_t n) {
  float v = e->value;
  UDSP_ADSR_STAGE st = e->stage;
  for (size_t i = 0; i < n; i++) {
    out[i] = v;
    switch (st) {
    case UDSP_ADSR_ATTACK:
      v += e->attack;
      /* within half a step: the sum's rounding doesn't cost a frame */
      if (v >= 1.0f - 0.5f * e->attack) {
        v = 1.0f;
        st = UDSP_ADSR_DECAY;
      }
      break;
    case UDSP_ADSR_DECAY: {
      float d = (v - e->sustain) * e->decay;
      if (d < UDSP_ADSR_FLOOR && d > -UDSP_ADSR_FLOOR) {
        v = e->sustain;
        st = UDSP_ADSR_SUSTAIN;
      } else {
        v = e->sustain + d;
      }
      break;
    }
    case UDSP_ADSR_SUSTAIN:
      v = e->sustain; /* a sustain moved while held */
      break;
    case UDSP_ADSR_RELEASE:
      v *= e->release;
      if (v < UDSP_ADSR_FLOOR) {
        v = 0.0f;
        st = UDSP_ADSR_IDLE;
      }
      break;
    case UDSP_ADSR_IDLE:
      v = 0.0f;
      break;
    }
  }
  e->value = v;
  e->stage = st;
}
