#include "synth/smooth.h"

void udsp_smooth_init(UDSP_SMOOTH *s, float value, uint32_t frames) {
  s->value = s->target = value;
  s->step = 0.0f;
  s->left = 0;
  s->frames = frames;
}

void udsp_smooth_to(UDSP_SMOOTH *s, float target) {
  s->target = target;
  if (s->frames == 0 || target == s->value) {
    s->value = target;
    s->left = 0;
    return;
  }
  s->step = (target - s->value) / (float)s->frames;
  s->left = s->frames;
}

void udsp_smooth_render(UDSP_SMOOTH *s, float *out, size_t n) {
  size_t i = 0;
  for (; i < n && s->left; i++) {
    out[i] = udsp_smooth_next(s);
  }
  for (; i < n; i++) {
    out[i] = s->value;
  }
}

void udsp_smooth_skip(UDSP_SMOOTH *s, size_t n) {
  if (n >= s->left) {
    s->value = s->target;
    s->left = 0;
    return;
  }
  s->left -= (uint32_t)n;
  s->value = s->target - s->step * (float)s->left;
}
