#include "wtab.h"

#include "osc.h"

#define TURNS (1.0f / 4294967296.0f)

void udsp_wtab_build(UDSP_WTAB *w, const float *amp, size_t n) {
  /* the sine in the last level meanwhile (the fundamental's own, made
   * from it at the end) */
  const int top = UDSP_WTAB_LEVELS - 1;
  const float *sine = w->t[top];
  for (uint32_t i = 0; i < UDSP_WTAB_SIZE; i++) {
    w->t[top][i] = udsp_osc_sine((float)i / (float)UDSP_WTAB_SIZE);
  }
  n = n > UDSP_WTAB_HARMONICS ? UDSP_WTAB_HARMONICS : n;
  /* from the fundamental's level down, each the one above plus the
   * harmonics it adds: harmonic h at point i is sine[h * i mod N] */
  size_t have = 0;
  for (int l = top - 1; l >= 0; l--) {
    float *t = w->t[l];
    const size_t want = (size_t)UDSP_WTAB_HARMONICS >> l;
    for (uint32_t i = 0; i < UDSP_WTAB_SIZE; i++) {
      t[i] = l == top - 1 ? 0.0f : w->t[l + 1][i];
    }
    for (size_t h = have + 1; h <= want && h <= n; h++) {
      if (amp[h - 1] == 0.0f) {
        continue;
      }
      for (uint32_t i = 0; i < UDSP_WTAB_SIZE; i++) {
        t[i] += amp[h - 1] * sine[(uint32_t)h * i & (UDSP_WTAB_SIZE - 1)];
      }
    }
    have = want;
  }
  const float a1 = n ? amp[0] : 0.0f;
  for (uint32_t i = 0; i < UDSP_WTAB_SIZE; i++) {
    w->t[top][i] = a1 == 0.0f ? 0.0f : a1 * w->t[top][i];
  }
  float peak = 0.0f;
  for (uint32_t i = 0; i < UDSP_WTAB_SIZE; i++) {
    const float a = w->t[0][i] < 0.0f ? -w->t[0][i] : w->t[0][i];
    peak = a > peak ? a : peak;
  }
  const float k = peak > 0.0f ? 1.0f / peak : 0.0f;
  for (int l = 0; l < UDSP_WTAB_LEVELS; l++) {
    for (uint32_t i = 0; i < UDSP_WTAB_SIZE; i++) {
      w->t[l][i] *= k;
    }
    w->t[l][UDSP_WTAB_SIZE] = w->t[l][0];
  }
}

void udsp_wtab_recipe(UDSP_WTAB_SHAPE s, size_t n, float *amp) {
  for (size_t h = 1; h <= UDSP_WTAB_HARMONICS; h++) {
    float a = 0.0f;
    if (h <= n) {
      switch (s) {
      case UDSP_WTAB_SAW:
        a = 1.0f / (float)h;
        break;
      case UDSP_WTAB_SQUARE:
        a = h & 1 ? 1.0f / (float)h : 0.0f;
        break;
      case UDSP_WTAB_TRIANGLE: /* sine phase: the signs alternate */
        a = h & 1 ? (h & 2 ? -1.0f : 1.0f) / (float)(h * h) : 0.0f;
        break;
      }
    }
    amp[h - 1] = a;
  }
}

int udsp_wtab_level(uint32_t inc) {
  /* level l keeps 64 >> l harmonics: all of them below Nyquist (half a
   * turn, 2^31) while (64 >> l) * inc < 2^31 */
  for (int l = 0; l < UDSP_WTAB_LEVELS - 1; l++) {
    if (((uint64_t)UDSP_WTAB_HARMONICS >> l) * inc < (1ull << 31)) {
      return l;
    }
  }
  return UDSP_WTAB_LEVELS - 1;
}

/* osc.c's wrap, and a turn just under 0 that rounds up to 1 back to 0
 * (a table read there would pass its end) */
static float wrap(float t) {
  float w = t - (float)(int32_t)t;
  w = w < 0.0f ? w + 1.0f : w;
  return w < 1.0f ? w : 0.0f;
}

void udsp_fmop_render(UDSP_FMOP *o, const float *pm, const float *amp,
                      float *out, size_t n) {
  if (!o->wave && o->feedback == 0.0f) {
    UDSP_OSC s = {UDSP_OSC_SINE, o->phase, o->inc, 1};
    if (pm) {
      udsp_osc_render_pm(&s, pm, out, n);
    } else {
      udsp_osc_render(&s, out, n);
    }
    for (size_t i = 0; i < n; i++) {
      out[i] *= amp[i];
    }
    o->phase = s.phase;
    if (n >= 2) {
      o->y1 = out[n - 1];
      o->y2 = out[n - 2];
    } else if (n == 1) {
      o->y2 = o->y1;
      o->y1 = out[0];
    }
    return;
  }
  const int level = o->wave ? udsp_wtab_level(o->inc) : 0;
  uint32_t phase = o->phase;
  float y1 = o->y1, y2 = o->y2;
  for (size_t i = 0; i < n; i++) {
    float t = (float)phase * TURNS + o->feedback * (0.5f * (y1 + y2));
    if (pm) {
      t += pm[i];
    }
    t = wrap(t);
    const float v =
        o->wave ? udsp_wtab_at(o->wave, level, t) : udsp_osc_sine(t);
    y2 = y1;
    y1 = amp[i] * v;
    out[i] = y1;
    phase += o->inc;
  }
  o->phase = phase;
  o->y1 = y1;
  o->y2 = y2;
}
