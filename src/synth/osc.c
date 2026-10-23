#include "synth/osc.h"

#include <string.h>

#include "base/simd.h"

#define TURNS (1.0f / 4294967296.0f)

static float poly_blep(float t, float dt) {
  if (dt <= 0.0f) {
    return 0.0f;
  }
  if (t < dt) {
    t /= dt;
    return t + t - t * t - 1.0f;
  }
  if (t > 1.0f - dt) {
    t = (t - 1.0f) / dt;
    return t * t + t + t + 1.0f;
  }
  return 0.0f;
}

static float wave_at(UDSP_OSC_WAVE w, float t, float dt) {
  switch (w) {
  case UDSP_OSC_SQUARE: {
    float v = t < 0.5f ? 1.0f : -1.0f;
    v += poly_blep(t, dt);
    float t2 = t + 0.5f;
    if (t2 >= 1.0f) {
      t2 -= 1.0f;
    }
    return v - poly_blep(t2, dt);
  }
  case UDSP_OSC_SAW:
    return 2.0f * t - 1.0f - poly_blep(t, dt);
  case UDSP_OSC_TRIANGLE:
    return t < 0.5f ? 4.0f * t - 1.0f : 3.0f - 4.0f * t;
  default:
    return udsp_osc_sine(t);
  }
}

#if UDSP_W > 1
static UDSP_VF vsine(UDSP_VF t) {
  const UDSP_VF half = udsp_vsplat(0.5f), quarter = udsp_vsplat(0.25f);
  UDSP_VF x = udsp_vsel(t >= half, t - 1.0f, t);
  x = udsp_vsel(x > quarter, half - x, udsp_vsel(x < -quarter, -half - x, x));
  const UDSP_VF u = 6.283185307179586f * x;
  const UDSP_VF u2 = u * u;
  UDSP_VF p = udsp_vsplat(-1.0f / 362880.0f);
  p = p * u2 + (1.0f / 5040.0f);
  p = p * u2 - (1.0f / 120.0f);
  p = p * u2 + (1.0f / 6.0f);
  p = p * u2 - 1.0f;
  return -u * p;
}

/* poly_blep() per lane: both branches computed, the lanes picked */
static UDSP_VF vblep(UDSP_VF t, float dt) {
  const UDSP_VF vdt = udsp_vsplat(dt), zero = udsp_vsplat(0.0f);
  UDSP_VF a = t / vdt;
  a = a + a - a * a - 1.0f;
  UDSP_VF b = (t - 1.0f) / vdt;
  b = b * b + b + b + 1.0f;
  return udsp_vsel(t < vdt, a, udsp_vsel(t > 1.0f - vdt, b, zero));
}

static UDSP_VF vwave(UDSP_OSC_WAVE w, UDSP_VF t, float dt) {
  const UDSP_VF half = udsp_vsplat(0.5f), one = udsp_vsplat(1.0f);
  switch (w) {
  case UDSP_OSC_SQUARE: {
    UDSP_VF v = udsp_vsel(t < half, one, -one);
    UDSP_VF t2 = t + 0.5f;
    t2 = udsp_vsel(t2 >= one, t2 - 1.0f, t2);
    if (dt > 0.0f) {
      v += vblep(t, dt);
      v -= vblep(t2, dt);
    }
    return v;
  }
  case UDSP_OSC_SAW: {
    UDSP_VF v = 2.0f * t - 1.0f;
    return dt > 0.0f ? v - vblep(t, dt) : v;
  }
  case UDSP_OSC_TRIANGLE:
    return udsp_vsel(t < half, 4.0f * t - 1.0f, 3.0f - 4.0f * t);
  default:
    return vsine(t);
  }
}
#endif

void udsp_osc_init(UDSP_OSC *o, UDSP_OSC_WAVE wave, uint32_t seed) {
  o->wave = wave;
  o->phase = 0;
  o->inc = 0;
  o->rng = seed ? seed : 0x9e3779b9u;
}

uint32_t udsp_osc_inc_millihz(uint32_t millihz, uint32_t rate) {
  if (rate == 0) {
    return 0;
  }
  return (uint32_t)(((uint64_t)millihz << 32) / ((uint64_t)rate * 1000ull));
}

uint32_t udsp_osc_inc_hz(double hz, uint32_t rate) {
  if (rate == 0 || hz <= 0.0) {
    return 0;
  }
  double inc = hz / rate * 4294967296.0 + 0.5;
  return inc >= 4294967295.0 ? 0xffffffffu : (uint32_t)inc;
}

void udsp_osc_render(UDSP_OSC *o, float *out, size_t n) {
  size_t i = 0;
  if (o->wave == UDSP_OSC_NOISE) {
    uint32_t x = o->rng;
    for (; i < n; i++) {
      x ^= x << 13;
      x ^= x >> 17;
      x ^= x << 5;
      out[i] = (float)(int32_t)x * (1.0f / 2147483648.0f);
    }
    o->rng = x;
    o->phase += (uint32_t)n * o->inc;
    return;
  }
  const float dt = (float)o->inc * TURNS;
  uint32_t phase = o->phase;
#if UDSP_W > 1
  UDSP_VU ph;
  for (int k = 0; k < UDSP_W; k++) {
    ph[k] = phase + (uint32_t)k * o->inc;
  }
  const uint32_t step = (uint32_t)UDSP_W * o->inc;
  for (; i + UDSP_W <= n; i += UDSP_W) {
    UDSP_VF t = __builtin_convertvector(ph, UDSP_VF) * TURNS;
    udsp_vstore(out + i, vwave(o->wave, t, dt));
    ph += step;
  }
  phase = ph[0];
#endif
  for (; i < n; i++) {
    out[i] = wave_at(o->wave, (float)phase * TURNS, dt);
    phase += o->inc;
  }
  o->phase = phase;
}

/* t wrapped into [0, 1): truncation, one down below 0 -- exact, the
 * same in a lane as in the scalar path */
static float wrap_turns(float t) {
  float w = t - (float)(int32_t)t;
  return w < 0.0f ? w + 1.0f : w;
}

void udsp_osc_render_pm(UDSP_OSC *o, const float *pm, float *out, size_t n) {
  if (o->wave == UDSP_OSC_NOISE) {
    udsp_osc_render(o, out, n);
    return;
  }
  size_t i = 0;
  const float dt = (float)o->inc * TURNS;
  uint32_t phase = o->phase;
#if UDSP_W > 1
  UDSP_VU ph;
  for (int k = 0; k < UDSP_W; k++) {
    ph[k] = phase + (uint32_t)k * o->inc;
  }
  const uint32_t step = (uint32_t)UDSP_W * o->inc;
  const UDSP_VF zero = udsp_vsplat(0.0f), one = udsp_vsplat(1.0f);
  for (; i + UDSP_W <= n; i += UDSP_W) {
    UDSP_VF t =
        __builtin_convertvector(ph, UDSP_VF) * TURNS + udsp_vload(pm + i);
    t = t -
        __builtin_convertvector(__builtin_convertvector(t, UDSP_VI), UDSP_VF);
    t = udsp_vsel(t < zero, t + one, t);
    udsp_vstore(out + i, vwave(o->wave, t, dt));
    ph += step;
  }
  phase = ph[0];
#endif
  for (; i < n; i++) {
    out[i] = wave_at(o->wave, wrap_turns((float)phase * TURNS + pm[i]), dt);
    phase += o->inc;
  }
  o->phase = phase;
}

static const char WAVE_NAMES[UDSP_OSC_WAVES][9] = {
    "sine", "square", "saw", "triangle", "noise",
};

UDSP_OSC_WAVE udsp_osc_wave_by_name(const char *name) {
  for (int i = 0; i < UDSP_OSC_WAVES; i++) {
    if (strcmp(name, WAVE_NAMES[i]) == 0) {
      return (UDSP_OSC_WAVE)i;
    }
  }
  return UDSP_OSC_WAVES;
}

const char *udsp_osc_wave_name(UDSP_OSC_WAVE w) {
  return w < UDSP_OSC_WAVES ? WAVE_NAMES[w] : "?";
}
