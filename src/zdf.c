#include "zdf.h"

#include "osc.h"

void udsp_zdf_reset(UDSP_ZDF *f) { f->ic1 = f->ic2 = 0.0f; }

void udsp_zdf_set(UDSP_ZDF *f, float cutoff_hz, float res, UDSP_ZDF_MODE mode,
                  unsigned rate) {
  float t = cutoff_hz / (float)rate;
  const float lo = 10.0f / (float)rate;
  t = t < lo ? lo : (t > 0.49f ? 0.49f : t);
  /* tan(pi t) = sin / cos, a half turn's angle as turns: t / 2 */
  const float g = udsp_osc_sine(0.5f * t) / udsp_osc_sine(0.5f * t + 0.25f);
  res = res < 0.0f ? 0.0f : (res > 1.0f ? 1.0f : res);
  f->k = 2.0f - 1.98f * res;
  f->a1 = 1.0f / (1.0f + g * (g + f->k));
  f->a2 = g * f->a1;
  f->a3 = g * f->a2;
  f->mode = mode;
}

void udsp_zdf_mono(UDSP_ZDF *f, float *io, size_t n) {
  float ic1 = f->ic1, ic2 = f->ic2;
  const float k = f->k, a1 = f->a1, a2 = f->a2, a3 = f->a3;
  const UDSP_ZDF_MODE mode = f->mode;
  for (size_t i = 0; i < n; i++) {
    const float v0 = io[i];
    const float v3 = v0 - ic2;
    const float v1 = a1 * ic1 + a2 * v3;
    const float v2 = ic2 + a2 * ic1 + a3 * v3;
    ic1 = 2.0f * v1 - ic1;
    ic2 = 2.0f * v2 - ic2;
    switch (mode) {
    case UDSP_ZDF_BAND:
      io[i] = v1;
      break;
    case UDSP_ZDF_HIGH:
      io[i] = v0 - k * v1 - v2;
      break;
    case UDSP_ZDF_NOTCH:
      io[i] = v0 - k * v1;
      break;
    default:
      io[i] = v2;
      break;
    }
  }
  f->ic1 = ic1;
  f->ic2 = ic2;
}
