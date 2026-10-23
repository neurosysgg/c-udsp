#include "fx/zdf.h"

#include "base/simd.h"
#include "synth/osc.h"

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

/* one filter, its loop as it always was but the mode a constant (in
 * locals, in a function of its own: inside zdf_group(), as zdf_rows()'s
 * first row, clang packs ic1 and ic2 into a vector, 1.5x slower) */
static inline __attribute__((always_inline)) void
zdf_one(UDSP_ZDF *f, float *io, size_t n, UDSP_ZDF_MODE mode) {
  float ic1 = f->ic1, ic2 = f->ic2;
  const float k = f->k, a1 = f->a1, a2 = f->a2, a3 = f->a3;
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

/* rows filters side by side, frame by frame: their chains (~20 cycles a
 * frame through ic2) overlap. Each filter's operations those of its own
 * loop, in order; the mode a constant, so its switch is gone too */
static inline __attribute__((always_inline)) void
zdf_rows(UDSP_ZDF *const f[], float *const io[], size_t rows, size_t n,
         UDSP_ZDF_MODE mode) {
  float *x[UDSP_ZDF_N];
  float c1[UDSP_ZDF_N], c2[UDSP_ZDF_N], k[UDSP_ZDF_N];
  float a1[UDSP_ZDF_N], a2[UDSP_ZDF_N], a3[UDSP_ZDF_N];
  for (size_t r = 0; r < rows; r++) {
    x[r] = io[r];
    c1[r] = f[r]->ic1;
    c2[r] = f[r]->ic2;
    k[r] = f[r]->k;
    a1[r] = f[r]->a1;
    a2[r] = f[r]->a2;
    a3[r] = f[r]->a3;
  }
  for (size_t i = 0; i < n; i++) {
    UDSP_UNROLL
    for (size_t r = 0; r < rows; r++) {
      const float v0 = x[r][i];
      const float v3 = v0 - c2[r];
      const float v1 = a1[r] * c1[r] + a2[r] * v3;
      const float v2 = c2[r] + a2[r] * c1[r] + a3[r] * v3;
      c1[r] = 2.0f * v1 - c1[r];
      c2[r] = 2.0f * v2 - c2[r];
      switch (mode) {
      case UDSP_ZDF_BAND:
        x[r][i] = v1;
        break;
      case UDSP_ZDF_HIGH:
        x[r][i] = v0 - k[r] * v1 - v2;
        break;
      case UDSP_ZDF_NOTCH:
        x[r][i] = v0 - k[r] * v1;
        break;
      default:
        x[r][i] = v2;
        break;
      }
    }
  }
  for (size_t r = 0; r < rows; r++) {
    f[r]->ic1 = c1[r];
    f[r]->ic2 = c2[r];
  }
}

#define ZDF_ROWS(mode)                                                         \
  switch (rows) {                                                              \
  case 1:                                                                      \
    udsp_zdf_mono(f[0], io[0], n);                                             \
    break;                                                                     \
  case 2:                                                                      \
    zdf_rows(f, io, 2, n, mode);                                               \
    break;                                                                     \
  case 3:                                                                      \
    zdf_rows(f, io, 3, n, mode);                                               \
    break;                                                                     \
  default:                                                                     \
    zdf_rows(f, io, 4, n, mode);                                               \
    break;                                                                     \
  }

/* 1..UDSP_ZDF_N filters of one mode (gcc's SLP would pack a filter's
 * ic1 and ic2 into a vector: shuffles on the chain, 1.6x slower) */
UDSP_NO_SLP static void zdf_group(UDSP_ZDF *const f[], float *const io[],
                                  size_t rows, size_t n, UDSP_ZDF_MODE mode) {
  switch (mode) {
  case UDSP_ZDF_BAND:
    ZDF_ROWS(UDSP_ZDF_BAND);
    break;
  case UDSP_ZDF_HIGH:
    ZDF_ROWS(UDSP_ZDF_HIGH);
    break;
  case UDSP_ZDF_NOTCH:
    ZDF_ROWS(UDSP_ZDF_NOTCH);
    break;
  default:
    ZDF_ROWS(UDSP_ZDF_LOW);
    break;
  }
}

/* (no SLP: gcc packs ic1 and ic2 as in zdf_group()) */
UDSP_NO_SLP void udsp_zdf_mono(UDSP_ZDF *f, float *io, size_t n) {
  switch (f->mode) {
  case UDSP_ZDF_BAND:
    zdf_one(f, io, n, UDSP_ZDF_BAND);
    break;
  case UDSP_ZDF_HIGH:
    zdf_one(f, io, n, UDSP_ZDF_HIGH);
    break;
  case UDSP_ZDF_NOTCH:
    zdf_one(f, io, n, UDSP_ZDF_NOTCH);
    break;
  default:
    zdf_one(f, io, n, UDSP_ZDF_LOW);
    break;
  }
}

void udsp_zdf_monos(UDSP_ZDF *const f[], float *const io[], size_t count,
                    size_t n) {
  for (size_t j = 0; j < count;) {
    /* a run of one mode, up to UDSP_ZDF_N */
    size_t rows = 1;
    while (j + rows < count && rows < UDSP_ZDF_N &&
           f[j + rows]->mode == f[j]->mode) {
      rows++;
    }
    zdf_group(f + j, io + j, rows, n, f[j]->mode);
    j += rows;
  }
}
