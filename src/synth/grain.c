#include "grain.h"

#include "base/trig.h"
#include "pitch.h"

static float hann(double x) {
  const double s = udsp_sin_pi(x);
  return (float)(s * s);
}

/* the perc's fall: from 1 at x = 0 to 2^-10 at x = 1, then to 0 over
 * the last sixteenth so it closes */
static double fall(double x) {
  double y = udsp_exp2(-10.0 * x);
  return x > 0.9375 ? y * (1.0 - x) * 16.0 : y;
}

void udsp_grain_window(float *table, UDSP_GRAIN_SHAPE shape) {
  for (int i = 0; i <= UDSP_GRAIN_TABLE; i++) {
    const double x = (double)i / UDSP_GRAIN_TABLE;
    float y = 0.0f;
    switch (shape) {
    case UDSP_GRAIN_TRIANGLE:
      y = (float)(x < 0.5 ? 2.0 * x : 2.0 - 2.0 * x);
      break;
    case UDSP_GRAIN_TUKEY:
      y = x < 0.25 ? hann(2.0 * x) : (x > 0.75 ? hann(2.0 * (1.0 - x)) : 1.0f);
      break;
    case UDSP_GRAIN_PERC:
    case UDSP_GRAIN_SWELL: {
      const double u = shape == UDSP_GRAIN_PERC ? x : 1.0 - x;
      /* the rise over the first 1/32 */
      y = (float)(u < 1.0 / 32.0 ? u * 32.0 * fall(u) : fall(u));
      break;
    }
    default:
      y = hann(x);
      break;
    }
    table[i] = y;
  }
  table[0] = table[UDSP_GRAIN_TABLE] = 0.0f;
}

void udsp_grain_start(UDSP_GRAIN *g, int64_t pos, int64_t step, uint32_t len,
                      float gl, float gr) {
  *g = (UDSP_GRAIN){pos, step, len, len, 0, gl, gr};
}

static inline float at(const float *x, uint64_t n, uint32_t ch, int64_t i) {
  return i < 0 || (uint64_t)i >= n ? 0.0f : x[(uint64_t)i * ch];
}

static inline float cubic(float a, float b, float c, float d, float t) {
  const float c1 = 0.5f * (c - a);
  const float c2 = a - 2.5f * b + 2.0f * c - 0.5f * d;
  const float c3 = 0.5f * (d - a) + 1.5f * (b - c);
  return ((c3 * t + c2) * t + c1) * t + b;
}

/* channel x (a stride of ch frames' floats) at Q32 position p */
static inline float read_at(const float *x, uint64_t n, uint32_t ch,
                            int64_t p) {
  const int64_t i = p >> 32;
  const float t =
      (float)(uint32_t)((uint64_t)p >> 8 & 0xffffffu) * (1.0f / 16777216.0f);
  if (i >= 1 && (uint64_t)i + 2 < n) {
    const float *f = x + (uint64_t)(i - 1) * ch;
    return cubic(f[0], f[ch], f[2 * ch], f[3 * ch], t);
  }
  return cubic(at(x, n, ch, i - 1), at(x, n, ch, i), at(x, n, ch, i + 1),
               at(x, n, ch, i + 2), t);
}

uint32_t udsp_grain_render(UDSP_GRAIN *g, const float *samples, uint32_t ch,
                           uint64_t n, const float *table, float *out,
                           uint32_t frames) {
  const uint32_t k = frames < g->left ? frames : g->left;
  /* the table's length over len frames, Q32: rounded down, so it never
   * reaches past the table's end */
  const uint64_t inc = ((uint64_t)UDSP_GRAIN_TABLE << 32) / g->len;
  uint64_t w = g->age * inc;
  int64_t p = g->pos;
  for (uint32_t i = 0; i < k; i++) {
    const uint32_t wi = (uint32_t)(w >> 32);
    const float wt =
        (float)(uint32_t)(w >> 8 & 0xffffffu) * (1.0f / 16777216.0f);
    const float a = table[wi], b = table[wi + 1];
    const float e = a + (b - a) * wt;
    const float l = read_at(samples, n, ch, p);
    const float r = ch > 1 ? read_at(samples + 1, n, ch, p) : l;
    out[2 * i] += l * e * g->gl;
    out[2 * i + 1] += r * e * g->gr;
    p += g->step;
    w += inc;
  }
  g->pos = p;
  g->age += k;
  g->left -= k;
  return k;
}
