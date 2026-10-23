#include "offline/sinc.h"

#include "base/simd.h"

#define ATTEN_DB 120.0
#define PI 3.14159265358979323846

static uint32_t gcd(uint32_t a, uint32_t b) {
  while (b) {
    uint32_t t = a % b;
    a = b;
    b = t;
  }
  return a;
}

typedef struct {
  uint32_t up, down, taps, rows;
  int interp;
} SHAPE;

static int shape(SHAPE *sh, uint32_t src, uint32_t dst, uint32_t zeros) {
  if (!src || !dst || src > UDSP_SINC_RATE_MAX || dst > UDSP_SINC_RATE_MAX ||
      !zeros) {
    return 0;
  }
  uint32_t g = gcd(src, dst);
  sh->up = dst / g;
  sh->down = src / g;
  if (src == dst) {
    sh->taps = 4;
    sh->rows = 1;
    sh->interp = 0;
    return 1;
  }
  /* downsampling stretches the kernel by the ratio: the same zero
   * crossings at the lower rate's spacing */
  double fc = dst < src ? (double)dst / src : 1.0;
  double half = zeros / fc;
  uint32_t taps = 2 * (uint32_t)half + (half > (uint32_t)half ? 2 : 0);
  sh->taps = (taps + 3) & ~3u;
  sh->interp = sh->up > UDSP_SINC_PHASES_MAX;
  sh->rows = sh->interp ? UDSP_SINC_PHASES_MAX + 1 : sh->up;
  return 1;
}

size_t udsp_sinc_floats(uint32_t src_rate, uint32_t dst_rate, uint32_t zeros) {
  SHAPE sh;
  return shape(&sh, src_rate, dst_rate, zeros) ? (size_t)sh.rows * sh.taps : 0;
}

/* sin(pi y) in plain double arithmetic (libm's sin may round differently
 * between libcs; this is the same bits everywhere, -ffp-contract=off) */
static double sin_pi(double y) {
  int64_t k = (int64_t)(y * 0.5);
  y -= 2.0 * (double)k; /* (-2, 2) */
  if (y > 1.0) {
    y -= 2.0;
  } else if (y < -1.0) {
    y += 2.0;
  }
  if (y > 0.5) {
    y = 1.0 - y;
  } else if (y < -0.5) {
    y = -1.0 - y;
  }
  double x = PI * y, x2 = x * x, term = x, s = x;
  for (int i = 1; i <= 12; i++) {
    term *= -x2 / (double)((2 * i) * (2 * i + 1));
    s += term;
  }
  return s;
}

/* the modified Bessel function I0 by its series, a fixed number of terms
 * (z <= the Kaiser beta here, ~12.3: done long before 64) */
static double bessel_i0(double z) {
  double q = z * z * 0.25, term = 1.0, s = 1.0;
  for (int k = 1; k <= 64; k++) {
    term *= q / (double)(k * k);
    s += term;
  }
  return s;
}

static double my_sqrt(double x) {
  return x > 0.0 ? sqrt(x) : 0.0; /* IEEE: correctly rounded everywhere */
}

/* tap k of a row for an output at fraction f (0..1) past an input */
static double tap(double f, uint32_t k, uint32_t taps, double fcut, double beta,
                  double i0_beta) {
  double half = taps * 0.5;
  double x = f + half - 1.0 - (double)k;
  double u = x / half;
  double w = u <= -1.0 || u >= 1.0
                 ? (u == -1.0 || u == 1.0 ? 1.0 / i0_beta : 0.0)
                 : bessel_i0(beta * my_sqrt(1.0 - u * u)) / i0_beta;
  double sinc = x == 0.0 ? 2.0 * fcut : sin_pi(2.0 * fcut * x) / (PI * x);
  return sinc * w;
}

/* a row, normalised in double to sum to 1 (computed twice rather than
 * held: a downsampling row can be thousands of taps) */
static void fill_row(float *row, double f, uint32_t taps, double fcut,
                     double beta, double i0_beta) {
  double sum = 0.0;
  for (uint32_t k = 0; k < taps; k++) {
    sum += tap(f, k, taps, fcut, beta, i0_beta);
  }
  for (uint32_t k = 0; k < taps; k++) {
    row[k] = (float)(tap(f, k, taps, fcut, beta, i0_beta) / sum);
  }
}

int udsp_sinc_init(UDSP_SINC *s, uint32_t src_rate, uint32_t dst_rate,
                   uint32_t zeros, float *table) {
  SHAPE sh;
  if (!shape(&sh, src_rate, dst_rate, zeros)) {
    return 0;
  }
  s->up = sh.up;
  s->down = sh.down;
  s->taps = sh.taps;
  s->rows = sh.rows;
  s->interp = sh.interp;
  s->table = table;
  if (src_rate == dst_rate) {
    table[0] = table[2] = table[3] = 0.0f;
    table[1] = 1.0f;
    return 1;
  }
  double fc = dst_rate < src_rate ? (double)dst_rate / src_rate : 1.0;
  /* Kaiser's: the transition band for ATTEN_DB over taps, centred so the
   * stopband starts at the lower rate's Nyquist */
  double delta = (ATTEN_DB - 8.0) / (2.285 * 2.0 * PI * sh.taps);
  double fcut = fc * 0.5 - delta * 0.5;
  double beta = 0.1102 * (ATTEN_DB - 8.7);
  double i0_beta = bessel_i0(beta);
  uint32_t phases = sh.interp ? UDSP_SINC_PHASES_MAX : sh.up;
  for (uint32_t r = 0; r < sh.rows; r++) {
    fill_row(table + (size_t)r * sh.taps, (double)r / phases, sh.taps, fcut,
             beta, i0_beta);
  }
  return 1;
}

uint64_t udsp_sinc_length(const UDSP_SINC *s, uint64_t in_frames) {
  return in_frames / s->down * s->up +
         ((in_frames % s->down) * s->up + s->down - 1) / s->down;
}

/* output o's input frame (the one at or before it) and its phase,
 * without o * down overflowing */
static void locate(const UDSP_SINC *s, uint64_t o, int64_t *i, uint32_t *ph) {
  uint64_t r = (o % s->up) * s->down;
  *i = (int64_t)(o / s->up * s->down + r / s->up);
  *ph = (uint32_t)(r % s->up);
}

void udsp_sinc_span(const UDSP_SINC *s, uint64_t first, size_t n, int64_t *lo,
                    int64_t *hi) {
  int64_t i0, i1;
  uint32_t ph;
  locate(s, first, &i0, &ph);
  locate(s, first + (n ? n - 1 : 0), &i1, &ph);
  *lo = i0 - (int64_t)(s->taps / 2) + 1;
  *hi = i1 + (int64_t)(s->taps / 2) + 1;
}

/* the dot product as four lanes, lane j summing taps j, j + 4, ... in
 * order, then (0 + 1) + (2 + 3): the vector and the plain build add the
 * same numbers in the same order. Row b is a second row to blend in by
 * mu (the interpolated table), NULL for none. x == NULL: read through
 * the window's bounds (the input's edges), silence outside. */
static float dot(const float *a, const float *b, float mu, const float *x,
                 const float *in, int64_t in_first, size_t in_n, int64_t j0,
                 uint32_t taps) {
  float acc[4] = {0.0f, 0.0f, 0.0f, 0.0f};
#if defined(UDSP_HAVE_V4)
  UDSP_V4 v = {0.0f, 0.0f, 0.0f, 0.0f};
  for (uint32_t k = 0; k < taps; k += 4) {
    UDSP_V4 c = udsp_v4load(a + k);
    if (b) {
      c = c + (UDSP_V4){mu, mu, mu, mu} * (udsp_v4load(b + k) - c);
    }
    UDSP_V4 xv;
    if (x) {
      xv = udsp_v4load(x + k);
    } else {
      float t[4];
      for (int j = 0; j < 4; j++) {
        int64_t at = j0 + (int64_t)k + j - in_first;
        t[j] = at >= 0 && (uint64_t)at < in_n ? in[at] : 0.0f;
      }
      xv = udsp_v4load(t);
    }
    v = v + c * xv;
  }
  udsp_v4store(acc, v);
#else
  for (uint32_t k = 0; k < taps; k += 4) {
    for (int j = 0; j < 4; j++) {
      float c = a[k + j];
      if (b) {
        c = c + mu * (b[k + j] - c);
      }
      float xv;
      if (x) {
        xv = x[k + j];
      } else {
        int64_t at = j0 + (int64_t)k + j - in_first;
        xv = at >= 0 && (uint64_t)at < in_n ? in[at] : 0.0f;
      }
      acc[j] = acc[j] + c * xv;
    }
  }
#endif
  return (acc[0] + acc[1]) + (acc[2] + acc[3]);
}

void udsp_sinc_run(const UDSP_SINC *s, const float *in, int64_t in_first,
                   size_t in_n, uint64_t first, float *out, size_t n) {
  const int64_t half = (int64_t)(s->taps / 2);
  for (size_t o = 0; o < n; o++) {
    int64_t i;
    uint32_t ph;
    locate(s, first + o, &i, &ph);
    const float *a, *b = NULL;
    float mu = 0.0f;
    if (s->interp) {
      /* the row below the exact fraction and the one above */
      double at = (double)ph / s->up * UDSP_SINC_PHASES_MAX;
      uint32_t r = (uint32_t)at;
      mu = (float)(at - r);
      a = s->table + (size_t)r * s->taps;
      b = a + s->taps;
    } else {
      a = s->table + (size_t)ph * s->taps;
    }
    int64_t j0 = i - half + 1;
    int inside =
        j0 >= in_first && j0 + (int64_t)s->taps <= in_first + (int64_t)in_n;
    out[o] = dot(a, b, mu, inside ? in + (j0 - in_first) : NULL, in, in_first,
                 in_n, j0, s->taps);
  }
}

/* --- varispeed ------------------------------------------------------------ */

#define VS_ATTEN_DB 100.0
/* the most taps an output reads: both sides at the widest, padded */
#define VS_TAPS_MAX(z) (2u * (z) * UDSP_VSINC_SPEED_MAX + 8u)

size_t udsp_vsinc_floats(uint32_t zeros) {
  return zeros ? (size_t)zeros * UDSP_VSINC_PHASES + 2 : 0;
}

void udsp_vsinc_init(UDSP_VSINC *k, uint32_t zeros, float *table) {
  k->zeros = zeros;
  k->table = table;
  const double taps = 2.0 * zeros;
  const double delta = (VS_ATTEN_DB - 8.0) / (2.285 * 2.0 * PI * taps);
  const double fcut = 0.5 - delta * 0.5;
  const double beta = 0.1102 * (VS_ATTEN_DB - 8.7);
  const double i0_beta = bessel_i0(beta);
  const uint32_t n = zeros * UDSP_VSINC_PHASES;
  for (uint32_t i = 0; i <= n + 1; i++) {
    const double x = (double)i / UDSP_VSINC_PHASES;
    const double u = x / zeros;
    const double w =
        u >= 1.0 ? 0.0 : bessel_i0(beta * my_sqrt(1.0 - u * u)) / i0_beta;
    const double sinc = i == 0 ? 2.0 * fcut : sin_pi(2.0 * fcut * x) / (PI * x);
    table[i] = (float)(sinc * w);
  }
}

/* the kernel's half width at step, in input frames */
static double vs_reach(const UDSP_VSINC *k, uint64_t step) {
  double s = (double)step / (double)UDSP_VSINC_ONE;
  s = s < 1.0 ? 1.0 : s > UDSP_VSINC_SPEED_MAX ? UDSP_VSINC_SPEED_MAX : s;
  return k->zeros * s;
}

void udsp_vsinc_span(const UDSP_VSINC *k, int64_t frame, uint32_t frac,
                     uint64_t step, size_t n, int64_t *lo, int64_t *hi) {
  const int64_t r = (int64_t)vs_reach(k, step) + 1;
  const uint64_t last = (uint64_t)frac + (n ? (uint64_t)(n - 1) * step : 0);
  *lo = frame - r;
  *hi = frame + (int64_t)(last >> 32) + r + 1;
}

void udsp_vsinc_run(const UDSP_VSINC *k, const float *in, uint32_t channels,
                    int64_t in_first, size_t in_n, int64_t frame, uint32_t frac,
                    uint64_t step, float *out, size_t out_stride, size_t n) {
  float w[VS_TAPS_MAX(UDSP_VSINC_ZEROS)];
  const double reach = vs_reach(k, step);
  /* table points an input frame of distance is */
  const double scale = UDSP_VSINC_PHASES * k->zeros / reach;
  const uint32_t limit = k->zeros * UDSP_VSINC_PHASES;
  const int64_t half = (int64_t)reach + 1;
  uint32_t taps = (uint32_t)(2 * half + 1);
  taps = (taps + 3) & ~3u;
  if (k->zeros > UDSP_VSINC_ZEROS || taps > VS_TAPS_MAX(UDSP_VSINC_ZEROS)) {
    return; /* a wider kernel than the scratch holds: not one made here */
  }
  for (size_t o = 0; o < n; o++) {
    const uint64_t at = (uint64_t)frac + (uint64_t)o * step;
    const int64_t f = frame + (int64_t)(at >> 32);
    const uint32_t fr = (uint32_t)at;
    float *y = out + o * out_stride;
    const int64_t rel = f - in_first;
    if (step == UDSP_VSINC_ONE && fr == 0) { /* the input as it is */
      for (uint32_t c = 0; c < channels; c++) {
        y[c] = rel >= 0 && (uint64_t)rel < in_n ? in[(size_t)rel * channels + c]
                                                : 0.0f;
      }
      continue;
    }
    /* the weights, for every channel: taps j0 .. j0 + taps - 1 */
    const int64_t j0 = f - half + 1;
    const double p = (double)fr / (double)UDSP_VSINC_ONE;
    float sum4[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    for (uint32_t t = 0; t < taps; t++) {
      double x = (double)(j0 + (int64_t)t - f) - p;
      x = (x < 0.0 ? -x : x) * scale;
      float v = 0.0f;
      if (x < (double)limit) {
        const uint32_t i = (uint32_t)x;
        const float mu = (float)(x - (double)i);
        v = k->table[i] + mu * (k->table[i + 1] - k->table[i]);
      }
      w[t] = v;
      sum4[t & 3] = sum4[t & 3] + v;
    }
    const float sum = (sum4[0] + sum4[1]) + (sum4[2] + sum4[3]);
    const float norm = sum != 0.0f ? 1.0f / sum : 0.0f;
    const int inside =
        j0 >= in_first && j0 + (int64_t)taps <= in_first + (int64_t)in_n;
    for (uint32_t c = 0; c < channels; c++) {
      float acc[4] = {0.0f, 0.0f, 0.0f, 0.0f};
      for (uint32_t t = 0; t < taps; t += 4) {
        for (uint32_t j = 0; j < 4; j++) {
          const int64_t r = j0 + (int64_t)(t + j) - in_first;
          const float xv = inside || (r >= 0 && (uint64_t)r < in_n)
                               ? in[(size_t)r * channels + c]
                               : 0.0f;
          acc[j] = acc[j] + w[t + j] * xv;
        }
      }
      y[c] = ((acc[0] + acc[1]) + (acc[2] + acc[3])) * norm;
    }
  }
}
