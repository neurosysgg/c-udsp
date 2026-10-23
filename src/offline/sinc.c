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

/* taps k .. k + 3 of row a (blended into row b by mu: the interpolated
 * table; b NULL for none) times the input's. x == NULL: read through the
 * window's bounds (the input's edges), silence outside. */
#if defined(UDSP_HAVE_V4)
static inline UDSP_V4 block(const float *a, const float *b, float mu,
                            const float *x, const float *in, int64_t in_first,
                            size_t in_n, int64_t j0, uint32_t k) {
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
  return c * xv;
}
#else
static inline void block(const float *a, const float *b, float mu,
                         const float *x, const float *in, int64_t in_first,
                         size_t in_n, int64_t j0, uint32_t k, float *acc) {
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

/* the dot product as sixteen lanes: the 4-tap blocks taken in turn by
 * four 4-lane sums (block n by sum n & 3, a sum's lane j summing its
 * blocks' tap j in order), then the sums (0 + 1) + (2 + 3) and their
 * lanes the same way. Four chains in flight instead of one (one was ~3
 * cycles a block of latency), and the vector and the plain build add
 * the same numbers in the same order. */
static inline __attribute__((always_inline)) float
dot_as(const float *a, const float *b, float mu, const float *x,
       const float *in, int64_t in_first, size_t in_n, int64_t j0,
       uint32_t taps, const int blend, const int direct) {
  /* constants where inlined: the tests leave the loop (gcc -O2 doesn't
   * unswitch it itself) */
  b = blend ? b : NULL;
  x = direct ? x : NULL;
  float acc[4];
  uint32_t k = 0;
#if defined(UDSP_HAVE_V4)
  UDSP_V4 v[4] = {{0}, {0}, {0}, {0}};
  for (; k + 16 <= taps; k += 16) {
    UDSP_UNROLL
    for (uint32_t q = 0; q < 4; q++) {
      v[q] = v[q] + block(a, b, mu, x, in, in_first, in_n, j0, k + 4 * q);
    }
  }
  for (uint32_t q = 0; k < taps; k += 4, q++) {
    v[q] = v[q] + block(a, b, mu, x, in, in_first, in_n, j0, k);
  }
  udsp_v4store(acc, (v[0] + v[1]) + (v[2] + v[3]));
#else
  float v[4][4] = {{0}};
  for (; k + 16 <= taps; k += 16) {
    for (uint32_t q = 0; q < 4; q++) {
      block(a, b, mu, x, in, in_first, in_n, j0, k + 4 * q, v[q]);
    }
  }
  for (uint32_t q = 0; k < taps; k += 4, q++) {
    block(a, b, mu, x, in, in_first, in_n, j0, k, v[q]);
  }
  for (int j = 0; j < 4; j++) {
    acc[j] = (v[0][j] + v[1][j]) + (v[2][j] + v[3][j]);
  }
#endif
  return (acc[0] + acc[1]) + (acc[2] + acc[3]);
}

static float dot(const float *a, const float *b, float mu, const float *x,
                 const float *in, int64_t in_first, size_t in_n, int64_t j0,
                 uint32_t taps) {
  if (b) {
    return x ? dot_as(a, b, mu, x, in, in_first, in_n, j0, taps, 1, 1)
             : dot_as(a, b, mu, x, in, in_first, in_n, j0, taps, 1, 0);
  }
  return x ? dot_as(a, b, mu, x, in, in_first, in_n, j0, taps, 0, 1)
           : dot_as(a, b, mu, x, in, in_first, in_n, j0, taps, 0, 0);
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

/* an output's weights, taps (a multiple of 4) of them from the tap at
 * d0 frames from its frame (1 - the half width), p its fraction: four
 * taps a step, each lane the scalar's operations (the distance in double,
 * the table read linearly) and sum4 the four lanes' sums */
static void vs_weights(const float *table, int64_t d0, double p, double scale,
                       uint32_t limit, uint32_t taps, float *w, float *sum4) {
#if defined(__SSE2__) && UDSP_W > 1
  __m128 sum = _mm_setzero_ps();
  const __m128i lim = _mm_set1_epi32((int)limit);
#if !defined(__AVX__)
  const __m128d sign = _mm_set1_pd(-0.0), sc = _mm_set1_pd(scale),
                pv = _mm_set1_pd(p);
#endif
  for (uint32_t t = 0; t < taps; t += 4) {
    const double d = (double)(d0 + (int64_t)t);
#if defined(__AVX__)
    const __m256d x0 = _mm256_sub_pd(
        _mm256_add_pd(_mm256_set1_pd(d), _mm256_set_pd(3.0, 2.0, 1.0, 0.0)),
        _mm256_set1_pd(p));
    const __m256d x = _mm256_mul_pd(_mm256_andnot_pd(_mm256_set1_pd(-0.0), x0),
                                    _mm256_set1_pd(scale));
    const __m128i i = _mm256_cvttpd_epi32(x);
    const __m128 mu = _mm256_cvtpd_ps(_mm256_sub_pd(x, _mm256_cvtepi32_pd(i)));
#else
    const __m128d xa = _mm_mul_pd(
        _mm_andnot_pd(
            sign,
            _mm_sub_pd(_mm_add_pd(_mm_set1_pd(d), _mm_set_pd(1.0, 0.0)), pv)),
        sc);
    const __m128d xb = _mm_mul_pd(
        _mm_andnot_pd(
            sign,
            _mm_sub_pd(_mm_add_pd(_mm_set1_pd(d), _mm_set_pd(3.0, 2.0)), pv)),
        sc);
    const __m128i ia = _mm_cvttpd_epi32(xa), ib = _mm_cvttpd_epi32(xb);
    const __m128i i = _mm_unpacklo_epi64(ia, ib);
    const __m128 mu =
        _mm_movelh_ps(_mm_cvtpd_ps(_mm_sub_pd(xa, _mm_cvtepi32_pd(ia))),
                      _mm_cvtpd_ps(_mm_sub_pd(xb, _mm_cvtepi32_pd(ib))));
#endif
    /* x < limit, limit whole: its whole part below it; the others read
     * the table's start and weigh 0 */
    const __m128i ok = _mm_cmplt_epi32(i, lim);
    const __m128i at = _mm_and_si128(i, ok);
#if defined(__AVX2__)
    const __m128 t0 = _mm_i32gather_ps(table, at, 4);
    const __m128 t1 = _mm_i32gather_ps(table + 1, at, 4);
#else
    int ix[4];
    _mm_storeu_si128((__m128i *)ix, at);
    const __m128 t0 =
        _mm_set_ps(table[ix[3]], table[ix[2]], table[ix[1]], table[ix[0]]);
    const __m128 t1 = _mm_set_ps(table[ix[3] + 1], table[ix[2] + 1],
                                 table[ix[1] + 1], table[ix[0] + 1]);
#endif
    const __m128 v =
        _mm_and_ps(_mm_add_ps(t0, _mm_mul_ps(mu, _mm_sub_ps(t1, t0))),
                   _mm_castsi128_ps(ok));
    _mm_storeu_ps(w + t, v);
    sum = _mm_add_ps(sum, v);
  }
  _mm_storeu_ps(sum4, sum);
#else
  for (int j = 0; j < 4; j++) {
    sum4[j] = 0.0f;
  }
  for (uint32_t t = 0; t < taps; t += 4) {
    UDSP_UNROLL
    for (uint32_t j = 0; j < 4; j++) {
      double x = (double)(d0 + (int64_t)(t + j)) - p;
      x = (x < 0.0 ? -x : x) * scale;
      float v = 0.0f;
      if (x < (double)limit) {
        const uint32_t i = (uint32_t)x;
        const float mu = (float)(x - (double)i);
        v = table[i] + mu * (table[i + 1] - table[i]);
      }
      w[t + j] = v;
      sum4[j] = sum4[j] + v;
    }
  }
#endif
}

/* both channels of a stereo input's taps (all inside it) weighed: the
 * scalar loop's four accumulators a channel, in lanes */
static void vs_dot2(const float *w, const float *x, uint32_t taps, float norm,
                    float *y) {
#if UDSP_HAVE_V4
  UDSP_V4 al = {0.0f, 0.0f, 0.0f, 0.0f}, ar = al;
  for (uint32_t t = 0; t < taps; t += 4) {
    const UDSP_V4 a = udsp_v4load(x + 2 * t), b = udsp_v4load(x + 2 * t + 4);
    const UDSP_V4 wv = udsp_v4load(w + t);
    al = al + wv * UDSP_SHUF4(a, b, 0, 2, 4, 6);
    ar = ar + wv * UDSP_SHUF4(a, b, 1, 3, 5, 7);
  }
  y[0] = ((al[0] + al[1]) + (al[2] + al[3])) * norm;
  y[1] = ((ar[0] + ar[1]) + (ar[2] + ar[3])) * norm;
#else
  for (uint32_t c = 0; c < 2; c++) {
    float acc[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    for (uint32_t t = 0; t < taps; t += 4) {
      for (uint32_t j = 0; j < 4; j++) {
        acc[j] = acc[j] + w[t + j] * x[(size_t)(t + j) * 2 + c];
      }
    }
    y[c] = ((acc[0] + acc[1]) + (acc[2] + acc[3])) * norm;
  }
#endif
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
    float sum4[4];
    vs_weights(k->table, 1 - half, p, scale, limit, taps, w, sum4);
    const float sum = (sum4[0] + sum4[1]) + (sum4[2] + sum4[3]);
    const float norm = sum != 0.0f ? 1.0f / sum : 0.0f;
    const int inside =
        j0 >= in_first && j0 + (int64_t)taps <= in_first + (int64_t)in_n;
    for (uint32_t c = 0; c < channels; c++) {
      float acc[4] = {0.0f, 0.0f, 0.0f, 0.0f};
      if (inside && channels == 2 && c == 0) {
        /* both channels at once: y[1]'s too */
        vs_dot2(w, in + (size_t)(j0 - in_first) * 2, taps, norm, y);
        break;
      }
      if (inside) {
        const float *x = in + (size_t)(j0 - in_first) * channels + c;
        for (uint32_t t = 0; t < taps; t += 4) {
          UDSP_UNROLL
          for (uint32_t j = 0; j < 4; j++) {
            acc[j] = acc[j] + w[t + j] * x[(size_t)(t + j) * channels];
          }
        }
      } else {
        for (uint32_t t = 0; t < taps; t += 4) {
          for (uint32_t j = 0; j < 4; j++) {
            const int64_t r = j0 + (int64_t)(t + j) - in_first;
            const float xv = r >= 0 && (uint64_t)r < in_n
                                 ? in[(size_t)r * channels + c]
                                 : 0.0f;
            acc[j] = acc[j] + w[t + j] * xv;
          }
        }
      }
      y[c] = ((acc[0] + acc[1]) + (acc[2] + acc[3])) * norm;
    }
  }
}
