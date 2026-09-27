#include "fft.h"

#include "simd.h"

#include <math.h>

#ifndef M_PI
/* M_PI is a POSIX/BSD math.h extension, not standard C11 -- define it
 * ourselves rather than relying on a feature-test macro being set. */
#define M_PI 3.14159265358979323846
#endif

static void swapf(float *a, float *b) {
  float t = *a;
  *a = *b;
  *b = t;
}

void udsp_fft_forward(float *re, float *im, size_t n) {
  /* bit-reversal permutation -- standard iterative in-place reordering so
   * the butterfly passes below can work purely on adjacent/strided pairs */
  for (size_t i = 1, j = 0; i < n; i++) {
    size_t bit = n >> 1;
    for (; j & bit; bit >>= 1) {
      j ^= bit;
    }
    j ^= bit;
    if (i < j) {
      swapf(&re[i], &re[j]);
      swapf(&im[i], &im[j]);
    }
  }

  for (size_t len = 2; len <= n; len <<= 1) {
    double ang = -2.0 * M_PI / (double)len;
    float wr = (float)cos(ang), wi = (float)sin(ang);
    for (size_t i = 0; i < n; i += len) {
      float cur_r = 1.0f, cur_i = 0.0f;
      for (size_t k = 0; k < len / 2; k++) {
        size_t a = i + k, b = i + k + len / 2;
        float ur = re[a], ui = im[a];
        float vr = re[b] * cur_r - im[b] * cur_i;
        float vi = re[b] * cur_i + im[b] * cur_r;

        re[a] = ur + vr;
        im[a] = ui + vi;
        re[b] = ur - vr;
        im[b] = ui - vi;

        float next_r = cur_r * wr - cur_i * wi;
        float next_i = cur_r * wi + cur_i * wr;
        cur_r = next_r;
        cur_i = next_i;
      }
    }
  }
}

void udsp_window_hann(const float *in, float *out, size_t n) {
  if (n < 2) {
    if (n == 1) {
      out[0] = in[0];
    }
    return;
  }
  for (size_t i = 0; i < n; i++) {
    double w = 0.5 - 0.5 * cos(2.0 * M_PI * (double)i / (double)(n - 1));
    out[i] = in[i] * (float)w;
  }
}

void udsp_fft_magnitude(const float *re, const float *im, size_t n, float *out) {
  size_t i = 0;
  for (; i + UDSP_W <= n; i += UDSP_W) {
    UDSP_VF r = udsp_vload(re + i), m = udsp_vload(im + i);
    udsp_vstore(out + i, udsp_vsqrt(r * r + m * m));
  }
  for (; i < n; i++) {
    out[i] = sqrtf(re[i] * re[i] + im[i] * im[i]);
  }
}

/* --- the table-driven transforms --- */

void udsp_fft_table(float *table, size_t n) {
  table[0] = table[n] = 0.0f;
  for (size_t h = 1; h < n; h <<= 1) {
    for (size_t k = 0; k < h; k++) {
      double ang = -M_PI * (double)k / (double)h;
      table[h + k] = (float)cos(ang);
      table[n + h + k] = (float)sin(ang);
    }
  }
}

/* the bits of i (below 2^16) reversed, as a bits-wide number */
static const unsigned char REV8[256] = {
#define R2(x) x, x + 128, x + 64, x + 192
#define R4(x) R2(x), R2(x + 32), R2(x + 16), R2(x + 48)
#define R6(x) R4(x), R4(x + 8), R4(x + 4), R4(x + 12)
    R6(0), R6(2), R6(1), R6(3),
#undef R2
#undef R4
#undef R6
};

static inline size_t bitrev(size_t i, int bits) {
  size_t r = ((size_t)REV8[i & 0xff] << 8) | REV8[(i >> 8) & 0xff];
  return r >> (16 - bits);
}

static int log2n(size_t n) {
  int b = 0;
  while ((1u << b) < n) {
    b++;
  }
  return b;
}

/* one radix-2 stage h: pairs (a, a + h), twiddle e^(-i pi k / h) */
static void stage2(float *restrict re, float *restrict im, size_t m, size_t h,
                   const float *restrict tw_re, const float *restrict tw_im) {
  for (size_t i = 0; i < m; i += 2 * h) {
    size_t k = 0;
    for (; k + UDSP_W <= h; k += UDSP_W) {
      size_t a = i + k, b = a + h;
      UDSP_VF wr = udsp_vload(tw_re + h + k), wi = udsp_vload(tw_im + h + k);
      UDSP_VF br = udsp_vload(re + b), bi = udsp_vload(im + b);
      UDSP_VF ar = udsp_vload(re + a), ai = udsp_vload(im + a);
      UDSP_VF vr = br * wr - bi * wi, vi = br * wi + bi * wr;
      udsp_vstore(re + a, ar + vr);
      udsp_vstore(im + a, ai + vi);
      udsp_vstore(re + b, ar - vr);
      udsp_vstore(im + b, ai - vi);
    }
    for (; k < h; k++) {
      size_t a = i + k, b = a + h;
      float wr = tw_re[h + k], wi = tw_im[h + k];
      float vr = re[b] * wr - im[b] * wi, vi = re[b] * wi + im[b] * wr;
      re[b] = re[a] - vr, im[b] = im[a] - vi;
      re[a] += vr, im[a] += vi;
    }
  }
}

/* stages h and 2h in one trip through the data: points a, b, c, d = a +
 * 0, h, 2h, 3h. Stage h pairs (a,b) and (c,d) with w1 = e^(-i pi k/h);
 * stage 2h pairs (a,c) with w2 = e^(-i pi k/2h) and (b,d) with
 * e^(-i pi (k+h)/2h) = -i w2, so its twiddle costs nothing. Half the
 * loads and stores of two radix-2 stages, the same multiplies. */
#define RADIX4_BODY(V, LOAD, STORE)                                          \
  /* in this order, loading each value right before its use: at most ~12   \
   * vectors live, inside the 16 registers x86-64 has (loaded all up front \
   * they spilled, most under AVX) */                                      \
  V wr = LOAD(tw_re + h + k), wi = LOAD(tw_im + h + k);                     \
  V xr = LOAD(re + a + h), xi = LOAD(im + a + h);                           \
  V tr = xr * wr - xi * wi, ti = xr * wi + xi * wr;                         \
  xr = LOAD(re + a), xi = LOAD(im + a);                                     \
  V a1r = xr + tr, a1i = xi + ti, b1r = xr - tr, b1i = xi - ti;             \
  xr = LOAD(re + a + 3 * h), xi = LOAD(im + a + 3 * h);                     \
  tr = xr * wr - xi * wi, ti = xr * wi + xi * wr;                           \
  xr = LOAD(re + a + 2 * h), xi = LOAD(im + a + 2 * h);                     \
  V c1r = xr + tr, c1i = xi + ti, d1r = xr - tr, d1i = xi - ti;             \
  wr = LOAD(tw_re + 2 * h + k), wi = LOAD(tw_im + 2 * h + k);               \
  tr = c1r * wr - c1i * wi, ti = c1r * wi + c1i * wr;                       \
  STORE(re + a, a1r + tr), STORE(im + a, a1i + ti);                         \
  STORE(re + a + 2 * h, a1r - tr), STORE(im + a + 2 * h, a1i - ti);         \
  /* d1 * w2 * -i: (x + iy) * -i = y - ix */                                \
  tr = d1r * wr - d1i * wi, ti = d1r * wi + d1i * wr;                       \
  STORE(re + a + h, b1r + ti), STORE(im + a + h, b1i - tr);                 \
  STORE(re + a + 3 * h, b1r - ti), STORE(im + a + 3 * h, b1i + tr);

static float sload(const float *p) { return *p; }
static void sstore(float *p, float v) { *p = v; }

static void stage4(float *restrict re, float *restrict im, size_t m, size_t h,
                   const float *restrict tw_re, const float *restrict tw_im) {
  for (size_t i = 0; i < m; i += 4 * h) {
    size_t k = 0;
    for (; k + UDSP_W <= h; k += UDSP_W) {
      size_t a = i + k;
      RADIX4_BODY(UDSP_VF, udsp_vload, udsp_vstore)
    }
#ifdef UDSP_HAVE_V4
    for (; k + 4 <= h; k += 4) {
      size_t a = i + k;
      RADIX4_BODY(UDSP_V4, udsp_v4load, udsp_v4store)
    }
#endif
    for (; k < h; k++) {
      size_t a = i + k;
      RADIX4_BODY(float, sload, sstore)
    }
  }
}

/* the butterflies of an m-point transform whose input is already in
 * bit-reversed order; tw_re/tw_im: a table's halves (any n >= m: stage h's
 * twiddles sit at [h + k] whatever the table's size) */
static void butterflies(float *restrict re, float *restrict im, size_t m,
                        const float *restrict tw_re,
                        const float *restrict tw_im, int first) {
  if (m == 2) {
    float r = re[0], i = im[0];
    re[0] = r + re[1], im[0] = i + im[1];
    re[1] = r - re[1], im[1] = i - im[1];
    return;
  }
  /* stages h = 1 and 2 as one radix-4 pass: their twiddles are 1 and -i,
   * so no multiply at all (udsp_fft_real() does this pass itself, fused
   * with its input's permutation: first = 0) */
  for (size_t i = 0; first && i < m; i += 4) {
    float r0 = re[i] + re[i + 1], i0 = im[i] + im[i + 1];
    float r1 = re[i] - re[i + 1], i1 = im[i] - im[i + 1];
    float r2 = re[i + 2] + re[i + 3], i2 = im[i + 2] + im[i + 3];
    float r3 = re[i + 2] - re[i + 3], i3 = im[i + 2] - im[i + 3];
    re[i] = r0 + r2, im[i] = i0 + i2;
    re[i + 2] = r0 - r2, im[i + 2] = i0 - i2;
    /* (r3 + i*i3) * -i = i3 - i*r3 */
    re[i + 1] = r1 + i3, im[i + 1] = i1 - r3;
    re[i + 3] = r1 - i3, im[i + 3] = i1 + r3;
  }
  size_t h = 4;
  for (; 4 * h <= m; h *= 4) {
    stage4(re, im, m, h, tw_re, tw_im);
  }
  if (h < m) {
    stage2(re, im, m, h, tw_re, tw_im);
  }
}

void udsp_fft_forward_table(float *re, float *im, size_t n, const float *table) {
  if (n < 2) {
    return;
  }
  int bits = log2n(n);
  for (size_t i = 1; i < n; i++) {
    size_t j = bits <= 16 ? bitrev(i, bits) : 0;
    if (bits > 16) {
      for (int b = 0; b < bits; b++) {
        j |= ((i >> b) & 1) << (bits - 1 - b);
      }
    }
    if (i < j) {
      swapf(&re[i], &re[j]);
      swapf(&im[i], &im[j]);
    }
  }
  butterflies(re, im, n, table, table + n, 1);
}

#ifdef UDSP_HAVE_V4
/* udsp_fft_real()'s input split into z, put in bit-reversed order and put
 * through the first radix-4 pass, in one trip. Group g of that pass (the
 * points 4g..4g+3 in bit-reversed order) holds z[r], z[r + m/2],
 * z[r + m/4], z[r + 3m/4] with r = the bits-2 reversal of g: so four lanes
 * of consecutive r come from four contiguous stretches of the input, the
 * butterflies run across lanes, and a 4x4 transpose turns the lanes into
 * one 16-byte store per group. The scalar path's scattered point-by-point
 * stores, a third of the whole transform's time, become a quarter as many
 * vector stores. */
static inline __attribute__((always_inline)) void
first_pass_real(const float *restrict in, const float *restrict w,
                float *restrict re, float *restrict im, size_t m, int bits) {
  size_t q = m / 4;
  for (size_t r = 0; r < q; r += 4) {
    UDSP_V4 xr[4], xi[4];
    static const size_t QUARTER[4] = {0, 2, 1, 3}; /* x_t = z[r + t'q] */
    UDSP_UNROLL
    for (int t = 0; t < 4; t++) {
      size_t at = 2 * (r + QUARTER[t] * q);
      UDSP_V4 a = udsp_v4load(in + at), b = udsp_v4load(in + at + 4);
      if (w) { /* the window applied on the way in: nothing stored for it */
        a *= udsp_v4load(w + at);
        b *= udsp_v4load(w + at + 4);
      }
      xr[t] = UDSP_SHUF4(a, b, 0, 2, 4, 6);
      xi[t] = UDSP_SHUF4(a, b, 1, 3, 5, 7);
    }
    UDSP_V4 r0 = xr[0] + xr[1], i0 = xi[0] + xi[1];
    UDSP_V4 r1 = xr[0] - xr[1], i1 = xi[0] - xi[1];
    UDSP_V4 r2 = xr[2] + xr[3], i2 = xi[2] + xi[3];
    UDSP_V4 r3 = xr[2] - xr[3], i3 = xi[2] - xi[3];
    /* y0..y3 per lane; (r3 + i*i3) * -i = i3 - i*r3 */
    UDSP_V4 yr[4] = {r0 + r2, r1 + i3, r0 - r2, r1 - i3};
    UDSP_V4 yi[4] = {i0 + i2, i1 - r3, i0 - i2, i1 + r3};
    UDSP_V4 *ys[2] = {yr, yi};
    float *outs[2] = {re, im};
    UDSP_UNROLL
    for (int c = 0; c < 2; c++) {
      UDSP_V4 *y = ys[c];
      UDSP_V4 t0 = UDSP_SHUF4(y[0], y[1], 0, 4, 1, 5);
      UDSP_V4 t1 = UDSP_SHUF4(y[0], y[1], 2, 6, 3, 7);
      UDSP_V4 t2 = UDSP_SHUF4(y[2], y[3], 0, 4, 1, 5);
      UDSP_V4 t3 = UDSP_SHUF4(y[2], y[3], 2, 6, 3, 7);
      UDSP_V4 g[4] = {UDSP_SHUF4(t0, t2, 0, 1, 4, 5),
                      UDSP_SHUF4(t0, t2, 2, 3, 6, 7),
                      UDSP_SHUF4(t1, t3, 0, 1, 4, 5),
                      UDSP_SHUF4(t1, t3, 2, 3, 6, 7)};
      UDSP_UNROLL
      for (int j = 0; j < 4; j++) {
        udsp_v4store(outs[c] + 4 * bitrev(r + (size_t)j, bits - 2), g[j]);
      }
    }
  }
}
#endif

/* udsp_fft_real(), with in times w (NULL: no window) */
static inline __attribute__((always_inline)) void
fft_real(const float *in, const float *w, float *re, float *im, size_t n,
         const float *table) {
  size_t m = n / 2;
  int bits = log2n(m);
  /* z[k] = in[2k] + i*in[2k+1], stored bit-reversed right away: the
   * permutation costs nothing on top of the copy that splits the input */
#ifdef UDSP_HAVE_V4
  if (m >= 16 && bits <= 18) {
    first_pass_real(in, w, re, im, m, bits);
    butterflies(re, im, m, table, table + n, 0);
  } else
#endif
  {
    if (bits <= 16) {
      for (size_t k = 0; k < m; k++) {
        size_t j = bitrev(k, bits);
        re[j] = w ? in[2 * k] * w[2 * k] : in[2 * k];
        im[j] = w ? in[2 * k + 1] * w[2 * k + 1] : in[2 * k + 1];
      }
    } else {
      for (size_t k = 0; k < m; k++) {
        re[k] = w ? in[2 * k] * w[2 * k] : in[2 * k];
        im[k] = w ? in[2 * k + 1] * w[2 * k + 1] : in[2 * k + 1];
      }
      for (size_t i = 1; i < m; i++) {
        size_t j = 0;
        for (int b = 0; b < bits; b++) {
          j |= ((i >> b) & 1) << (bits - 1 - b);
        }
        if (i < j) {
          swapf(&re[i], &re[j]);
          swapf(&im[i], &im[j]);
        }
      }
    }
    butterflies(re, im, m, table, table + n, 1);
  }

  /* untangle the two real transforms z carries (even samples in its real
   * part, odd in its imaginary): with Z the transform of z,
   *   X[k]   = E + W^k O
   *   X[m-k] = conj(E - W^k O)
   *   E = (Z[k] + conj(Z[m-k])) / 2,  O = -i (Z[k] - conj(Z[m-k])) / 2,
   * W^k = e^(-2 pi i k / n): the table's last stage (h = m) */
  const float *wr = table + m, *wi = table + n + m;
  float r0 = re[0], i0 = im[0];
  re[0] = r0 + i0, im[0] = 0.0f;
  re[m] = r0 - i0, im[m] = 0.0f;
  size_t k = 1;
  /* W bins from k up and W from m-k down per step, the latter's lanes
   * reversed so lane j of both is one pair (k+j, m-k-j); while the two
   * blocks don't meet */
  for (; 2 * (k + UDSP_W - 1) < m; k += UDSP_W) {
    size_t d = m - k - (UDSP_W - 1);
    UDSP_VF zr = udsp_vload(re + k), zi = udsp_vload(im + k);
    UDSP_VF yr = udsp_vrev(udsp_vload(re + d)), yi = udsp_vrev(udsp_vload(im + d));
    UDSP_VF er = (zr + yr) * 0.5f, ei = (zi - yi) * 0.5f;
    UDSP_VF or_ = (zi + yi) * 0.5f, oi = (yr - zr) * 0.5f;
    UDSP_VF c = udsp_vload(wr + k), sn = udsp_vload(wi + k);
    UDSP_VF tr = c * or_ - sn * oi, ti = c * oi + sn * or_;
    udsp_vstore(re + k, er + tr);
    udsp_vstore(im + k, ei + ti);
    udsp_vstore(re + d, udsp_vrev(er - tr));
    udsp_vstore(im + d, udsp_vrev(ti - ei));
  }
  for (; k <= m / 2; k++) {
    float zr = re[k], zi = im[k], yr = re[m - k], yi = im[m - k];
    float er = 0.5f * (zr + yr), ei = 0.5f * (zi - yi);
    float or_ = 0.5f * (zi + yi), oi = 0.5f * (yr - zr);
    float tr = wr[k] * or_ - wi[k] * oi, ti = wr[k] * oi + wi[k] * or_;
    re[k] = er + tr, im[k] = ei + ti;
    re[m - k] = er - tr, im[m - k] = ti - ei;
  }
}

void udsp_fft_real(const float *in, float *re, float *im, size_t n,
                   const float *table) {
  fft_real(in, NULL, re, im, n, table);
}

void udsp_fft_real_windowed(const float *in, const float *w, float *re,
                            float *im, size_t n, const float *table) {
  fft_real(in, w, re, im, n, table);
}

void udsp_fft_power(const float *re, const float *im, size_t n, float *out) {
  size_t i = 0;
  for (; i + UDSP_W <= n; i += UDSP_W) {
    UDSP_VF r = udsp_vload(re + i), m = udsp_vload(im + i);
    udsp_vstore(out + i, r * r + m * m);
  }
  for (; i < n; i++) {
    out[i] = re[i] * re[i] + im[i] * im[i];
  }
}

void udsp_window_hann_table(float *w, size_t n) {
  if (n == 1) {
    w[0] = 1.0f;
    return;
  }
  for (size_t i = 0; i < n; i++) {
    w[i] = (float)(0.5 - 0.5 * cos(2.0 * M_PI * (double)i / (double)(n - 1)));
  }
}

void udsp_window_apply(const float *in, const float *w, float *out, size_t n) {
  size_t i = 0;
  for (; i + UDSP_W <= n; i += UDSP_W) {
    udsp_vstore(out + i, udsp_vload(in + i) * udsp_vload(w + i));
  }
  for (; i < n; i++) {
    out[i] = in[i] * w[i];
  }
}
