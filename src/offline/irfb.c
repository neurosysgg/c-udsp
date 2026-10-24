#include "offline/irfb.h"

#include <math.h>
#include <string.h>

#include "base/fft.h"

size_t udsp_irfb_frames(size_t len, size_t d, float g, size_t cap) {
  const double a = fabs((double)g);
  if (len >= cap) {
    return cap;
  }
  if (a == 0.0) {
    return len;
  }
  if (a >= 1.0) {
    return cap;
  }
  /* each trip round the loop is len + d frames longer and a times as
   * loud: -90 dB after log(10^-4.5) / log(a) of them */
  const double trips = ceil(-4.5 * 2.302585092994046 / log(a));
  const double n = (double)len + trips * (double)(len + d);
  return n < (double)cap ? (size_t)n : cap;
}

size_t udsp_irfb_size(size_t frames) {
  size_t n = 4;
  while (n < 2 * frames) {
    n *= 2;
  }
  return n;
}

/* h padded to n, its bins 0..n/2 into re / im */
static void spectrum(const float *h, size_t len, size_t stride, size_t n,
                     const float *table, float *re, float *im, float *pad) {
  for (size_t i = 0; i < len; i++) {
    pad[i] = h[i * stride];
  }
  memset(pad + len, 0, (n - len) * sizeof(float));
  udsp_fft_real(pad, re, im, n, table);
}

float udsp_irfb_peak(const float *h, size_t len, size_t stride, size_t n,
                     const float *table, float *scratch) {
  float *re = scratch, *im = re + n / 2 + 1, *pad = im + n / 2 + 1;
  spectrum(h, len, stride, n, table, re, im, pad);
  float top = 0.0f;
  for (size_t k = 0; k <= n / 2; k++) {
    const float m = re[k] * re[k] + im[k] * im[k];
    top = m > top ? m : top;
  }
  return sqrtf(top);
}

void udsp_irfb(const float *h, size_t len, size_t stride, float g, size_t d,
               float *y, size_t frames, size_t ystride, size_t n,
               const float *table, float *scratch) {
  float *re = scratch, *im = re + n / 2 + 1, *pad = im + n / 2 + 1;
  spectrum(h, len, stride, n, table, re, im, pad);
  const size_t half = n / 2;
  for (size_t k = 0; k <= half; k++) {
    /* e^(-i 2 pi k d / n): the table's last stage holds it for j < n / 2,
     * the other half its negative */
    size_t j = (size_t)(((unsigned long long)k * d) % n);
    float wr = 1.0f, wi = 0.0f;
    if (j) {
      const int neg = j >= half;
      j -= neg ? half : 0;
      wr = table[half + j];
      wi = table[n + half + j];
      wr = neg ? -wr : wr;
      wi = neg ? -wi : wi;
    }
    /* Y = H / (1 - g W H) */
    const float hr = re[k], hi = im[k];
    const float lr = g * (wr * hr - wi * hi), li = g * (wr * hi + wi * hr);
    const float dr = 1.0f - lr, di = -li;
    const float den = dr * dr + di * di;
    const float s = den > 0.0f ? 1.0f / (den * (float)n) : 0.0f;
    re[k] = (hr * dr + hi * di) * s;
    im[k] = (hi * dr - hr * di) * s;
  }
  udsp_fft_real_inverse(re, im, pad, n, table);
  for (size_t i = 0; i < frames; i++) {
    y[i * ystride] = pad[i];
  }
}
