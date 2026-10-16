#include "r128.h"

#include "analyze.h"

#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define R128_BIN_LOW -70.0
#define R128_BIN_STEP 0.1

/* the prototypes' parameters as libebur128 and ffmpeg's ebur128 fit them
 * to the spec's 48 kHz table: a high shelf (+4 dB above ~1.7 kHz, the
 * head) and a high-pass at ~38 Hz (the RLB curve) */
void udsp_r128_kweight(int sample_rate, double k[2][5]) {
  double fs = (double)sample_rate;

  double f0 = 1681.974450955533, g = 3.999843853973347, q = 0.7071752369554196;
  double kk = tan(M_PI * f0 / fs);
  double vh = pow(10.0, g / 20.0);
  double vb = pow(vh, 0.4996667741545416);
  double a0 = 1.0 + kk / q + kk * kk;
  k[0][0] = (vh + vb * kk / q + kk * kk) / a0;
  k[0][1] = 2.0 * (kk * kk - vh) / a0;
  k[0][2] = (vh - vb * kk / q + kk * kk) / a0;
  k[0][3] = 2.0 * (kk * kk - 1.0) / a0;
  k[0][4] = (1.0 - kk / q + kk * kk) / a0;

  f0 = 38.13547087602444;
  q = 0.5003270373238773;
  kk = tan(M_PI * f0 / fs);
  a0 = 1.0 + kk / q + kk * kk;
  k[1][0] = 1.0;
  k[1][1] = -2.0;
  k[1][2] = 1.0;
  k[1][3] = 2.0 * (kk * kk - 1.0) / a0;
  k[1][4] = (1.0 - kk / q + kk * kk) / a0;
}

static void default_weights(float *w, int channels) {
  static const float surround[UDSP_R128_MAX_CHANNELS + 1][UDSP_R128_MAX_CHANNELS] = {
      [4] = {1.0f, 1.0f, 1.41f, 1.41f},
      [5] = {1.0f, 1.0f, 1.0f, 1.41f, 1.41f},
      [6] = {1.0f, 1.0f, 1.0f, 0.0f, 1.41f, 1.41f},
      [7] = {1.0f, 1.0f, 1.0f, 0.0f, 1.41f, 1.41f, 1.41f},
      [8] = {1.0f, 1.0f, 1.0f, 0.0f, 1.41f, 1.41f, 1.41f, 1.41f},
  };
  for (int c = 0; c < channels; c++) {
    w[c] = channels < 4 ? 1.0f : surround[channels][c];
  }
}

void udsp_r128_init(UDSP_R128 *m, int sample_rate, int channels) {
  memset(m, 0, sizeof(*m));
  udsp_r128_kweight(sample_rate, m->k);
  m->channels = channels;
  default_weights(m->weight, channels);
  /* exact at every common rate; 11025 Hz rounds to 1103 */
  m->hop = (size_t)((double)sample_rate / 10.0 + 0.5);
}

/* one sample through one channel's K-weighting, both stages in
 * transposed Direct Form II as loudness.c's; doubles, since stage 2's
 * poles sit at ~0.995 where float state drifts. z is the four states */
static inline double kweight_step(const double k[2][5], double z[4],
                                  double v) {
  const double *f = k[0], *g = k[1];
  double y = f[0] * v + z[0];
  z[0] = f[1] * v - f[3] * y + z[1];
  z[1] = f[2] * v - f[4] * y;
  double o = g[0] * y + z[2];
  z[2] = g[1] * y - g[3] * o + z[3];
  z[3] = g[2] * y - g[4] * o;
  return o;
}

/* channels a and b's K-weighted sums of squares over n frames (b < 0:
 * a alone), in one loop: each channel is one serial chain of
 * multiply-adds, two of them overlap (~1.8x over a channel at a time),
 * the same operations per channel, so the same bits */
static void channel_sumsq(const double k[2][5], double (*z)[4],
                          const float *x, size_t n, int stride, int a, int b,
                          double *sa, double *sb) {
  double za[4], zb[4] = {0}, qa = 0.0, qb = 0.0;
  memcpy(za, z[a], sizeof(za));
  if (b < 0) {
    for (size_t i = 0; i < n; i++) {
      double o = kweight_step(k, za, x[i * (size_t)stride + (size_t)a]);
      qa += o * o;
    }
  } else {
    memcpy(zb, z[b], sizeof(zb));
    for (size_t i = 0; i < n; i++) {
      const float *frame = x + i * (size_t)stride;
      double oa = kweight_step(k, za, frame[a]);
      double ob = kweight_step(k, zb, frame[b]);
      qa += oa * oa;
      qb += ob * ob;
    }
    memcpy(z[b], zb, sizeof(zb));
  }
  memcpy(z[a], za, sizeof(za));
  *sa = qa;
  *sb = qb;
}

static void add_block(UDSP_R128_HIST *h, double energy) {
  if (energy <= 0.0) {
    return;
  }
  double l = -0.691 + 10.0 * log10(energy);
  if (l <= R128_BIN_LOW) {
    return;
  }
  int bin = (int)((l - R128_BIN_LOW) / R128_BIN_STEP);
  if (bin >= UDSP_R128_BINS) {
    bin = UDSP_R128_BINS - 1;
  }
  h->count[bin]++;
  h->energy[bin] += energy;
}

void udsp_r128_add(UDSP_R128 *m, const float *interleaved, size_t frames) {
  float peak = udsp_analyze_peak(interleaved, frames * (size_t)m->channels);
  if (peak > m->hist.peak) {
    m->hist.peak = peak;
  }
  while (frames > 0) {
    size_t n = m->hop - m->acc_n;
    if (n > frames) {
      n = frames;
    }
    /* the weighted channels two at a time, added in channel order */
    int c = 0;
    for (;;) {
      while (c < m->channels && m->weight[c] == 0.0f) {
        c++;
      }
      int d = c + 1;
      while (d < m->channels && m->weight[d] == 0.0f) {
        d++;
      }
      if (c >= m->channels) {
        break;
      }
      double sc, sd;
      channel_sumsq(m->k, m->z, interleaved, n, m->channels, c,
                    d < m->channels ? d : -1, &sc, &sd);
      m->acc += m->weight[c] * sc;
      if (d >= m->channels) {
        break;
      }
      m->acc += m->weight[d] * sd;
      c = d + 1;
    }
    m->acc_n += n;
    interleaved += n * (size_t)m->channels;
    frames -= n;
    if (m->acc_n == m->hop) {
      m->sub[m->hops % 4] = m->acc;
      m->hops++;
      m->acc = 0.0;
      m->acc_n = 0;
      if (m->hops >= 4) {
        double sum = m->sub[0] + m->sub[1] + m->sub[2] + m->sub[3];
        add_block(&m->hist, sum / (double)(4 * m->hop));
      }
    }
  }
}

static double to_lufs(double energy) { return -0.691 + 10.0 * log10(energy); }

double udsp_r128_integrated(const UDSP_R128_HIST *h) {
  uint64_t count = 0;
  double energy = 0.0;
  for (int i = 0; i < UDSP_R128_BINS; i++) {
    count += h->count[i];
    energy += h->energy[i];
  }
  if (count == 0) {
    return -HUGE_VAL;
  }
  double gate = to_lufs(energy / (double)count) - 10.0;

  count = 0;
  energy = 0.0;
  for (int i = 0; i < UDSP_R128_BINS; i++) {
    if (h->count[i] == 0) {
      continue;
    }
    double low = R128_BIN_LOW + i * R128_BIN_STEP;
    /* the bin the gate falls in: in or out by its own blocks' loudness */
    if (low + R128_BIN_STEP <= gate ||
        (low < gate && to_lufs(h->energy[i] / (double)h->count[i]) <= gate)) {
      continue;
    }
    count += h->count[i];
    energy += h->energy[i];
  }
  /* the relative gate sits 10 LU under a mean the loudest block is at or
   * above, so something always passes it */
  return to_lufs(energy / (double)count);
}

void udsp_r128_hist_add(UDSP_R128_HIST *dst, const UDSP_R128_HIST *src) {
  for (int i = 0; i < UDSP_R128_BINS; i++) {
    dst->count[i] += src->count[i];
    dst->energy[i] += src->energy[i];
  }
  if (src->peak > dst->peak) {
    dst->peak = src->peak;
  }
}
