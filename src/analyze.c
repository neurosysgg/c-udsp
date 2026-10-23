#include "analyze.h"

#include "simd.h"

#include <math.h>

/* floats summed in float lanes before the chunk's total goes into the
 * double: a float lane stays precise over a chunk this long, and the
 * double keeps a whole take's sum (millions of samples) precise */
#ifndef RMS_CHUNK /* overridable for bench/sweep.py */
#define RMS_CHUNK 4096
#endif

static double sumsq_chunk(const float *s, size_t n) {
  UDSP_VF acc[UDSP_ACC];
  for (int a = 0; a < UDSP_ACC; a++) {
    acc[a] = udsp_vsplat(0.0f);
  }
  size_t i = 0;
  for (; i + UDSP_ACC * UDSP_W <= n; i += UDSP_ACC * UDSP_W) {
    UDSP_UNROLL
    for (int a = 0; a < UDSP_ACC; a++) {
      UDSP_VF v = udsp_vload(s + i + (size_t)a * UDSP_W);
      acc[a] += v * v;
    }
  }
  for (int a = 1; a < UDSP_ACC; a++) {
    acc[0] += acc[a];
  }
  double sum = udsp_vsum(acc[0]);
  for (; i < n; i++) {
    sum += (double)s[i] * s[i];
  }
  return sum;
}

static double sum_chunk(const float *s, size_t n) {
  UDSP_VF acc[UDSP_ACC];
  for (int a = 0; a < UDSP_ACC; a++) {
    acc[a] = udsp_vsplat(0.0f);
  }
  size_t i = 0;
  for (; i + UDSP_ACC * UDSP_W <= n; i += UDSP_ACC * UDSP_W) {
    UDSP_UNROLL
    for (int a = 0; a < UDSP_ACC; a++) {
      acc[a] += udsp_vload(s + i + (size_t)a * UDSP_W);
    }
  }
  for (int a = 1; a < UDSP_ACC; a++) {
    acc[0] += acc[a];
  }
  double sum = udsp_vsum(acc[0]);
  for (; i < n; i++) {
    sum += s[i];
  }
  return sum;
}

double udsp_analyze_sum(const float *x, size_t n) {
  double sum = 0.0;
  for (size_t i = 0; i < n; i += RMS_CHUNK) {
    sum += sum_chunk(x + i, n - i < RMS_CHUNK ? n - i : RMS_CHUNK);
  }
  return sum;
}

float udsp_analyze_rms(const float *samples, size_t n) {
  if (n == 0) {
    return 0.0f;
  }
  double sumsq = 0.0;
  for (size_t i = 0; i < n; i += RMS_CHUNK) {
    sumsq += sumsq_chunk(samples + i, n - i < RMS_CHUNK ? n - i : RMS_CHUNK);
  }
  return (float)sqrt(sumsq / (double)n);
}

float udsp_analyze_peak(const float *samples, size_t n) {
  UDSP_VF acc[UDSP_ACC];
  for (int a = 0; a < UDSP_ACC; a++) {
    acc[a] = udsp_vsplat(0.0f);
  }
  size_t i = 0;
  for (; i + UDSP_ACC * UDSP_W <= n; i += UDSP_ACC * UDSP_W) {
    UDSP_UNROLL
    for (int a = 0; a < UDSP_ACC; a++) {
      acc[a] = udsp_vmax(
          udsp_vabs(udsp_vload(samples + i + (size_t)a * UDSP_W)), acc[a]);
    }
  }
  for (int a = 1; a < UDSP_ACC; a++) {
    acc[0] = udsp_vmax(acc[a], acc[0]);
  }
  float peak = udsp_vhmax(acc[0]);
  for (; i < n; i++) {
    float a = fabsf(samples[i]);
    peak = a > peak ? a : peak;
  }
  return peak;
}

void udsp_analyze_minmax(const float *samples, size_t n, float *min_out,
                         float *max_out) {
  if (n == 0) {
    *min_out = *max_out = 0.0f;
    return;
  }
  /* min and max are two chains already: half the accumulators each */
  enum { ACC = UDSP_ACC / 2 > 0 ? UDSP_ACC / 2 : 1 };
  UDSP_VF lo[ACC], hi[ACC];
  for (int a = 0; a < ACC; a++) {
    lo[a] = hi[a] = udsp_vsplat(samples[0]);
  }
  size_t i = 0;
  for (; i + ACC * UDSP_W <= n; i += ACC * UDSP_W) {
    UDSP_UNROLL
    for (int a = 0; a < ACC; a++) {
      UDSP_VF v = udsp_vload(samples + i + (size_t)a * UDSP_W);
      lo[a] = udsp_vmin(v, lo[a]);
      hi[a] = udsp_vmax(v, hi[a]);
    }
  }
  for (int a = 1; a < ACC; a++) {
    lo[0] = udsp_vmin(lo[a], lo[0]);
    hi[0] = udsp_vmax(hi[a], hi[0]);
  }
  float mn = udsp_vhmin(lo[0]), mx = udsp_vhmax(hi[0]);
  for (; i < n; i++) {
    mn = samples[i] < mn ? samples[i] : mn;
    mx = samples[i] > mx ? samples[i] : mx;
  }
  *min_out = mn;
  *max_out = mx;
}

void udsp_mono_downmix(const float *interleaved, size_t frames, int channels,
                       float *mono_out) {
  if (channels <= 1) {
    memcpy(mono_out, interleaved, frames * sizeof(float));
    return;
  }
  if (channels == 2) {
    /* (l + r) * 0.5 in float is the double average rounded: the halving
     * is exact */
    size_t i = 0;
#ifdef UDSP_HAVE_V4
    for (; i + 4 <= frames; i += 4) {
      UDSP_V4 a = udsp_v4load(interleaved + 2 * i);
      UDSP_V4 b = udsp_v4load(interleaved + 2 * i + 4);
      udsp_v4store(
          mono_out + i,
          (UDSP_SHUF4(a, b, 0, 2, 4, 6) + UDSP_SHUF4(a, b, 1, 3, 5, 7)) * 0.5f);
    }
#endif
    for (; i < frames; i++) {
      mono_out[i] = (interleaved[2 * i] + interleaved[2 * i + 1]) * 0.5f;
    }
    return;
  }
  for (size_t i = 0; i < frames; i++) {
    double sum = 0.0;
    const float *frame = interleaved + i * (size_t)channels;
    for (int c = 0; c < channels; c++) {
      sum += frame[c];
    }
    mono_out[i] = (float)(sum / (double)channels);
  }
}

float udsp_analyze_dbfs(float amplitude) {
  if (amplitude <= 0.0f) {
    return UDSP_DBFS_FLOOR;
  }
  float db = 20.0f * log10f(amplitude);
  return db < UDSP_DBFS_FLOOR ? UDSP_DBFS_FLOOR : db;
}
