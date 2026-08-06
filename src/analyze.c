#include "analyze.h"

#include <math.h>

float udsp_analyze_rms(const float *samples, size_t n) {
  if (n == 0) {
    return 0.0f;
  }
  double sumsq = 0.0;
  for (size_t i = 0; i < n; i++) {
    double v = samples[i];
    sumsq += v * v;
  }
  return (float)sqrt(sumsq / (double)n);
}

float udsp_analyze_peak(const float *samples, size_t n) {
  float peak = 0.0f;
  for (size_t i = 0; i < n; i++) {
    float a = fabsf(samples[i]);
    if (a > peak) {
      peak = a;
    }
  }
  return peak;
}

void udsp_mono_downmix(const float *interleaved, size_t frames, int channels,
                       float *mono_out) {
  if (channels <= 1) {
    for (size_t i = 0; i < frames; i++) {
      mono_out[i] = interleaved[i];
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
