#include "udsp_test.h"

#include "analysis/analyze.h"

#include <math.h>
#include <stdlib.h>

static void test_rms_known_values(void) {
  float zeros[8] = {0};
  UDSP_TEST_ASSERT(udsp_analyze_rms(zeros, 8) == 0.0f, "RMS of silence is 0");

  float ones[8];
  for (int i = 0; i < 8; i++) {
    ones[i] = 1.0f;
  }
  float rms = udsp_analyze_rms(ones, 8);
  UDSP_TEST_ASSERT(fabsf(rms - 1.0f) < 0.0001f,
                   "RMS of a constant 1.0 signal is 1.0, unscaled");

  UDSP_TEST_ASSERT(udsp_analyze_rms(ones, 0) == 0.0f,
                   "RMS of an empty buffer is 0, not a NaN from 0/0");
}

static void test_rms_is_channel_agnostic(void) {
  float interleaved[4] = {1.0f, -1.0f, 1.0f, -1.0f};
  float rms = udsp_analyze_rms(interleaved, 4);
  UDSP_TEST_ASSERT(
      fabsf(rms - 1.0f) < 0.0001f,
      "RMS treats interleaved multi-channel samples as one flat stream");
}

static void test_peak(void) {
  float samples[5] = {0.1f, -0.9f, 0.5f, -0.3f, 0.2f};
  UDSP_TEST_ASSERT(fabsf(udsp_analyze_peak(samples, 5) - 0.9f) < 0.0001f,
                   "peak is the max absolute value, sign-independent");
  UDSP_TEST_ASSERT(udsp_analyze_peak(samples, 0) == 0.0f,
                   "peak of an empty buffer is 0");
}

static void test_mono_downmix(void) {
  float stereo[4] = {1.0f, -1.0f, 1.0f, 0.0f};
  float mono[2];
  udsp_mono_downmix(stereo, 2, 2, mono);
  UDSP_TEST_ASSERT(mono[0] == 0.0f,
                   "downmix averages a fully out-of-phase stereo frame to 0");
  UDSP_TEST_ASSERT(fabsf(mono[1] - 0.5f) < 0.0001f,
                   "downmix averages an in-phase stereo frame correctly");

  float mono_src[3] = {0.1f, 0.2f, 0.3f};
  float mono_out[3];
  udsp_mono_downmix(mono_src, 3, 1, mono_out);
  UDSP_TEST_ASSERT(mono_out[0] == 0.1f && mono_out[1] == 0.2f &&
                       mono_out[2] == 0.3f,
                   "downmix with channels<=1 is a plain copy");
}

static void test_dbfs(void) {
  UDSP_TEST_ASSERT(udsp_analyze_dbfs(1.0f) < 0.001f &&
                       udsp_analyze_dbfs(1.0f) > -0.001f,
                   "full-scale amplitude is ~0 dBFS");
  UDSP_TEST_ASSERT(udsp_analyze_dbfs(0.0f) == UDSP_DBFS_FLOOR,
                   "silence clamps to the dBFS floor, not -inf");
  UDSP_TEST_ASSERT(
      udsp_analyze_dbfs(-1.0f) == UDSP_DBFS_FLOOR,
      "a non-positive amplitude (invalid, but shouldn't crash) also clamps to the floor");
}

static unsigned int g_seed = 5;
static float noise(void) {
  g_seed = g_seed * 1103515245u + 12345u;
  return (float)((g_seed >> 8) & 0xffff) / 32768.0f - 1.0f;
}

/* every length 0..100 (all the vector paths' remainders) against plain
 * loops, the extreme placed at each position in turn */
static void test_vector_paths_every_length(void) {
  float x[101] = {0};
  int bad_peak = 0, bad_minmax = 0, bad_rms = 0, bad_sum = 0;
  for (size_t n = 0; n <= 100; n++) {
    for (size_t at = 0; at < (n ? n : 1); at++) {
      for (size_t i = 0; i < n; i++) {
        x[i] = 0.5f * noise();
      }
      if (n) {
        x[at] = at % 2 ? -0.9f : 0.9f;
      }
      float peak = 0, lo = n ? x[0] : 0, hi = n ? x[0] : 0;
      double ss = 0;
      for (size_t i = 0; i < n; i++) {
        peak = fabsf(x[i]) > peak ? fabsf(x[i]) : peak;
        lo = x[i] < lo ? x[i] : lo;
        hi = x[i] > hi ? x[i] : hi;
        ss += (double)x[i] * x[i];
      }
      float glo, ghi;
      udsp_analyze_minmax(x, n, &glo, &ghi);
      bad_peak += udsp_analyze_peak(x, n) != peak;
      bad_minmax += glo != lo || ghi != hi;
      float rms = n ? (float)sqrt(ss / (double)n) : 0.0f;
      bad_rms += fabsf(udsp_analyze_rms(x, n) - rms) > 1e-6f * (rms + 1e-6f);
      double sum = 0;
      for (size_t i = 0; i < n; i++) {
        sum += x[i];
      }
      bad_sum += fabs(udsp_analyze_sum(x, n) - sum) > 1e-5;
    }
  }
  UDSP_TEST_ASSERT(
      !bad_peak,
      "peak = the plain loop's at every length and position (%d off)",
      bad_peak);
  UDSP_TEST_ASSERT(
      !bad_minmax,
      "minmax = the plain loop's at every length and position (%d off)",
      bad_minmax);
  UDSP_TEST_ASSERT(!bad_rms,
                   "rms = the double loop's to 1e-6 at every length (%d off)",
                   bad_rms);
  UDSP_TEST_ASSERT(!bad_sum,
                   "sum = the double loop's to 1e-5 at every length (%d off)",
                   bad_sum);
  float lo = 1, hi = 1;
  udsp_analyze_minmax(x, 0, &lo, &hi);
  UDSP_TEST_ASSERT(lo == 0.0f && hi == 0.0f, "minmax of nothing is 0, 0");
  float neg[3] = {-0.5f, -0.2f, -0.7f};
  udsp_analyze_minmax(neg, 3, &lo, &hi);
  UDSP_TEST_ASSERT(lo == -0.7f && hi == -0.2f,
                   "minmax of all-negative samples keeps their signs");
}

/* a whole take: float lanes flushed into a double per chunk stay precise
 * where one float sum wouldn't */
static void test_rms_long_take(void) {
  size_t n = 20u * 60 * 48000; /* 20 min at 48 kHz */
  float *x = malloc(n * sizeof(float));
  double ss = 0;
  for (size_t i = 0; i < n; i++) {
    x[i] = 0.25f * noise();
    ss += (double)x[i] * x[i];
  }
  double want = sqrt(ss / (double)n);
  double got = udsp_analyze_rms(x, n);
  UDSP_TEST_ASSERT(fabs(got - want) / want < 1e-6,
                   "rms over 20 min is precise to 1e-6 "
                   "(%.9f vs %.9f)",
                   got, want);
  free(x);
}

static void test_downmix_stereo_exact(void) {
  float st[2 * 37], mono[37];
  for (int i = 0; i < 2 * 37; i++) {
    st[i] = noise();
  }
  udsp_mono_downmix(st, 37, 2, mono);
  int same = 1;
  for (int i = 0; i < 37; i++) {
    same &= mono[i] == (float)(((double)st[2 * i] + st[2 * i + 1]) / 2.0);
  }
  UDSP_TEST_ASSERT(
      same, "the stereo downmix is the double average, rounded, exactly");
}

int main(void) {
  test_vector_paths_every_length();
  test_rms_long_take();
  test_downmix_stereo_exact();
  test_rms_known_values();
  test_rms_is_channel_agnostic();
  test_peak();
  test_mono_downmix();
  test_dbfs();
  return udsp_test_summary();
}
