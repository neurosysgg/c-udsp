#include "udsp_test.h"

#include "analyze.h"

#include <math.h>

static void test_rms_known_values(void) {
  float zeros[8] = {0};
  UDSP_TEST_ASSERT(udsp_analyze_rms(zeros, 8) == 0.0f, "RMS of silence is 0");

  float ones[8];
  for (int i = 0; i < 8; i++) {
    ones[i] = 1.0f;
  }
  float rms = udsp_analyze_rms(ones, 8);
  UDSP_TEST_ASSERT(fabsf(rms - 1.0f) < 0.0001f, "RMS of a constant 1.0 signal is 1.0, unscaled");

  UDSP_TEST_ASSERT(udsp_analyze_rms(ones, 0) == 0.0f, "RMS of an empty buffer is 0, not a NaN from 0/0");
}

static void test_rms_is_channel_agnostic(void) {
  float interleaved[4] = {1.0f, -1.0f, 1.0f, -1.0f};
  float rms = udsp_analyze_rms(interleaved, 4);
  UDSP_TEST_ASSERT(fabsf(rms - 1.0f) < 0.0001f,
                  "RMS treats interleaved multi-channel samples as one flat stream");
}

static void test_peak(void) {
  float samples[5] = {0.1f, -0.9f, 0.5f, -0.3f, 0.2f};
  UDSP_TEST_ASSERT(fabsf(udsp_analyze_peak(samples, 5) - 0.9f) < 0.0001f,
                  "peak is the max absolute value, sign-independent");
  UDSP_TEST_ASSERT(udsp_analyze_peak(samples, 0) == 0.0f, "peak of an empty buffer is 0");
}

static void test_mono_downmix(void) {
  float stereo[4] = {1.0f, -1.0f, 1.0f, 0.0f};
  float mono[2];
  udsp_mono_downmix(stereo, 2, 2, mono);
  UDSP_TEST_ASSERT(mono[0] == 0.0f, "downmix averages a fully out-of-phase stereo frame to 0");
  UDSP_TEST_ASSERT(fabsf(mono[1] - 0.5f) < 0.0001f, "downmix averages an in-phase stereo frame correctly");

  float mono_src[3] = {0.1f, 0.2f, 0.3f};
  float mono_out[3];
  udsp_mono_downmix(mono_src, 3, 1, mono_out);
  UDSP_TEST_ASSERT(mono_out[0] == 0.1f && mono_out[1] == 0.2f && mono_out[2] == 0.3f,
                  "downmix with channels<=1 is a plain copy");
}

static void test_dbfs(void) {
  UDSP_TEST_ASSERT(udsp_analyze_dbfs(1.0f) < 0.001f && udsp_analyze_dbfs(1.0f) > -0.001f,
                  "full-scale amplitude is ~0 dBFS");
  UDSP_TEST_ASSERT(udsp_analyze_dbfs(0.0f) == UDSP_DBFS_FLOOR, "silence clamps to the dBFS floor, not -inf");
  UDSP_TEST_ASSERT(udsp_analyze_dbfs(-1.0f) == UDSP_DBFS_FLOOR,
                  "a non-positive amplitude (invalid, but shouldn't crash) also clamps to the floor");
}

int main(void) {
  test_rms_known_values();
  test_rms_is_channel_agnostic();
  test_peak();
  test_mono_downmix();
  test_dbfs();
  return udsp_test_summary();
}
