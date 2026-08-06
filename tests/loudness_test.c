#include "udsp_test.h"

#include "loudness.h"

#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define UDSP_LOUDNESS_TEST_N 8192

static void test_silence_is_floor(void) {
  float silence[UDSP_LOUDNESS_TEST_N] = {0};
  float scratch[UDSP_LOUDNESS_TEST_N];
  float lufs = udsp_loudness_lufs_momentary(silence, UDSP_LOUDNESS_TEST_N, scratch);
  UDSP_TEST_ASSERT(lufs == UDSP_LUFS_FLOOR, "pure silence reads exactly the LUFS floor");
}

static void test_empty_buffer_is_floor(void) {
  float scratch[1];
  UDSP_TEST_ASSERT(udsp_loudness_lufs_momentary(NULL, 0, scratch) == UDSP_LUFS_FLOOR,
                  "n=0 reads the floor rather than dividing by zero");
}

static void test_tone_reads_above_floor(void) {
  static float sig[UDSP_LOUDNESS_TEST_N];
  static float scratch[UDSP_LOUDNESS_TEST_N];
  for (int i = 0; i < UDSP_LOUDNESS_TEST_N; i++) {
    sig[i] = sinf(2.0f * (float)M_PI * 1000.0f * (float)i / 44100.0f);
  }
  float lufs = udsp_loudness_lufs_momentary(sig, UDSP_LOUDNESS_TEST_N, scratch);
  UDSP_TEST_ASSERT(lufs > UDSP_LUFS_FLOOR + 10.0f, "a full-scale 1kHz tone reads well above the floor");
}

static void test_louder_signal_reads_higher(void) {
  static float quiet[UDSP_LOUDNESS_TEST_N];
  static float loud[UDSP_LOUDNESS_TEST_N];
  static float scratch[UDSP_LOUDNESS_TEST_N];
  for (int i = 0; i < UDSP_LOUDNESS_TEST_N; i++) {
    float s = sinf(2.0f * (float)M_PI * 1000.0f * (float)i / 44100.0f);
    quiet[i] = s * 0.1f;
    loud[i] = s * 1.0f;
  }
  float lufs_quiet = udsp_loudness_lufs_momentary(quiet, UDSP_LOUDNESS_TEST_N, scratch);
  float lufs_loud = udsp_loudness_lufs_momentary(loud, UDSP_LOUDNESS_TEST_N, scratch);
  UDSP_TEST_ASSERT(lufs_loud > lufs_quiet, "a 10x louder signal reads a higher LUFS value");
}

int main(void) {
  test_silence_is_floor();
  test_empty_buffer_is_floor();
  test_tone_reads_above_floor();
  test_louder_signal_reads_higher();
  return udsp_test_summary();
}
