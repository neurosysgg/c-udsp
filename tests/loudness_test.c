#include "udsp_test.h"

#include "analysis/loudness.h"

#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define UDSP_LOUDNESS_TEST_N 8192

static void test_silence_is_floor(void) {
  float silence[UDSP_LOUDNESS_TEST_N] = {0};
  float scratch[UDSP_LOUDNESS_TEST_N];
  float lufs =
      udsp_loudness_lufs_momentary(silence, UDSP_LOUDNESS_TEST_N, scratch);
  UDSP_TEST_ASSERT(lufs == UDSP_LUFS_FLOOR,
                   "pure silence reads exactly the LUFS floor");
}

static void test_empty_buffer_is_floor(void) {
  float scratch[1];
  UDSP_TEST_ASSERT(udsp_loudness_lufs_momentary(NULL, 0, scratch) ==
                       UDSP_LUFS_FLOOR,
                   "n=0 reads the floor rather than dividing by zero");
}

static void test_tone_reads_above_floor(void) {
  static float sig[UDSP_LOUDNESS_TEST_N];
  static float scratch[UDSP_LOUDNESS_TEST_N];
  for (int i = 0; i < UDSP_LOUDNESS_TEST_N; i++) {
    sig[i] = sinf(2.0f * (float)M_PI * 1000.0f * (float)i / 44100.0f);
  }
  float lufs = udsp_loudness_lufs_momentary(sig, UDSP_LOUDNESS_TEST_N, scratch);
  UDSP_TEST_ASSERT(lufs > UDSP_LUFS_FLOOR + 10.0f,
                   "a full-scale 1kHz tone reads well above the floor");
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
  float lufs_quiet =
      udsp_loudness_lufs_momentary(quiet, UDSP_LOUDNESS_TEST_N, scratch);
  float lufs_loud =
      udsp_loudness_lufs_momentary(loud, UDSP_LOUDNESS_TEST_N, scratch);
  UDSP_TEST_ASSERT(lufs_loud > lufs_quiet,
                   "a 10x louder signal reads a higher LUFS value");
}

/* the old two-pass Direct Form I cascade, as the reference */
static double ref_lufs(const float *x, size_t n) {
  static const double B[2][5] = {
      {1.53512485958697, -2.69169618940638, 1.19839281085285, -1.69065929318241,
       0.73248077421585},
      {1.0, -2.0, 1.0, -1.99004745483398, 0.99007225036621}};
  static float buf[48000];
  for (size_t i = 0; i < n; i++) {
    buf[i] = x[i];
  }
  for (int st = 0; st < 2; st++) {
    double x1 = 0, x2 = 0, y1 = 0, y2 = 0;
    for (size_t i = 0; i < n; i++) {
      double x0 = buf[i];
      double y0 = B[st][0] * x0 + B[st][1] * x1 + B[st][2] * x2 -
                  B[st][3] * y1 - B[st][4] * y2;
      buf[i] = (float)y0;
      x2 = x1, x1 = x0, y2 = y1, y1 = y0;
    }
  }
  double ss = 0;
  for (size_t i = 0; i < n; i++) {
    ss += (double)buf[i] * buf[i];
  }
  return -0.691 + 10.0 * log10(ss / (double)n);
}

static void test_calibration_and_reference(void) {
  static float x[19200], scratch[19200];
  /* BS.1770's calibration: a 997 Hz sine at 0 dBFS in one channel reads
   * -3.01 LKFS */
  for (int i = 0; i < 19200; i++) {
    x[i] = (float)sin(2 * M_PI * 997.0 * i / 48000.0);
  }
  float lufs = udsp_loudness_lufs_momentary(x, 19200, scratch);
  UDSP_TEST_ASSERT(fabsf(lufs + 3.01f) < 0.05f,
                   "a 0 dBFS 997 Hz sine reads -3.01 LUFS (%.3f)", lufs);
  unsigned int seed = 11;
  for (int i = 0; i < 19200; i++) {
    seed = seed * 1103515245u + 12345u;
    x[i] = 0.3f * ((float)((seed >> 8) & 0xffff) / 32768.0f - 1.0f);
  }
  lufs = udsp_loudness_lufs_momentary(x, 19200, scratch);
  double want = ref_lufs(x, 19200);
  UDSP_TEST_ASSERT(fabs(lufs - want) < 0.001,
                   "the one-pass cascade reads the two-pass one's "
                   "LUFS on noise (%.4f vs %.4f)",
                   lufs, want);
}

int main(void) {
  test_calibration_and_reference();
  test_silence_is_floor();
  test_empty_buffer_is_floor();
  test_tone_reads_above_floor();
  test_louder_signal_reads_higher();
  return udsp_test_summary();
}
