#include "udsp_test.h"

#include "spectrogram.h"
#include "spectrum.h"

#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void test_init_is_silent(void) {
  UDSP_SPECTROGRAM sg;
  udsp_spectrogram_init(&sg);
  for (int c = 0; c < UDSP_SPECTROGRAM_COLS; c++) {
    for (int b = 0; b < UDSP_SPECTROGRAM_BINS; b++) {
      UDSP_TEST_ASSERT(
          sg.intensity[c][b] == 0,
          "a freshly-initialized spectrogram is all-zero (col %d bin %d)", c,
          b);
    }
  }
}

static void test_push_scrolls_and_appends(void) {
  UDSP_SPECTROGRAM sg;
  udsp_spectrogram_init(&sg);

  static float tone[UDSP_FFT_N];
  static float re[UDSP_FFT_N], im[UDSP_FFT_N];
  for (int i = 0; i < UDSP_FFT_N; i++) {
    tone[i] = sinf(2.0f * (float)M_PI * 1000.0f * (float)i / 44100.0f);
  }

  udsp_spectrogram_push(&sg, tone, 44100, re, im);

  int any_nonzero_in_last_col = 0;
  for (int b = 0; b < UDSP_SPECTROGRAM_BINS; b++) {
    if (sg.intensity[UDSP_SPECTROGRAM_COLS - 1][b] != 0) {
      any_nonzero_in_last_col = 1;
    }
  }
  UDSP_TEST_ASSERT(
      any_nonzero_in_last_col,
      "the newest column (rightmost) holds the just-pushed spectrum");

  unsigned char last_col_snapshot[UDSP_SPECTROGRAM_BINS];
  memcpy(last_col_snapshot, sg.intensity[UDSP_SPECTROGRAM_COLS - 1],
         sizeof(last_col_snapshot));

  udsp_spectrogram_push(&sg, tone, 44100, re, im);
  int scrolled = memcmp(sg.intensity[UDSP_SPECTROGRAM_COLS - 2],
                        last_col_snapshot, sizeof(last_col_snapshot)) == 0;
  UDSP_TEST_ASSERT(
      scrolled,
      "a second push scrolls the first column's data one step toward index 0");
}

int main(void) {
  test_init_is_silent();
  test_push_scrolls_and_appends();
  return udsp_test_summary();
}
