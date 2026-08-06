#include "udsp_test.h"

#include "analyze.h"
#include "spectrogram.h"
#include "spectrum.h"

#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define UDSP_SPECTRUM_TEST_BARS 28

static void test_silence_is_floor_every_bar(void) {
  static float silence[UDSP_FFT_N];
  static float re[UDSP_FFT_N], im[UDSP_FFT_N];
  float bars[UDSP_SPECTRUM_TEST_BARS];
  udsp_spectrum_bars(silence, 44100, re, im, bars, UDSP_SPECTRUM_TEST_BARS);
  for (int i = 0; i < UDSP_SPECTRUM_TEST_BARS; i++) {
    UDSP_TEST_ASSERT(bars[i] == UDSP_DBFS_FLOOR, "silence reads the dBFS floor in bar %d", i);
  }
}

static void test_tone_peaks_in_expected_bar(void) {
  static float tone[UDSP_FFT_N];
  static float re[UDSP_FFT_N], im[UDSP_FFT_N];
  int sample_rate = 44100;
  float freq = 1000.0f;
  for (int i = 0; i < UDSP_FFT_N; i++) {
    tone[i] = sinf(2.0f * (float)M_PI * freq * (float)i / (float)sample_rate);
  }

  float bars[UDSP_SPECTRUM_TEST_BARS];
  udsp_spectrum_bars(tone, sample_rate, re, im, bars, UDSP_SPECTRUM_TEST_BARS);

  int loudest = 0;
  for (int i = 1; i < UDSP_SPECTRUM_TEST_BARS; i++) {
    if (bars[i] > bars[loudest]) {
      loudest = i;
    }
  }

  /* log-spaced from 20Hz to Nyquist (22050Hz): 1kHz lands at
   * log(1000/20)/log(22050/20) * 28 ~= bar 15.6, so 15 -- this pins the
   * exact placement so a future change to the log-spacing math (or the
   * bin/frequency conversion) that shifts it gets caught immediately. */
  UDSP_TEST_ASSERT(loudest == 15, "a 1kHz tone's loudest bar lands at bar 15/28 (got %d)", loudest);
  UDSP_TEST_ASSERT(bars[loudest] > -20.0f, "the loudest bar reads well above the floor for a full-scale tone");
}

static void test_bin_count_is_flexible(void) {
  static float tone[UDSP_FFT_N];
  static float re[UDSP_FFT_N], im[UDSP_FFT_N];
  for (int i = 0; i < UDSP_FFT_N; i++) {
    tone[i] = sinf(2.0f * (float)M_PI * 1000.0f * (float)i / 44100.0f);
  }
  float bars[UDSP_SPECTROGRAM_BINS];
  udsp_spectrum_bars(tone, 44100, re, im, bars, UDSP_SPECTROGRAM_BINS);
  int any_above_floor = 0;
  for (int i = 0; i < UDSP_SPECTROGRAM_BINS; i++) {
    if (bars[i] > UDSP_DBFS_FLOOR) {
      any_above_floor = 1;
    }
  }
  UDSP_TEST_ASSERT(any_above_floor, "udsp_spectrum_bars() works with a different n_bars (spectrogram's bin count)");
}

int main(void) {
  test_silence_is_floor_every_bar();
  test_tone_peaks_in_expected_bar();
  test_bin_count_is_flexible();
  return udsp_test_summary();
}
