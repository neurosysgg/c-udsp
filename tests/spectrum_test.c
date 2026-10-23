#include "udsp_test.h"

#include "analysis/analyze.h"
#include "analysis/spectrogram.h"
#include "analysis/spectrum.h"

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
    UDSP_TEST_ASSERT(bars[i] == UDSP_DBFS_FLOOR,
                     "silence reads the dBFS floor in bar %d", i);
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
  UDSP_TEST_ASSERT(loudest == 15,
                   "a 1kHz tone's loudest bar lands at bar 15/28 (got %d)",
                   loudest);
  UDSP_TEST_ASSERT(
      bars[loudest] > -20.0f,
      "the loudest bar reads well above the floor for a full-scale tone");
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
  UDSP_TEST_ASSERT(
      any_above_floor,
      "udsp_spectrum_bars() works with a different n_bars (spectrogram's bin count)");
}

/* the planned path gives the old path's bars: noise and tones, several
 * rates and bar counts */
static void test_planned_matches_bars(void) {
  static float sig[UDSP_FFT_N], re[UDSP_FFT_N], im[UDSP_FFT_N];
  static float scratch[UDSP_SPECTRUM_SCRATCH];
  static UDSP_SPECTRUM_PLAN plan;
  static float want[UDSP_SPECTRUM_MAX_BARS], got[UDSP_SPECTRUM_MAX_BARS];
  static const int RATES[] = {8000, 44100, 48000, 96000};
  static const int BARS[] = {1, 28, 48, 128, 1024};
  unsigned int seed = 3;
  for (int kind = 0; kind < 2; kind++) {
    for (int i = 0; i < UDSP_FFT_N; i++) {
      seed = seed * 1103515245u + 12345u;
      float noise = (float)((seed >> 8) & 0xffff) / 32768.0f - 1.0f;
      sig[i] = kind ? 0.5f * (float)sin(2 * M_PI * 3000.0 * i / 48000.0) +
                          0.01f * noise
                    : noise;
    }
    for (int r = 0; r < 4; r++) {
      for (int b = 0; b < 5; b++) {
        udsp_spectrum_plan(&plan, UDSP_FFT_N, RATES[r], BARS[b]);
        udsp_spectrum_bars(sig, RATES[r], re, im, want, BARS[b]);
        udsp_spectrum_bars_planned(&plan, sig, scratch, got);
        float worst = 0.0f;
        for (int i = 0; i < BARS[b]; i++) {
          float d = fabsf(want[i] - got[i]);
          worst = d > worst ? d : worst;
        }
        UDSP_TEST_ASSERT(
            worst < 0.01f,
            "planned bars = udsp_spectrum_bars(): %s, %d Hz, %d bars "
            "(worst %.4f dB)",
            kind ? "tone" : "noise", RATES[r], BARS[b], worst);
      }
    }
  }
}

static void test_plan_edges(void) {
  static UDSP_SPECTRUM_PLAN plan;
  static float sig[UDSP_FFT_N], scratch[UDSP_SPECTRUM_SCRATCH];
  float bars[4];
  udsp_spectrum_plan(&plan, UDSP_FFT_N, 44100, 0);
  UDSP_TEST_ASSERT(plan.n_bars == 1, "a plan for 0 bars makes 1");
  udsp_spectrum_plan(&plan, UDSP_FFT_N, 44100, 1 << 20);
  UDSP_TEST_ASSERT(plan.n_bars == UDSP_SPECTRUM_MAX_BARS,
                   "a plan's bars are capped");
  for (int i = 0; i < UDSP_FFT_N; i++) {
    sig[i] = 1.0f;
  }
  udsp_spectrum_plan(&plan, UDSP_FFT_N, 30, 4); /* Nyquist below MIN_HZ */
  udsp_spectrum_bars_planned(&plan, sig, scratch, bars);
  UDSP_TEST_ASSERT(bars[0] == UDSP_DBFS_FLOOR && bars[3] == UDSP_DBFS_FLOOR,
                   "a rate with nothing above MIN_HZ reads the floor");
  /* power: a full-scale sine at an exact bin reads about -6 dB (the Hann
   * window's coherent gain, the documented simplification) in that bin */
  for (int i = 0; i < UDSP_FFT_N; i++) {
    sig[i] = (float)sin(2 * M_PI * 64.0 * i / UDSP_FFT_N);
  }
  udsp_spectrum_plan(&plan, UDSP_FFT_N, 48000, 28);
  const float *p = udsp_spectrum_power(&plan, sig, scratch);
  float db = 10.0f * log10f(p[64]);
  UDSP_TEST_ASSERT(db > -6.5f && db < -5.5f,
                   "a full-scale sine's bin reads ~-6 dB (%.2f)", db);

  /* the FFT size is the plan's */
  udsp_spectrum_plan(&plan, 1000, 48000, 28);
  UDSP_TEST_ASSERT(plan.n == UDSP_FFT_N,
                   "a size that isn't a power of two gets UDSP_FFT_N");
  udsp_spectrum_plan(&plan, 1 << 20, 48000, 28);
  UDSP_TEST_ASSERT(plan.n == UDSP_FFT_N,
                   "a size over the maximum gets UDSP_FFT_N");
}

/* short and long plans agree on where a tone is; each one's bin reads
 * the same level (the window's gain, normalised by its own size) */
static void test_plan_sizes(void) {
  static UDSP_SPECTRUM_PLAN small, big;
  static float sig[UDSP_FFT_N], scratch[UDSP_SPECTRUM_SCRATCH];
  float a[28], b[28];
  for (int i = 0; i < UDSP_FFT_N; i++) {
    sig[i] = (float)sin(2 * M_PI * 1000.0 * i / 48000.0);
  }
  int sizes[] = {16, 256, 2048, 4096};
  udsp_spectrum_plan(&big, UDSP_FFT_N, 48000, 28);
  udsp_spectrum_bars_planned(&big, sig, scratch, b);
  int want = 0;
  for (int i = 1; i < 28; i++) {
    want = b[i] > b[want] ? i : want;
  }
  for (int s = 0; s < 4; s++) {
    udsp_spectrum_plan(&small, sizes[s], 48000, 28);
    udsp_spectrum_bars_planned(&small, sig, scratch, a);
    int got = 0;
    for (int i = 1; i < 28; i++) {
      got = a[i] > a[got] ? i : got;
    }
    UDSP_TEST_ASSERT(
        sizes[s] < 256 || (got == want && fabsf(a[got] - b[want]) < 1.5f),
        "n=%d: a 1 kHz tone's loudest bar is the 8192-point plan's "
        "(bar %d, %.1f vs %.1f dB)",
        sizes[s], got, a[got], b[want]);
  }
}

int main(void) {
  test_planned_matches_bars();
  test_plan_edges();
  test_plan_sizes();
  test_silence_is_floor_every_bar();
  test_tone_peaks_in_expected_bar();
  test_bin_count_is_flexible();
  return udsp_test_summary();
}
