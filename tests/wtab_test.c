#include "udsp_test.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "synth/osc.h"
#include "synth/wtab.h"

#define RATE 48000
#define PI 3.14159265358979323846

/* the sine-phase amplitude of harmonic h over a period of n points */
static double harmonic(const float *x, size_t n, size_t h) {
  double s = 0.0, c = 0.0;
  for (size_t i = 0; i < n; i++) {
    s += x[i] * sin(2.0 * PI * (double)(h * i) / (double)n);
    c += x[i] * cos(2.0 * PI * (double)(h * i) / (double)n);
  }
  return sqrt(s * s + c * c) * 2.0 / (double)n;
}

static void test_recipes(void) {
  float a[UDSP_WTAB_HARMONICS];
  udsp_wtab_recipe(UDSP_WTAB_SAW, 3, a);
  UDSP_TEST_ASSERT(a[0] == 1.0f && a[1] == 0.5f && a[2] == 1.0f / 3.0f &&
                       a[3] == 0.0f && a[63] == 0.0f,
                   "Saw3: 1, 1/2, 1/3, nothing past");
  udsp_wtab_recipe(UDSP_WTAB_SQUARE, 64, a);
  UDSP_TEST_ASSERT(a[0] == 1.0f && a[1] == 0.0f && a[2] == 1.0f / 3.0f &&
                       a[62] == 1.0f / 63.0f,
                   "a square's odd harmonics");
  udsp_wtab_recipe(UDSP_WTAB_TRIANGLE, 64, a);
  UDSP_TEST_ASSERT(a[0] == 1.0f && a[2] == -1.0f / 9.0f && a[4] == 1.0f / 25.0f,
                   "a triangle's in sine phase: 1/h^2, the signs alternating");
}

static void test_levels(void) {
  UDSP_WTAB *w = malloc(sizeof *w);
  float a[UDSP_WTAB_HARMONICS];
  udsp_wtab_recipe(UDSP_WTAB_SAW, 64, a);
  udsp_wtab_build(w, a, 64);
  float peak = 0.0f;
  for (size_t i = 0; i < UDSP_WTAB_SIZE; i++) {
    peak = fmaxf(peak, fabsf(w->t[0][i]));
  }
  UDSP_TEST_ASSERT(peak == 1.0f, "the full table's peak at 1");
  int kept = 1, gone = 1, same = 1;
  const double h1 = harmonic(w->t[0], UDSP_WTAB_SIZE, 1);
  for (int l = 0; l < UDSP_WTAB_LEVELS; l++) {
    const size_t top = (size_t)UDSP_WTAB_HARMONICS >> l;
    same &= fabs(harmonic(w->t[l], UDSP_WTAB_SIZE, 1) - h1) < 1e-6;
    kept &= fabs(harmonic(w->t[l], UDSP_WTAB_SIZE, top) - h1 / top) < 1e-6;
    for (size_t h = top + 1; h <= 2 * top && h <= 80; h++) {
      gone &= harmonic(w->t[l], UDSP_WTAB_SIZE, h) < 1e-6;
    }
    same &= w->t[l][UDSP_WTAB_SIZE] == w->t[l][0];
  }
  UDSP_TEST_ASSERT(same, "the fundamental as loud at every level; each "
                         "period's first point again at its end");
  UDSP_TEST_ASSERT(kept && gone, "level l: harmonics to 64 >> l at 1/h, "
                                 "none past");
  udsp_wtab_recipe(UDSP_WTAB_TRIANGLE, 64, a);
  udsp_wtab_build(w, a, 64);
  UDSP_TEST_ASSERT(w->t[0][UDSP_WTAB_SIZE / 4] == 1.0f &&
                       fabsf(w->t[0][3 * UDSP_WTAB_SIZE / 4] + 1.0f) < 1e-6f &&
                       fabsf(w->t[0][UDSP_WTAB_SIZE / 8] - 0.5f) < 0.01f,
                   "a triangle: 1 at a quarter, -1 at three, 0.5 between");
  memset(a, 0, sizeof a);
  udsp_wtab_build(w, a, 64);
  UDSP_TEST_ASSERT(w->t[0][100] == 0.0f && w->t[6][100] == 0.0f,
                   "no harmonics: silence");
  free(w);
}

static void test_level_of(void) {
  UDSP_TEST_ASSERT(udsp_wtab_level(udsp_osc_inc_hz(100.0, RATE)) == 0 &&
                       udsp_wtab_level(udsp_osc_inc_hz(440.0, RATE)) == 1 &&
                       udsp_wtab_level(udsp_osc_inc_hz(750.0, RATE)) == 2 &&
                       udsp_wtab_level(udsp_osc_inc_hz(20000.0, RATE)) == 6 &&
                       udsp_wtab_level(udsp_osc_inc_hz(30000.0, RATE)) == 6,
                   "100 Hz all 64, 440 32 (54 fit), 750 16 (Nyquist is "
                   "32's), 20 kHz the fundamental, past Nyquist too");
}

static void test_fmop(void) {
  enum { N = 512 };
  float amp[N], pm[N], got[N], want[N];
  for (size_t i = 0; i < N; i++) {
    amp[i] = 0.5f + 0.5f * (float)i / N;
    pm[i] = 0.3f * udsp_osc_sine((float)i / 37.0f - (float)(i / 37));
  }
  UDSP_OSC s;
  udsp_osc_init(&s, UDSP_OSC_SINE, 1);
  s.inc = udsp_osc_inc_hz(440.0, RATE);
  udsp_osc_render_pm(&s, pm, want, N);
  UDSP_FMOP o = {NULL, 0, udsp_osc_inc_hz(440.0, RATE), 0.0f, 0.0f, 0.0f};
  udsp_fmop_render(&o, pm, amp, got, N);
  int eq = 1;
  for (size_t i = 0; i < N; i++) {
    eq &= got[i] == want[i] * amp[i];
  }
  UDSP_TEST_ASSERT(eq && o.phase == s.phase,
                   "a sine without feedback: the oscillator's PM times amp");

  /* a saw at 750 Hz: a period 64 frames, each on a table point, level 2
   * (16 harmonics); nothing past them */
  UDSP_WTAB *w = malloc(sizeof *w);
  float a[UDSP_WTAB_HARMONICS];
  udsp_wtab_recipe(UDSP_WTAB_SAW, 64, a);
  udsp_wtab_build(w, a, 64);
  for (size_t i = 0; i < N; i++) {
    amp[i] = 1.0f;
  }
  o = (UDSP_FMOP){w, 0, 1u << 26, 0.0f, 0.0f, 0.0f};
  udsp_fmop_render(&o, NULL, amp, got, N);
  double top = 0.0;
  for (size_t h = 17; h <= 32; h++) {
    top = fmax(top, harmonic(got, 64, h));
  }
  UDSP_TEST_ASSERT(top < 1e-6 &&
                       fabs(harmonic(got, 64, 16) -
                            harmonic(w->t[0], UDSP_WTAB_SIZE, 16)) < 1e-6,
                   "a saw at 750 Hz: its 16th harmonic as the full table's, "
                   "none past it (%g)",
                   top);

  /* feedback and a table, cut into odd blocks: the same bits */
  for (int k = 0; k < 2; k++) {
    UDSP_FMOP x = {k ? w : NULL, 7, udsp_osc_inc_hz(311.1, RATE), 0.4f, 0, 0};
    UDSP_FMOP y = x;
    udsp_fmop_render(&x, pm, amp, want, N);
    for (size_t at = 0, len = 1; at < N; at += len, len = len * 2 + 1) {
      const size_t m = at + len > N ? N - at : len;
      udsp_fmop_render(&y, pm + at, amp + at, got + at, m);
    }
    int same = memcmp(got, want, sizeof got) == 0, bounded = 1, moved = 0;
    for (size_t i = 0; i < N; i++) {
      bounded &= fabsf(want[i]) <= 1.0f;
    }
    UDSP_FMOP z = {k ? w : NULL, 7, udsp_osc_inc_hz(311.1, RATE), 0.0f, 0, 0};
    udsp_fmop_render(&z, pm, amp, got, N);
    for (size_t i = 0; i < N; i++) {
      moved |= got[i] != want[i];
    }
    UDSP_TEST_ASSERT(same && bounded && moved,
                     "%s with feedback: the same bits in any blocks, within "
                     "its level, not what it is without",
                     k ? "a table" : "the sine");
  }
  free(w);
}

int main(void) {
  test_recipes();
  test_levels();
  test_level_of();
  test_fmop();
  return udsp_test_summary();
}
