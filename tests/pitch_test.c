#include "udsp_test.h"

#include <math.h>

#include "osc.h"
#include "pitch.h"

static void test_octaves_are_exact_shifts(void) {
  UDSP_TEST_ASSERT(udsp_pitch_ratio_q32(0) == 1ull << 32, "unison is 1.0");
  UDSP_TEST_ASSERT(udsp_pitch_ratio_q32(12) == 1ull << 33,
                   "an octave up is 2.0");
  UDSP_TEST_ASSERT(udsp_pitch_ratio_q32(-12) == 1ull << 31,
                   "an octave down is 0.5");
  UDSP_TEST_ASSERT(udsp_pitch_ratio_q32(-24) == 1ull << 30, "two down is 0.25");
  UDSP_TEST_ASSERT(udsp_pitch_ratio_q32(19) == udsp_pitch_ratio_q32(7) << 1,
                   "a twelfth is a fifth an octave up");
}

static void test_downward_uses_floor_division(void) {
  /* -1 is octave -1 / pitch class 11, not octave 0 / class -1 */
  UDSP_TEST_ASSERT(udsp_pitch_ratio_q32(-1) == (1979448ull << 12) >> 1,
                   "a semitone down is B of the octave below");
  UDSP_TEST_ASSERT(udsp_pitch_ratio_q32(-1) < udsp_pitch_ratio_q32(0),
                   "and lower than unison");
}

static void test_note_frequencies(void) {
  UDSP_TEST_ASSERT(udsp_note_millihz(69) == 440000, "A4 is 440 Hz");
  UDSP_TEST_ASSERT(udsp_note_millihz(60) == 261625, "middle C is 261.625 Hz");
  UDSP_TEST_ASSERT(udsp_note_millihz(12) == 16351, "the lowest note");
  UDSP_TEST_ASSERT(udsp_note_millihz(119) == 7902133, "the highest note");
  UDSP_TEST_ASSERT(udsp_note_millihz(11) == 0 && udsp_note_millihz(120) == 0,
                   "outside the range is 0");
}

static void test_exp2(void) {
  UDSP_TEST_ASSERT(udsp_exp2(0.0) == 1.0 && udsp_exp2(3.0) == 8.0 &&
                       udsp_exp2(-1.0) == 0.5 && udsp_exp2(-10.0) == 1.0 / 1024,
                   "whole powers are exact");
  double worst = 0.0;
  for (int i = -4000; i <= 4000; i++) {
    double x = i * 0.00731;
    double e = fabs(udsp_exp2(x) / exp2(x) - 1.0);
    worst = e > worst ? e : worst;
  }
  UDSP_TEST_ASSERT(worst < 1e-13, "within 1e-13 of libm over +-29 (%.2g)",
                   worst);
  UDSP_TEST_ASSERT(udsp_exp2(100.0) == udsp_exp2(60.0), "clamped at 60");
}

static void test_note_inc(void) {
  const uint32_t a4 = udsp_osc_inc_hz(440.0, 48000);
  UDSP_TEST_ASSERT(udsp_note_inc(69.0, 48000) == a4, "69.0 is 440 Hz");
  uint32_t a5 = udsp_note_inc(81.0, 48000);
  UDSP_TEST_ASSERT(a5 >= 2 * a4 - 1 && a5 <= 2 * a4 + 1,
                   "an octave up doubles it");
  double cents = 1200.0 * log2((double)udsp_note_inc(69.5, 48000) / a4);
  UDSP_TEST_ASSERT(fabs(cents - 50.0) < 0.01, "69.5 is a quarter tone up");
  UDSP_TEST_ASSERT(udsp_note_inc(-30.0, 48000) > 0,
                   "below MIDI 0 still sounds");
  UDSP_TEST_ASSERT(udsp_note_inc(200.0, 48000) ==
                       udsp_osc_inc_hz(23952.0, 48000),
                   "capped just under half the rate");
}

int main(void) {
  test_exp2();
  test_note_inc();
  test_octaves_are_exact_shifts();
  test_downward_uses_floor_division();
  test_note_frequencies();
  return udsp_test_summary();
}
