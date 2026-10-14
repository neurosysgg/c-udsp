#include "udsp_test.h"

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

int main(void) {
  test_octaves_are_exact_shifts();
  test_downward_uses_floor_division();
  test_note_frequencies();
  return udsp_test_summary();
}
