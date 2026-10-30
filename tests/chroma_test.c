#include "udsp_test.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "analysis/chroma.h"
#include "analysis/spectrum.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define N 8192
#define RATE 12000 /* a low band's: a semitone a few bins at 55 Hz */

static UDSP_SPECTRUM_PLAN plan;
static float x[N], scratch[UDSP_SPECTRUM_SCRATCH];

static double note_hz(int midi) { return 440.0 * pow(2.0, (midi - 69) / 12.0); }

/* notes (MIDI, n of them), harmonics of each at 1/h */
static size_t peaks_of(const int *notes, int n, int harmonics, UDSP_PEAK *out,
                       size_t cap) {
  memset(x, 0, sizeof x);
  for (int k = 0; k < n; k++) {
    for (int h = 1; h <= harmonics; h++) {
      for (int i = 0; i < N; i++) {
        x[i] +=
            (float)(0.2 / h * sin(2 * M_PI * h * note_hz(notes[k]) * i / RATE));
      }
    }
  }
  udsp_spectrum_plan(&plan, N, RATE, 16);
  const float *p = udsp_spectrum_power(&plan, x, scratch);
  return udsp_chroma_peaks(p, N, RATE, 25.0f, 5000.0f, 1e-5f, out, cap);
}

static void test_a_peak_to_the_cent(void) {
  UDSP_PEAK p[8];
  const int a2[] = {45};
  const size_t n = peaks_of(a2, 1, 1, p, 8);
  UDSP_TEST_ASSERT(n == 1, "n == 1");
  float cents;
  UDSP_TEST_ASSERT(udsp_chroma_note(p[0].hz, &cents) == 45,
                   "udsp_chroma_note(p[0].hz, &cents) == 45");
  UDSP_TEST_ASSERT(fabsf(cents) < 5.0f, "fabsf(cents) < 5.0f");
  UDSP_TEST_ASSERT(udsp_chroma_note(-1.0f, NULL) == -1,
                   "udsp_chroma_note(-1.0f, NULL) == -1");
}

static void test_the_strongest_kept_in_order(void) {
  UDSP_PEAK p[3];
  const int notes[] = {45};
  const size_t n = peaks_of(notes, 1, 6, p, 3);
  UDSP_TEST_ASSERT(n == 3, "n == 3");
  UDSP_TEST_ASSERT(p[0].hz < p[1].hz && p[1].hz < p[2].hz,
                   "p[0].hz < p[1].hz && p[1].hz < p[2].hz");
  UDSP_TEST_ASSERT(udsp_chroma_note(p[0].hz, NULL) == 45,
                   "udsp_chroma_note(p[0].hz, NULL) == 45");
}

static void test_a_root_under_its_harmonics(void) {
  UDSP_PEAK p[32];
  const int notes[] = {33}; /* A1 */
  const size_t n = peaks_of(notes, 1, 8, p, 32);
  UDSP_TEST_ASSERT(udsp_chroma_note(udsp_chroma_root(p, n), NULL) == 33,
                   "udsp_chroma_note(udsp_chroma_root(p, n), NULL) == 33");
  UDSP_TEST_ASSERT(udsp_chroma_root(NULL, 0) == 0.0f,
                   "udsp_chroma_root(NULL, 0) == 0.0f");
}

static void test_chords(void) {
  UDSP_PEAK p[32];
  float chroma[12] = {0};
  const int am[] = {57, 60, 64};
  size_t n = peaks_of(am, 3, 1, p, 32);
  udsp_chroma_fold(p, n, chroma);
  int root, chord;
  UDSP_TEST_ASSERT(udsp_chroma_chord(chroma, &root, &chord) > 0.9f,
                   "udsp_chroma_chord(chroma, &root, &chord) > 0.9f");
  UDSP_TEST_ASSERT(root == 9 && chord == UDSP_CHORD_MINOR,
                   "root == 9 && chord == UDSP_CHORD_MINOR");
  UDSP_TEST_ASSERT(!strcmp(udsp_chord_name(chord), "m"),
                   "!strcmp(udsp_chord_name(chord), 'm')");

  memset(chroma, 0, sizeof chroma);
  const int g7[] = {55, 59, 62, 65};
  n = peaks_of(g7, 4, 1, p, 32);
  udsp_chroma_fold(p, n, chroma);
  udsp_chroma_chord(chroma, &root, &chord);
  UDSP_TEST_ASSERT(root == 7 && chord == UDSP_CHORD_DOM7,
                   "root == 7 && chord == UDSP_CHORD_DOM7");

  const float none[12] = {0};
  UDSP_TEST_ASSERT(udsp_chroma_chord(none, &root, &chord) == 0.0f,
                   "udsp_chroma_chord(none, &root, &chord) == 0.0f");
}

static void test_keys(void) {
  /* C major's and A minor's notes, each scale's tonic chord tones heavier */
  float c_major[12] = {0}, a_minor[12] = {0};
  const int scale[] = {0, 2, 4, 5, 7, 9, 11};
  for (int i = 0; i < 7; i++) {
    c_major[scale[i]] = a_minor[scale[i]] = 1.0f;
  }
  c_major[0] += 2;
  c_major[4] += 1;
  c_major[7] += 1.5f;
  a_minor[9] += 2;
  a_minor[0] += 1;
  a_minor[4] += 1.5f;
  int root, minor;
  UDSP_TEST_ASSERT(udsp_chroma_key(c_major, &root, &minor) > 0.8f,
                   "udsp_chroma_key(c_major, &root, &minor) > 0.8f");
  UDSP_TEST_ASSERT(root == 0 && minor == 0, "root == 0 && minor == 0");
  UDSP_TEST_ASSERT(udsp_chroma_key(a_minor, &root, &minor) > 0.8f,
                   "udsp_chroma_key(a_minor, &root, &minor) > 0.8f");
  UDSP_TEST_ASSERT(root == 9 && minor == 1, "root == 9 && minor == 1");
  const float flat[12] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  UDSP_TEST_ASSERT(udsp_chroma_key(flat, &root, &minor) == 0.0f,
                   "udsp_chroma_key(flat, &root, &minor) == 0.0f");
}

int main(void) {
  test_a_peak_to_the_cent();
  test_the_strongest_kept_in_order();
  test_a_root_under_its_harmonics();
  test_chords();
  test_keys();
  return udsp_test_summary();
}
