#include "udsp_test.h"

#include "synth/grain.h"

#include <math.h>
#include <string.h>

#define N 4096

static float TABLE[UDSP_GRAIN_SHAPES][UDSP_GRAIN_TABLE + 1];
static float X[N], ST[N * 2];
static float A[N * 2], B[N * 2];

static void test_windows(void) {
  for (int s = 0; s < UDSP_GRAIN_SHAPES; s++) {
    udsp_grain_window(TABLE[s], (UDSP_GRAIN_SHAPE)s);
    float top = 0.0f;
    for (int i = 0; i <= UDSP_GRAIN_TABLE; i++) {
      top = TABLE[s][i] > top ? TABLE[s][i] : top;
      UDSP_TEST_ASSERT(TABLE[s][i] >= 0.0f && TABLE[s][i] <= 1.0f, "in 0..1");
    }
    UDSP_TEST_ASSERT(TABLE[s][0] == 0.0f && TABLE[s][UDSP_GRAIN_TABLE] == 0.0f,
                     "closed at both ends");
    UDSP_TEST_ASSERT(top > 0.8f, "reaches near 1");
  }
  UDSP_TEST_ASSERT(TABLE[UDSP_GRAIN_HANN][512] == 1.0f, "hann's middle");
  UDSP_TEST_ASSERT(fabsf(TABLE[UDSP_GRAIN_HANN][256] - 0.5f) < 1e-6f,
                   "hann's quarter");
  UDSP_TEST_ASSERT(TABLE[UDSP_GRAIN_TRIANGLE][256] == 0.5f, "a triangle");
  UDSP_TEST_ASSERT(TABLE[UDSP_GRAIN_TUKEY][400] == 1.0f, "tukey's flat");
  for (int i = 0; i <= UDSP_GRAIN_TABLE; i++) {
    UDSP_TEST_ASSERT(TABLE[UDSP_GRAIN_PERC][i] ==
                         TABLE[UDSP_GRAIN_SWELL][UDSP_GRAIN_TABLE - i],
                     "the swell is the perc backwards");
  }
  UDSP_TEST_ASSERT(TABLE[UDSP_GRAIN_PERC][64] >
                       4.0f * TABLE[UDSP_GRAIN_PERC][512],
                   "the perc falls");
}

/* a grain over a constant is its window, sample for sample */
static void test_window_shapes_a_constant(void) {
  for (int i = 0; i < N; i++) {
    X[i] = 1.0f;
  }
  UDSP_GRAIN g;
  udsp_grain_start(&g, (int64_t)100 << 32, 1ll << 32, 1024, 0.5f, 0.25f);
  memset(A, 0, sizeof A);
  uint32_t k = udsp_grain_render(&g, X, 1, N, TABLE[UDSP_GRAIN_HANN], A, 2000);
  UDSP_TEST_ASSERT(k == 1024 && g.left == 0, "plays its length and stops");
  UDSP_TEST_ASSERT(A[0] == 0.0f, "from the window's zero");
  for (int i = 0; i < 1024; i++) {
    const float w = TABLE[UDSP_GRAIN_HANN][i];
    UDSP_TEST_ASSERT(fabsf(A[2 * i] - w * 0.5f) < 1e-6f, "left: the window");
    UDSP_TEST_ASSERT(fabsf(A[2 * i + 1] - w * 0.25f) < 1e-6f,
                     "right: its gain");
  }
  UDSP_TEST_ASSERT(A[2 * 1024] == 0.0f, "nothing past it");
}

/* a ramp read at twice the speed, and backwards, is exact */
static void test_pitch_and_direction(void) {
  for (int i = 0; i < N; i++) {
    X[i] = (float)i;
  }
  const int64_t one = 1ll << 32;
  float flat[UDSP_GRAIN_TABLE + 1];
  udsp_grain_window(flat, UDSP_GRAIN_TUKEY);
  UDSP_GRAIN g;
  udsp_grain_start(&g, 100 * one, 2 * one, 1024, 1.0f, 1.0f);
  memset(A, 0, sizeof A);
  udsp_grain_render(&g, X, 1, N, flat, A, 1024);
  UDSP_TEST_ASSERT(A[2 * 512] == 1124.0f, "two frames a frame");
  udsp_grain_start(&g, 3000 * one, -one / 2, 1024, 1.0f, 1.0f);
  memset(A, 0, sizeof A);
  udsp_grain_render(&g, X, 1, N, flat, A, 1024);
  UDSP_TEST_ASSERT(A[2 * 512] == 2744.0f, "backwards at half speed");
  udsp_grain_start(&g, 100 * one + one / 2, one, 1024, 1.0f, 1.0f);
  memset(A, 0, sizeof A);
  udsp_grain_render(&g, X, 1, N, flat, A, 1024);
  UDSP_TEST_ASSERT(A[2 * 512] == 612.5f, "between frames on a ramp");
}

static void test_stereo_and_stride(void) {
  for (int i = 0; i < N; i++) {
    ST[2 * i] = 1.0f;
    ST[2 * i + 1] = -1.0f;
  }
  float flat[UDSP_GRAIN_TABLE + 1];
  udsp_grain_window(flat, UDSP_GRAIN_TUKEY);
  UDSP_GRAIN g;
  udsp_grain_start(&g, 10ll << 32, 1ll << 32, 1024, 1.0f, 1.0f);
  memset(A, 0, sizeof A);
  udsp_grain_render(&g, ST, 2, N, flat, A, 1024);
  UDSP_TEST_ASSERT(A[2 * 512] == 1.0f && A[2 * 512 + 1] == -1.0f,
                   "each side its own");
  /* adds to what's there */
  udsp_grain_start(&g, 10ll << 32, 1ll << 32, 1024, 1.0f, 1.0f);
  udsp_grain_render(&g, ST, 2, N, flat, A, 1024);
  UDSP_TEST_ASSERT(A[2 * 512] == 2.0f, "mixed in");
}

/* past either end it's silence, never a read outside */
static void test_ends(void) {
  float s[8] = {1, 1, 1, 1, 1, 1, 1, 1};
  float flat[UDSP_GRAIN_TABLE + 1];
  udsp_grain_window(flat, UDSP_GRAIN_TUKEY);
  UDSP_GRAIN g;
  udsp_grain_start(&g, -(40ll << 32), 1ll << 32, 1024, 1.0f, 1.0f);
  memset(A, 0, sizeof A);
  udsp_grain_render(&g, s, 1, 8, flat, A, 1024);
  UDSP_TEST_ASSERT(A[2 * 512] == 0.0f, "far past the end: silence");
  udsp_grain_start(&g, -(500ll << 32), 1ll << 32, 1024, 1.0f, 1.0f);
  memset(A, 0, sizeof A);
  udsp_grain_render(&g, s, 1, 8, flat, A, 1024);
  UDSP_TEST_ASSERT(A[2 * 503] == 1.0f, "the sample as it passes");
  UDSP_TEST_ASSERT(A[2 * 600] == 0.0f, "then nothing");
}

/* the same bits however the frames are cut */
static void test_any_block(void) {
  for (int i = 0; i < N; i++) {
    X[i] = sinf((float)i * 0.03f);
  }
  const int cuts[] = {1, 7, 32, 333};
  UDSP_GRAIN g;
  udsp_grain_start(&g, 77ll << 32, 0x12345678ll, 1000, 0.7f, 0.3f);
  memset(A, 0, sizeof A);
  udsp_grain_render(&g, X, 1, N, TABLE[UDSP_GRAIN_PERC], A, 1000);
  for (size_t c = 0; c < sizeof cuts / sizeof cuts[0]; c++) {
    udsp_grain_start(&g, 77ll << 32, 0x12345678ll, 1000, 0.7f, 0.3f);
    memset(B, 0, sizeof B);
    for (uint32_t at = 0; at < 1000;) {
      at += udsp_grain_render(&g, X, 1, N, TABLE[UDSP_GRAIN_PERC], B + 2 * at,
                              (uint32_t)cuts[c]);
    }
    int same = 1;
    for (int i = 0; i < 2 * N; i++) {
      same = same && A[i] == B[i];
    }
    UDSP_TEST_ASSERT(same, "the same bits at any cut");
  }
}

int main(void) {
  test_windows();
  test_window_shapes_a_constant();
  test_pitch_and_direction();
  test_stereo_and_stride();
  test_ends();
  test_any_block();
  return udsp_test_summary();
}
