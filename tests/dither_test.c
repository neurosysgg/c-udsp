#include "udsp_test.h"

#include <math.h>
#include <string.h>

#include "dither.h"

#define N 200000

static float buf[N];

/* the noise alone: within one step either side, centred, triangular
 * (the middle half holds 3/4 of it), the same for the same seed */
static void test_noise(void) {
  UDSP_DITHER d;
  udsp_dither_init(&d, 12345);
  for (int i = 0; i < N; i++) {
    buf[i] = 0.0f;
  }
  udsp_dither_tpdf(&d, buf, N, 16);
  const float step = 1.0f / 32768.0f;
  double sum = 0.0;
  int out = 0, middle = 0;
  for (int i = 0; i < N; i++) {
    sum += buf[i];
    out += fabsf(buf[i]) > step;
    middle += fabsf(buf[i]) < step * 0.5f;
  }
  const double share = (double)middle / N;
  UDSP_TEST_ASSERT(out == 0 && fabs(sum / N) < step * 0.01 && share > 0.74 &&
                       share < 0.76,
                   "±1 step, centred, triangular (%.4f in the middle half)",
                   share);
  UDSP_DITHER e;
  udsp_dither_init(&e, 12345);
  float again[64] = {0};
  udsp_dither_tpdf(&e, again, 64, 16);
  UDSP_TEST_ASSERT(memcmp(again, buf, sizeof again) == 0,
                   "the same seed, the same noise");
  UDSP_DITHER z;
  udsp_dither_init(&z, 0);
  UDSP_TEST_ASSERT(z.state == 1, "seed 0 taken as 1");
}

/* a level a quarter step up, dithered and quantized: its average kept
 * (without, it rounds to 0); 24 bits' step 256 times finer */
static void test_linear(void) {
  UDSP_DITHER d;
  udsp_dither_init(&d, 7);
  const float step = 1.0f / 32768.0f;
  for (int i = 0; i < N; i++) {
    buf[i] = 0.25f * step;
  }
  udsp_dither_tpdf(&d, buf, N, 16);
  double sum = 0.0;
  for (int i = 0; i < N; i++) {
    sum += roundf(buf[i] / step);
  }
  UDSP_TEST_ASSERT(fabs(sum / N - 0.25) < 0.01,
                   "a quarter step kept on average (%.4f)", sum / N);
  for (int i = 0; i < 1000; i++) {
    buf[i] = 0.0f;
  }
  udsp_dither_tpdf(&d, buf, 1000, 24);
  float most = 0.0f;
  for (int i = 0; i < 1000; i++) {
    most = fabsf(buf[i]) > most ? fabsf(buf[i]) : most;
  }
  UDSP_TEST_ASSERT(most <= step / 256.0f && most > step / 1024.0f,
                   "24 bits: within its own step");
}

int main(void) {
  test_noise();
  test_linear();
  return udsp_test_summary();
}
