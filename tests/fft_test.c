#include "udsp_test.h"

#include "fft.h"

#include <math.h>

static void test_hann_window(void) {
  float in[8], out[8];
  for (int i = 0; i < 8; i++) {
    in[i] = 1.0f;
  }
  udsp_window_hann(in, out, 8);
  UDSP_TEST_ASSERT(out[0] < 0.01f, "Hann window is ~0 at the first sample");
  UDSP_TEST_ASSERT(out[7] < 0.01f, "Hann window is ~0 at the last sample");
  UDSP_TEST_ASSERT(out[3] > 0.9f || out[4] > 0.9f, "Hann window peaks near the center");
}

static void test_hann_window_in_place(void) {
  float buf[4] = {1.0f, 1.0f, 1.0f, 1.0f};
  udsp_window_hann(buf, buf, 4);
  UDSP_TEST_ASSERT(buf[0] < 0.01f, "Hann window works in-place (in == out)");
}

static void test_fft_dc_signal(void) {
  size_t n = 8;
  float re[8], im[8];
  for (size_t i = 0; i < n; i++) {
    re[i] = 1.0f;
    im[i] = 0.0f;
  }
  udsp_fft_forward(re, im, n);
  float mag[8];
  udsp_fft_magnitude(re, im, n, mag);
  UDSP_TEST_ASSERT(mag[0] > 7.9f, "a constant (DC) signal's entire energy lands in bin 0");
  for (size_t i = 1; i < n; i++) {
    UDSP_TEST_ASSERT(mag[i] < 0.01f, "every non-DC bin is ~0 for a pure DC signal (bin %zu)", i);
  }
}

static void test_fft_nyquist_signal(void) {
  /* alternating +1/-1 is the highest frequency an 8-sample-per-cycle
   * discrete signal can represent -- all its energy belongs at bin n/2. */
  size_t n = 8;
  float re[8], im[8];
  for (size_t i = 0; i < n; i++) {
    re[i] = (i % 2 == 0) ? 1.0f : -1.0f;
    im[i] = 0.0f;
  }
  udsp_fft_forward(re, im, n);
  float mag[8];
  udsp_fft_magnitude(re, im, n, mag);
  UDSP_TEST_ASSERT(mag[4] > 7.9f, "an alternating +-1 signal's energy lands at the Nyquist bin (n/2)");
}

int main(void) {
  test_hann_window();
  test_hann_window_in_place();
  test_fft_dc_signal();
  test_fft_nyquist_signal();
  return udsp_test_summary();
}
