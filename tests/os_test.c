#include "udsp_test.h"

#include "fx/os.h"
#include "fx/sos.h"

#include <math.h>

#define R 48000
#define N 9600 /* base frames */
#define PI 3.14159265358979323846

static float in[2 * N], hi[2 * 4 * N], back[2 * N];

/* a sine of amp at hz in both channels, n frames */
static void sine(float *x, size_t n, double hz, double rate, double amp) {
  for (size_t i = 0; i < n; i++) {
    x[2 * i] = x[2 * i + 1] = (float)(amp * sin(2 * PI * hz * i / rate));
  }
}

/* n floats equal one by one */
static int same(const float *a, const float *b, size_t n) {
  int ok = 1;
  for (size_t i = 0; i < n; i++) {
    ok &= a[i] == b[i];
  }
  return ok;
}

/* the left's amplitude at hz over frames [from, to), and its phase */
static double amp_at(const float *x, size_t from, size_t to, double hz,
                     double rate, double *phase) {
  double s = 0, c = 0;
  for (size_t i = from; i < to; i++) {
    s += x[2 * i] * sin(2 * PI * hz * i / rate);
    c += x[2 * i] * cos(2 * PI * hz * i / rate);
  }
  if (phase) {
    *phase = atan2(c, s);
  }
  return 2 * sqrt(s * s + c * c) / (double)(to - from);
}

static void test_dc(void) {
  UDSP_OS o;
  udsp_os_reset(&o);
  for (size_t i = 0; i < 2 * N; i++) {
    in[i] = 0.5f;
  }
  udsp_os_up(&o, in, hi, N, 4);
  udsp_os_down(&o, hi, back, N, 4);
  UDSP_TEST_ASSERT(fabsf(hi[2 * 4 * N - 1] - 0.5f) < 1e-6f &&
                       fabsf(hi[2 * 4 * N - 4] - 0.5f) < 1e-6f &&
                       fabsf(back[2 * N - 2] - 0.5f) < 1e-6f,
                   "DC through 4x up and down: unchanged (%g, %g)",
                   (double)hi[2 * 4 * N - 1], (double)back[2 * N - 2]);
}

static void test_round_trip(void) {
  UDSP_OS o;
  udsp_os_reset(&o);
  sine(in, N, 1000, R, 0.5);
  udsp_os_up(&o, in, hi, N, 4);
  udsp_os_down(&o, hi, back, N, 4);
  double p_in, p_out;
  const double a_in = amp_at(in, 4800, N, 1000, R, &p_in);
  const double a_out = amp_at(back, 4800, N, 1000, R, &p_out);
  const double late = (p_in - p_out) / (2 * PI * 1000) * R;
  UDSP_TEST_ASSERT(fabs(20 * log10(a_out / a_in)) < 1e-4 && late > 4.5 &&
                       late < 4.75,
                   "1 kHz up and down at 4x: its level (%.6f dB), 4.6 "
                   "frames late (%.3f)",
                   20 * log10(a_out / a_in), late);
  udsp_os_reset(&o);
  sine(in, N, 20000, R, 0.5);
  udsp_os_up(&o, in, hi, N, 4);
  udsp_os_down(&o, hi, back, N, 4);
  const double a20 = amp_at(back, 4800, N, 20000, R, NULL);
  UDSP_TEST_ASSERT(fabs(20 * log10(a20 / 0.5)) < 1e-3, "20 kHz too (%.6f dB)",
                   20 * log10(a20 / 0.5));
}

static void test_images_and_aliases(void) {
  UDSP_OS o;
  udsp_os_reset(&o);
  sine(in, N, 15000, R, 1.0);
  udsp_os_up(&o, in, hi, N, 2);
  const double img = amp_at(hi, 4800, 2 * N, 33000, 2 * R, NULL);
  UDSP_TEST_ASSERT(20 * log10(img) < -88,
                   "up 2x: 15 kHz's image at 33 kHz "
                   "under -88 dB (%.1f)",
                   20 * log10(img));
  udsp_os_reset(&o);
  sine(hi, 2 * N, 30000, 2 * R, 1.0);
  udsp_os_down(&o, hi, back, N, 2);
  const double al = amp_at(back, 4800, N, 18000, R, NULL);
  UDSP_TEST_ASSERT(20 * log10(al) < -88,
                   "down 2x: 30 kHz folding to 18 kHz "
                   "under -88 dB (%.1f)",
                   20 * log10(al));
  udsp_os_reset(&o);
  sine(hi, 4 * N, 80000, 4 * R, 1.0);
  udsp_os_down(&o, hi, back, N, 4);
  const double al4 = amp_at(back, 4800, N, 16000, R, NULL);
  UDSP_TEST_ASSERT(20 * log10(al4) < -78,
                   "down 4x: 80 kHz folding to 16 kHz "
                   "under -78 dB (%.1f)",
                   20 * log10(al4));
}

static void test_blocks(void) {
  static float hi2[2 * 4 * N], back2[2 * N];
  UDSP_OS a, b;
  udsp_os_reset(&a);
  udsp_os_reset(&b);
  for (size_t i = 0; i < 2 * N; i++) {
    in[i] = (float)sin(i * 0.37) * (i & 1 ? 0.3f : 0.9f);
  }
  udsp_os_up(&a, in, hi, N, 4);
  udsp_os_down(&a, hi, back, N, 4);
  static const size_t step[] = {1, 7, 333};
  for (size_t i = 0, k = 0; i < N; k++) {
    size_t n = step[k % 3];
    n = n > N - i ? N - i : n;
    udsp_os_up(&b, in + 2 * i, hi2 + 8 * i, n, 4);
    udsp_os_down(&b, hi2 + 8 * i, back2 + 2 * i, n, 4);
    i += n;
  }
  UDSP_TEST_ASSERT(same(hi, hi2, 2 * 4 * N) && same(back, back2, 2 * N),
                   "blocks of 1, 7 and 333 frames: the same bits");
}

static void test_sos(void) {
  static const UDSP_SOS_SECTION half = {1, 0, 0, -0.5f, 0};
  UDSP_SOS f;
  udsp_sos_init(&f, &half, 1);
  float x[2 * 16] = {0};
  x[0] = 1;
  x[1] = 2;
  udsp_sos_stereo(&f, x, 16);
  int ok = 1;
  for (int i = 0; i < 16; i++) {
    ok &= x[2 * i] == ldexpf(1, -i) && x[2 * i + 1] == ldexpf(2, -i);
  }
  UDSP_TEST_ASSERT(ok, "a section y = x + y[-1] / 2: halves, each channel");
  static const UDSP_SOS_SECTION two[2] = {{0.2f, 0.3f, 0.1f, -0.4f, 0.2f},
                                          {1, -1.5f, 0.6f, 0.1f, -0.3f}};
  UDSP_SOS g, g0, g1;
  udsp_sos_init(&g, two, 2);
  udsp_sos_init(&g0, two, 1);
  udsp_sos_init(&g1, two + 1, 1);
  float y[2 * 64], z[2 * 64];
  for (int i = 0; i < 2 * 64; i++) {
    y[i] = z[i] = (float)sin(i * 1.3);
  }
  udsp_sos_stereo(&g, y, 20);
  udsp_sos_stereo(&g, y + 40, 44);
  udsp_sos_stereo(&g0, z, 64);
  udsp_sos_stereo(&g1, z, 64);
  UDSP_TEST_ASSERT(same(y, z, 2 * 64),
                   "two sections: each in turn, in any blocks");
}

int main(void) {
  test_dc();
  test_round_trip();
  test_images_and_aliases();
  test_blocks();
  test_sos();
  return udsp_test_summary();
}
