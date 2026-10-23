#include "udsp_test.h"

#include "offline/sinc.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static float *g_table;

static int make(UDSP_SINC *s, uint32_t src, uint32_t dst) {
  size_t n = udsp_sinc_floats(src, dst, UDSP_SINC_ZEROS);
  free(g_table);
  g_table = malloc(n * sizeof *g_table);
  return g_table && udsp_sinc_init(s, src, dst, UDSP_SINC_ZEROS, g_table);
}

static float *sine(double hz, uint32_t rate, size_t n) {
  float *x = malloc(n * sizeof *x);
  for (size_t i = 0; i < n; i++) {
    x[i] = (float)(0.5 * sin(2.0 * M_PI * hz * (double)i / rate));
  }
  return x;
}

static void test_shape(void) {
  UDSP_TEST_ASSERT(udsp_sinc_floats(44100, 48000, 128) == 160u * 256,
                   "44.1 -> 48 kHz: 160 phases of 256 taps (%zu)",
                   udsp_sinc_floats(44100, 48000, 128));
  UDSP_TEST_ASSERT(udsp_sinc_floats(48000, 24000, 128) == 512,
                   "48 -> 24 kHz: one phase, the kernel twice as long");
  UDSP_TEST_ASSERT(udsp_sinc_floats(0, 48000, 128) == 0 &&
                       udsp_sinc_floats(48000, 2000000, 128) == 0 &&
                       udsp_sinc_floats(48000, 44100, 0) == 0,
                   "a rate of 0, past the limit, no zeros: refused");
  UDSP_SINC s;
  UDSP_TEST_ASSERT(make(&s, 44100, 48000) && s.up == 160 && s.down == 147 &&
                       !s.interp,
                   "the ratio reduced");
  UDSP_TEST_ASSERT(
      udsp_sinc_length(&s, 44100) == 48000 && udsp_sinc_length(&s, 1) == 2 &&
          udsp_sinc_length(&s, 147) == 160 && udsp_sinc_length(&s, 148) == 162,
      "the output's length: ceil(in * up / down)");
  UDSP_TEST_ASSERT(make(&s, 44100, 47999) && s.interp &&
                       s.rows == UDSP_SINC_PHASES_MAX + 1,
                   "47999 Hz: too many phases, interpolated rows");
}

static void test_dc(void) {
  static const uint32_t rates[][2] = {
      {44100, 48000}, {48000, 44100}, {96000, 48000}, {22050, 48000}};
  for (size_t r = 0; r < sizeof rates / sizeof *rates; r++) {
    UDSP_SINC s;
    make(&s, rates[r][0], rates[r][1]);
    float in[4096], out[512];
    for (size_t i = 0; i < 4096; i++) {
      in[i] = 0.25f;
    }
    /* far enough from both ends for every tap to land in the input */
    udsp_sinc_run(&s, in, 0, 4096, 1000, out, 512);
    double worst = 0.0;
    for (size_t i = 0; i < 512; i++) {
      double e = fabs(out[i] - 0.25);
      worst = e > worst ? e : worst;
    }
    UDSP_TEST_ASSERT(worst < 1e-6, "%u -> %u: a constant passes (off by %.2g)",
                     rates[r][0], rates[r][1], worst);
  }
}

static void test_same_rate(void) {
  UDSP_SINC s;
  make(&s, 48000, 48000);
  float in[64], out[64];
  for (int i = 0; i < 64; i++) {
    in[i] = (float)i - 20.5f;
  }
  udsp_sinc_run(&s, in, 0, 64, 0, out, 64);
  UDSP_TEST_ASSERT(!memcmp(in, out, sizeof in), "the same rate: a copy");
}

/* a tone resampled is the tone at the new rate, in the passband */
static void sine_matches(uint32_t src, uint32_t dst, double hz, double tol) {
  UDSP_SINC s;
  make(&s, src, dst);
  const size_t n = 8192;
  float *in = sine(hz, src, n);
  uint64_t out_n = udsp_sinc_length(&s, n);
  float *out = malloc(out_n * sizeof *out);
  udsp_sinc_run(&s, in, 0, n, 0, out, out_n);
  double worst = 0.0;
  for (uint64_t i = 1000; i + 1000 < out_n; i++) {
    double want = 0.5 * sin(2.0 * M_PI * hz * (double)i / dst);
    double e = fabs(out[i] - want);
    worst = e > worst ? e : worst;
  }
  UDSP_TEST_ASSERT(worst < tol, "%g Hz, %u -> %u: off by %.2g (< %.2g)", hz,
                   src, dst, worst, tol);
  free(in);
  free(out);
}

static void test_sine(void) {
  /* float's own rounding is ~3e-8 on 0.5; the kernel's ripple adds a
   * little near the band's edge */
  sine_matches(44100, 48000, 1000.0, 2e-6);
  sine_matches(44100, 48000, 15000.0, 2e-6);
  sine_matches(48000, 44100, 1000.0, 2e-6);
  sine_matches(96000, 44100, 10000.0, 2e-6);
  sine_matches(44100, 47999, 1000.0, 2e-6);
  sine_matches(44100, 47999, 15000.0, 1e-5);
}

/* past the lower rate's Nyquist: gone (120 dB) */
static void test_stopband(void) {
  UDSP_SINC s;
  make(&s, 48000, 44100);
  const size_t n = 8192;
  float *in = sine(23000.0, 48000, n);
  float out[4096];
  udsp_sinc_run(&s, in, 0, n, 1000, out, 4096);
  double peak = 0.0;
  for (size_t i = 0; i < 4096; i++) {
    peak = fabs(out[i]) > peak ? fabs(out[i]) : peak;
  }
  double db = 20.0 * log10(peak / 0.5 + 1e-30);
  UDSP_TEST_ASSERT(db < -110.0, "23 kHz into 44.1: %.1f dB", db);
  free(in);
}

/* outputs cut into blocks, each from a window of just its span, are
 * the whole run's bits -- the edges (silence outside the input) too */
static void cut(uint32_t src, uint32_t dst) {
  UDSP_SINC s;
  make(&s, src, dst);
  const size_t n = 3000;
  float *in = sine(3000.0, src, n);
  for (size_t i = 0; i < n; i += 7) {
    in[i] += 0.3f; /* not just a tone */
  }
  uint64_t out_n = udsp_sinc_length(&s, n);
  float *whole = malloc(out_n * sizeof *whole);
  float *parts = malloc(out_n * sizeof *parts);
  udsp_sinc_run(&s, in, 0, n, 0, whole, out_n);
  for (uint64_t o = 0; o < out_n;) {
    size_t len = 1 + (size_t)(o * 7 % 97);
    len = o + len > out_n ? (size_t)(out_n - o) : len;
    int64_t lo, hi;
    udsp_sinc_span(&s, o, len, &lo, &hi);
    lo = lo < 0 ? 0 : lo;
    hi = hi > (int64_t)n ? (int64_t)n : hi;
    udsp_sinc_run(&s, in + lo, lo, (size_t)(hi - lo), o, parts + o, len);
    o += len;
  }
  UDSP_TEST_ASSERT(!memcmp(whole, parts, out_n * sizeof *whole),
                   "%u -> %u: blocks from their spans are the whole run", src,
                   dst);
  free(in);
  free(whole);
  free(parts);
}

static void test_cut(void) {
  cut(44100, 48000);
  cut(48000, 44100);
  cut(44100, 47999);
}

/* an impulse's response is the row itself, scaled */
static void test_impulse(void) {
  UDSP_SINC s;
  make(&s, 44100, 48000);
  float in[1] = {1.0f}, out[400];
  udsp_sinc_run(&s, in, 500, 1, 0, out, 400);
  double sum = 0.0, peak = 0.0;
  for (size_t i = 0; i < 400; i++) {
    sum += out[i];
    peak = fabs(out[i]) > peak ? fabs(out[i]) : peak;
  }
  /* an impulse alone at frame 500 of the input: 160/147 outputs a frame,
   * each tap row summing to 1 -> the outputs sum to ~up / down */
  UDSP_TEST_ASSERT(sum == 0.0 && peak == 0.0,
                   "an impulse out of every output's reach: silence");
  udsp_sinc_run(&s, in, 200, 1, 0, out, 400);
  sum = 0.0;
  for (size_t i = 0; i < 400; i++) {
    sum += out[i];
  }
  UDSP_TEST_ASSERT(fabs(sum - 160.0 / 147.0) < 1e-4,
                   "an impulse's response sums to up / down (%.6f)", sum);
}

/* --- varispeed ----------------------------------------------------------- */

static float *g_vtable;
static UDSP_VSINC g_vs;

static uint64_t speed_q(double s) { return (uint64_t)(s * UDSP_VSINC_ONE); }

static void test_vs_identity(void) {
  g_vtable = malloc(udsp_vsinc_floats(UDSP_VSINC_ZEROS) * sizeof *g_vtable);
  udsp_vsinc_init(&g_vs, UDSP_VSINC_ZEROS, g_vtable);
  UDSP_TEST_ASSERT(udsp_vsinc_floats(32) == 2 * (32u * UDSP_VSINC_PHASES + 2) &&
                       udsp_vsinc_floats(0) == 0,
                   "the table's size");
  float in[2 * 300], out[2 * 200];
  for (size_t i = 0; i < 600; i++) {
    in[i] = (float)((i * 7919) % 1000) / 1000.0f - 0.5f;
  }
  udsp_vsinc_run(&g_vs, in, 2, 0, 300, 50, 0, UDSP_VSINC_ONE, out, 2, 200);
  UDSP_TEST_ASSERT(!memcmp(out, in + 100, sizeof out),
                   "speed 1 at a whole frame: the input itself, stereo");
}

static void test_vs_dc(void) {
  static float in[4000], out[500];
  for (size_t i = 0; i < 4000; i++) {
    in[i] = 0.5f;
  }
  static const double speeds[] = {0.25, 0.7, 1.0, 1.37, 3.1, 6.0};
  for (size_t k = 0; k < sizeof speeds / sizeof *speeds; k++) {
    udsp_vsinc_run(&g_vs, in, 1, 0, 4000, 300, 0x12345678u, speed_q(speeds[k]),
                   out, 1, 500);
    double worst = 0.0;
    for (size_t i = 0; i < 500; i++) {
      worst = fabs(out[i] - 0.5) > worst ? fabs(out[i] - 0.5) : worst;
    }
    UDSP_TEST_ASSERT(worst < 1e-6, "a constant through at speed %.2f (%.2g)",
                     speeds[k], worst);
  }
}

/* a sine at speed s against the sine it should become */
static double vs_sine_error(double hz, double s) {
  const uint32_t rate = 48000;
  const size_t n_in = 24000, n = 4000;
  float *in = sine(hz, rate, n_in);
  float *out = malloc(n * sizeof *out);
  const int64_t frame = 1000;
  udsp_vsinc_run(&g_vs, in, 1, 0, n_in, frame, 0, speed_q(s), out, 1, n);
  const double step = (double)speed_q(s) / UDSP_VSINC_ONE;
  double worst = 0.0;
  for (size_t i = 0; i < n; i++) {
    double ideal = 0.5 * sin(2.0 * M_PI * hz * (frame + i * step) / rate);
    worst = fabs(out[i] - ideal) > worst ? fabs(out[i] - ideal) : worst;
  }
  free(in);
  free(out);
  return worst;
}

static void test_vs_sine(void) {
  double e1 = vs_sine_error(1000.0, 1.5), e2 = vs_sine_error(5000.0, 0.7),
         e3 = vs_sine_error(440.0, 1.0001);
  UDSP_TEST_ASSERT(e1 < 1e-4 && e2 < 1e-4 && e3 < 1e-4,
                   "a sine read faster, slower, a hair off: the sine it "
                   "should be (%.2g, %.2g, %.2g; -80 dB)",
                   e1, e2, e3);
}

/* what would fold back when read faster is filtered out instead */
static void test_vs_alias(void) {
  const size_t n_in = 48000, n = 8000;
  float *in = sine(18000.0, 48000, n_in);
  float *out = malloc(n * sizeof *out);
  udsp_vsinc_run(&g_vs, in, 1, 0, n_in, 2000, 0, speed_q(2.0), out, 1, n);
  double sq = 0.0;
  for (size_t i = 0; i < n; i++) {
    sq += (double)out[i] * out[i];
  }
  const double db = 10.0 * log10(sq / n / 0.125 + 1e-30);
  UDSP_TEST_ASSERT(db < -90.0,
                   "18 kHz read twice as fast (36 kHz: past Nyquist): gone, "
                   "not aliased to 12 kHz (%.1f dB)",
                   db);
  free(in);
  free(out);
}

/* however the outputs are cut into blocks, the same bits; a window of
 * the span's frames gives them all */
static void test_vs_cut(void) {
  const size_t n_in = 20000, n = 3000;
  float *in = sine(3000.0, 48000, n_in);
  float *whole = malloc(n * sizeof *whole), *cut = malloc(n * sizeof *cut);
  const uint64_t step = speed_q(1.2345);
  udsp_vsinc_run(&g_vs, in, 1, 0, n_in, 500, 77, step, whole, 1, n);
  int64_t frame = 500;
  uint32_t frac = 77;
  int same_span = 1;
  for (size_t done = 0; done < n;) {
    size_t len = n - done < 37 ? n - done : 37;
    int64_t lo, hi;
    udsp_vsinc_span(&g_vs, frame, frac, step, len, &lo, &hi);
    /* only the span's frames given: still every output exactly */
    udsp_vsinc_run(&g_vs, in + lo, 1, lo, (size_t)(hi - lo), frame, frac, step,
                   cut + done, 1, len);
    const uint64_t at = (uint64_t)frac + len * step;
    frame += (int64_t)(at >> 32);
    frac = (uint32_t)at;
    done += len;
    same_span &= lo >= 0 && hi <= (int64_t)n_in;
  }
  UDSP_TEST_ASSERT(same_span && !memcmp(whole, cut, n * sizeof *whole),
                   "cut into blocks of 37, each given its span alone: the "
                   "same bits as at once");
  free(in);
  free(whole);
  free(cut);
}

/* udsp_vsinc_run as it was a tap at a time (before the weights went
 * into lanes, session 257): the bits the run is held to */
static void ref_vsinc(const UDSP_VSINC *k, const float *in, uint32_t channels,
                      int64_t in_first, size_t in_n, int64_t frame,
                      uint32_t frac, uint64_t step, float *out, size_t n) {
  float w[2 * UDSP_VSINC_ZEROS * UDSP_VSINC_SPEED_MAX + 8];
  double s = (double)step / (double)UDSP_VSINC_ONE;
  s = s < 1.0 ? 1.0 : s > UDSP_VSINC_SPEED_MAX ? UDSP_VSINC_SPEED_MAX : s;
  const double reach = k->zeros * s;
  const double scale = UDSP_VSINC_PHASES * k->zeros / reach;
  const uint32_t limit = k->zeros * UDSP_VSINC_PHASES;
  const int64_t half = (int64_t)reach + 1;
  uint32_t taps = (uint32_t)(2 * half + 1);
  taps = (taps + 3) & ~3u;
  for (size_t o = 0; o < n; o++) {
    const uint64_t at = (uint64_t)frac + (uint64_t)o * step;
    const int64_t f = frame + (int64_t)(at >> 32);
    const uint32_t fr = (uint32_t)at;
    float *y = out + o * channels;
    const int64_t rel = f - in_first;
    if (step == UDSP_VSINC_ONE && fr == 0) {
      for (uint32_t c = 0; c < channels; c++) {
        y[c] = rel >= 0 && (uint64_t)rel < in_n ? in[(size_t)rel * channels + c]
                                                : 0.0f;
      }
      continue;
    }
    const int64_t j0 = f - half + 1;
    const double p = (double)fr / (double)UDSP_VSINC_ONE;
    float sum4[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    for (uint32_t t = 0; t < taps; t += 4) {
      for (uint32_t j = 0; j < 4; j++) {
        double x = (double)(j0 + (int64_t)(t + j) - f) - p;
        x = (x < 0.0 ? -x : x) * scale;
        float v = 0.0f;
        if (x < (double)limit) {
          const uint32_t i = (uint32_t)x;
          const float mu = (float)(x - (double)i);
          v = k->table[2 * i] + mu * (k->table[2 * i + 2] - k->table[2 * i]);
        }
        w[t + j] = v;
        sum4[j] = sum4[j] + v;
      }
    }
    const float sum = (sum4[0] + sum4[1]) + (sum4[2] + sum4[3]);
    const float norm = sum != 0.0f ? 1.0f / sum : 0.0f;
    for (uint32_t c = 0; c < channels; c++) {
      float acc[4] = {0.0f, 0.0f, 0.0f, 0.0f};
      for (uint32_t t = 0; t < taps; t += 4) {
        for (uint32_t j = 0; j < 4; j++) {
          const int64_t r = j0 + (int64_t)(t + j) - in_first;
          const float xv = r >= 0 && (uint64_t)r < in_n
                               ? in[(size_t)r * channels + c]
                               : 0.0f;
          acc[j] = acc[j] + w[t + j] * xv;
        }
      }
      y[c] = ((acc[0] + acc[1]) + (acc[2] + acc[3])) * norm;
    }
  }
}

/* slower, about 1, faster, past the widest; mono, stereo and three
 * channels; windows the kernel runs off at both ends: the old bits */
static void test_vs_reference(void) {
  enum { IN = 3000, OUT = 600 };
  static float in[IN * 3], got[OUT * 3], want[OUT * 3];
  uint32_t rng = 5;
  for (int i = 0; i < IN * 3; i++) {
    rng ^= rng << 13;
    rng ^= rng >> 17;
    rng ^= rng << 5;
    in[i] = (float)(int32_t)rng * (1.0f / 2147483648.0f);
  }
  static const double SPEED[] = {0.37, 0.9999, 1.0, 1.0001, 1.5, 3.7, 6.0};
  int same = 1;
  float worst = 0.0f;
  for (uint32_t ch = 1; ch <= 3; ch++) {
    for (int s = 0; s < 7; s++) {
      const uint64_t step = (uint64_t)(SPEED[s] * (double)UDSP_VSINC_ONE);
      /* from before the window's start, from its middle */
      static const int64_t FROM[] = {-40, 900};
      for (int f = 0; f < 2; f++) {
        const size_t n = SPEED[s] > 4.0 ? OUT / 2 : OUT;
        udsp_vsinc_run(&g_vs, in, ch, 0, IN / 2, FROM[f], 0x9e3779b9u, step,
                       got, ch, n);
        ref_vsinc(&g_vs, in, ch, 0, IN / 2, FROM[f], 0x9e3779b9u, step, want,
                  n);
        same &= !memcmp(got, want, n * ch * sizeof(float));
        for (size_t i = 0; i < n * ch; i++) {
          const float d = fabsf(got[i] - want[i]);
          worst = d > worst ? d : worst;
        }
      }
    }
  }
#if defined(__clang__) && defined(__FMA__)
  /* clang fuses a * b + c by default (DESIGN.md "Bits"): the scalar
   * reference's taps, not the run's lanes */
  UDSP_TEST_ASSERT(worst < 1e-5f,
                   "speeds 0.37 .. 6, 1-3 channels, both "
                   "edges: the tap-at-a-time answers (%g, fused)",
                   (double)worst);
#else
  UDSP_TEST_ASSERT(same,
                   "speeds 0.37 .. 6, 1-3 channels, both edges: the "
                   "tap-at-a-time bits (%g)",
                   (double)worst);
#endif
}

int main(void) {
  test_shape();
  test_dc();
  test_same_rate();
  test_sine();
  test_stopband();
  test_cut();
  test_impulse();
  test_vs_identity();
  test_vs_dc();
  test_vs_sine();
  test_vs_alias();
  test_vs_cut();
  test_vs_reference();
  free(g_table);
  free(g_vtable);
  return udsp_test_summary();
}
