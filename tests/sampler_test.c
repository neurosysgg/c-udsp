#include "udsp_test.h"

#include "sampler.h"

#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <string.h>

static float RAMP[16];

static UDSP_SAMPLE ramp_sample(uint32_t rate) {
  for (int i = 0; i < 16; i++) {
    RAMP[i] = (float)i;
  }
  UDSP_SAMPLE s = {RAMP, 16, rate, 1, 0, 0, 60};
  return s;
}

static void test_root_note_plays_unchanged(void) {
  UDSP_SAMPLE s = ramp_sample(48000);
  UDSP_SAMPLER_VOICE v;
  udsp_sampler_note_on(&v, &s, 60, 48000);
  float out[8];
  udsp_sampler_render(&v, out, 4);
  UDSP_TEST_ASSERT(out[0] == 0.0f && out[2] == 1.0f && out[6] == 3.0f,
                   "at the root, frame for frame");
  UDSP_TEST_ASSERT(out[1] == out[0] && out[7] == out[6], "mono on both sides");
}

static void test_octave_up_skips_frames(void) {
  UDSP_SAMPLE s = ramp_sample(48000);
  UDSP_SAMPLER_VOICE v;
  udsp_sampler_note_on(&v, &s, 72, 48000);
  float out[6];
  udsp_sampler_render(&v, out, 3);
  UDSP_TEST_ASSERT(out[0] == 0.0f && out[2] == 2.0f && out[4] == 4.0f,
                   "an octave up reads every second frame");
}

static void test_rate_conversion_interpolates(void) {
  UDSP_SAMPLE s = ramp_sample(24000);
  UDSP_SAMPLER_VOICE v;
  udsp_sampler_note_on(&v, &s, 60, 48000);
  float out[8];
  udsp_sampler_render(&v, out, 4);
  UDSP_TEST_ASSERT(out[0] == 0.0f && out[2] == 0.5f && out[4] == 1.0f &&
                       out[6] == 1.5f,
                   "a 24 kHz sample at 48 kHz: halfway values between frames");
}

static void test_one_shot_ends(void) {
  UDSP_SAMPLE s = ramp_sample(48000);
  UDSP_SAMPLER_VOICE v;
  udsp_sampler_note_on(&v, &s, 60, 48000);
  float out[40];
  udsp_sampler_render(&v, out, 20);
  UDSP_TEST_ASSERT(out[28] == 14.0f, "the last interpolable frame plays");
  UDSP_TEST_ASSERT(out[30] == 0.0f && out[38] == 0.0f && !v.playing,
                   "then silence, and the voice stops");
}

static void test_loop_wraps_to_its_start(void) {
  UDSP_SAMPLE s = ramp_sample(48000);
  s.loop_start = 4;
  s.loop_end = 8;
  UDSP_SAMPLER_VOICE v;
  udsp_sampler_note_on(&v, &s, 60, 48000);
  float out[24];
  udsp_sampler_render(&v, out, 12);
  UDSP_TEST_ASSERT(out[14] == 7.0f && out[16] == 4.0f && out[22] == 7.0f,
                   "frames 4..7 repeat after the first pass");
  UDSP_TEST_ASSERT(v.playing, "a loop keeps playing");
}

static void test_silent_voice(void) {
  UDSP_SAMPLER_VOICE v = {0};
  float out[4] = {1, 1, 1, 1};
  udsp_sampler_render(&v, out, 2);
  UDSP_TEST_ASSERT(out[0] == 0.0f && out[3] == 0.0f, "no note: silence");
  UDSP_SAMPLE s = ramp_sample(0);
  udsp_sampler_note_on(&v, &s, 60, 0);
  UDSP_TEST_ASSERT(!v.playing, "a zero output rate doesn't start");
}

static void test_step(void) {
  UDSP_TEST_ASSERT(udsp_sampler_step_q32(0.0, 48000, 48000) == 1ull << 32 &&
                       udsp_sampler_step_q32(12.0, 48000, 48000) == 2ull
                                                                        << 32 &&
                       udsp_sampler_step_q32(-12.0, 48000, 48000) == 1ull << 31,
                   "steps: unison 1, an octave up 2, down 1/2, exactly");
  uint64_t s = udsp_sampler_step_q32(0.0, 44100, 48000);
  UDSP_TEST_ASSERT(s == (uint64_t)(44100.0 / 48000.0 * 4294967296.0),
                   "44.1 kHz played at 48: 0.91875");
  double r = (double)udsp_sampler_step_q32(0.5, 48000, 48000) / 4294967296.0;
  UDSP_TEST_ASSERT(fabs(r - 1.0293022366) < 1e-9, "a quarter tone: 2^(1/24)");
  UDSP_TEST_ASSERT(udsp_sampler_step_q32(120.0, 48000, 48000) == 64ull << 32 &&
                       udsp_sampler_step_q32(0.0, 48000, 0) == 0,
                   "capped at 64; no rate, no step");
  UDSP_TEST_ASSERT(udsp_resample_window(1ull << 32, 1ull << 32, 4) == 7 &&
                       udsp_resample_window(3ull << 31, 2ull << 32, 3) == 8 &&
                       udsp_resample_window(1ull << 32, 1ull << 32, 0) == 0,
                   "the window: frames to floor(last) + 2");
}

static void test_cubic_exact(void) {
  float w[40 * 2], out[16 * 2];
  for (int i = 0; i < 40; i++) {
    w[i] = 0.25f; /* a constant, mono */
  }
  udsp_resample_cubic(w, 1, (1ull << 32) + 12345, 0x9abcdef0ull, out, 16);
  int ok = 1;
  for (int i = 0; i < 32; i++) {
    ok &= out[i] == 0.25f;
  }
  UDSP_TEST_ASSERT(ok, "a constant stays it at any fraction");
  for (int i = 0; i < 40; i++) {
    w[i * 2] = (float)i;
    w[i * 2 + 1] = 100.0f - (float)i * 2.0f;
  }
  uint64_t step = 3ull << 30; /* 0.75 */
  uint64_t end = udsp_resample_cubic(w, 2, 1ull << 32, step, out, 16);
  ok = end == (1ull << 32) + 16 * step;
  for (int i = 0; i < 16; i++) {
    double pos = 1.0 + 0.75 * i;
    ok &= fabs(out[i * 2] - pos) < 1e-5 &&
          fabs(out[i * 2 + 1] - (100.0 - 2.0 * pos)) < 1e-5;
  }
  UDSP_TEST_ASSERT(ok, "a ramp, stereo: on the line between frames");
  udsp_resample_cubic(w, 2, 5ull << 32, 2ull << 32, out, 8);
  ok = 1;
  for (int i = 0; i < 8; i++) {
    ok &= out[i * 2] == w[(5 + 2 * i) * 2] &&
          out[i * 2 + 1] == w[(5 + 2 * i) * 2 + 1];
  }
  UDSP_TEST_ASSERT(ok, "whole positions: the frames themselves");
}

/* a 1 kHz sine at 44.1 kHz played at 48 kHz: within -80 dB of the true
 * one, where linear interpolation is ~-50 dB */
static void test_cubic_sine(void) {
  enum { N = 2048 };
  static float w[N], out[N * 2];
  for (int i = 0; i < N; i++) {
    w[i] = (float)sin(2.0 * M_PI * 1000.0 * i / 44100.0);
  }
  uint64_t step = udsp_sampler_step_q32(0.0, 44100, 48000);
  udsp_resample_cubic(w, 1, 1ull << 32, step, out, 1500);
  double worst = 0, worst_lin = 0;
  for (int i = 0; i < 1500; i++) {
    double pos = 1.0 + (double)step * i / 4294967296.0;
    double want = sin(2.0 * M_PI * 1000.0 * pos / 44100.0);
    double e = fabs(out[i * 2] - want);
    worst = e > worst ? e : worst;
    int k = (int)pos;
    double lin = w[k] + (w[k + 1] - w[k]) * (pos - k);
    e = fabs(lin - want);
    worst_lin = e > worst_lin ? e : worst_lin;
  }
  UDSP_TEST_ASSERT(worst < 1e-4 && worst < worst_lin / 20,
                   "a sine resampled: off by %.2g (linear %.2g)", worst,
                   worst_lin);
}

int main(void) {
  test_root_note_plays_unchanged();
  test_octave_up_skips_frames();
  test_rate_conversion_interpolates();
  test_one_shot_ends();
  test_loop_wraps_to_its_start();
  test_silent_voice();
  test_step();
  test_cubic_exact();
  test_cubic_sine();
  return udsp_test_summary();
}
