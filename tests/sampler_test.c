#include "udsp_test.h"

#include "sampler.h"

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

int main(void) {
  test_root_note_plays_unchanged();
  test_octave_up_skips_frames();
  test_rate_conversion_interpolates();
  test_one_shot_ends();
  test_loop_wraps_to_its_start();
  test_silent_voice();
  return udsp_test_summary();
}
