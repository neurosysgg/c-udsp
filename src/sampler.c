#include "sampler.h"

#include <string.h>

#include "pitch.h"

void udsp_sampler_note_on(UDSP_SAMPLER_VOICE *v, const UDSP_SAMPLE *s, int note,
                          uint32_t out_rate) {
  v->smp = s;
  v->pos_q32 = 0;
  v->playing = 0;
  if (!s || !s->frames || s->frame_count == 0 || out_rate == 0) {
    return;
  }
  /* the transposition, then the rate as a separate multiply and divide:
   * folded into one shifted expression it overflows 64 bits and wraps
   * silently (utracker's fifth above middle C came out at 0.13) */
  uint64_t ratio = udsp_pitch_ratio_q32(note - (int)s->root_note);
  uint64_t step = ratio * (uint64_t)s->sample_rate / (uint64_t)out_rate;
  const uint64_t step_max = 64ull << 32;
  v->step_q32 = step > step_max ? step_max : step;
  v->playing = 1;
}

void udsp_sampler_stop(UDSP_SAMPLER_VOICE *v) { v->playing = 0; }

void udsp_sampler_render(UDSP_SAMPLER_VOICE *v, float *out, size_t frames) {
  const UDSP_SAMPLE *s = v->smp;
  if (!v->playing || !s || !s->frames) {
    memset(out, 0, frames * 2 * sizeof(float));
    return;
  }
  const uint32_t ch = s->channels;
  const int looping =
      s->loop_end > s->loop_start && s->loop_end <= s->frame_count;
  const uint64_t loop_end_q32 = (uint64_t)s->loop_end << 32;
  const uint64_t loop_len_q32 = loop_end_q32 - ((uint64_t)s->loop_start << 32);

  for (size_t f = 0; f < frames; f++) {
    if (looping) {
      /* minus the loop's length, not back to its start: the fraction is
       * kept, or a loop that isn't a whole number of output frames sheds
       * one each time round (a slow detune) */
      while (v->pos_q32 >= loop_end_q32) {
        v->pos_q32 -= loop_len_q32;
      }
    } else if ((v->pos_q32 >> 32) + 1 >= s->frame_count) {
      v->playing = 0;
      memset(out + f * 2, 0, (frames - f) * 2 * sizeof(float));
      return;
    }
    const uint64_t idx = v->pos_q32 >> 32;
    const float frac =
        (float)(uint32_t)(v->pos_q32 & 0xffffffffull) * (1.0f / 4294967296.0f);
    const float *a = s->frames + idx * ch;
    /* past a loop's last frame comes its first, not the file's next */
    const float *b = looping && idx + 1 >= s->loop_end
                         ? s->frames + (uint64_t)s->loop_start * ch
                         : a + ch;
    const float l = a[0] + (b[0] - a[0]) * frac;
    out[f * 2] = l;
    out[f * 2 + 1] = ch == 2 ? a[1] + (b[1] - a[1]) * frac : l;
    v->pos_q32 += v->step_q32;
  }
}

uint64_t udsp_sampler_step_q32(double semitones, uint32_t src_rate,
                               uint32_t out_rate) {
  if (out_rate == 0) {
    return 0;
  }
  double step = udsp_exp2(semitones / 12.0) * (double)src_rate /
                (double)out_rate * 4294967296.0;
  const double step_max = 64.0 * 4294967296.0;
  return step >= step_max ? (uint64_t)step_max : (uint64_t)step;
}

uint64_t udsp_resample_window(uint64_t pos_q32, uint64_t step_q32,
                              size_t frames) {
  if (frames == 0) {
    return 0;
  }
  uint64_t last = pos_q32 + step_q32 * (uint64_t)(frames - 1);
  return (last >> 32) + 3;
}

static inline float hermite(float ym1, float y0, float y1, float y2, float x) {
  const float c1 = 0.5f * (y1 - ym1);
  const float c2 = ym1 - 2.5f * y0 + 2.0f * y1 - 0.5f * y2;
  const float c3 = 0.5f * (y2 - ym1) + 1.5f * (y0 - y1);
  return ((c3 * x + c2) * x + c1) * x + y0;
}

uint64_t udsp_resample_cubic(const float *window, uint32_t ch, uint64_t pos_q32,
                             uint64_t step_q32, float *out, size_t frames) {
  for (size_t f = 0; f < frames; f++) {
    const uint64_t idx = pos_q32 >> 32;
    const float x =
        (float)(uint32_t)(pos_q32 & 0xffffffffull) * (1.0f / 4294967296.0f);
    const float *p = window + (idx - 1) * ch;
    const float l = hermite(p[0], p[ch], p[2 * ch], p[3 * ch], x);
    out[f * 2] = l;
    out[f * 2 + 1] =
        ch == 2 ? hermite(p[1], p[ch + 1], p[2 * ch + 1], p[3 * ch + 1], x) : l;
    pos_q32 += step_q32;
  }
  return pos_q32;
}
