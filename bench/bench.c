/* bench: what each c-µdsp kernel costs, on the machine and flags it was
 * built with. Each kernel is timed on synthetic audio: its iteration count
 * calibrated to ~10 ms per run, then REPS runs, reporting the median and
 * the fastest per call and per sample. One line per kernel, tab-separated
 * (bench/matrix.py collects the builds into one table):
 *
 *   kernel  n  median_ns  min_ns  ns_per_sample
 *
 * Run pinned to one core (taskset) on an otherwise quiet machine; the
 * matrix script does. -k NAME runs only the kernels whose name has NAME
 * in it; -r REPS sets the runs (default 21). */
#define _DEFAULT_SOURCE

#include "adsr.h"
#include "analyze.h"
#include "ballistics.h"
#include "delay.h"
#include "dyn.h"
#include "eq.h"
#include "fdn.h"
#include "fft.h"
#include "gain.h"
#include "lfo.h"
#include "loudness.h"
#include "openv.h"
#include "osc.h"
#include "r128.h"
#include "ramp.h"
#include "sampler.h"
#include "scope.h"
#include "shape.h"
#include "sinc.h"
#include "spectrogram.h"
#include "spectrum.h"
#include "svf.h"
#include "wtab.h"
#include "yin.h"
#include "zdf.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* the transform size of every kernel that takes one: 2048, what the
 * baseline was measured at (UDSP_FFT_N was 2048 then) and what live
 * displays use; the old spectrum API is fixed at UDSP_FFT_N (8192 since
 * ctui-mus's c16af1d), timed as *_old8k */
#define N 2048
#define BIG (1u << 20) /* a long stretch: ~22 s at 48 kHz, 4 MB of floats */
#define LUFS_N 19200   /* 400 ms at 48 kHz: one momentary-loudness window */

static double now_ns(void) {
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return (double)t.tv_sec * 1e9 + (double)t.tv_nsec;
}

/* results land here so no kernel is optimised away */
static volatile float g_sink;

static float *g_in, *g_stereo, *g_re, *g_im, *g_out, *g_big, *g_table,
    *g_hann;
static float g_level[64], g_target[64], g_peak[64], g_age[64];
static UDSP_SPECTROGRAM g_sg;
static UDSP_R128 g_r128;

static unsigned int g_seed = 12345;
static float noise(void) {
  g_seed = g_seed * 1103515245u + 12345u;
  return (float)((g_seed >> 8) & 0xffff) / 32768.0f - 1.0f;
}

/* music-ish: a bass tone, a wobbling melody, some noise */
static void fill(float *s, size_t n, int ch) {
  for (size_t i = 0; i < n; i++) {
    double t = (double)(i / (size_t)ch) / 48000.0;
    s[i] = (float)(0.4 * sin(2 * M_PI * 55.0 * t) +
                   0.3 * sin(2 * M_PI * (440.0 + 20.0 * sin(t * 3)) * t)) +
           0.05f * noise();
  }
}

typedef struct {
  const char *name;
  size_t n; /* samples (frames for downmix) one call covers */
  void (*run)(void);
} KERNEL;

/* --- the kernels, as their callers use them --- */

static void k_hann(void) {
  udsp_window_hann(g_in, g_re, N);
  g_sink = g_re[100];
}

static void k_fft(void) {
  /* a copy in, then the transform: the copy (~100 ns) is part of every
   * real caller too, and keeps the input from blowing up run after run */
  memcpy(g_re, g_in, N * sizeof(float));
  memset(g_im, 0, N * sizeof(float));
  udsp_fft_forward(g_re, g_im, N);
  g_sink = g_re[7];
}

static void k_hann_table(void) {
  udsp_window_apply(g_in, g_hann, g_re, N);
  g_sink = g_re[100];
}

static void k_fft_table(void) {
  memcpy(g_re, g_in, N * sizeof(float));
  memset(g_im, 0, N * sizeof(float));
  udsp_fft_forward_table(g_re, g_im, N, g_table);
  g_sink = g_re[7];
}

static void k_fft_real(void) {
  udsp_fft_real(g_in, g_re, g_im, N, g_table);
  g_sink = g_re[7];
}

static void k_power(void) {
  udsp_fft_power(g_in, g_in + N / 2, N / 2 + 1, g_out);
  g_sink = g_out[3];
}

/* what a spectrum takes now: window table, real FFT, power */
static void k_block(void) {
  udsp_window_apply(g_in, g_hann, g_out, N);
  udsp_fft_real(g_out, g_re, g_im, N, g_table);
  udsp_fft_power(g_re, g_im, N / 2 + 1, g_re);
  g_sink = g_re[9];
}

/* the same with the window applied inside the transform */
static void k_block_fused(void) {
  udsp_fft_real_windowed(g_in, g_hann, g_re, g_im, N, g_table);
  udsp_fft_power(g_re, g_im, N / 2 + 1, g_re);
  g_sink = g_re[9];
}

static void k_magnitude(void) {
  udsp_fft_magnitude(g_in, g_in + N / 2, N / 2, g_out);
  g_sink = g_out[3];
}

static void k_bars(void) {
  float bars[28];
  udsp_spectrum_bars(g_in, 48000, g_re, g_im, bars, 28);
  g_sink = bars[5];
}

static void k_bars128(void) {
  float bars[128];
  udsp_spectrum_bars(g_in, 48000, g_re, g_im, bars, 128);
  g_sink = bars[5];
}

static UDSP_SPECTRUM_PLAN g_plan28 __attribute__((aligned(64))),
    g_plan128 __attribute__((aligned(64)));
static float g_scratch[UDSP_SPECTRUM_SCRATCH] __attribute__((aligned(64)));

static void k_bars_planned(void) {
  float bars[28];
  udsp_spectrum_bars_planned(&g_plan28, g_in, g_scratch, bars);
  g_sink = bars[5];
}

static void k_bars128_planned(void) {
  float bars[128];
  udsp_spectrum_bars_planned(&g_plan128, g_in, g_scratch, bars);
  g_sink = bars[5];
}

static UDSP_SPECTRUM_PLAN g_plan8k __attribute__((aligned(64)));
static void k_bars128_planned8k(void) {
  float bars[128];
  udsp_spectrum_bars_planned(&g_plan8k, g_in, g_scratch, bars);
  g_sink = bars[5];
}

static void k_rms(void) { g_sink = udsp_analyze_rms(g_in, N); }
static void k_rms_big(void) { g_sink = udsp_analyze_rms(g_big, BIG); }
static void k_peak(void) { g_sink = udsp_analyze_peak(g_in, N); }
static void k_peak_big(void) { g_sink = udsp_analyze_peak(g_big, BIG); }

static void k_sum(void) { g_sink = (float)udsp_analyze_sum(g_in, N / 2); }

static void k_minmax(void) {
  float lo, hi;
  udsp_analyze_minmax(g_in, N, &lo, &hi);
  g_sink = lo + hi;
}
static void k_minmax_big(void) {
  float lo, hi;
  udsp_analyze_minmax(g_big, BIG, &lo, &hi);
  g_sink = lo + hi;
}

static void k_downmix(void) {
  udsp_mono_downmix(g_stereo, N, 2, g_out);
  g_sink = g_out[9];
}

static void k_dbfs(void) {
  float acc = 0;
  for (int i = 0; i < 28; i++) {
    acc += udsp_analyze_dbfs(g_in[i] * g_in[i]);
  }
  g_sink = acc;
}

static void k_scope(void) {
  udsp_scope_extract(g_in, N, g_out);
  g_sink = g_out[1];
}

static void k_ballistics(void) {
  udsp_ballistics(g_level, g_target, 28, 1.0f / 30, 0.015f, 0.25f);
  udsp_peak_hold(g_peak, g_age, g_level, 28, 1.0f / 30, 0.7f, 30.0f);
  g_sink = g_peak[3];
}

static void k_lufs(void) {
  g_sink = udsp_loudness_lufs_momentary(g_big, LUFS_N, g_out);
}

/* a program's loudness: the long stretch as stereo, ~11 s */
static void k_r128(void) {
  udsp_r128_init(&g_r128, 48000, 2);
  udsp_r128_add(&g_r128, g_big, BIG / 2);
  g_sink = (float)udsp_r128_integrated(&g_r128.hist);
}

static void k_spectrogram(void) {
  udsp_spectrogram_push(&g_sg, g_in, 48000, g_re, g_im);
  g_sink = g_sg.intensity[0][0];
}

/* --- the synthesis half: N mono samples or N / 2 stereo frames --- */

static UDSP_OSC g_osc;
static void osc_wave(UDSP_OSC_WAVE w) {
  g_osc.wave = w;
  g_osc.inc = udsp_osc_inc_hz(440.0, 48000);
  udsp_osc_render(&g_osc, g_out, N);
  g_sink = g_out[7];
}
static void k_osc_sine(void) { osc_wave(UDSP_OSC_SINE); }
static void k_osc_saw(void) { osc_wave(UDSP_OSC_SAW); }
static void k_osc_square(void) { osc_wave(UDSP_OSC_SQUARE); }
static void k_osc_triangle(void) { osc_wave(UDSP_OSC_TRIANGLE); }
static void k_osc_noise(void) { osc_wave(UDSP_OSC_NOISE); }

static UDSP_SAMPLE g_smp;
static UDSP_SAMPLER_VOICE g_voice;
static void k_sampler(void) {
  if (!g_voice.playing || g_voice.pos_q32 > ((uint64_t)BIG / 4 << 32)) {
    udsp_sampler_note_on(&g_voice, &g_smp, 67, 48000);
  }
  udsp_sampler_render(&g_voice, g_out, N / 2);
  g_sink = g_out[5];
}

static uint64_t g_rpos = 1ull << 32;
static void k_resample(void) {
  uint64_t step = udsp_sampler_step_q32(7.0, 44100, 48000);
  if (udsp_resample_window(g_rpos, step, N / 2) + 4 > BIG / 2) {
    g_rpos = 1ull << 32;
  }
  g_rpos = udsp_resample_cubic(g_smp.frames, 2, g_rpos, step, g_out, N / 2);
  g_sink = g_out[5];
}

/* the sinc resampler, 44.1 -> 48 kHz, one channel (a clip imported) */
static UDSP_SINC g_sinc;
static float *g_sinc_table;
static uint64_t g_sinc_pos;
static void k_sinc(void) {
  if (!g_sinc_table) {
    g_sinc_table = malloc(udsp_sinc_floats(44100, 48000, UDSP_SINC_ZEROS) *
                          sizeof(float));
    udsp_sinc_init(&g_sinc, 44100, 48000, UDSP_SINC_ZEROS, g_sinc_table);
  }
  if (g_sinc_pos + N > udsp_sinc_length(&g_sinc, BIG / 2)) {
    g_sinc_pos = 0;
  }
  udsp_sinc_run(&g_sinc, g_smp.frames, 0, BIG / 2, g_sinc_pos, g_out, N);
  g_sinc_pos += N;
  g_sink = g_out[5];
}

/* the varispeed sinc at 1.0595 (a semitone up), stereo: a warped clip
 * (N outputs: N / 2 frames of two channels) */
static UDSP_VSINC g_vs;
static float *g_vs_table;
static uint64_t g_vs_pos;
static void k_vsinc(void) {
  if (!g_vs_table) {
    g_vs_table = malloc(udsp_vsinc_floats(UDSP_VSINC_ZEROS) * sizeof(float));
    udsp_vsinc_init(&g_vs, UDSP_VSINC_ZEROS, g_vs_table);
  }
  const uint64_t step = (uint64_t)(1.0595 * UDSP_VSINC_ONE);
  if ((g_vs_pos >> 32) + N > BIG / 4) {
    g_vs_pos = 0;
  }
  udsp_vsinc_run(&g_vs, g_smp.frames, 2, 0, BIG / 4,
                 (int64_t)(g_vs_pos >> 32) + 64, (uint32_t)g_vs_pos, step,
                 g_out, 2, N / 2);
  g_vs_pos += (uint64_t)(N / 2) * step;
  g_sink = g_out[5];
}

static float g_yin_scratch[44100 / 40 + 2];
static void k_yin(void) {
  g_sink = udsp_yin(g_smp.frames, 4096, 44100, 40.0f, 4000.0f, 0.15f,
                    g_yin_scratch, NULL);
}

static UDSP_SVF g_svf;
static void k_svf(void) {
  memcpy(g_out, g_stereo, N * sizeof(float));
  udsp_svf_stereo(&g_svf, g_out, N / 2);
  g_sink = g_out[3];
}

/* FM: a sine carrier, its phase moved by a modulator's swing */
static UDSP_OSC g_carrier;
static void k_osc_pm(void) {
  udsp_osc_render_pm(&g_carrier, g_in, g_out, N);
  g_sink = g_out[5];
}

/* Operator's operators (d8): a sine with feedback (serial), a saw's
 * table under PM, and its envelope */
static UDSP_WTAB g_wtab;
static UDSP_FMOP g_op_fb, g_op_tab;
static float g_amp[N];
static void k_fmop_fb(void) {
  udsp_fmop_render(&g_op_fb, NULL, g_amp, g_out, N);
  g_sink = g_out[5];
}
static void k_fmop_table(void) {
  udsp_fmop_render(&g_op_tab, g_in, g_amp, g_out, N);
  g_sink = g_out[5];
}
static UDSP_OPENV g_openv;
static void k_openv(void) {
  udsp_openv_gate(&g_openv, g_openv.stage != UDSP_OPENV_RETURN);
  udsp_openv_render(&g_openv, g_out, N);
  g_sink = g_out[5];
}

/* a bell and a 48 dB cut: one section and four */
static UDSP_EQ_BAND g_eq_bell, g_eq_cut;
static void k_eq(void) {
  memcpy(g_out, g_in, N * sizeof(float));
  udsp_eq_stereo(&g_eq_bell, g_out, N / 2);
  g_sink = g_out[3];
}
static void k_eq48(void) {
  memcpy(g_out, g_in, N * sizeof(float));
  udsp_eq_stereo(&g_eq_cut, g_out, N / 2);
  g_sink = g_out[3];
}

/* a compressor's gain from a stereo key (peak), and a limiter over 1 ms
 * of lookahead */
static UDSP_COMP g_comp;
static void k_comp(void) {
  udsp_comp_gain(&g_comp, g_in, g_out, N / 2);
  g_sink = g_out[3];
}
static UDSP_LIMITER g_lim;
static float g_lim_mem[UDSP_LIMITER_FLOATS(48)];
static void k_limiter(void) {
  memcpy(g_out, g_in, N * sizeof(float));
  udsp_limiter_stereo(&g_lim, g_out, N / 2);
  g_sink = g_out[3];
}

/* the reverb: eight modulated lines, the diffusers, a predelay */
static UDSP_FDN g_fdn;
static float *g_fdn_mem;
static void k_fdn(void) {
  udsp_fdn_stereo(&g_fdn, g_in, g_out, N / 2);
  g_sink = g_out[3];
}

static UDSP_ZDF g_zdf;
static void k_zdf(void) {
  memcpy(g_out, g_in, N * sizeof(float));
  udsp_zdf_mono(&g_zdf, g_out, N);
  g_sink = g_out[3];
}

/* an envelope that never settles: attack, decay and release in turn */
static UDSP_ADSR g_adsr;
static void k_adsr(void) {
  udsp_adsr_gate(&g_adsr, g_adsr.stage != UDSP_ADSR_SUSTAIN);
  udsp_adsr_render(&g_adsr, g_out, N);
  g_sink = g_out[7];
}

static UDSP_DELAY g_delay;
static float *g_ring;
static void k_delay(void) {
  memcpy(g_out, g_stereo, N * sizeof(float));
  udsp_delay_stereo(&g_delay, g_out, N / 2, 12000, 0.4f, 0.35f);
  g_sink = g_out[3];
}

static void k_softclip(void) {
  memcpy(g_out, g_in, N * sizeof(float));
  udsp_softclip(g_out, N, 3.0f);
  g_sink = g_out[9];
}

static void k_mix_add(void) {
  udsp_mix_add(g_out, g_in, N, 0.5f);
  g_sink = g_out[9];
}

static void k_gain_stereo(void) {
  memcpy(g_out, g_stereo, N * sizeof(float));
  udsp_gain_stereo(g_out, N / 2, 0.7f, 0.3f);
  g_sink = g_out[9];
}

static UDSP_RAMP g_ramp = {0.8f, 0.8f, 1.0f / 144};
static void k_ramp_held(void) {
  memcpy(g_out, g_stereo, N * sizeof(float));
  udsp_ramp_stereo(&g_ramp, g_out, N / 2);
  g_sink = g_out[9];
}

static const KERNEL KERNELS[] = {
    {"hann", N, k_hann},
    {"hann_table", N, k_hann_table},
    {"fft", N, k_fft},
    {"fft_table", N, k_fft_table},
    {"fft_real", N, k_fft_real},
    {"power", N / 2 + 1, k_power},
    {"block", N, k_block},
    {"block_fused", N, k_block_fused},
    {"magnitude", N / 2, k_magnitude},
    {"bars28_old8k", UDSP_FFT_N, k_bars},
    {"bars28_planned", N, k_bars_planned},
    {"bars128_old8k", UDSP_FFT_N, k_bars128},
    {"bars128_planned", N, k_bars128_planned},
    {"bars128_planned8k", UDSP_FFT_N, k_bars128_planned8k},
    {"spectrogram_old8k", UDSP_FFT_N, k_spectrogram},
    {"rms", N, k_rms},
    {"rms_big", BIG, k_rms_big},
    {"peak", N, k_peak},
    {"peak_big", BIG, k_peak_big},
    {"sum", N / 2, k_sum},
    {"minmax", N, k_minmax},
    {"minmax_big", BIG, k_minmax_big},
    {"downmix2", N, k_downmix},
    {"dbfs28", 28, k_dbfs},
    {"scope", N, k_scope},
    {"ballistics28", 28, k_ballistics},
    {"lufs", LUFS_N, k_lufs},
    {"r128", BIG, k_r128},
    {"osc_sine", N, k_osc_sine},
    {"osc_saw", N, k_osc_saw},
    {"osc_square", N, k_osc_square},
    {"osc_triangle", N, k_osc_triangle},
    {"osc_noise", N, k_osc_noise},
    {"sampler", N, k_sampler},
    {"resample", N, k_resample},
    {"sinc", N, k_sinc},
    {"vsinc", N, k_vsinc},
    {"yin4096", 4096, k_yin},
    {"svf", N, k_svf},
    {"eq", N, k_eq},
    {"eq48", N, k_eq48},
    {"comp", N / 2, k_comp},
    {"limiter", N / 2, k_limiter},
    {"fdn", N / 2, k_fdn},
    {"zdf", N, k_zdf},
    {"osc_pm", N, k_osc_pm},
    {"adsr", N, k_adsr},
    {"fmop_fb", N, k_fmop_fb},
    {"fmop_table", N, k_fmop_table},
    {"openv", N, k_openv},
    {"delay", N, k_delay},
    {"softclip", N, k_softclip},
    {"mix_add", N, k_mix_add},
    {"gain_stereo", N, k_gain_stereo},
    {"ramp_held", N, k_ramp_held},
};

static int cmp(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return x < y ? -1 : x > y;
}

static void time_kernel(const KERNEL *k, int reps) {
  /* calibrate: double the iterations until one run takes >= 10 ms */
  long iters = 1;
  for (;;) {
    double t0 = now_ns();
    for (long i = 0; i < iters; i++) {
      k->run();
    }
    if (now_ns() - t0 >= 1e7 || iters >= (1L << 30)) {
      break;
    }
    iters *= 2;
  }
  double *per = malloc((size_t)reps * sizeof(double));
  for (int r = 0; r < reps; r++) {
    double t0 = now_ns();
    for (long i = 0; i < iters; i++) {
      k->run();
    }
    per[r] = (now_ns() - t0) / (double)iters;
  }
  qsort(per, (size_t)reps, sizeof(double), cmp);
  double med = per[reps / 2];
  printf("%s\t%zu\t%.1f\t%.1f\t%.4f\n", k->name, k->n, med, per[0],
         med / (double)k->n);
  free(per);
}

int main(int argc, char **argv) {
  const char *only = NULL;
  int reps = 21;
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-k") && i + 1 < argc) {
      only = argv[++i];
    } else if (!strcmp(argv[i], "-r") && i + 1 < argc) {
      reps = atoi(argv[++i]);
    } else {
      fprintf(stderr, "usage: %s [-k NAME] [-r REPS]\n", argv[0]);
      return 2;
    }
  }
  if (reps < 1) {
    reps = 1;
  }
  size_t fft = UDSP_FFT_N * sizeof(float);
  g_in = aligned_alloc(64, fft);
  g_stereo = aligned_alloc(64, 2 * fft);
  g_re = aligned_alloc(64, fft);
  g_im = aligned_alloc(64, fft);
  g_out = aligned_alloc(64, LUFS_N * sizeof(float));
  g_big = aligned_alloc(64, BIG * sizeof(float));
  g_table = aligned_alloc(64, UDSP_FFT_TABLE_FLOATS(N) * sizeof(float));
  g_hann = aligned_alloc(64, fft);
  udsp_fft_table(g_table, N);
  udsp_window_hann_table(g_hann, N);
  fill(g_in, UDSP_FFT_N, 1);
  fill(g_stereo, 2 * UDSP_FFT_N, 2);
  fill(g_big, BIG, 1);
  for (int i = 0; i < 64; i++) {
    g_target[i] = -60.0f + (float)(i % 7) * 9.0f;
    g_level[i] = g_peak[i] = -30.0f;
    g_age[i] = 0.0f;
  }
  udsp_spectrogram_init(&g_sg);
  udsp_osc_init(&g_osc, UDSP_OSC_SINE, 1);
  g_smp = (UDSP_SAMPLE){g_big, BIG / 2, 44100, 2, 0, 0, 60};
  g_fdn_mem = calloc(udsp_fdn_floats(48000), sizeof(float));
  udsp_fdn_init(&g_fdn, g_fdn_mem, 48000);
  udsp_fdn_set(&g_fdn, 1.0f, 2.0f, 6000.0f, 0.6f, 0.5f, 10.0f, 1.0f);
  udsp_comp_reset(&g_comp);
  udsp_comp_set(&g_comp, -18.0f, 4.0f, 6.0f, 5.0f, 100.0f, 0, 48000);
  udsp_limiter_init(&g_lim, g_lim_mem, 48);
  udsp_limiter_set(&g_lim, 48, 0.5f, 50.0f, 48000);
  udsp_eq_reset(&g_eq_bell);
  udsp_eq_set(&g_eq_bell, UDSP_EQ_BELL, 1200.0f, 1.0f, 6.0f, 48000);
  udsp_eq_reset(&g_eq_cut);
  udsp_eq_set(&g_eq_cut, UDSP_EQ_LOW_CUT_48, 120.0f, 0.71f, 0.0f, 48000);
  udsp_svf_reset(&g_svf);
  udsp_svf_set(&g_svf, 1200.0f, 0.5f, UDSP_SVF_LOW, 48000);
  udsp_zdf_reset(&g_zdf);
  udsp_zdf_set(&g_zdf, 1200.0f, 0.5f, UDSP_ZDF_LOW, 48000);
  udsp_osc_init(&g_carrier, UDSP_OSC_SINE, 1);
  g_carrier.inc = udsp_osc_inc_hz(440.0, 48000);
  float harm[UDSP_WTAB_HARMONICS];
  udsp_wtab_recipe(UDSP_WTAB_SAW, 64, harm);
  udsp_wtab_build(&g_wtab, harm, 64);
  for (size_t i = 0; i < N; i++) {
    g_amp[i] = 0.8f;
  }
  g_op_fb = (UDSP_FMOP){NULL, 0, udsp_osc_inc_hz(440.0, 48000), 0.3f, 0, 0};
  g_op_tab = (UDSP_FMOP){&g_wtab, 0, udsp_osc_inc_hz(110.0, 48000), 0, 0, 0};
  UDSP_OPENV_SHAPE sh = {0.1f,  1.0f,  0.5f, 0.0f,  5.0f, 30.0f,
                         40.0f, 20.0f, 0.5f, -0.5f, 0.3f, UDSP_OPENV_LOOP};
  udsp_openv_init(&g_openv, &sh, 48000);
  memset(&g_adsr, 0, sizeof(g_adsr));
  udsp_adsr_set(&g_adsr, 5.0f, 30.0f, 0.5f, 40.0f, 48000);
  g_ring = aligned_alloc(64, UDSP_DELAY_RING_FLOATS(24000) * sizeof(float));
  udsp_delay_init(&g_delay, g_ring, 24000);
  memset(g_out, 0, LUFS_N * sizeof(float));
  udsp_spectrum_plan(&g_plan28, N, 48000, 28);
  udsp_spectrum_plan(&g_plan128, N, 48000, 128);
  udsp_spectrum_plan(&g_plan8k, UDSP_FFT_N, 48000, 128);

  for (size_t i = 0; i < sizeof(KERNELS) / sizeof(KERNELS[0]); i++) {
    if (!only || strstr(KERNELS[i].name, only)) {
      time_kernel(&KERNELS[i], reps);
    }
  }
  return 0;
}
