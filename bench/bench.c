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

#include "analyze.h"
#include "ballistics.h"
#include "fft.h"
#include "loudness.h"
#include "scope.h"
#include "spectrogram.h"
#include "spectrum.h"

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

static void k_spectrogram(void) {
  udsp_spectrogram_push(&g_sg, g_in, 48000, g_re, g_im);
  g_sink = g_sg.intensity[0][0];
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
