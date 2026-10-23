/* reduce: how to write a float reduction so it runs at the machine's
 * speed without -ffast-math -- the study behind c-µdsp's kernels
 * (bench/README.md "Reductions").
 *
 * Without -ffast-math a compiler may not reassociate float adds, so
 * `for (i) sum += x[i]*x[i]` is one dependency chain: one add per add
 * latency (3-4 cycles), vectorised or not. The fix is independent
 * accumulators, written so the compiler keeps them apart: an array of
 * lanes (plain C, the SLP vectoriser packs them into registers) or GCC
 * vector extensions (explicit, every compiler level). This times each
 * shape on the same data, for a sum of squares and a max of |x|:
 *
 *   reduce [-r REPS]   ->  shape  n  median_ns  ns_per_sample */
#define _DEFAULT_SOURCE

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef float v8f __attribute__((vector_size(32), aligned(4)));
typedef int v8i __attribute__((vector_size(32), aligned(4)));

/* a > b ? a : b per lane (vector ?: is C++ only in gcc): the compiler
 * turns the blend back into maxps */
static inline v8f vmax(v8f a, v8f b) {
  v8i m = a > b;
  return (v8f)((m & (v8i)a) | (~m & (v8i)b));
}

static double now_ns(void) {
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return (double)t.tv_sec * 1e9 + (double)t.tv_nsec;
}

static volatile float g_sink;

/* --- sums of squares --- */

__attribute__((noinline)) static float ss_double(const float *x, size_t n) {
  double s = 0;
  for (size_t i = 0; i < n; i++) {
    s += (double)x[i] * x[i];
  }
  return (float)s;
}

__attribute__((noinline)) static float ss_float(const float *x, size_t n) {
  float s = 0;
  for (size_t i = 0; i < n; i++) {
    s += x[i] * x[i];
  }
  return s;
}

#define SS_LANES(L)                                                            \
  __attribute__((noinline)) static float ss_lanes##L(const float *x,           \
                                                     size_t n) {               \
    float acc[L] = {0};                                                        \
    size_t i = 0;                                                              \
    for (; i + L <= n; i += L) {                                               \
      for (int j = 0; j < L; j++) {                                            \
        acc[j] += x[i + j] * x[i + j];                                         \
      }                                                                        \
    }                                                                          \
    float s = 0;                                                               \
    for (int j = 0; j < L; j++) {                                              \
      s += acc[j];                                                             \
    }                                                                          \
    for (; i < n; i++) {                                                       \
      s += x[i] * x[i];                                                        \
    }                                                                          \
    return s;                                                                  \
  }
SS_LANES(4)
SS_LANES(8)
SS_LANES(16)
SS_LANES(32)

#define SS_VEC(A)                                                              \
  __attribute__((noinline)) static float ss_vec##A(const float *x, size_t n) { \
    v8f acc[A];                                                                \
    memset(acc, 0, sizeof(acc));                                               \
    size_t i = 0;                                                              \
    for (; i + 8 * A <= n; i += 8 * A) {                                       \
      for (int a = 0; a < A; a++) {                                            \
        v8f v = *(const v8f *)(x + i + 8 * a);                                 \
        acc[a] += v * v;                                                       \
      }                                                                        \
    }                                                                          \
    for (int a = 1; a < A; a++) {                                              \
      acc[0] += acc[a];                                                        \
    }                                                                          \
    float s = 0;                                                               \
    for (int j = 0; j < 8; j++) {                                              \
      s += acc[0][j];                                                          \
    }                                                                          \
    for (; i < n; i++) {                                                       \
      s += x[i] * x[i];                                                        \
    }                                                                          \
    return s;                                                                  \
  }
SS_VEC(1)
SS_VEC(2)
SS_VEC(4)

/* --- max |x| --- */

__attribute__((noinline)) static float mx_plain(const float *x, size_t n) {
  float m = 0;
  for (size_t i = 0; i < n; i++) {
    float a = fabsf(x[i]);
    if (a > m) {
      m = a;
    }
  }
  return m;
}

#define MX_LANES(L)                                                            \
  __attribute__((noinline)) static float mx_lanes##L(const float *x,           \
                                                     size_t n) {               \
    float acc[L] = {0};                                                        \
    size_t i = 0;                                                              \
    for (; i + L <= n; i += L) {                                               \
      for (int j = 0; j < L; j++) {                                            \
        float a = fabsf(x[i + j]);                                             \
        acc[j] = a > acc[j] ? a : acc[j];                                      \
      }                                                                        \
    }                                                                          \
    float m = 0;                                                               \
    for (int j = 0; j < L; j++) {                                              \
      m = acc[j] > m ? acc[j] : m;                                             \
    }                                                                          \
    for (; i < n; i++) {                                                       \
      float a = fabsf(x[i]);                                                   \
      m = a > m ? a : m;                                                       \
    }                                                                          \
    return m;                                                                  \
  }
MX_LANES(8)
MX_LANES(16)
MX_LANES(32)

#define MX_VEC(A)                                                              \
  __attribute__((noinline)) static float mx_vec##A(const float *x, size_t n) { \
    v8f acc[A];                                                                \
    memset(acc, 0, sizeof(acc));                                               \
    const v8i abs = {~0u >> 1, ~0u >> 1, ~0u >> 1, ~0u >> 1,                   \
                     ~0u >> 1, ~0u >> 1, ~0u >> 1, ~0u >> 1};                  \
    size_t i = 0;                                                              \
    for (; i + 8 * A <= n; i += 8 * A) {                                       \
      for (int a = 0; a < A; a++) {                                            \
        v8f v = *(const v8f *)(x + i + 8 * a);                                 \
        acc[a] = vmax((v8f)((v8i)v & abs), acc[a]);                            \
      }                                                                        \
    }                                                                          \
    for (int a = 1; a < A; a++) {                                              \
      acc[0] = vmax(acc[a], acc[0]);                                           \
    }                                                                          \
    float m = 0;                                                               \
    for (int j = 0; j < 8; j++) {                                              \
      m = acc[0][j] > m ? acc[0][j] : m;                                       \
    }                                                                          \
    for (; i < n; i++) {                                                       \
      float a = fabsf(x[i]);                                                   \
      m = a > m ? a : m;                                                       \
    }                                                                          \
    return m;                                                                  \
  }
MX_VEC(1)
MX_VEC(2)
MX_VEC(4)

typedef struct {
  const char *name;
  float (*fn)(const float *, size_t);
} SHAPE;

static const SHAPE SHAPES[] = {
    {"ss_double", ss_double},   {"ss_float", ss_float},
    {"ss_lanes4", ss_lanes4},   {"ss_lanes8", ss_lanes8},
    {"ss_lanes16", ss_lanes16}, {"ss_lanes32", ss_lanes32},
    {"ss_vec1", ss_vec1},       {"ss_vec2", ss_vec2},
    {"ss_vec4", ss_vec4},       {"mx_plain", mx_plain},
    {"mx_lanes8", mx_lanes8},   {"mx_lanes16", mx_lanes16},
    {"mx_lanes32", mx_lanes32}, {"mx_vec1", mx_vec1},
    {"mx_vec2", mx_vec2},       {"mx_vec4", mx_vec4},
};

static int cmp(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return x < y ? -1 : x > y;
}

int main(int argc, char **argv) {
  int reps = 21;
  if (argc == 3 && !strcmp(argv[1], "-r")) {
    reps = atoi(argv[2]);
  }
  static const size_t SIZES[] = {2048, 1u << 20};
  float *x = aligned_alloc(64, (1u << 20) * sizeof(float));
  unsigned int seed = 1;
  for (size_t i = 0; i < (1u << 20); i++) {
    seed = seed * 1103515245u + 12345u;
    x[i] = (float)((seed >> 8) & 0xffff) / 32768.0f - 1.0f;
  }
  double *per = malloc((size_t)reps * sizeof(double));
  for (size_t s = 0; s < 2; s++) {
    size_t n = SIZES[s];
    for (size_t k = 0; k < sizeof(SHAPES) / sizeof(SHAPES[0]); k++) {
      long iters = 1;
      for (;;) {
        double t0 = now_ns();
        for (long i = 0; i < iters; i++) {
          g_sink = SHAPES[k].fn(x, n);
        }
        if (now_ns() - t0 >= 1e7) {
          break;
        }
        iters *= 2;
      }
      for (int r = 0; r < reps; r++) {
        double t0 = now_ns();
        for (long i = 0; i < iters; i++) {
          g_sink = SHAPES[k].fn(x, n);
        }
        per[r] = (now_ns() - t0) / (double)iters;
      }
      qsort(per, (size_t)reps, sizeof(double), cmp);
      printf("%s/%s\t%zu\t%.1f\t%.1f\t%.4f\n", SHAPES[k].name,
             n > 4096 ? "big" : "l1", n, per[reps / 2], per[0],
             per[reps / 2] / (double)n);
    }
  }
  return 0;
}
