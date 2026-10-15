#ifndef UDSP_SIMD_H
#define UDSP_SIMD_H

/* internal: the vector type every hot loop here is written against
 * (bench/README.md "Reductions" has the measurements behind it).
 *
 * GCC vector extensions (gcc and clang) at the target's native width: 16
 * bytes on baseline x86-64/SSE2 or NEON, 32 with AVX. Explicit vectors
 * because both compilers vectorise these at -O2 already and the same way
 * every time, where auto-vectorisation of the same loops comes and goes
 * with the optimisation level. Never wider than the target: a 32-byte
 * vector on an SSE2 target is split and, for anything but + and *,
 * scalarised (6x slower than a plain loop in the study).
 *
 * UDSP_SCALAR (or a compiler without the extension) makes a vector one
 * float: the same kernels as plain C, the reference the tests hold the
 * vector build to. */

#include <math.h>
#include <string.h>

#if defined(__GNUC__) && !defined(UDSP_SCALAR)
#if defined(__AVX__)
#define UDSP_W 8
#else
#define UDSP_W 4
#endif
typedef float UDSP_VF __attribute__((vector_size(UDSP_W * 4)));
typedef int UDSP_VI __attribute__((vector_size(UDSP_W * 4)));
/* phase accumulators (osc.c): uint32 lanes that wrap like the scalar
 * one, converted with __builtin_convertvector (correctly rounded, as a
 * scalar (float) cast is) */
typedef unsigned UDSP_VU __attribute__((vector_size(UDSP_W * 4)));
/* a where m is set (all ones), else b: the blend vector ?: would be */
static inline UDSP_VF udsp_vsel(UDSP_VI m, UDSP_VF a, UDSP_VF b) {
  return (UDSP_VF)((m & (UDSP_VI)a) | (~m & (UDSP_VI)b));
}
/* four lanes whatever the target: for loops shorter than UDSP_W (an FFT's
 * narrow stages under AVX) */
#define UDSP_HAVE_V4 1
typedef float UDSP_V4 __attribute__((vector_size(16)));

static inline UDSP_V4 udsp_v4load(const float *p) {
  UDSP_V4 v;
  memcpy(&v, p, sizeof(v));
  return v;
}
static inline void udsp_v4store(float *p, UDSP_V4 v) {
  memcpy(p, &v, sizeof(v));
}
/* lanes picked from a (0-3) and b (4-7): __builtin_shufflevector (clang,
 * gcc 12+), else gcc's own __builtin_shuffle */
#if defined(__clang__) || __GNUC__ >= 12
#define UDSP_SHUF4(a, b, i, j, k, l) __builtin_shufflevector(a, b, i, j, k, l)
#else
typedef int UDSP_V4I __attribute__((vector_size(16)));
#define UDSP_SHUF4(a, b, i, j, k, l)                                          \
  __builtin_shuffle(a, b, (UDSP_V4I){i, j, k, l})
#endif

static inline UDSP_VF udsp_vload(const float *p) {
  UDSP_VF v;
  memcpy(&v, p, sizeof(v));
  return v;
}
static inline void udsp_vstore(float *p, UDSP_VF v) {
  memcpy(p, &v, sizeof(v));
}
static inline UDSP_VF udsp_vsplat(float x) {
  return (UDSP_VF){0} + x;
}
/* a > b ? a : b per lane, and |a|: the target's own instructions where it
 * has them. The generic form (a compare and a blend: vector ?: is C++-only
 * in gcc) is right everywhere, but gcc doesn't fold it back into maxps and
 * runs it in the integer domain: 3x slower in the peak kernel. maxps is
 * exactly a > b ? a : b, so every build gives the scalar build's answer */
#if defined(__AVX__) || defined(__SSE__)
#include <immintrin.h>
#endif
static inline UDSP_VF udsp_vmax(UDSP_VF a, UDSP_VF b) {
#if defined(__AVX__)
  return (UDSP_VF)_mm256_max_ps((__m256)a, (__m256)b);
#elif defined(__SSE__)
  return (UDSP_VF)_mm_max_ps((__m128)a, (__m128)b);
#else
  UDSP_VI m = a > b;
  return (UDSP_VF)((m & (UDSP_VI)a) | (~m & (UDSP_VI)b));
#endif
}
static inline UDSP_VF udsp_vmin(UDSP_VF a, UDSP_VF b) {
#if defined(__AVX__)
  return (UDSP_VF)_mm256_min_ps((__m256)a, (__m256)b);
#elif defined(__SSE__)
  return (UDSP_VF)_mm_min_ps((__m128)a, (__m128)b);
#else
  UDSP_VI m = a < b;
  return (UDSP_VF)((m & (UDSP_VI)a) | (~m & (UDSP_VI)b));
#endif
}
static inline UDSP_VF udsp_vabs(UDSP_VF a) {
#if defined(__AVX__)
  return (UDSP_VF)_mm256_andnot_ps(_mm256_set1_ps(-0.0f), (__m256)a);
#elif defined(__SSE__)
  return (UDSP_VF)_mm_andnot_ps(_mm_set1_ps(-0.0f), (__m128)a);
#else
  return (UDSP_VF)((UDSP_VI)a & ((UDSP_VI){0} + 0x7fffffff));
#endif
}
/* sqrt per lane: sqrtps (a scalar sqrtf() loop doesn't vectorise without
 * -fno-math-errno: errno for negative input) */
static inline UDSP_VF udsp_vsqrt(UDSP_VF v) {
#if defined(__AVX__)
  return (UDSP_VF)_mm256_sqrt_ps((__m256)v);
#elif defined(__SSE__)
  return (UDSP_VF)_mm_sqrt_ps((__m128)v);
#else
  for (int j = 0; j < UDSP_W; j++) {
    v[j] = sqrtf(v[j]);
  }
  return v;
#endif
}

/* the lanes in reverse order */
static inline UDSP_VF udsp_vrev(UDSP_VF v) {
#if defined(__AVX__)
  __m256 r = _mm256_permute_ps((__m256)v, 0x1b); /* 3 2 1 0 in each half */
  return (UDSP_VF)_mm256_permute2f128_ps(r, r, 0x01);
#elif defined(__SSE__)
  return (UDSP_VF)_mm_shuffle_ps((__m128)v, (__m128)v, 0x1b);
#else
  UDSP_VF r;
  for (int j = 0; j < UDSP_W; j++) {
    r[j] = v[UDSP_W - 1 - j];
  }
  return r;
#endif
}

static inline float udsp_vsum(UDSP_VF v) {
  float s = 0.0f;
  for (int j = 0; j < UDSP_W; j++) {
    s += v[j];
  }
  return s;
}
static inline float udsp_vhmax(UDSP_VF v) {
  float m = v[0];
  for (int j = 1; j < UDSP_W; j++) {
    m = v[j] > m ? v[j] : m;
  }
  return m;
}
static inline float udsp_vhmin(UDSP_VF v) {
  float m = v[0];
  for (int j = 1; j < UDSP_W; j++) {
    m = v[j] < m ? v[j] : m;
  }
  return m;
}
#else
#define UDSP_W 1
typedef float UDSP_VF;

static inline UDSP_VF udsp_vload(const float *p) { return *p; }
static inline void udsp_vstore(float *p, UDSP_VF v) { *p = v; }
static inline UDSP_VF udsp_vsplat(float x) { return x; }
static inline UDSP_VF udsp_vmax(UDSP_VF a, UDSP_VF b) { return a > b ? a : b; }
static inline UDSP_VF udsp_vmin(UDSP_VF a, UDSP_VF b) { return a < b ? a : b; }
static inline UDSP_VF udsp_vabs(UDSP_VF a) { return fabsf(a); }
static inline UDSP_VF udsp_vrev(UDSP_VF v) { return v; }
static inline UDSP_VF udsp_vsqrt(UDSP_VF v) { return sqrtf(v); }
static inline float udsp_vsum(UDSP_VF v) { return v; }
static inline float udsp_vhmax(UDSP_VF v) { return v; }
static inline float udsp_vhmin(UDSP_VF v) { return v; }
#endif

/* independent accumulators a reduction keeps: FP add/FMA latency (3-4
 * cycles on Zen 3) x the two pipes that can take one each cycle. Fewer
 * leaves the loop waiting on its own chain -- the one thing unrolling
 * (-funroll-loops) doesn't fix, since without -ffast-math the compiler
 * may not split a float sum itself. Overridable (-DUDSP_ACC=N) for
 * bench/sweep.py, which measures the choice per kernel and build. */
#ifndef UDSP_ACC
#define UDSP_ACC 8
#endif

/* on each loop over the accumulators: gcc's -O2 doesn't fully unroll
 * them by itself, and then the array stays in memory instead of
 * registers (a store-to-load trip per step, 2-3x slower). Written as a
 * _Pragma so a macro can carry it. UDSP_NO_UNROLL drops it (bench/
 * sweep.py's knob; never for a real build) */
#ifdef UDSP_NO_UNROLL
#define UDSP_UNROLL
#else
#define UDSP_UNROLL _Pragma("GCC unroll 16")
#endif

#endif
