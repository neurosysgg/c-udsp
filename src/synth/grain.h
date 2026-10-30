#ifndef UDSP_GRAIN_H
#define UDSP_GRAIN_H

#include <stddef.h>
#include <stdint.h>

/* A grain: a short stretch of a sample read at a pitch, forwards or
 * backwards, under a window -- what a granular instrument starts by the
 * dozen (its scheduling, positions and sprays are the caller's). The
 * sample is whole in memory and stays the caller's; a frame outside it
 * reads as silence, so a grain may run off either end.
 *
 * Fixed point as the sampler's voice (Q32, a step per output frame,
 * negative backwards) and the window a table read by a 32-bit phase, so
 * a grain is the same bits however its frames are cut into calls. The
 * read is Catmull-Rom, as udsp_resample_cubic's. */

#define UDSP_GRAIN_TABLE 1024

typedef enum {
  UDSP_GRAIN_HANN,
  UDSP_GRAIN_TRIANGLE,
  UDSP_GRAIN_TUKEY, /* flat, with raised cosines over its outer quarters */
  UDSP_GRAIN_PERC,  /* a fast rise, an exponential fall: 0 at both ends */
  UDSP_GRAIN_SWELL, /* the perc backwards */
  UDSP_GRAIN_SHAPES,
} UDSP_GRAIN_SHAPE;

typedef struct {
  int64_t pos;   /* Q32, in the sample's frames */
  int64_t step;  /* Q32 source frames per frame, negative backwards */
  uint32_t left; /* frames to go; 0: done */
  uint32_t len;  /* its whole length */
  uint32_t age;  /* frames played */
  float gl, gr;  /* its level on each side */
} UDSP_GRAIN;

/* shape's window into table (UDSP_GRAIN_TABLE + 1 values, 0 at both
 * ends, at most 1), from udsp_sin_pi: the same in every libc */
void udsp_grain_window(float *table, UDSP_GRAIN_SHAPE shape);

/* g over len >= 1 frames from Q32 position pos, stepping by step */
void udsp_grain_start(UDSP_GRAIN *g, int64_t pos, int64_t step, uint32_t len,
                      float gl, float gr);

/* adds up to frames STEREO frames of g into out, from n frames of ch
 * (>= 1) interleaved channels -- the first two; mono on both sides --
 * under table; returns the frames it played (fewer: it's done) */
uint32_t udsp_grain_render(UDSP_GRAIN *g, const float *samples, uint32_t ch,
                           uint64_t n, const float *table, float *out,
                           uint32_t frames);

#endif
