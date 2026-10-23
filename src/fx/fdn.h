#ifndef UDSP_FDN_H
#define UDSP_FDN_H

#include <stddef.h>
#include <stdint.h>

#include "fx/line.h"

/* A reverb: a feedback delay network of eight lines (Jot's) mixed by a
 * Hadamard matrix (orthogonal: the loop loses only what the lines'
 * gains take), each line's gain set so it falls 60 dB in the decay time
 * whatever its length, and a one-pole low-pass in each for the highs
 * dying first (damping). Ahead of it a predelay and four Schroeder
 * allpasses (diffusion: the attack smeared into a wash); the lines'
 * reads move slowly around their lengths (modulation: no metallic
 * ring); two orthogonal sums of the lines out, left and right, a width
 * between them. Its memory is the caller's (udsp_fdn_floats). */

#define UDSP_FDN_N 8
#define UDSP_FDN_DIFFUSERS 4
/* the largest size, and the longest predelay */
#define UDSP_FDN_SIZE_MAX 2.0f
#define UDSP_FDN_PREDELAY_MAX_MS 250.0f

typedef struct {
  unsigned rate;
  UDSP_LINE line[UDSP_FDN_N];
  float len[UDSP_FDN_N];      /* frames, now */
  float g[UDSP_FDN_N];        /* per pass */
  float lp[UDSP_FDN_N], damp; /* the low-passes' states, their pole */
  float depth;                /* the reads' swing, frames */
  float rot_c, rot_s, mc, ms; /* the modulation: a rotation a frame */
  uint32_t t;                 /* frames since init */
  UDSP_LINE pre;
  uint32_t predelay;
  UDSP_LINE ap[UDSP_FDN_DIFFUSERS];
  uint32_t ap_len[UDSP_FDN_DIFFUSERS];
  float ap_g;
  float width;
} UDSP_FDN;

/* the floats udsp_fdn_init() wants at rate */
size_t udsp_fdn_floats(unsigned rate);

void udsp_fdn_init(UDSP_FDN *f, float *mem, unsigned rate);
/* silence in every line (the settings kept) */
void udsp_fdn_clear(UDSP_FDN *f);

/* size 0.2 .. UDSP_FDN_SIZE_MAX (1: lines of 30-75 ms), decay seconds
 * (60 dB), damping's corner in Hz, diffusion and modulation 0..1,
 * predelay in ms, width 0 (mono) .. 1 */
void udsp_fdn_set(UDSP_FDN *f, float size, float decay_s, float damp_hz,
                  float diffusion, float modulation, float predelay_ms,
                  float width);

/* frames of in (interleaved stereo) into wet (the reverb alone; may be
 * in) */
void udsp_fdn_stereo(UDSP_FDN *f, const float *in, float *wet, size_t frames);

/* frames gone by unheard (a host's sleep in silence): the modulation
 * turned on exactly as udsp_fdn_stereo() would have, the lines left as
 * they are */
void udsp_fdn_idle(UDSP_FDN *f, size_t frames);

#endif
