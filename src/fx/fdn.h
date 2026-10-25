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
 * between them. Its memory is the caller's (udsp_fdn_floats). The
 * damping can be a first-order high shelf instead (udsp_fdn_set_shelf:
 * the highs losing so many dB a second whatever the decay, Kilohearts'
 * Reverb's). */

#define UDSP_FDN_N 8
#define UDSP_FDN_DIFFUSERS 4
/* the largest size, and the longest predelay */
#define UDSP_FDN_SIZE_MAX 2.0f
#define UDSP_FDN_PREDELAY_MAX_MS 250.0f
/* the early taps': how many, the latest at size 1 */
#define UDSP_FDN_TAPS 16
#define UDSP_FDN_TAP_MAX_MS 125.0f

/* an early reflection: so many ms back at size 1 (scaled by the size),
 * from an input channel to an output's, at a gain */
typedef struct {
  float ms;
  int in, out; /* 0 left, 1 right */
  float gain;
} UDSP_FDN_TAP;

typedef struct {
  unsigned rate;
  UDSP_LINE line[UDSP_FDN_N];
  float len[UDSP_FDN_N];      /* frames, now */
  float g[UDSP_FDN_N];        /* per pass */
  float lp[UDSP_FDN_N], damp; /* the low-passes' states, their pole */
  float hf[UDSP_FDN_N];       /* the shelf's highs kept a pass (0: none) */
  float lp2[UDSP_FDN_N], damp2, hf2[UDSP_FDN_N]; /* a second shelf's */
  float depth;                                   /* the reads' swing, frames */
  float rot_c, rot_s, mc, ms; /* the modulation: a rotation a frame */
  uint32_t t;                 /* frames since init */
  UDSP_LINE pre;
  uint32_t predelay;
  UDSP_LINE ap[UDSP_FDN_DIFFUSERS];
  uint32_t ap_len[UDSP_FDN_DIFFUSERS];
  float ap_g;
  float width;
  UDSP_LINE in[2]; /* the input's, for the taps */
  UDSP_FDN_TAP tap[UDSP_FDN_TAPS];
  uint32_t tap_at[UDSP_FDN_TAPS]; /* frames back, now */
  int taps;
  float early, keep, feed; /* the taps' share; what's left of the input
                            * then, the taps into the lines */
  float size;
  float diffuse;                /* the tail's gain */
  int shelf;                    /* the damping a shelf, not the low-pass */
  float shelf_db_s;             /* its loss a second */
  float shelf2_hz, shelf2_db_s; /* the second's (0 dB: none) */
} UDSP_FDN;

/* the floats udsp_fdn_init() wants at rate */
size_t udsp_fdn_floats(unsigned rate);

void udsp_fdn_init(UDSP_FDN *f, float *mem, unsigned rate);
/* silence in every line (the settings kept) */
void udsp_fdn_clear(UDSP_FDN *f);

/* size 0.2 .. UDSP_FDN_SIZE_MAX (1: lines of 30-75 ms), decay seconds
 * (60 dB), damping's corner in Hz, diffusion and modulation 0..1,
 * predelay in ms, width 0 (mono) .. 1 (.. 2: the side over the mid) */
void udsp_fdn_set(UDSP_FDN *f, float size, float decay_s, float damp_hz,
                  float diffusion, float modulation, float predelay_ms,
                  float width);

/* on: the damping a first-order high shelf at damp_hz (udsp_fdn_set's),
 * the highs losing db_s dB a second over the decay (0: none), then a
 * second at hz2 losing db2_s; off the one-pole low-pass. Before
 * udsp_fdn_set (which keeps it): the lines' lengths set their share a
 * pass, and on they're whole frames (a fraction's cubic dulls the highs a
 * little every pass) */
void udsp_fdn_set_shelf(UDSP_FDN *f, int on, float db_s, float hz2,
                        float db2_s);

/* early reflections (Kilohearts' Reverb's): n taps (up to
 * UDSP_FDN_TAPS, UDSP_FDN_TAP_MAX_MS), crossfaded in by early (0..1):
 * the wet gets early times the taps (before the width), the lines the
 * input at 1 - early + early keep and early feed times the taps' mono
 * sum along a row of their own (a tail apart from the input's). n 0:
 * none, the bits as before. The size scales the taps' times */
void udsp_fdn_set_early(UDSP_FDN *f, const UDSP_FDN_TAP *taps, int n,
                        float early, float keep, float feed);

/* the tail's gain (the taps' as they are); 1 at init */
void udsp_fdn_set_diffuse(UDSP_FDN *f, float gain);

/* frames of in (interleaved stereo) into wet (the reverb alone; may be
 * in) */
void udsp_fdn_stereo(UDSP_FDN *f, const float *in, float *wet, size_t frames);

/* frames gone by unheard (a host's sleep in silence): the modulation
 * turned on exactly as udsp_fdn_stereo() would have, the lines left as
 * they are */
void udsp_fdn_idle(UDSP_FDN *f, size_t frames);

#endif
