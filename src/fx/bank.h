#ifndef UDSP_BANK_H
#define UDSP_BANK_H

#include <stddef.h>
#include <stdint.h>

/* A bank of band-pass filters side by side, and their envelopes: a
 * vocoder's analysis and its carrier's bands. Each band is 1-4 Simper
 * state-variable band-passes in a row (zdf.h's), each normalized to
 * unity at its centre; the bands are kept as arrays across the bank, so
 * a frame runs UDSP_W bands at a time (100 bands of order 2, a frame:
 * ~25 vector steps). Stateful as zdf.h is: the caller owns the struct. */

#define UDSP_BANK_MAX 104 /* bands, padded to a whole number of vectors */
#define UDSP_BANK_ORDER 4

typedef struct {
  uint32_t n, lanes, order; /* bands; n rounded up to vectors; stages */
  float a1[UDSP_BANK_MAX], a2[UDSP_BANK_MAX], a3[UDSP_BANK_MAX];
  float k[UDSP_BANK_MAX];    /* a stage's 1 / Q: its centre's unity gain */
  float gain[UDSP_BANK_MAX]; /* the band's, on its last stage */
  float ic1[UDSP_BANK_ORDER][UDSP_BANK_MAX],
      ic2[UDSP_BANK_ORDER][UDSP_BANK_MAX];
} UDSP_BANK;

/* n bands (1 .. UDSP_BANK_MAX) centred on hz[], each q[] its whole
 * band's Q (centre / its -3 dB width: the stages' own Q widened so
 * their product keeps that width) and gain[] (NULL: 1), order stages
 * (1 .. 4); centres clamped to 10 Hz .. 0.49 of the rate. The state is
 * kept where n and order don't change (a retune glides), else zeroed. */
void udsp_bank_set(UDSP_BANK *b, uint32_t n, const float *hz, const float *q,
                   const float *gain, uint32_t order, unsigned rate);

/* a band's equivalent noise bandwidth in Hz (centre hz, Q q, order
 * stages; the narrow-band form, |H|^2 ~ 1 / (1 + x^2) a stage): what
 * white noise of variance v gives v times 2 of over the rate */
float udsp_bank_noise_width(float hz, float q, uint32_t order);

/* band i's response at hz, complex (its gain and every stage): what a
 * caller sums over the bands to know how they overlap */
void udsp_bank_response(const UDSP_BANK *b, uint32_t i, float hz, unsigned rate,
                        float *re, float *im);

/* zeroes the state */
void udsp_bank_reset(UDSP_BANK *b);

/* frames of x (every stride-th float) through every band: out[f * lanes
 * + band], the lanes past n zero */
void udsp_bank_run(UDSP_BANK *b, const float *x, size_t stride, float *out,
                   size_t frames);

/* The bands' envelopes: each band's |x| through a one-pole that rises
 * by attack and falls by release, held hold frames after its last rise
 * before it falls. A sine's envelope is its amplitude with a fast attack
 * and a slow release, 2 / pi of it (its mean) with the two alike. */
typedef struct {
  uint32_t lanes;
  float attack, release; /* the one-pole's step: 1 - e^(-1 / frames) */
  float hold;            /* frames */
  float env[UDSP_BANK_MAX], held[UDSP_BANK_MAX];
} UDSP_FOLLOW;

/* times in frames (0: at once); env, held kept */
void udsp_follow_set(UDSP_FOLLOW *f, uint32_t lanes, float attack,
                     float release, float hold);

/* io as udsp_bank_run's out, its bands replaced by their envelopes */
void udsp_follow_run(UDSP_FOLLOW *f, float *io, size_t frames);

#endif
