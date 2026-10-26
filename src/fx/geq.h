#ifndef UDSP_GEQ_H
#define UDSP_GEQ_H

#include "fx/eq.h"

/* A graphic EQ's bells: Kilohearts' Carve EQ's, measured through the
 * plugin (ctui-wm's notes/vst-replace.md "Carve EQ, measured"). 31 bands
 * a third of an octave apart, 1000 * 2^((i - 17) / 3) Hz (19.7 Hz ..
 * 20.2 kHz), each a matched bell (eq.h's UDSP_EQ_BELL_MATCHED) at Q
 * 2.1445; the curve drawn as a gain a band, the bells' gains solved so
 * the cascade passes through it at the centres -- the interaction's
 * inverse at a 6 dB prototype, then six Newton steps on the cascade's
 * own magnitudes. Targets +-16 dB (Carve's), bells +-30 (eq.c's). A
 * band whose centre is past 0.49 of the rate is left out. */

#define UDSP_GEQ_BANDS 31
#define UDSP_GEQ_Q 2.1445f
#define UDSP_GEQ_RANGE 16.0f

typedef struct {
  unsigned rate;
  int n; /* the bands in use: those under 0.49 of the rate */
  double hz[UDSP_GEQ_BANDS];
  double phi[UDSP_GEQ_BANDS]; /* sin^2(pi hz / rate) at each centre */
  double inv[UDSP_GEQ_BANDS][UDSP_GEQ_BANDS];
} UDSP_GEQ;

/* band i's centre in Hz */
double udsp_geq_hz(int i);

void udsp_geq_init(UDSP_GEQ *q, unsigned rate);

/* the bells' gains (dB) for a curve (dB at the centres, clamped to the
 * range); the bands left out get 0 */
void udsp_geq_solve(const UDSP_GEQ *q, const float *target_db, float *gain_db);

/* band i as a bell of gain_db */
void udsp_geq_band(const UDSP_GEQ *q, int i, float gain_db, UDSP_EQ_BAND *b);

/* the cascade's gain (dB) at hz for the bells' gains: its curve */
float udsp_geq_response(const UDSP_GEQ *q, const float *gain_db, double hz);

#endif
