#ifndef UDSP_CHROMA_H
#define UDSP_CHROMA_H

#include <stddef.h>

/* What's being played, read off a power spectrum (udsp_spectrum_power()):
 * its peaks, how they fold onto the twelve pitch classes, the root under
 * them, the chord and the key they fit. Every function stateless: a
 * caller wanting a steadier answer averages the chroma it passes in.
 *
 * Peaks and not bins: a Hann window spreads one partial over four bins,
 * into the neighbouring semitones of everything under a few hundred Hz.
 * Only a local maximum counts, its frequency refined by a parabola
 * through the log power. */

typedef struct {
  float hz;
  float power; /* as the spectrum's: a full-scale sine near 1 */
} UDSP_PEAK;

/* the peaks of the bins 0..n/2 of an n-point spectrum between fmin and
 * fmax over floor (power), the strongest cap of them, by rising frequency;
 * how many */
size_t udsp_chroma_peaks(const float *power, int n, int rate, float fmin,
                         float fmax, float floor, UDSP_PEAK *out, size_t cap);

/* the MIDI note (69 = 440 Hz) nearest hz, its pitch class its % 12; cents
 * (-50..50) unless NULL */
int udsp_chroma_note(float hz, float *cents);

/* peaks folded into chroma[12] (0 C .. 11 B), their amplitudes added in;
 * chroma is added to, not cleared */
void udsp_chroma_fold(const UDSP_PEAK *p, size_t n, float chroma[12]);

/* the fundamental the peaks best stand on, in Hz, 0 with none: each peak
 * scored by the harmonics of it that are there (1/h of each), the lowest
 * within reach of the best -- so a bass note under its overtones, not its
 * loudest partial */
float udsp_chroma_root(const UDSP_PEAK *p, size_t n);

typedef enum {
  UDSP_CHORD_MAJOR,
  UDSP_CHORD_MINOR,
  UDSP_CHORD_DOM7,
  UDSP_CHORD_MAJ7,
  UDSP_CHORD_MIN7,
  UDSP_CHORD_DIM,
  UDSP_CHORD_SUS4,
  UDSP_CHORD_SUS2,
  UDSP_CHORD_POWER,
  UDSP_CHORDS,
} UDSP_CHORD;

/* its suffix: "", "m", "7", "maj7", "m7", "dim", "sus4", "sus2", "5" */
const char *udsp_chord_name(int chord);

/* the chord template chroma fits best (cosine similarity, the score,
 * 0..1), *root its pitch class, *chord its UDSP_CHORD; 0 for silence */
float udsp_chroma_chord(const float chroma[12], int *root, int *chord);

/* the key chroma fits best: Krumhansl and Kessler's profiles correlated
 * at all twelve roots, *root its pitch class, *minor 1 for minor; the
 * correlation (-1..1), 0 for silence */
float udsp_chroma_key(const float chroma[12], int *root, int *minor);

#endif
