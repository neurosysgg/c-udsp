#include "analysis/chroma.h"

#include <math.h>
#include <string.h>

size_t udsp_chroma_peaks(const float *power, int n, int rate, float fmin,
                         float fmax, float floor, UDSP_PEAK *out, size_t cap) {
  const float bin = (float)rate / (float)n;
  int k0 = (int)(fmin / bin);
  int k1 = (int)(fmax / bin) + 1;
  k0 = k0 < 2 ? 2 : k0;
  k1 = k1 > n / 2 - 1 ? n / 2 - 1 : k1;
  size_t got = 0;
  for (int k = k0; k < k1; k++) {
    const float m = power[k];
    if (m < floor || m < power[k - 1] || m <= power[k + 1]) {
      continue;
    }
    const float a = logf(power[k - 1] + 1e-20f), b = logf(m + 1e-20f),
                c = logf(power[k + 1] + 1e-20f);
    const float d = a - 2.0f * b + c;
    const float off = d < 0.0f ? 0.5f * (a - c) / d : 0.0f;
    if (got == cap) {
      /* full: the weakest goes, if this one is stronger */
      size_t low = 0;
      for (size_t i = 1; i < got; i++) {
        low = out[i].power < out[low].power ? i : low;
      }
      if (out[low].power >= m) {
        continue;
      }
      memmove(&out[low], &out[low + 1], (got - low - 1) * sizeof(*out));
      got--;
    }
    out[got].hz = ((float)k + off) * bin;
    out[got++].power = m;
  }
  return got;
}

int udsp_chroma_note(float hz, float *cents) {
  if (hz <= 0.0f) {
    return -1;
  }
  const float x = 69.0f + 12.0f * log2f(hz / 440.0f);
  const float note = roundf(x);
  if (cents) {
    *cents = 100.0f * (x - note);
  }
  return (int)note;
}

void udsp_chroma_fold(const UDSP_PEAK *p, size_t n, float chroma[12]) {
  for (size_t i = 0; i < n; i++) {
    const int note = udsp_chroma_note(p[i].hz, NULL);
    if (note >= 0) {
      chroma[note % 12] += sqrtf(p[i].power);
    }
  }
}

float udsp_chroma_root(const UDSP_PEAK *p, size_t n) {
  enum { HARMONICS = 8 };
  const float tol = 0.03f, reach = 0.7f;
  float scores[64], best = 0.0f;
  n = n > 64 ? 64 : n;
  for (size_t i = 0; i < n; i++) {
    float sum = 0.0f;
    for (int h = 1; h <= HARMONICS; h++) {
      const float want = p[i].hz * (float)h;
      float hit = 0.0f;
      for (size_t j = i; j < n; j++) {
        if (fabsf(p[j].hz - want) < tol * want) {
          const float a = sqrtf(p[j].power);
          hit = a > hit ? a : hit;
        }
      }
      sum += hit / (float)h;
    }
    scores[i] = sum;
    best = sum > best ? sum : best;
  }
  for (size_t i = 0; i < n; i++) { /* by rising frequency */
    if (best > 0.0f && scores[i] >= reach * best) {
      return p[i].hz;
    }
  }
  return 0.0f;
}

static const unsigned short CHORD_MASK[UDSP_CHORDS] = {
    [UDSP_CHORD_MAJOR] = 1u << 0 | 1u << 4 | 1u << 7,
    [UDSP_CHORD_MINOR] = 1u << 0 | 1u << 3 | 1u << 7,
    [UDSP_CHORD_DOM7] = 1u << 0 | 1u << 4 | 1u << 7 | 1u << 10,
    [UDSP_CHORD_MAJ7] = 1u << 0 | 1u << 4 | 1u << 7 | 1u << 11,
    [UDSP_CHORD_MIN7] = 1u << 0 | 1u << 3 | 1u << 7 | 1u << 10,
    [UDSP_CHORD_DIM] = 1u << 0 | 1u << 3 | 1u << 6,
    [UDSP_CHORD_SUS4] = 1u << 0 | 1u << 5 | 1u << 7,
    [UDSP_CHORD_SUS2] = 1u << 0 | 1u << 2 | 1u << 7,
    [UDSP_CHORD_POWER] = 1u << 0 | 1u << 7,
};

const char *udsp_chord_name(int chord) {
  static const char *const NAMES[UDSP_CHORDS] = {
      "", "m", "7", "maj7", "m7", "dim", "sus4", "sus2", "5"};
  return chord >= 0 && chord < UDSP_CHORDS ? NAMES[chord] : "";
}

float udsp_chroma_chord(const float chroma[12], int *root, int *chord) {
  float norm = 0.0f;
  for (int i = 0; i < 12; i++) {
    norm += chroma[i] * chroma[i];
  }
  *root = 0;
  *chord = UDSP_CHORD_MAJOR;
  if (norm <= 0.0f) {
    return 0.0f;
  }
  float best = -1.0f;
  for (int c = 0; c < UDSP_CHORDS; c++) {
    int tones = 0;
    for (int i = 0; i < 12; i++) {
      tones += (CHORD_MASK[c] >> i) & 1;
    }
    for (int r = 0; r < 12; r++) {
      float dot = 0.0f;
      for (int i = 0; i < 12; i++) {
        if ((CHORD_MASK[c] >> i) & 1) {
          dot += chroma[(r + i) % 12];
        }
      }
      /* a dyad fits anything with those two in it: only with no third */
      const float weight = c == UDSP_CHORD_POWER ? 0.85f : 1.0f;
      const float score = weight * dot / sqrtf(norm * (float)tones);
      if (score > best) {
        best = score;
        *root = r;
        *chord = c;
      }
    }
  }
  return best;
}

float udsp_chroma_key(const float chroma[12], int *root, int *minor) {
  static const float PROFILE[2][12] = {
      {6.35f, 2.23f, 3.48f, 2.33f, 4.38f, 4.09f, 2.52f, 5.19f, 2.39f, 3.66f,
       2.29f, 2.88f},
      {6.33f, 2.68f, 3.52f, 5.38f, 2.60f, 3.53f, 2.54f, 4.75f, 3.98f, 2.69f,
       3.34f, 3.17f},
  };
  *root = 0;
  *minor = 0;
  float mean = 0.0f, var = 0.0f;
  for (int i = 0; i < 12; i++) {
    mean += chroma[i] / 12.0f;
  }
  for (int i = 0; i < 12; i++) {
    var += (chroma[i] - mean) * (chroma[i] - mean);
  }
  /* all alike (down to rounding) says nothing */
  if (var <= 1e-6f * mean * mean) {
    return 0.0f;
  }
  float best = -2.0f;
  for (int m = 0; m < 2; m++) {
    float pmean = 0.0f, pvar = 0.0f;
    for (int i = 0; i < 12; i++) {
      pmean += PROFILE[m][i] / 12.0f;
    }
    for (int i = 0; i < 12; i++) {
      pvar += (PROFILE[m][i] - pmean) * (PROFILE[m][i] - pmean);
    }
    for (int r = 0; r < 12; r++) {
      float dot = 0.0f;
      for (int i = 0; i < 12; i++) {
        dot += (chroma[(r + i) % 12] - mean) * (PROFILE[m][i] - pmean);
      }
      const float corr = dot / sqrtf(var * pvar);
      if (corr > best) {
        best = corr;
        *root = r;
        *minor = m;
      }
    }
  }
  return best;
}
