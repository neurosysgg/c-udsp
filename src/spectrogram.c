#include "spectrogram.h"

#include "analyze.h"
#include "spectrum.h"

#include <string.h>

static unsigned char quantize_dbfs(float dbfs) {
  float t = (dbfs - UDSP_DBFS_FLOOR) / (0.0f - UDSP_DBFS_FLOOR);
  t = t < 0.0f ? 0.0f : (t > 1.0f ? 1.0f : t);
  return (unsigned char)(t * 255.0f + 0.5f);
}

void udsp_spectrogram_init(UDSP_SPECTROGRAM *sg) {
  memset(sg->intensity, 0, sizeof(sg->intensity));
}

void udsp_spectrogram_push(UDSP_SPECTROGRAM *sg, const float *mono_tail,
                           int sample_rate, float *scratch_re, float *scratch_im) {
  float bars_dbfs[UDSP_SPECTROGRAM_BINS];
  udsp_spectrum_bars(mono_tail, sample_rate, scratch_re, scratch_im, bars_dbfs,
                     UDSP_SPECTROGRAM_BINS);

  /* scroll every column one step toward index 0, then fill the freed last
   * column -- a plain memmove over the whole grid; simple and correct,
   * not the fastest possible (a true ring index would avoid the copy),
   * a deliberate "simple first" tradeoff for a payload this small
   * (COLS*BINS bytes). */
  memmove(&sg->intensity[0], &sg->intensity[1],
         sizeof(sg->intensity) - sizeof(sg->intensity[0]));

  for (int bin = 0; bin < UDSP_SPECTROGRAM_BINS; bin++) {
    sg->intensity[UDSP_SPECTROGRAM_COLS - 1][bin] = quantize_dbfs(bars_dbfs[bin]);
  }
}
