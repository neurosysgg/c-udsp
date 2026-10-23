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

void udsp_spectrogram_push_mag(UDSP_SPECTROGRAM *sg, const float *mag,
                               int sample_rate) {
  float bars_dbfs[UDSP_SPECTROGRAM_BINS];
  udsp_spectrum_bin(mag, sample_rate, bars_dbfs, UDSP_SPECTROGRAM_BINS);

  /* scroll every column one step toward index 0, then fill the freed last
   * column -- a plain memmove over the whole grid; simple and correct,
   * not the fastest possible (a true ring index would avoid the copy),
   * a deliberate "simple first" tradeoff for a payload this small
   * (COLS*BINS bytes). Through the grid's own bytes, not a row pointer:
   * a fortified memmove sizes &intensity[0] (or [0][0]) as one row. */
  unsigned char *grid = (unsigned char *)sg;
  memmove(grid, grid + sizeof(sg->intensity[0]),
          sizeof(sg->intensity) - sizeof(sg->intensity[0]));

  for (int bin = 0; bin < UDSP_SPECTROGRAM_BINS; bin++) {
    sg->intensity[UDSP_SPECTROGRAM_COLS - 1][bin] =
        quantize_dbfs(bars_dbfs[bin]);
  }
}

void udsp_spectrogram_push(UDSP_SPECTROGRAM *sg, const float *mono_tail,
                           int sample_rate, float *scratch_re,
                           float *scratch_im) {
  udsp_spectrum_magnitude(mono_tail, scratch_re, scratch_im, scratch_re);
  udsp_spectrogram_push_mag(sg, scratch_re, sample_rate);
}
