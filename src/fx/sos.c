#include "sos.h"

#include <string.h>

void udsp_sos_init(UDSP_SOS *f, const UDSP_SOS_SECTION *sec, int n) {
  memcpy(f->sec, sec, (size_t)n * sizeof *sec);
  f->n = n;
  udsp_sos_reset(f);
}

void udsp_sos_reset(UDSP_SOS *f) { memset(f->s, 0, sizeof f->s); }

void udsp_sos_stereo(UDSP_SOS *f, float *io, size_t frames) {
  for (int k = 0; k < f->n; k++) {
    const UDSP_SOS_SECTION *q = &f->sec[k];
    for (int c = 0; c < 2; c++) {
      float s1 = f->s[k][c][0], s2 = f->s[k][c][1];
      for (size_t i = 0; i < frames; i++) {
        const float x = io[2 * i + c];
        const float y = q->b0 * x + s1;
        s1 = q->b1 * x - q->a1 * y + s2;
        s2 = q->b2 * x - q->a2 * y;
        io[2 * i + c] = y;
      }
      f->s[k][c][0] = s1;
      f->s[k][c][1] = s2;
    }
  }
}
