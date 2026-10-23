#include "fx/delay.h"

#include <string.h>

void udsp_delay_init(UDSP_DELAY *d, float *ring, uint32_t cap) {
  d->ring = ring;
  d->cap = cap;
  udsp_delay_reset(d);
}

void udsp_delay_reset(UDSP_DELAY *d) {
  memset(d->ring, 0, UDSP_DELAY_RING_FLOATS(d->cap) * sizeof(float));
  d->write = 0;
}

void udsp_delay_stereo(UDSP_DELAY *d, float *io, size_t frames, uint32_t len,
                       float feedback, float mix) {
  if (d->cap == 0) {
    return;
  }
  len = len < 1 ? 1 : (len > d->cap ? d->cap : len);
  uint32_t w = d->write, r = (w + d->cap - len) % d->cap;
  for (size_t i = 0; i < frames; i++) {
    for (int c = 0; c < 2; c++) {
      const float dry = io[i * 2 + c];
      const float wet = d->ring[r * 2 + c];
      d->ring[w * 2 + c] = dry + wet * feedback;
      io[i * 2 + c] = dry + wet * mix;
    }
    w = w + 1 == d->cap ? 0 : w + 1;
    r = r + 1 == d->cap ? 0 : r + 1;
  }
  d->write = w;
}
