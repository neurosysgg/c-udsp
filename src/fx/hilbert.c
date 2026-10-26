#include "fx/hilbert.h"

#include "base/trig.h"

/* Niemitalo's coefficients, squared */
static const float CHAIN[2][4] = {
    {0.6923878f * 0.6923878f, 0.9360654322959f * 0.9360654322959f,
     0.9882295226860f * 0.9882295226860f, 0.9987488452737f * 0.9987488452737f},
    {0.4021921162426f * 0.4021921162426f, 0.8561710882420f * 0.8561710882420f,
     0.9722909545651f * 0.9722909545651f, 0.9952884791278f * 0.9952884791278f}};

static float chain(UDSP_HILBERT_SECTION *c, const float *a, float x) {
  for (int k = 0; k < 4; k++) {
    UDSP_HILBERT_SECTION *s = &c[k];
    const float y = a[k] * (x + s->y2) - s->x2;
    s->x2 = s->x1;
    s->x1 = x;
    s->y2 = s->y1;
    s->y1 = y;
    x = y;
  }
  return x;
}

void udsp_hilbert_step(UDSP_HILBERT *h, float x, float *re, float *im) {
  *re = h->late;
  h->late = chain(h->s[0], CHAIN[0], x);
  *im = chain(h->s[1], CHAIN[1], x);
}

void udsp_hilbert_fir(float *h, size_t taps) {
  const long c = (long)(taps - 1) / 2;
  for (long i = 0; i < (long)taps; i++) {
    const long k = i - c;
    const double s = udsp_sin_pi((double)i / (double)(taps - 1));
    h[i] =
        k % 2 ? (float)(2.0 / (3.141592653589793 * (double)k) * s * s) : 0.0f;
  }
}
