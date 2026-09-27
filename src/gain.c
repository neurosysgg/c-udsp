#include "gain.h"

void udsp_gain_apply(float *buf, size_t n, float gain) {
  for (size_t i = 0; i < n; i++) {
    buf[i] *= gain;
  }
}
