#include "analysis/ballistics.h"

#include <math.h>

static float coeff(float dt, float tau) {
  return tau <= 0.0f ? 1.0f : 1.0f - expf(-dt / tau);
}

void udsp_ballistics(float *level, const float *target, size_t n, float dt,
                     float attack_s, float release_s) {
  float up = coeff(dt, attack_s), down = coeff(dt, release_s);
  for (size_t i = 0; i < n; i++) {
    float k = target[i] > level[i] ? up : down;
    level[i] += (target[i] - level[i]) * k;
  }
}

void udsp_peak_hold(float *peak, float *age, const float *level, size_t n,
                    float dt, float hold_s, float fall_per_s) {
  for (size_t i = 0; i < n; i++) {
    if (level[i] >= peak[i]) {
      peak[i] = level[i];
      age[i] = 0.0f;
      continue;
    }
    /* the part of this step past the hold time is spent falling */
    float falling = age[i] + dt - (age[i] > hold_s ? age[i] : hold_s);
    age[i] += dt;
    if (falling > 0.0f) {
      peak[i] -= fall_per_s * falling;
      if (peak[i] < level[i]) {
        peak[i] = level[i];
      }
    }
  }
}
