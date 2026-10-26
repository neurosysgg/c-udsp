#include "base/trig.h"

#include <stdint.h>

#define PI 3.14159265358979323846

double udsp_sin_pi(double y) {
  int64_t k = (int64_t)(y * 0.5);
  y -= 2.0 * (double)k; /* (-2, 2) */
  if (y > 1.0) {
    y -= 2.0;
  } else if (y < -1.0) {
    y += 2.0;
  }
  if (y > 0.5) {
    y = 1.0 - y;
  } else if (y < -0.5) {
    y = -1.0 - y;
  }
  double x = PI * y, x2 = x * x, term = x, s = x;
  for (int i = 1; i <= 12; i++) {
    term *= -x2 / (double)((2 * i) * (2 * i + 1));
    s += term;
  }
  return s;
}
