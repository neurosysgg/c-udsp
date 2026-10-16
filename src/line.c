#include "line.h"

#include <string.h>

uint32_t udsp_line_size(uint32_t max) {
  uint32_t n = 8;
  while (n < max + 4) {
    n <<= 1;
  }
  return n;
}

void udsp_line_init(UDSP_LINE *l, float *buf, uint32_t size) {
  l->buf = buf;
  l->mask = size - 1;
  udsp_line_clear(l);
}

void udsp_line_clear(UDSP_LINE *l) {
  memset(l->buf, 0, (size_t)(l->mask + 1) * sizeof(float));
  l->w = 0;
}
