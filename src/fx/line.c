#include "fx/line.h"

#include <string.h>

#include "base/simd.h"

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

void udsp_line_reads(const UDSP_LINE *l, const float *d, float *out, size_t n) {
  size_t j = 0;
#if UDSP_W > 1
  UDSP_VI lane;
  for (int k = 0; k < UDSP_W; k++) {
    lane[k] = k;
  }
  const UDSP_VI mask = (UDSP_VI){0} + (int)l->mask;
  for (; j + UDSP_W <= n; j += UDSP_W) {
    const UDSP_VF dv = udsp_vload(d + j);
    const UDSP_VI i = __builtin_convertvector(dv, UDSP_VI);
    const UDSP_VF f = dv - __builtin_convertvector(i, UDSP_VF);
    /* frame j's write position, i frames back: its tap i */
    const UDSP_VI at = (UDSP_VI){0} + (int)(l->w + (uint32_t)j) + lane - i;
    UDSP_VF y0, y1, y2;
    const UDSP_VF ym1 =
        udsp_vtaps(l->buf, (UDSP_VI){0}, at & mask, mask, &y0, &y1, &y2);
    udsp_vstore(out + j, udsp_vhermite(ym1, y0, y1, y2, f));
  }
#endif
  for (; j < n; j++) {
    const UDSP_LINE at = {l->buf, l->mask, (l->w + (uint32_t)j) & l->mask};
    out[j] = udsp_line_read(&at, d[j]);
  }
}

void udsp_line_write(UDSP_LINE *l, const float *x, size_t n) {
  for (size_t j = 0; j < n;) {
    size_t k = l->mask + 1 - l->w;
    k = k < n - j ? k : n - j;
    memcpy(l->buf + l->w, x + j, k * sizeof(float));
    l->w = (l->w + (uint32_t)k) & l->mask;
    j += k;
  }
}
