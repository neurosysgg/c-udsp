#include "os.h"

#include <string.h>

/* hiir's PolyphaseIir2Designer: 8 coefficients at a transition of 0.0313
 * of the doubled rate, 4 at 0.14 */
static const float COEF[UDSP_OS_STAGES][UDSP_OS_COEFS] = {
    {0.0462761761f, 0.169424325f, 0.332916371f, 0.500246586f, 0.648729628f,
     0.771148382f, 0.871007787f, 0.957684798f},
    {0.0636761356f, 0.236978866f, 0.486229081f, 0.802361914f},
};
static const int N_COEF[UDSP_OS_STAGES] = {8, 4};

void udsp_os_reset(UDSP_OS *o) { memset(o, 0, sizeof *o); }

/* the allpasses of one path (every second, from path) over v */
static float path(UDSP_OS_CHAIN *ch, int stage, int p, int c, float v) {
  const float *k = COEF[stage];
  for (int i = p; i < N_COEF[stage]; i += 2) {
    const float t = (v - ch->y[i][c]) * k[i] + ch->x[i][c];
    ch->x[i][c] = v;
    ch->y[i][c] = t;
    v = t;
  }
  return v;
}

/* one frame doubled at stage s into out[0..3] */
static void up1(UDSP_OS *o, int s, const float *in, float *out) {
  for (int c = 0; c < 2; c++) {
    out[c] = path(&o->up[s], s, 0, c, in[c]);
    out[2 + c] = path(&o->up[s], s, 1, c, in[c]);
  }
}

/* two frames (in[0..3]) halved at stage s */
static void down1(UDSP_OS *o, int s, const float *in, float *out) {
  for (int c = 0; c < 2; c++) {
    const float a = path(&o->down[s], s, 0, c, in[c]);
    const float b = path(&o->down[s], s, 1, c, o->odd[s][c]);
    o->odd[s][c] = in[2 + c];
    out[c] = 0.5f * (a + b);
  }
}

void udsp_os_up(UDSP_OS *o, const float *in, float *out, size_t frames,
                unsigned factor) {
  for (size_t i = 0; i < frames; i++) {
    if (factor == 2) {
      up1(o, 0, in + 2 * i, out + 4 * i);
    } else {
      float two[4];
      up1(o, 0, in + 2 * i, two);
      up1(o, 1, two, out + 8 * i);
      up1(o, 1, two + 2, out + 8 * i + 4);
    }
  }
}

void udsp_os_down(UDSP_OS *o, const float *in, float *out, size_t frames,
                  unsigned factor) {
  for (size_t i = 0; i < frames; i++) {
    if (factor == 2) {
      down1(o, 0, in + 4 * i, out + 2 * i);
    } else {
      float two[4];
      down1(o, 1, in + 8 * i, two);
      down1(o, 1, in + 8 * i + 4, two + 2);
      down1(o, 0, two, out + 2 * i);
    }
  }
}
