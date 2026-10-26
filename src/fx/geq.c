#include "fx/geq.h"

#include "base/trig.h"
#include "fx/dyn.h"
#include "synth/pitch.h"

#define PROTO_DB 6.0
#define STEPS 6

double udsp_geq_hz(int i) { return 1000.0 * udsp_exp2((i - 17) / 3.0); }

/* band j's bell at gain_db: its power at centre i, in dB */
typedef struct {
  double b[3], a[3];
} BIQUAD;

static void bell(const UDSP_GEQ *q, int j, double gain_db, BIQUAD *f) {
  udsp_eq_bell_matched(q->hz[j], UDSP_GEQ_Q, gain_db, q->rate, f->b, f->a);
}

static double power_db(const BIQUAD *f, double phi) {
  return 0.5 * (double)udsp_db_from_gain(
                   (float)udsp_eq_biquad_power(f->b, f->a, phi));
}

/* a's inverse into inv by Gauss-Jordan with partial pivoting, n x n */
static void invert(double a[UDSP_GEQ_BANDS][UDSP_GEQ_BANDS],
                   double inv[UDSP_GEQ_BANDS][UDSP_GEQ_BANDS], int n) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      inv[i][j] = i == j;
    }
  }
  for (int c = 0; c < n; c++) {
    int p = c;
    for (int r = c + 1; r < n; r++) {
      if ((a[r][c] < 0 ? -a[r][c] : a[r][c]) >
          (a[p][c] < 0 ? -a[p][c] : a[p][c])) {
        p = r;
      }
    }
    for (int j = 0; j < n; j++) {
      double t = a[c][j];
      a[c][j] = a[p][j];
      a[p][j] = t;
      t = inv[c][j];
      inv[c][j] = inv[p][j];
      inv[p][j] = t;
    }
    const double d = 1.0 / a[c][c];
    for (int j = 0; j < n; j++) {
      a[c][j] *= d;
      inv[c][j] *= d;
    }
    for (int r = 0; r < n; r++) {
      if (r != c && a[r][c] != 0.0) {
        const double m = a[r][c];
        for (int j = 0; j < n; j++) {
          a[r][j] -= m * a[c][j];
          inv[r][j] -= m * inv[c][j];
        }
      }
    }
  }
}

void udsp_geq_init(UDSP_GEQ *q, unsigned rate) {
  q->rate = rate;
  q->n = 0;
  for (int i = 0; i < UDSP_GEQ_BANDS; i++) {
    q->hz[i] = udsp_geq_hz(i);
    const double s = udsp_sin_pi(q->hz[i] / rate);
    q->phi[i] = s * s;
    q->n += q->hz[i] < 0.49 * rate;
  }
  double b[UDSP_GEQ_BANDS][UDSP_GEQ_BANDS];
  for (int j = 0; j < q->n; j++) {
    BIQUAD f;
    bell(q, j, PROTO_DB, &f);
    for (int i = 0; i < q->n; i++) {
      b[i][j] = power_db(&f, q->phi[i]) / PROTO_DB;
    }
  }
  invert(b, q->inv, q->n);
}

static float clampf(float x, float r) { return x < -r ? -r : (x > r ? r : x); }

void udsp_geq_solve(const UDSP_GEQ *q, const float *target_db, float *gain_db) {
  const int n = q->n;
  double t[UDSP_GEQ_BANDS], g[UDSP_GEQ_BANDS], e[UDSP_GEQ_BANDS];
  for (int i = 0; i < n; i++) {
    t[i] = clampf(target_db[i], UDSP_GEQ_RANGE);
    e[i] = t[i];
    g[i] = 0.0;
  }
  /* the first step from flat is the interaction's inverse alone */
  for (int s = 0; s <= STEPS; s++) {
    for (int j = 0; j < n; j++) {
      double d = 0.0;
      for (int i = 0; i < n; i++) {
        d += q->inv[j][i] * e[i];
      }
      g[j] = clampf((float)(g[j] + d), 30.0f);
    }
    if (s == STEPS) {
      break;
    }
    for (int i = 0; i < n; i++) {
      e[i] = t[i];
    }
    for (int j = 0; j < n; j++) {
      if (g[j] == 0.0) {
        continue;
      }
      BIQUAD f;
      bell(q, j, g[j], &f);
      for (int i = 0; i < n; i++) {
        e[i] -= power_db(&f, q->phi[i]);
      }
    }
  }
  for (int j = 0; j < UDSP_GEQ_BANDS; j++) {
    gain_db[j] = j < n ? (float)g[j] : 0.0f;
  }
}

void udsp_geq_band(const UDSP_GEQ *q, int i, float gain_db, UDSP_EQ_BAND *b) {
  udsp_eq_set(b, UDSP_EQ_BELL_MATCHED, (float)q->hz[i], UDSP_GEQ_Q, gain_db,
              q->rate);
}

float udsp_geq_response(const UDSP_GEQ *q, const float *gain_db, double hz) {
  const double s = udsp_sin_pi(hz / q->rate), phi = s * s;
  double db = 0.0;
  for (int j = 0; j < q->n; j++) {
    if (gain_db[j] != 0.0f) {
      BIQUAD f;
      bell(q, j, gain_db[j], &f);
      db += power_db(&f, phi);
    }
  }
  return (float)db;
}
