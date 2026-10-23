#include "udsp_test.h"

#include "analysis/r128.h"

#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* fed in odd-sized chunks, so a hop's end falls mid-chunk */
#define CHUNK 1237

typedef struct {
  double dbfs; /* per channel; -HUGE_VAL silent */
  double seconds;
} SEGMENT;

/* a 1 kHz sine, the same phase in every channel, each channel at its own
 * level for each segment (EBU Tech 3341's test signals) */
static void feed(UDSP_R128 *m, int rate, int channels,
                 const double *dbfs_per_ch, double seconds, double *phase) {
  static float buf[CHUNK * UDSP_R128_MAX_CHANNELS];
  size_t frames = (size_t)(seconds * rate + 0.5);
  while (frames > 0) {
    size_t n = frames < CHUNK ? frames : CHUNK;
    for (size_t i = 0; i < n; i++) {
      double s = sin(*phase);
      *phase += 2.0 * M_PI * 1000.0 / rate;
      for (int c = 0; c < channels; c++) {
        buf[i * channels + c] = (float)(s * pow(10.0, dbfs_per_ch[c] / 20.0));
      }
    }
    udsp_r128_add(m, buf, n);
    frames -= n;
  }
}

static double stereo_program(int rate, const SEGMENT *seg, int nseg) {
  static UDSP_R128 m;
  udsp_r128_init(&m, rate, 2);
  double phase = 0.0;
  for (int i = 0; i < nseg; i++) {
    double lv[2] = {seg[i].dbfs, seg[i].dbfs};
    feed(&m, rate, 2, lv, seg[i].seconds, &phase);
  }
  return udsp_r128_integrated(&m.hist);
}

static void test_kweight_at_48k_is_the_spec_table(void) {
  static const double spec[2][5] = {
      {1.53512485958697, -2.69169618940638, 1.19839281085285, -1.69065929318241,
       0.73248077421585},
      {1.0, -2.0, 1.0, -1.99004745483398, 0.99007225036621}};
  double k[2][5];
  udsp_r128_kweight(48000, k);
  double worst = 0.0;
  for (int s = 0; s < 2; s++) {
    for (int j = 0; j < 5; j++) {
      worst = fmax(worst, fabs(k[s][j] - spec[s][j]));
    }
  }
  UDSP_TEST_ASSERT(worst < 1e-8,
                   "48 kHz K-weighting matches BS.1770's table (%.2e)", worst);
}

static void test_tech3341(void) {
  static const SEGMENT c1[] = {{-23, 20}};
  static const SEGMENT c2[] = {{-33, 20}};
  static const SEGMENT c3[] = {{-36, 10}, {-23, 60}, {-36, 10}};
  static const SEGMENT c4[] = {
      {-72, 10}, {-36, 10}, {-23, 20}, {-36, 10}, {-72, 10}};
  static const SEGMENT c5[] = {{-26, 20}, {-20, 20.1}, {-26, 20}};
  struct {
    const char *name;
    const SEGMENT *seg;
    int n;
    double want;
  } cases[] = {
      {"1 (-23 dBFS)", c1, 1, -23},      {"2 (-33 dBFS)", c2, 1, -33},
      {"3 (relative gate)", c3, 3, -23}, {"4 (both gates)", c4, 5, -23},
      {"5 (-26/-20/-26)", c5, 3, -23},
  };
  for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
    double l = stereo_program(48000, cases[i].seg, cases[i].n);
    UDSP_TEST_ASSERT(fabs(l - cases[i].want) < 0.1,
                     "Tech 3341 case %s reads %.1f LUFS (%.4f)", cases[i].name,
                     cases[i].want, l);
  }

  static UDSP_R128 m;
  udsp_r128_init(&m, 48000, 5);
  double lv[5] = {-28, -28, -24, -30, -30}, phase = 0.0;
  feed(&m, 48000, 5, lv, 20, &phase);
  double l = udsp_r128_integrated(&m.hist);
  UDSP_TEST_ASSERT(fabs(l + 23.0) < 0.1,
                   "Tech 3341 case 6 (5.0, surround weights) reads "
                   "-23.0 LUFS (%.4f)",
                   l);
}

static void test_other_rates(void) {
  static const SEGMENT c1[] = {{-23, 20}};
  static const int rates[] = {44100, 96000, 192000, 22050};
  for (size_t i = 0; i < sizeof(rates) / sizeof(rates[0]); i++) {
    double l = stereo_program(rates[i], c1, 1);
    UDSP_TEST_ASSERT(fabs(l + 23.0) < 0.1,
                     "-23 dBFS at %d Hz reads -23.0 LUFS (%.4f)", rates[i], l);
  }
}

static void test_nothing_to_measure(void) {
  static const SEGMENT silent[] = {{-HUGE_VAL, 5}};
  static const SEGMENT short_tone[] = {{-10, 0.39}};
  static const SEGMENT under_gate[] = {{-80, 5}};
  UDSP_TEST_ASSERT(stereo_program(44100, silent, 1) == -HUGE_VAL,
                   "silence measures nothing");
  UDSP_TEST_ASSERT(stereo_program(44100, short_tone, 1) == -HUGE_VAL,
                   "under 400 ms measures nothing (no whole block)");
  UDSP_TEST_ASSERT(stereo_program(44100, under_gate, 1) == -HUGE_VAL,
                   "-80 dBFS is under the absolute gate");
}

static void test_chunking_does_not_matter(void) {
  static UDSP_R128 a, b;
  static float buf[44100 * 3 * 2];
  unsigned int seed = 7;
  for (size_t i = 0; i < sizeof(buf) / sizeof(buf[0]); i++) {
    seed = seed * 1103515245u + 12345u;
    buf[i] = 0.4f * ((float)((seed >> 8) & 0xffff) / 32768.0f - 1.0f);
  }
  size_t frames = sizeof(buf) / sizeof(buf[0]) / 2;
  udsp_r128_init(&a, 44100, 2);
  udsp_r128_add(&a, buf, frames);
  udsp_r128_init(&b, 44100, 2);
  for (size_t at = 0, n = 1; at < frames; at += n, n = n * 3 + 1) {
    udsp_r128_add(&b, buf + at * 2, at + n > frames ? frames - at : n);
  }
  UDSP_TEST_ASSERT(
      memcmp(a.hist.count, b.hist.count, sizeof(a.hist.count)) == 0 &&
          udsp_r128_integrated(&a.hist) == udsp_r128_integrated(&b.hist) &&
          a.hist.peak == b.hist.peak,
      "one chunk or many, the same blocks and loudness (%.4f)",
      udsp_r128_integrated(&a.hist));
}

static void test_album_and_peak(void) {
  static UDSP_R128 t1, t2;
  static UDSP_R128_HIST album;
  double phase = 0.0, lv1[2] = {-23, -23}, lv2[2] = {-33, -33};
  udsp_r128_init(&t1, 48000, 2);
  feed(&t1, 48000, 2, lv1, 20, &phase);
  udsp_r128_init(&t2, 48000, 2);
  feed(&t2, 48000, 2, lv2, 20, &phase);
  memset(&album, 0, sizeof(album));
  udsp_r128_hist_add(&album, &t1.hist);
  udsp_r128_hist_add(&album, &t2.hist);
  /* both pass the gate (-33 is above -25.6 - 10): the energies' mean */
  double want = -23.0 + 10.0 * log10((1.0 + 0.1) / 2.0);
  double l = udsp_r128_integrated(&album);
  UDSP_TEST_ASSERT(fabs(l - want) < 0.02,
                   "an album of -23 and -33 reads their energies' "
                   "mean, %.2f LUFS (%.4f)",
                   want, l);
  float p = (float)pow(10.0, -23.0 / 20.0);
  UDSP_TEST_ASSERT(fabsf(album.peak - p) < 1e-4f,
                   "the album's peak is its loudest track's, "
                   "%.4f (%.4f)",
                   p, album.peak);
}

static void test_lfe_is_not_counted(void) {
  static UDSP_R128 a, b;
  double phase = 0.0, with[6] = {-28, -28, -24, -10, -30, -30},
         without[6] = {-28, -28, -24, -HUGE_VAL, -30, -30};
  udsp_r128_init(&a, 48000, 6);
  feed(&a, 48000, 6, with, 5, &phase);
  phase = 0.0;
  udsp_r128_init(&b, 48000, 6);
  feed(&b, 48000, 6, without, 5, &phase);
  UDSP_TEST_ASSERT(
      udsp_r128_integrated(&a.hist) == udsp_r128_integrated(&b.hist),
      "5.1's LFE adds nothing (%.4f)", udsp_r128_integrated(&a.hist));
}

int main(void) {
  test_kweight_at_48k_is_the_spec_table();
  test_tech3341();
  test_other_rates();
  test_nothing_to_measure();
  test_chunking_does_not_matter();
  test_album_and_peak();
  test_lfe_is_not_counted();
  return udsp_test_summary();
}
