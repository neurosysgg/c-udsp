#ifndef UDSP_TEST_H
#define UDSP_TEST_H

/* header-only, no dependency beyond libc -- plain C unit tests for this
 * library, mirroring the assert/summary pair vendor projects (ctui,
 * ctui-mus) already use for their own headless tests. */

#include <stdio.h>

static int udsp_test_pass = 0;
static int udsp_test_fail = 0;

#define UDSP_TEST_ASSERT(cond, ...)                                            \
  do {                                                                         \
    if (cond) {                                                                \
      udsp_test_pass++;                                                        \
      printf("  ok - " __VA_ARGS__);                                           \
      printf("\n");                                                            \
    } else {                                                                   \
      udsp_test_fail++;                                                        \
      printf("  FAIL - " __VA_ARGS__);                                         \
      printf(" (%s:%d)\n", __FILE__, __LINE__);                                \
    }                                                                          \
  } while (0)

static int udsp_test_summary(void) {
  printf("-- %d passed, %d failed --\n", udsp_test_pass, udsp_test_fail);
  return udsp_test_fail == 0 ? 0 : 1;
}

#endif
