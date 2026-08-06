#include "scope.h"

void udsp_scope_extract(const float *mono, size_t n, float *trace_out) {
  size_t offset = 0;

  if (n > UDSP_SCOPE_POINTS) {
    size_t search_end = n - UDSP_SCOPE_POINTS;
    size_t search_start =
        search_end > UDSP_SCOPE_SEARCH_WINDOW ? search_end - UDSP_SCOPE_SEARCH_WINDOW : 0;

    int found = 0;
    for (size_t i = search_start; i < search_end; i++) {
      if (mono[i] <= 0.0f && mono[i + 1] > 0.0f) {
        offset = i;
        found = 1;
        break;
      }
    }
    if (!found) {
      offset = search_end; /* untriggered fallback: most recent points */
    }
  }

  size_t available = n > offset ? n - offset : 0;
  size_t copy = available < UDSP_SCOPE_POINTS ? available : UDSP_SCOPE_POINTS;
  for (size_t i = 0; i < copy; i++) {
    trace_out[i] = mono[offset + i];
  }
  for (size_t i = copy; i < UDSP_SCOPE_POINTS; i++) {
    trace_out[i] = 0.0f;
  }
}
