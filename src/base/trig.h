#ifndef UDSP_TRIG_H
#define UDSP_TRIG_H

/* sin(pi y) in plain double arithmetic: libm's sin may round differently
 * between libcs, this is the same bits everywhere (-ffp-contract=off) --
 * a resampler's kernel, a filter's design. Within 1e-15 of sin. */
double udsp_sin_pi(double y);

#endif
