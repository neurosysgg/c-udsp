#ifndef UDSP_CONV_H
#define UDSP_CONV_H

#include <stddef.h>

/* Uniformly partitioned convolution (overlap-save): an impulse response
 * cut into partitions of `block` frames, each one's spectrum kept, and
 * a signal convolved with it a block at a time -- each block's spectrum
 * (of it and the block before) into a delay line of spectra, multiplied
 * with the partitions' and summed, the sum's inverse the block's output.
 * Its cost per frame is a transform's plus one complex multiply-add per
 * partition and bin, whatever the response's length; its latency is one
 * block (a block's output is ready once the block is in).
 *
 * The pieces are separate so a caller wires its channels: a delay line
 * per input, one accumulator per output, as many multiply-adds into it
 * as inputs reach it (a true-stereo response: four, two inputs' lines
 * each into both outputs). Everything is caller-owned floats; a table is
 * udsp_fft_table()'s for 2 * block. */

/* floats one spectrum (a partition's, a line's slot, an accumulator)
 * holds: block + 1 bins, the real parts then the imaginary ones */
#define UDSP_CONV_SPECTRUM(block) (2 * ((block) + 1))

/* the partitions a response of frames needs (at least 1) */
size_t udsp_conv_parts(size_t frames, size_t block);

/* the response's partitions' spectra into spectra (udsp_conv_parts()
 * times UDSP_CONV_SPECTRUM floats): frames samples read at ir[i *
 * stride], the 1 / (2 * block) of the inverse folded in. scratch holds
 * 2 * block floats. */
void udsp_conv_ir(const float *ir, size_t frames, size_t stride, size_t block,
                  float *spectra, const float *table, float *scratch);

/* an input's delay line: its last parts blocks' spectra, newest at pos */
typedef struct {
  size_t block, parts, pos;
  float *line;  /* parts * UDSP_CONV_SPECTRUM(block) */
  float *frame; /* 2 * block: the block before and this one */
} UDSP_CONV_IN;

/* line and frame given and zeroed */
void udsp_conv_in_init(UDSP_CONV_IN *in, size_t block, size_t parts,
                       float *line, float *frame);

/* the next block (block samples at x[i * stride]) in: its spectrum (of
 * it and the block before) the line's newest */
void udsp_conv_in_push(UDSP_CONV_IN *in, const float *x, size_t stride,
                       const float *table);

/* acc (one spectrum) += the sum over p < parts of in's spectrum p
 * blocks old times the response's partition p (parts no more than the
 * line's) */
void udsp_conv_mac(const UDSP_CONV_IN *in, const float *spectra, size_t parts,
                   float *acc);

/* acc's block of output, block samples to out[i * stride]; acc is left
 * changed (zero it before the next block's multiply-adds). scratch holds
 * 2 * block floats. */
void udsp_conv_out(float *acc, size_t block, float *out, size_t stride,
                   const float *table, float *scratch);

/* --- the head, direct ------------------------------------------------
 *
 * A response's first frames convolved sample by sample, no latency: the
 * head of a non-uniform scheme (the first taps direct, partitions of a
 * small block next, of a big one after; each partition's block no
 * longer than the frames ahead of it, so its latency is hidden). Each
 * output sums its taps in the same order at any vector width (the lanes
 * run across outputs), so it's the same bits on every build. */

/* the head's taps reversed into hr (taps floats: the response's first
 * taps frames at ir[i * stride], zeros past frames) */
void udsp_conv_head(const float *ir, size_t frames, size_t stride, size_t taps,
                    float *hr);

/* out[i] += the sum over j < taps of h[j] x[i - j], for i < n: x points
 * at the first of n new samples, the taps - 1 before it the history the
 * caller keeps in front of them */
void udsp_conv_direct(const float *hr, size_t taps, const float *x, float *out,
                      size_t n);

#endif
