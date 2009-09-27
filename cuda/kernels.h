/*
 * Module    : Kernels
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 *
 * Functions which are internally implemented using CUDA, public interface.
 */

#ifndef __KERNELS_H__
#define __KERNELS_H__

#ifdef __cplusplus
extern "C" {
#endif

// -----------------------------------------------------------------------------
// Ion Series
// -----------------------------------------------------------------------------
void
addIons
(
    int                 charge,
    float               residual,
    float               *y_ions,
    int                 *out,
    unsigned int        len_ions,
    unsigned int        len_spec,
    unsigned int	offset
);


// -----------------------------------------------------------------------------
// Prelude
// -----------------------------------------------------------------------------
float fold_plusf(float *xs, int N);

void  zipWith_timesif(int *xs, float *ys, float *zs, int N);

void  permute_ui(unsigned int *in, unsigned int *out, unsigned int *indices, int length);
void  bpermute_f(float *in, float *out, unsigned int *indices, int length);
int   compact_f(float *in, float *out, unsigned int *flags, int length);

void  scanl_plusui(unsigned int *in, unsigned int *out, int N);
void  scanr_plusui(unsigned int *in, unsigned int *out, int N);

void  scanl1Seg_plusf(float *in, unsigned int *flags, float *out, int N);
void  scanr1Seg_plusf(float *in, unsigned int *flags, float *out, int N);

void  map_getAAMass(char *ions, float *masses, int N);

void  replicate_ui(unsigned int *out, unsigned int x, int N);

void  sort_f(float *vals, int N);
void  sortPairs_f(float *keys, float *vals, int N);

void  group_f(float *xs, unsigned int *flags, int N);

#ifdef __cplusplus
}
#endif
#endif

