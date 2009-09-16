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
    float               *b_ions,
    float               *y_ions,
    int                 *out,
    unsigned int        len_ions,
    unsigned int        len_spec
);


// -----------------------------------------------------------------------------
// Prelude
// -----------------------------------------------------------------------------
float fold_plusf(float *xs, int N);
void  zipWith_timesif(int *xs, float *ys, float *zs, int N);


#ifdef __cplusplus
}
#endif
#endif
/*
 * vim: syn=cuda
 */
