/* -----------------------------------------------------------------------------
 *
 * Module    : Algorithms
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/


#ifndef __ALGORITHMS_H__
#define __ALGORITHMS_H__

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Sparse-matrix dense-vector multiplication
 */
void
smvm_f(float *d_y, float *d_x, float *d_data, unsigned int *d_rowPtr, unsigned int *d_colIdx, unsigned int num_rows);

#ifdef __cplusplus
}
#endif
#endif
