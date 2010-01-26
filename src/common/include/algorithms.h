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
 * Scan
 */
void prescanl_plusui(const unsigned int *in, unsigned int *out, const unsigned int N);

/*
 * Sparse-matrix dense-vector multiplication
 */
void smvm_f(float *d_y, float *d_x, float *d_data, unsigned int *d_rowPtr, unsigned int *d_colIdx, unsigned int num_rows);


/*
 * Radix sort (in-place)
 */
void radixsort_f(float *d_keys, void *d_vals, unsigned int N);


/*
 * Enum
 */
void enumFromTo_i(int *out, int from, int to);

#ifdef __cplusplus
}
#endif
#endif
