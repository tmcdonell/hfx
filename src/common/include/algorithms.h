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
 * Scan and segmented-scan
 * Unlike Haskell-style scans, the size of the output array will not change.
 *   Pre-scans are exclusive:   prescanl  == init . scanl
 *   Post-scans are inclusive:  postscanl == scanl1
 */
void prescanl_plusui(const unsigned int *d_in, unsigned int *d_out, const unsigned int N);
void prescanr_plusui(const unsigned int *d_in, unsigned int *d_out, const unsigned int N);

void postsegscanr_plusf(const float *d_in, const unsigned int *d_flags, float *d_out, const unsigned int N);

/*
 * Sparse-matrix dense-vector multiplication
 */
void smvm_f(float *d_y, float *d_x, float *d_data, unsigned int *d_rowPtr, unsigned int *d_colIdx, unsigned int num_rows);

/*
 * Radix sort (in-place)
 */
void radixsort_f(float *d_keys, void *d_vals, unsigned int N);

/*
 * Permute (32-bit payload)
 */
void permute(const void *d_in, void *d_out, const unsigned int *d_indices, const unsigned int length);
void bpermute(const void *d_in, void *d_out, const unsigned int *d_indices, const unsigned int length);
unsigned int compact(const void *d_in, void *d_out, const unsigned int *d_flags, const unsigned int length);
unsigned int compactIndices(unsigned int *d_out, const unsigned int *d_flags, const unsigned int length);

/*
 * Replicate (32-bit symbol)
 */
void replicate(void *d_out, const unsigned int symbol, const unsigned int N);

/*
 * Filter
 */
unsigned int filterInRange_f(const float *d_in, float *d_out, const unsigned int length, const float min, const float max);
unsigned int findIndicesInRange_f(const float *d_in, unsigned int *d_out, const unsigned int length, const float min, const float max);

/*
 * Generate theoretical spectra
 */
void addIons(unsigned int *d_spec, const float *d_residuals, const float *d_ladder, const unsigned int *d_rowPtr, const unsigned int *d_inRangeIdx, const unsigned int num_inRange, const unsigned int max_charge, const unsigned int len_spec);

/*
 * Dense matrix-vector multiplication
 */
void mvm_if(float *d_y, const unsigned int *d_A, const float *d_x, const unsigned int m, const unsigned int n);


#ifdef __cplusplus
}
#endif
#endif
