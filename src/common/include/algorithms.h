/* -----------------------------------------------------------------------------
 *
 * Module    : Algorithms
 * Copyright : (c) [2009..2010] Trevor L. McDonell
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/


#ifndef __ALGORITHMS_H__
#define __ALGORITHMS_H__

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif


#if 0
/*
 * Scan and segmented-scan
 * Unlike Haskell-style scans, the size of the output array will not change.
 *   Pre-scans are exclusive:   prescanl  == init . scanl
 *   Post-scans are inclusive:  postscanl == scanl1
 */
void prescanl_plusui(const uint32_t *d_in, uint32_t *d_out, const uint32_t N);
void prescanr_plusui(const uint32_t *d_in, uint32_t *d_out, const uint32_t N);
#endif

/*
 * key-value sort (in-place)
 */
void sort_rf(float *d_keys, uint32_t *d_vals, uint32_t N);

#if 0
/*
 * Permute (32-bit payload)
 */
void permute(const void *d_in, void *d_out, const uint32_t *d_indices, const uint32_t length);
void bpermute(const void *d_in, void *d_out, const uint32_t *d_indices, const uint32_t length);
uint32_t compact(const void *d_in, void *d_out, const uint32_t *d_flags, const uint32_t length);
uint32_t compactIndices(uint32_t *d_out, const uint32_t *d_flags, const uint32_t length);
#endif
#if 0
/*
 * Filter
 */
uint32_t
filterInRange_f
(
    const float         *d_in,
    float               *d_out,
    const uint32_t      length,
    const float         min,
    const float         max
);
#endif

uint32_t
findIndicesInRange_f
(
    const float         *d_in,
    uint32_t            *d_out,
    const uint32_t      length,
    const float         min,
    const float         max
);

/*
 * Generate theoretical spectra
 */
void
addIons
(
    uint32_t            *d_spec,
    const float         *d_residual,
    const float         *d_mass,
    const uint8_t       *d_ions,
    const uint32_t      *d_tc,
    const uint32_t      *d_tn,
    const uint32_t      *d_idx,
    const uint32_t      num_idx,
    const uint32_t      max_charge,
    const uint32_t      len_spec
);

#if 0
void
addIons_inplace
(
    float               *d_score,
    const float         *d_spec,
    const float         *d_residual,
    const float         *d_ions,
    const uint32_t      *d_tc,
    const uint32_t      *d_tn,
    const uint32_t      *d_idx,
    const uint32_t      num_idx,
    const uint32_t      max_charge,
    const uint32_t      len_spec
);
#endif

/*
 * Dense matrix-vector multiplication
 */
void mvm_if(float *d_y, const uint32_t *d_A, const float *d_x, const uint32_t m, const uint32_t n);


#ifdef __cplusplus
}
#endif
#endif
