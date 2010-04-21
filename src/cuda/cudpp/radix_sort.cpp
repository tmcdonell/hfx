/* -----------------------------------------------------------------------------
 *
 * Module    : Sort
 * Copyright : (c) [2009..2010] Trevor L. McDonell
 * License   : BSD
 *
 *----------------------------------------------------------------------------*/

#include "algorithms.h"
#include "cudpp/utils.h"

#include <cudpp.h>


/*
 * In-place radix sort of values or key-value pairs. Values can be any 32-bit
 * type, as their payload is never inspected or manipulated.
 */
template <typename T>
static void
radix_sort
(
    T                   *d_keys,
    void                *d_vals,
    unsigned int        length,
    int                 bits    = 8 * sizeof(T)
)
{
    CUDPPHandle         plan;
    CUDPPConfiguration  cp;

    cp.datatype  = getType<T>();
    cp.algorithm = CUDPP_SORT_RADIX;
    cp.options   = (d_vals != NULL) ? CUDPP_OPTION_KEY_VALUE_PAIRS
                                    : CUDPP_OPTION_KEYS_ONLY;

    cudppPlan(&plan, cp, length, 1, 0);
    cudppSort(plan, d_keys, d_vals, bits, length);

    cudppDestroyPlan(plan);
}


/* -----------------------------------------------------------------------------
 * Instances
 * ---------------------------------------------------------------------------*/

void radixsort_f(float *d_keys, void *d_vals, unsigned int N)
{
    radix_sort<float>(d_keys, d_vals, N);
}

