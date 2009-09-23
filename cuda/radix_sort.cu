/*
 * Module    : Prelude
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 */

#include "kernels.h"

#include <cudpp.h>


template <typename T> CUDPPDatatype getType();
template <> CUDPPDatatype getType<float>() { return CUDPP_FLOAT; }
template <> CUDPPDatatype getType<unsigned int>() { return CUDPP_UINT; }


/*
 * In place sort of values or key-value pairs.
 */
template <typename T>
static void
sort
(
    int                 length,
    T                   *keys,
    T                   *vals   = NULL,
    int                 bits    = 8 * sizeof(T)
)
{
    CUDPPHandle         plan;
    CUDPPConfiguration  cp;

    cp.datatype  = getType<T>();
    cp.algorithm = CUDPP_SORT_RADIX;
    cp.options   = (vals != NULL) ? CUDPP_OPTION_KEY_VALUE_PAIRS
                                  : CUDPP_OPTION_KEYS_ONLY;

    cudppPlan(&plan, cp, length, 1, 0);
    cudppSort(plan, keys, vals, bits, length);

    cudppDestroyPlan(plan);
}


// -----------------------------------------------------------------------------
// Instances
// -----------------------------------------------------------------------------

void sort_f(float *vals, int length)
{
    sort<float>(length, vals);
}

void sortPairs_f(float *keys, float *vals, int length)
{
    sort<float>(length, keys, vals);
}

