/* -----------------------------------------------------------------------------
 *
 * Module    : Sort
 * Copyright : (c) [2009..2011] Trevor L. McDonell
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <stdint.h>

#include "algorithms.h"


void sort_rf(float *d_keys_raw, uint32_t *d_vals_raw, uint32_t N)
{
    thrust::device_ptr<float>    d_keys(d_keys_raw);
    thrust::device_ptr<uint32_t> d_vals(d_vals_raw);

    thrust::sort_by_key(d_keys, d_keys+N, d_vals, thrust::greater<float>());
}

