/*
 * Module    : Prelude
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 */

#ifndef __MAP_KERNEL__
#define __MAP_KERNEL__

#include "operator.h"

/*
 * Apply a function each element of an array. A single thread is used to compute
 * each result. The input and output array may be coincident.
 */
template <class fn, typename Ta, typename Tb>
__global__ void
map
(
    Ta *xs,
    Tb *out
)
{
    unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
    out[idx]         = fn::apply(xs[idx]);
}

#endif

