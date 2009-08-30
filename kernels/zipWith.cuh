/*
 * Module    : Prelude
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 */

#ifndef __ZIPWITH_KERNEL__
#define __ZIPWITH_KERNEL__

#include "operator.h"

/*
 * Combine two arrays using the given binary operator function. A single thread
 * is used to compute each result pair.
 */
template <class op, typename Ta, typename Tb, typename Tc>
__global__ void
zipWith
(
    Ta *xs,
    Tb *ys,
    Tc *out
)
{
    unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
    out[idx]         = op::apply(xs[idx], ys[idx]);
}

#endif

