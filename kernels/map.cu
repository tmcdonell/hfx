/*
 * Module    : Prelude
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 */

#include "utils.h"
#include "kernels.h"
#include "operator.h"

/*
 * Apply a function each element of an array. A single thread is used to compute
 * each result. The input and output array may be coincident.
 */
template <class fn, bool lengthIsPow2, typename Ta, typename Tb>
__global__ static void
map_core
(
    Ta  *xs,
    Tb  *out,
    int length
)
{
    unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;

    if (lengthIsPow2 || idx < length)
        out[idx] = fn::apply(xs[idx]);
}


template <class fn, typename Ta, typename Tb>
void
map
(
    Ta  *xs,
    Tb  *out,
    int length
)
{
    unsigned int threads = min(ceilPow2(length), 512);
    unsigned int blocks  = (length + threads - 1) / threads;

    if (isPow2(length))
        map_core< fn,true  ><<<blocks,threads>>>(xs, out, length);
    else
        map_core< fn,false ><<<blocks,threads>>>(xs, out, length);
}


// -----------------------------------------------------------------------------
// Instances
// -----------------------------------------------------------------------------

void map_fromIntegralf(int *xs, float *out, int N)
{
    map< fromIntegral<int,float> >(xs, out, N);
}

