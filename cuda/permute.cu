/*
 * Module    : Prelude
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 */

#include "utils.h"
#include "kernels.h"


/*
 * Permute an array according to the permutation indices. Assumes that
 * length(indices) << length(in), so makes no attempt to coalesce memory access.
 */
template <typename T, bool lengthIsPow2>
__global__ static void
permute_core
(
    const T             *in,
    T                   *out,
    int                 *indices,
    int                 length
)
{
    unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;

    if (lengthIsPow2 || idx < length)
        out[idx] = in[indices[idx]];
}


template <typename T>
static void
permute
(
    const T             *in,
    T                   *out,
    int                 *indices,
    int                 length
)
{
    unsigned int threads = min(ceilPow2(length), 512);
    unsigned int blocks  = (length + threads - 1) / threads;

    if (isPow2(length))
        permute_core< T,true  ><<<blocks,threads>>>(in, out, indices, length);
    else
        permute_core< T,false ><<<blocks,threads>>>(in, out, indices, length);
}


// -----------------------------------------------------------------------------
// Instances
// -----------------------------------------------------------------------------

void permute_f(float *in, float *out, int *indices, int length)
{
    permute<float>(in, out, indices, length);
}

