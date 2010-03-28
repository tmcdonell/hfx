/* -----------------------------------------------------------------------------
 *
 * Module    : Permute
 * Copyright : (c) [2009..2010] Trevor L. McDonell
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#include "utils.h"
#include "device.h"
#include "permute.h"
#include "algorithms.h"

#include <stdint.h>


static void
permute_control(uint32_t n, uint32_t &blocks, uint32_t &threads)
{
    threads = ceilPow2(n + PERMUTE_ELT_PER_THREAD - 1) / PERMUTE_ELT_PER_THREAD;
    threads = min(threads, MAX_THREADS);

    blocks  = (n + threads * PERMUTE_ELT_PER_THREAD - 1) / (threads * PERMUTE_ELT_PER_THREAD);
    blocks  = min(blocks, MAX_BLOCKS);
}


/*
 * Permute an array according to the permutation indices. This handles both
 * forward (scatter) and backward (gather) permutation, where:
 *
 *   bpermute :: [a] -> [Int] -> [a]
 *   bpermute v is = [ v!i | i <- is ]
 *
 * In this case, `length' specifies the number of elements in the `indices' and
 * `out' arrays.
 */
template <typename Ta, typename Tb, bool backward, bool lengthIsPow2>
__global__ static void
permute_core
(
    const Ta            *d_in,
    Tb                  *d_out,
    const uint32_t      *indices,
    const uint32_t      length
)
{
    uint32_t       idx;
    const uint32_t gridSize = __umul24(blockDim.x, gridDim.x);

    for (idx = __umul24(blockDim.x, blockIdx.x) + threadIdx.x; idx < length; idx += gridSize)
    {
        if (backward) d_out[idx]          = d_in[indices[idx]];
        else          d_out[indices[idx]] = d_in[idx];
    }
}


template <typename Ta, typename Tb, bool backward>
static void
permute
(
    const Ta            *d_in,
    Tb                  *d_out,
    const uint32_t      *indices,
    const uint32_t      length
)
{
    uint32_t threads;
    uint32_t blocks;

    permute_control(length, blocks, threads);

    if (isPow2(length)) permute_core< Ta,Tb,backward,true  ><<<blocks,threads>>>(d_in, d_out, indices, length);
    else                permute_core< Ta,Tb,backward,false ><<<blocks,threads>>>(d_in, d_out, indices, length);
}


// -----------------------------------------------------------------------------
// Instances
// -----------------------------------------------------------------------------

void permute(const void *d_in, void *d_out, const uint32_t *d_indices, const uint32_t length)
{
    permute<uint32_t,uint32_t,false>((const uint32_t*) d_in, (uint32_t*) d_out, d_indices, length);
}

void bpermute(const void *d_in, void *d_out, const uint32_t *d_indices, const uint32_t length)
{
    permute<uint32_t,uint32_t,true>((const uint32_t*) d_in, (uint32_t*) d_out, d_indices, length);
}

