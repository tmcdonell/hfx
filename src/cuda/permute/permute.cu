/* -----------------------------------------------------------------------------
 *
 * Module    : Permute
 * Copyright : (c) 2009 Trevor L. McDonell
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
    threads = min(ceilPow2(n), MAX_THREADS);
    blocks  = min((n + threads - 1) / threads, MAX_BLOCKS);
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
 *
 * An alternative to back permute, where we do not explicitly know the offset
 * indices, is to use an array of [0,1] flags specifying valid elements, which
 * we can exclusive-sum-scan to get the offsets. In this case, `length'
 * specifies the number of elements in the `in' array. The template parameter
 * `backward' specifies whether the offsets were calculated with a left or right
 * scan. A further alternative to the `compact' method, is to store the array
 * indices of the valid elements, rather than the values themselves.
 *
 * We return the number of valid elements found via `num_valid', but blunder on
 * ahead regardless of whether the `out' array is large enough or not.
 *
 * NOTE:
 *   Under most circumstances we should have Ta == Tb. The only exception is
 *   when compact and compactIdx are active, in which case Tb == uint32_t.
 */
template <typename Ta, typename Tb, bool backward, bool compact, bool compactIdx, bool lengthIsPow2>
__global__ static void
permute_core
(
    const Ta            *d_in,
    Tb                  *d_out,
    const uint32_t      *indices,
    const uint32_t      length,
    const uint32_t      *valid     = NULL,
    uint32_t            *num_valid = NULL
)
{
    uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;

    /*
     * Return the number of valid entries found
     */
    if (compact && idx == 0)
    {
        if (backward) num_valid[0] = valid[0] + indices[0];
        else          num_valid[0] = valid[length-1] + indices[length-1];
    }

    for (; idx < length; idx += gridDim.x)
    {
        if (compact)
        {
            if (valid[idx]) d_out[indices[idx]] = compactIdx ? idx : d_in[idx];
        }
        else
        {
            if (backward)   d_out[idx]          = d_in[indices[idx]];
            else            d_out[indices[idx]] = d_in[idx];
        }
    }
}


template <typename Ta, typename Tb, bool backward, bool compact, bool compactIdx>
static void
permute
(
    const Ta            *d_in,
    Tb                  *d_out,
    const uint32_t      *indices,
    const uint32_t      length,
    const uint32_t      *valid     = NULL,
    uint32_t            *num_valid = NULL
)
{
    uint32_t threads;
    uint32_t blocks;

    permute_control(length, blocks, threads);

    if (isPow2(length)) permute_core< Ta,Tb,backward,compact,compactIdx,true  ><<<blocks,threads>>>(d_in, d_out, indices, length, valid, num_valid);
    else                permute_core< Ta,Tb,backward,compact,compactIdx,false ><<<blocks,threads>>>(d_in, d_out, indices, length, valid, num_valid);
}


template <typename Ta, typename Tb, bool backward, bool compactIdx>
static unsigned int
compact
(
    const Ta            *d_in,
    Ta                  *d_out,
    const uint32_t      *d_flags,
    const uint32_t      length
)
{
    uint32_t N;
    uint32_t *num;
    uint32_t *indices;

    cudaMalloc((void**) &num, sizeof(uint32_t));
    cudaMalloc((void**) &indices, length * sizeof(uint32_t));

    if (backward) prescanr_plusui(d_flags, indices, length);
    else          prescanl_plusui(d_flags, indices, length);

    /*
     * At this point, we know exactly how many elements will be required for the
     * `out' array. Maybe we should allocate that array here?
     */
    permute<Ta,Tb,backward,true,compactIdx>(d_in, d_out, indices, length, d_flags, num);

    cudaMemcpy(&N, num, sizeof(uint32_t), cudaMemcpyDeviceToHost);
    cudaFree(num);
    cudaFree(indices);

    return N;
}


// -----------------------------------------------------------------------------
// Instances
// -----------------------------------------------------------------------------

void permute(const void *d_in, void *d_out, const uint32_t *d_indices, const uint32_t length)
{
    permute<uint32_t,uint32_t,false,false,false>((const uint32_t*) d_in, (uint32_t*) d_out, d_indices, length);
}

void bpermute(const void *d_in, void *d_out, const uint32_t *d_indices, const uint32_t length)
{
    permute<uint32_t,uint32_t,true,false,false>((const uint32_t*) d_in, (uint32_t*) d_out, d_indices, length);
}

unsigned int compact(const void *d_in, void *d_out, const uint32_t *d_flags, const uint32_t length)
{
    unsigned int N = compact<uint32_t,uint32_t,false,false>((const uint32_t*) d_in, (uint32_t*) d_out, d_flags, length);
    return N;
}

unsigned int compactIndices(uint32_t *d_out, const uint32_t *d_flags, const uint32_t length)
{
    unsigned int N = compact<uint32_t,uint32_t,false,true>(NULL, d_out, d_flags, length);
    return N;
}

