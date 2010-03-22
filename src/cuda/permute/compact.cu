/* -----------------------------------------------------------------------------
 *
 * Module    : Compact
 * Copyright : (c) 2010 Trevor L. McDonell
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#include "utils.h"
#include "device.h"
#include "permute.h"
#include "algorithms.h"

#include <stdint.h>


static void
compact_control(uint32_t n, uint32_t &blocks, uint32_t &threads)
{
    threads = ceilPow2(n + COMPACT_ELT_PER_THREAD - 1) / COMPACT_ELT_PER_THREAD;
    threads = min(threads, MAX_THREADS);

    blocks  = (n + threads * COMPACT_ELT_PER_THREAD - 1) / (threads * COMPACT_ELT_PER_THREAD);
    blocks  = min(blocks, MAX_BLOCKS);
}


/*
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
 *   We should have Ta == Tb for compact, and Tb == uint32_t when compacting
 *   indices, instead of the elements.
 */
template <typename Ta, typename Tb, bool backward, bool compactIdx, bool lengthIsPow2>
__global__ static void
compact_core
(
    const Ta            *d_in,
    Tb                  *d_out,
    const uint32_t      *d_indices,
    const uint32_t      *d_valid,
    uint32_t            *num_valid,
    const uint32_t      length
)
{
    const uint32_t      len2      = length / COMPACT_ELT_PER_THREAD;
    const uint32_t      gridSize  = __umul24(blockDim.x, gridDim.x);
    uint32_t            idx       = __umul24(blockDim.x, blockIdx.x) + threadIdx.x;
    uint2               *d_valid2 = (uint2*) d_valid;
    uint2               tmp;
    uint32_t            gid;

    /*
     * Return the number of valid entries found
     */
    if (idx == 0)
    {
        if (backward) num_valid[0] = d_valid[0] + d_indices[0];
        else          num_valid[0] = d_valid[length-1] + d_indices[length-1];
    }

    for (; idx < len2; idx += gridSize)
    {
        gid = idx * COMPACT_ELT_PER_THREAD;
        tmp = d_valid2[idx];

        /*
         * This is optimised for the case where the majority of elements are
         * invalid, and hence the write bandwidth is very low. Otherwise we
         * should also read the index and input arrays with multiple elements
         * per thread as well.
         */
        if (tmp.x) d_out[d_indices[gid  ]] = compactIdx ? gid   : d_in[gid  ];
        if (tmp.y) d_out[d_indices[gid+1]] = compactIdx ? gid+1 : d_in[gid+1];
//      if (tmp.z) d_out[d_indices[gid+2]] = compactIdx ? gid+2 : d_in[gid+2];
//      if (tmp.w) d_out[d_indices[gid+3]] = compactIdx ? gid+3 : d_in[gid+3];
    }

    idx += (COMPACT_ELT_PER_THREAD - 1) * len2;
    if (idx < length)
    {
        if (d_valid[idx]) d_out[d_indices[idx]] = compactIdx ? idx : d_in[idx];
    }
}


template <typename Ta, typename Tb, bool backward, bool compactIdx>
static unsigned int
compact
(
    const Ta            *d_in,
    Tb                  *d_out,
    const uint32_t      *d_flags,
    const uint32_t      length
)
{
    uint32_t N;
    uint32_t *num;
    uint32_t *indices;
    uint32_t threads;
    uint32_t blocks;

    cudaMalloc((void**) &num, sizeof(uint32_t));
    cudaMalloc((void**) &indices, length * sizeof(uint32_t));

    /*
     * Scan the [0,1] flags to determine the corresponding output array index
     * for valid elements
     */
    if (backward) prescanr_plusui(d_flags, indices, length);
    else          prescanl_plusui(d_flags, indices, length);

    compact_control(length, blocks, threads);

    if (isPow2(length)) compact_core< Ta,Tb,backward,compactIdx,true  ><<<blocks,threads>>>(d_in, d_out, indices, d_flags, num, length);
    else                compact_core< Ta,Tb,backward,compactIdx,false ><<<blocks,threads>>>(d_in, d_out, indices, d_flags, num, length);

    /*
     * Retrieve the number of elements actually stored in the output array
     */
    cudaMemcpy(&N, num, sizeof(uint32_t), cudaMemcpyDeviceToHost);
    cudaFree(indices);
    cudaFree(num);

    return N;
}


// -----------------------------------------------------------------------------
// Instances
// -----------------------------------------------------------------------------

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

