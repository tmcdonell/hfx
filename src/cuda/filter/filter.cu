/* -----------------------------------------------------------------------------
 *
 * Module    : Filter
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#include "utils.h"
#include "device.h"
#include "filter.h"
#include "algorithms.h"
#include "cudpp/type_vector.h"

#include <stdint.h>


static void
filter_control(uint32_t n, uint32_t &blocks, uint32_t &threads)
{
    threads = ceilPow2(n + FILTER_ELT_PER_THREAD - 1) / FILTER_ELT_PER_THREAD;
    threads = min(threads, MAX_THREADS);

    blocks  = (n + threads * FILTER_ELT_PER_THREAD - 1) / (threads * FILTER_ELT_PER_THREAD);
    blocks  = min(blocks, MAX_BLOCKS);
}


/*
 * Filter an array of elements predicated by lying inside some (inclusive)
 * range. The output array will be filled with zeros for elements that do not
 * pass, or the array index of the elements that do.
 *
 * The filtering result will be compacted to a dense array, with the total
 * number of elements fulfilling the predicate returned. The valid array
 * elements themselves or the indices of said elements may be returned.
 */
template <typename T, bool lengthIsPow2>
__global__ static void
filterInRange_core
(
    const T             *d_in,
    uint32_t            *d_valid,
    const uint32_t      length,
    const T             m,
    const T             n
)
{
    uint32_t            idx;
    uint4               tmp;
    uint4               *d_valid4 = (uint4*) d_valid;
    const uint32_t      len4      = length / FILTER_ELT_PER_THREAD;
    const uint32_t      gridSize  = __umul24(blockDim.x, gridDim.x);

    typename typeToVector<T,4>::Result  val;
    typename typeToVector<T,4>::Result* d_in4 = (typename typeToVector<T,4>::Result*) d_in;

    /*
     * Mark elements that pass the predicate with a [0,1] head flag. Can not
     * store the index directly, as `compact' needs to scan this result. Boo.
     */
    for (idx = __umul24(blockIdx.x, blockDim.x) + threadIdx.x; idx < len4; idx += gridSize)
    {
        val   = d_in4[idx];
        tmp.x = (m <= val.x && val.x <= n);
        tmp.y = (m <= val.y && val.y <= n);
        tmp.z = (m <= val.z && val.z <= n);
        tmp.w = (m <= val.w && val.w <= n);

        d_valid4[idx] = tmp;
    }

    /*
     * Translate the indices back into single-element coordinates. This method
     * ensures that that successive threads retain sequential indices.
     */
    idx += (FILTER_ELT_PER_THREAD - 1) * len4;
    if (idx < length)
    {
        T x          = d_in[idx];
        d_valid[idx] = (m <= x && x <= n);
    }
}


template <typename T, bool findIndices>
static unsigned int
filterInRange
(
    const T             *d_in,
    void                *d_out,
    const uint32_t      length,
    const T             m,
    const T             n
)
{
    uint32_t            threads;
    uint32_t            blocks;
    uint32_t            *d_valid;

    cudaMalloc((void**) &d_valid, length * sizeof(uint32_t));

    /*
     * Determine which elements satisfy the predicate
     */
    filter_control(length, blocks, threads);
    if (isPow2(length)) filterInRange_core< T,true ><<<blocks,threads>>>(d_in, d_valid, length, m, n);
    else                filterInRange_core< T,false><<<blocks,threads>>>(d_in, d_valid, length, m, n);

    /*
     * Compact these elements into a dense array. The output array should be at
     * least as large as the input, and the actual number of elements will be
     * returned.
     *
     * Doesn't like `d_out' having different types based on a template
     * parameter, so make untyped (compact casts internally anyway).
     */
    uint32_t N;
    if (findIndices) N = compactIndices((uint32_t*)  d_out, d_valid, length);
    else             N = compact((const void*) d_in, d_out, d_valid, length);

    cudaFree(d_valid);
    return N;
}


// -----------------------------------------------------------------------------
// Instances
// -----------------------------------------------------------------------------

unsigned int filterInRange_f(const float *d_in, float *d_out, const uint32_t length, const float m, const float n)
{
    unsigned int N = filterInRange<float,false>(d_in, d_out, length, m, n);
    return N;
}

unsigned int findIndicesInRange_f(const float *d_in, uint32_t *d_out, const uint32_t length, const float m, const float n)
{
    unsigned int N = filterInRange<float,true>(d_in, d_out, length, m, n);
    return N;
}

