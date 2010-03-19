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

#include <stdint.h>


static void
filter_control(uint32_t n, uint32_t &blocks, uint32_t &threads)
{
    threads = min(ceilPow2(n), MAX_THREADS);
    blocks  = min((n + threads - 1) / threads, MAX_BLOCKS);
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
    uint32_t idx;

    /*
     * Mark elements that pass the predicate with a [0,1] head flag. Can not
     * store the index directly, as `compact' needs to scan this result. Boo.
     */
    for (idx = blockIdx.x * blockDim.x + threadIdx.x; idx < length; idx += gridDim.x)
    {
        T val        = d_in[idx];
        d_valid[idx] = (m <= val && val <= n);
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

