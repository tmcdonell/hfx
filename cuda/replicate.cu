/*
 * Module    : Prelude
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 */

#include "kernels.h"
#include "operator.h"
#include "cudpp/cudpp_globals.h"

static void
replicate_control
(
    unsigned int        n,
    unsigned int        &blocks,
    unsigned int        &threads
)
{
    blocks  = max(1.0, (double)n / (SCAN_ELTS_PER_THREAD * CTA_SIZE));
    threads = blocks > 1 ? CTA_SIZE : ceil((double)n / SCAN_ELTS_PER_THREAD);
}


/*
 * Apply a function each element of an array. A single thread is used to compute
 * each result. The input and output array may be coincident.
 */
template <typename T>
__global__ static void
replicate_core
(
    T           *out,
    const T     symbol,
    const int   length
)
{
    /*
     * Each thread processes eight elements, this is the first index
     */
    unsigned int idx = blockIdx.x * (blockDim.x << 3) + threadIdx.x;

#pragma unroll
    for (unsigned int i = 0; i < SCAN_ELTS_PER_THREAD; ++i)
    {
	if (idx < length)
	    out[idx] = symbol;

	idx += blockDim.x;
    }
}


template <typename T>
void
replicate
(
    T           *out,
    const T     symbol,
    const int   length
)
{
    unsigned int threads;
    unsigned int blocks;

    replicate_control(length, blocks, threads);
    replicate_core< T ><<<blocks,threads>>>(out, symbol, length);
}


// -----------------------------------------------------------------------------
// Instances
// -----------------------------------------------------------------------------

void replicate_ui(unsigned int *out, unsigned int x, int N)
{
    replicate<unsigned int>(out, x, N);
}

