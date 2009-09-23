/*
 * Module    : Prelude
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 */

#include "kernels.h"
#include "operator.h"
#include "shared_mem.h"
#include "cudpp/cudpp_globals.h"


template <typename T>
static void
groupBy_control
(
    unsigned int        n,
    unsigned int        &blocks,
    unsigned int        &threads,
    unsigned int        &smem
)
{
    blocks  = max(1.0, ceil((double)n / (SCAN_ELTS_PER_THREAD * CTA_SIZE)));
    threads = blocks > 1 ? CTA_SIZE : ceil((double)n / SCAN_ELTS_PER_THREAD);
    smem    = threads * sizeof(T);
}


/*
 * Apply a function each element of an array. A single thread is used to compute
 * each result. The input and output array may be coincident.
 */
template <class eq, typename T>
__global__ static void
groupBy_core
(
    T                   *xs,
    unsigned int        *flags,
    int                 length
)
{
    /*
     * Store a section of the array, to quickly check neighbouring elements.
     * This has enough for (blockDim.x + 1) elements.
     */
    SharedMemory<T> smem;
    T* sdata = smem.getPointer();

    /*
     * Each thread processes eight elements, this is the first index
     */
    unsigned int idx = blockIdx.x * (blockDim.x << 3) + threadIdx.x;

    /*
     * Fist pass, special handling for the neighbouring element of thread zero
     */
    sdata[threadIdx.x + 1] = xs[idx];
    __syncthreads();

    if (threadIdx.x == 0)
    {
        if (idx == 0) flags[idx] = 1;
        else          flags[idx] = !eq::apply(xs[idx-1], sdata[threadIdx.x + 1]);
    }
    else
    {
        flags[idx] = !eq::apply(sdata[threadIdx.x], sdata[threadIdx.x + 1]);
    }
    idx += blockDim.x;

    /*
     * Each thread processes an additional seven elements
     */
#pragma unroll
    for (unsigned int i = 0; i < SCAN_ELTS_PER_THREAD - 1; ++i)
    {
        if (threadIdx.x == 0)
            sdata[0] = sdata[blockDim.x];

	if (idx < length)
        {
            __syncthreads();
            sdata[threadIdx.x + 1] = xs[idx];

            __syncthreads();
            flags[idx] = !eq::apply(sdata[threadIdx.x], sdata[threadIdx.x + 1]);
        }

	idx += blockDim.x;
    }
}


template <class eq, typename T>
void
groupBy
(
    T                   *xs,
    unsigned int        *flags,
    int                 length
)
{
    unsigned int threads;
    unsigned int blocks;
    unsigned int smem;

    groupBy_control<T>(length, blocks, threads, smem);
    groupBy_core< eq,T ><<<blocks,threads,smem>>>(xs, flags, length);
}


// -----------------------------------------------------------------------------
// Instances
// -----------------------------------------------------------------------------

void group_f(float *xs, unsigned int *flags, int length)
{
    groupBy< Eq<float>,float >(xs, flags, length);
}

