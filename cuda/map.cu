/*
 * Module    : Prelude
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 */

#include "mass.h"
#include "kernels.h"
#include "operator.h"
#include "cudpp/cudpp_globals.h"

static void
map_control
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
template <class fn, typename Ta, typename Tb>
__global__ static void
map_core
(
    Ta  *xs,
    Tb  *out,
    int length
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
	    out[idx] = fn::apply(xs[idx]);

	idx += blockDim.x;
    }
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
    unsigned int threads;
    unsigned int blocks;

    map_control(length, blocks, threads);
    map_core< fn,Ta,Tb ><<<blocks,threads>>>(xs, out, length);
}


// -----------------------------------------------------------------------------
// Instances
// -----------------------------------------------------------------------------

void map_getAAMass(char *ions, float *masses, int N)
{
    map< getAAMass<char,float> >(ions, masses, N);
}

