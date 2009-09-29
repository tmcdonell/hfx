/*
 * Module    : Prelude
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 */

#include "utils.h"
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
    threads = min(ceilPow2(n), CTA_SIZE);
    blocks  = (n + threads - 1) / threads;
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
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < length)
	out[idx] = symbol;
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

