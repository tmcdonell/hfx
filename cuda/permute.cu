/*
 * Module    : Prelude
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 */

#include "utils.h"
#include "kernels.h"
#include "cudpp/cudpp_globals.h"


static void
permute_control
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
 * Permute an array according to the permutation indices. This handles both
 * forward and backward permutation, where:
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
 * scan.
 *
 * We return the number of valid elements found via `num_valid', but blunder on
 * ahead regardless of whether the `out' array is large enough or not.
 */
template <typename T, bool backward, bool compact>
__global__ static void
permute_core
(
    const T             *in,
    T                   *out,
    unsigned int        *indices,
    unsigned int        length,
    unsigned int        *valid     = NULL,
    unsigned int        *num_valid = NULL
)
{
    /*
     * Each thread processes eight elements, this is the index of the first
     */
    unsigned int idx = blockIdx.x * (blockDim.x << 3) + threadIdx.x;

    /*
     * Return the number of valid entries found
     */
    if (compact && threadIdx.x == 0)
    {
        if (backward)
            num_valid[0] = valid[0] + indices[0];
        else
            num_valid[0] = valid[length-1] + indices[length-1];
    }

#pragma unroll
    for (unsigned int i = 0; i < SCAN_ELTS_PER_THREAD; ++i)
    {
        if (idx < length)
        {
            if (compact && valid[idx])
                out[indices[idx]] = in[idx];
            else if (backward)
                out[idx] = in[indices[idx]];
            else
                out[indices[idx]] = in[idx];
        }
        idx += blockDim.x;
    }
}


template <typename T, bool backward, bool compact>
static void
permute
(
    const T             *in,
    T                   *out,
    unsigned int        *indices,
    unsigned int        length,
    unsigned int        *valid     = NULL,
    unsigned int        *num_valid = NULL
)
{
    unsigned int threads;
    unsigned int blocks;

    permute_control(length, blocks, threads);
    permute_core< T,backward,compact >
        <<<blocks,threads>>>(in, out, indices, length, valid, num_valid);
}


template <typename T, bool backward>
static int
compact
(
    const T             *in,
    T                   *out,
    unsigned int        *flags,
    unsigned int        length
)
{
    unsigned int N;
    unsigned int *num;
    unsigned int *indices;

    cudaMalloc(&num, sizeof(unsigned int));
    cudaMalloc(&indices, length * sizeof(unsigned int));

    if (backward) scanr_plusui(flags, indices, length);
    else          scanl_plusui(flags, indices, length);

    /*
     * At this point, we know exactly how many elements will be required for the
     * `out' array. Maybe we should allocate that array here?
     */
    permute<T,backward,true>(in, out, indices, length, flags, num);

    cudaMemcpy(&N, num, sizeof(unsigned int), cudaMemcpyDeviceToHost);
    cudaFree(num);
    cudaFree(indices);

    return N;
}


// -----------------------------------------------------------------------------
// Instances
// -----------------------------------------------------------------------------

void permute_i(int *in, int *out, unsigned int *indices, int length)
{
    permute<int,false,false>(in, out, indices, length);
}


void bpermute_f(float *in, float *out, unsigned int *indices, int length)
{
    permute<float,true,false>(in, out, indices, length);
}


int  compact_f(float *in, float *out, unsigned int *flags, int length)
{
    int N = compact<float,false>(in, out, flags, length);
    return N;
}

