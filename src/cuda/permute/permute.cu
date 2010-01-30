/* -----------------------------------------------------------------------------
 *
 * Module    : Permute
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#include "utils.h"
#include "permute.h"
#include "algorithms.h"

#include <stdint.h>


static void
permute_control(uint32_t n, uint32_t &blocks, uint32_t &threads)
{
    threads = min(ceilPow2(n), MAX_THREADS);
    blocks  = (n + threads - 1) / threads;
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
template <typename T, bool backward, bool compact, bool lengthIsPow2>
__global__ static void
permute_core
(
    const T             *d_in,
    T                   *d_out,
    const uint32_t      *indices,
    const uint32_t      length,
    const uint32_t      *valid     = NULL,
    uint32_t            *num_valid = NULL
)
{
    const uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;

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

    if (lengthIsPow2 || idx < length)
    {
	if (compact && valid[idx])
	    d_out[indices[idx]] = d_in[idx];
	else if (backward)
	    d_out[idx] = d_in[indices[idx]];
	else
	    d_out[indices[idx]] = d_in[idx];
    }
}


template <typename T, bool backward, bool compact>
static void
permute
(
    const T             *d_in,
    T                   *d_out,
    const uint32_t      *indices,
    const uint32_t      length,
    const uint32_t      *valid     = NULL,
    uint32_t            *num_valid = NULL
)
{
    uint32_t threads;
    uint32_t blocks;

    permute_control(length, blocks, threads);

    if (isPow2(length)) permute_core< T,backward,compact,true  ><<<blocks,threads>>>(d_in, d_out, indices, length, valid, num_valid);
    else                permute_core< T,backward,compact,false ><<<blocks,threads>>>(d_in, d_out, indices, length, valid, num_valid);
}


template <typename T, bool backward>
static unsigned int
compact
(
    const T             *d_in,
    T                   *d_out,
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
    permute<T,backward,true>(d_in, d_out, indices, length, d_flags, num);

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
    permute<uint32_t,false,false>((const uint32_t*) d_in, (uint32_t*) d_out, d_indices, length);
}

void bpermute(const void *d_in, void *d_out, const uint32_t *d_indices, const uint32_t length)
{
    permute<uint32_t,true,false>((const uint32_t*) d_in, (uint32_t*) d_out, d_indices, length);
}

unsigned int compact(const void *d_in, void *d_out, const uint32_t *d_flags, const uint32_t length)
{
    unsigned int N = compact<uint32_t,false>((const uint32_t*) d_in, (uint32_t*) d_out, d_flags, length);
    return N;
}

