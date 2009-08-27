/*
 * Module    : Prelude
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 */

#ifndef __FOLD_KERNEL__
#define __FOLD_KERNEL__

#define MAX_THREADS     128
#define MAX_BLOCKS      64

#include "utils.h"
#include "operator.h"
#include "shared_mem.h"


/*
 * Compute multiple elements per thread sequentially. This reduces the overall
 * cost of the algorithm while keeping the work complexity O(n) and the step
 * complexity O(log n). c.f. Brent's Theorem optimisation.
 *
 * Stolen from the CUDA SDK examples
 */
template <unsigned int blockSize, bool lengthIsPow2, class BinaryOp, typename T>
__global__ void
fold_recursive
(
    const T     *d_xs,
    T           *d_ys,
    int         length
)
{
    SharedMemory<T> smem;
    T *scratch = smem.getPointer();

    /*
     * Calculate first level of reduction reading into shared memory
     */
    unsigned int i;
    unsigned int tid      = threadIdx.x;
    unsigned int gridSize = blockSize * 2 * gridDim.x;

    scratch[tid] = BinaryOp::identity();

    /*
     * Reduce multiple elements per thread. The number is determined by the
     * number of active thread blocks (via gridDim). More blocks will result in
     * a larger `gridSize', and hence fewer elements per thread
     *
     * The loop stride of `gridSize' is used to maintain coalescing.
     */
    for (i =  blockIdx.x * blockSize * 2 + tid; i <  length; i += gridSize)
    {
        scratch[tid] = BinaryOp::apply(scratch[tid], d_xs[i]);

        /*
         * Ensure we don't read out of bounds. This is optimised away if the
         * input length is a power of two
         */
        if (lengthIsPow2 || i + blockSize < length)
            scratch[tid] = BinaryOp::apply(scratch[tid], d_xs[i+blockSize]);
    }
    __syncthreads();

    /*
     * Now, calculate the reduction in shared memory
     */
    if (blockSize >= 512) { if (tid < 256) { scratch[tid] = BinaryOp::apply(scratch[tid], scratch[tid+256]); } __syncthreads(); }
    if (blockSize >= 256) { if (tid < 128) { scratch[tid] = BinaryOp::apply(scratch[tid], scratch[tid+128]); } __syncthreads(); }
    if (blockSize >= 128) { if (tid <  64) { scratch[tid] = BinaryOp::apply(scratch[tid], scratch[tid+ 64]); } __syncthreads(); }

#ifndef __DEVICE_EMULATION__
    if (tid < 32)
#endif
    {
        if (blockSize >= 64) { scratch[tid] = BinaryOp::apply(scratch[tid], scratch[tid+32]);  __EMUSYNC; }
        if (blockSize >= 32) { scratch[tid] = BinaryOp::apply(scratch[tid], scratch[tid+16]);  __EMUSYNC; }
        if (blockSize >= 16) { scratch[tid] = BinaryOp::apply(scratch[tid], scratch[tid+ 8]);  __EMUSYNC; }
        if (blockSize >=  8) { scratch[tid] = BinaryOp::apply(scratch[tid], scratch[tid+ 4]);  __EMUSYNC; }
        if (blockSize >=  4) { scratch[tid] = BinaryOp::apply(scratch[tid], scratch[tid+ 2]);  __EMUSYNC; }
        if (blockSize >=  2) { scratch[tid] = BinaryOp::apply(scratch[tid], scratch[tid+ 1]);  __EMUSYNC; }
    }

    /*
     * Write the results of this block back to global memory
     */
    if (tid == 0)
        d_ys[blockIdx.x] = scratch[0];
}


/*
 * Wrapper function for kernel launch
 */
template <class BinaryOp, typename T>
void
fold_dispatch
(
    const T     *d_xs,
    T           *d_ys,
    int         length,
    int         blocks,
    int         threads
)
{
    unsigned int smem = threads * sizeof(T);

    if (isPow2(length))
    {
        switch (threads)
        {
        case 512: fold_recursive<512,true,BinaryOp,T><<<blocks,threads,smem>>>(d_xs, d_ys, length); break;
        case 256: fold_recursive<256,true,BinaryOp,T><<<blocks,threads,smem>>>(d_xs, d_ys, length); break;
        case 128: fold_recursive<128,true,BinaryOp,T><<<blocks,threads,smem>>>(d_xs, d_ys, length); break;
        case  64: fold_recursive< 64,true,BinaryOp,T><<<blocks,threads,smem>>>(d_xs, d_ys, length); break;
        case  32: fold_recursive< 32,true,BinaryOp,T><<<blocks,threads,smem>>>(d_xs, d_ys, length); break;
        case  16: fold_recursive< 16,true,BinaryOp,T><<<blocks,threads,smem>>>(d_xs, d_ys, length); break;
        case   8: fold_recursive<  8,true,BinaryOp,T><<<blocks,threads,smem>>>(d_xs, d_ys, length); break;
        case   4: fold_recursive<  4,true,BinaryOp,T><<<blocks,threads,smem>>>(d_xs, d_ys, length); break;
        case   2: fold_recursive<  2,true,BinaryOp,T><<<blocks,threads,smem>>>(d_xs, d_ys, length); break;
        case   1: fold_recursive<  1,true,BinaryOp,T><<<blocks,threads,smem>>>(d_xs, d_ys, length); break;
        default:
            assert(!"Non-exhaustive patterns in match");
        }
    }
    else
    {
        switch (threads)
        {
        case 512: fold_recursive<512,false,BinaryOp,T><<<blocks,threads,smem>>>(d_xs, d_ys, length); break;
        case 256: fold_recursive<256,false,BinaryOp,T><<<blocks,threads,smem>>>(d_xs, d_ys, length); break;
        case 128: fold_recursive<128,false,BinaryOp,T><<<blocks,threads,smem>>>(d_xs, d_ys, length); break;
        case  64: fold_recursive< 64,false,BinaryOp,T><<<blocks,threads,smem>>>(d_xs, d_ys, length); break;
        case  32: fold_recursive< 32,false,BinaryOp,T><<<blocks,threads,smem>>>(d_xs, d_ys, length); break;
        case  16: fold_recursive< 16,false,BinaryOp,T><<<blocks,threads,smem>>>(d_xs, d_ys, length); break;
        case   8: fold_recursive<  8,false,BinaryOp,T><<<blocks,threads,smem>>>(d_xs, d_ys, length); break;
        case   4: fold_recursive<  4,false,BinaryOp,T><<<blocks,threads,smem>>>(d_xs, d_ys, length); break;
        case   2: fold_recursive<  2,false,BinaryOp,T><<<blocks,threads,smem>>>(d_xs, d_ys, length); break;
        case   1: fold_recursive<  1,false,BinaryOp,T><<<blocks,threads,smem>>>(d_xs, d_ys, length); break;
        default:
            assert(!"Non-exhaustive patterns in match");
        }
    }
}


/*
 * Compute the number of blocks and threads to use for the reduction kernel
 */
void
fold_control
(
    int         n,
    int         &blocks,
    int         &threads,
    int         maxThreads = MAX_THREADS,
    int         maxBlocks  = MAX_BLOCKS
)
{
    threads = (n < maxThreads*2) ? ceilPow2((n+1)/2) : maxThreads;
    blocks  = (n + threads * 2 - 1) / (threads * 2);
    blocks  = min(blocks, maxBlocks);
}


/*
 * Apply a binary operator to an array, reducing the array to a single value.
 * The reduction will take place in parallel, so the operator must be
 * associative.
 */
template <class BinaryOp, typename T>
T
fold
(
    const T     *d_xs,
    int         n
)
{
    int blocks;
    int threads;
    T   gpu_result;
    T*  d_data          = NULL;

    /*
     * Allocate temporary storage for the block-level reduction
     */
    fold_control(n, blocks, threads);
    cudaMalloc((void **) &d_data, sizeof(T) * blocks);

    /*
     * Recursively reduce the partial block sums to a single value
     */
    fold_dispatch<BinaryOp,T>(d_xs, d_data, n, blocks, threads);

    n = blocks;
    while (n > 1)
    {
        fold_control(n, blocks, threads);
        fold_dispatch<BinaryOp,T>(d_data, d_data, n, blocks, threads);

        n = (n + threads * 2 - 1) / (threads * 2);
    }
    assert(n == 1);

    /*
     * Read back the final result
     */
    cudaMemcpy(&gpu_result, d_data, sizeof(T), cudaMemcpyDeviceToHost);
    cudaFree(d_data);

    return gpu_result;
}

#endif

