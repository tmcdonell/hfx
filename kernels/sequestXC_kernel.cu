/*
 * Module    : Sequest
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 */


template <unsigned int blockSize, bool lengthIsPow2>
__global__ void
xcorr_dotP
(
    float        *specExp,
    int          *specThry,
    float        *dotP,
    unsigned int length
)
{
    __shared__ float scratch[];

    /*
     * Perform first level of reduction, reading from global memory and writing
     * to shared memory
     */
    unsigned int i        = blockIdx.x * blockSize * 2 + threadIdx.x;
    unsigned int tid      = threadIdx.x;
    unsigned int gridSize = blockSize * 2 * gridDim.x;

    scratch[tid]          = 0;

    /*
     * Reduce multiple elements per thread. The number is determined by the
     * number of active thread blocks (via gridDim). More blocks will result in
     * a larger `gridSize', and hence fewer elements per thread
     */
    while (i < length)
    {
        scratch[tid] += specExp[i] * specThry[i];

        /*
         * Ensure we don't read out of bounds. This is optimised away if the
         * input length is a power of two
         */
        if (lengthIsPow2 || i + blockSize < n)
        {
            scratch[tid] += specExp[i+blockSize] * specThry[i+blockSize];
        }

        /*
         * Loop stride to maintain coalescing
         */
        i += gridSize;
    }
    __syncthreads();

    /*
     * Now, calculate the reduction in shared memory
     */
    if (blockSize >= 512) { if (tid < 256) { scratch[tid] += scratch[tid+256]; } __syncthreads(); }
    if (blockSize >= 256) { if (tid < 128) { scratch[tid] += scratch[tid+128]; } __syncthreads(); }
    if (blockSize >= 128) { if (tid <  64) { scratch[tid] += scratch[tid+ 64]; } __syncthreads(); }

    if (tid < 32)
    {
        if (blockSize >= 64) { scratch[tid] += scratch[tid+32];  __syncthreads(); }
        if (blockSize >= 32) { scratch[tid] += scratch[tid+16];  __syncthreads(); }
        if (blockSize >= 16) { scratch[tid] += scratch[tid+ 8];  __syncthreads(); }
        if (blockSize >=  8) { scratch[tid] += scratch[tid+ 4];  __syncthreads(); }
        if (blockSize >=  4) { scratch[tid] += scratch[tid+ 2];  __syncthreads(); }
        if (blockSize >=  2) { scratch[tid] += scratch[tid+ 1];  __syncthreads(); }
    }

    /*
     * Write the results of this block back to global memory
     */
    if (tid == 0)
        dotP[blockIdx.x] = scratch[0];
}

