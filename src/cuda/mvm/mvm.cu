/* -----------------------------------------------------------------------------
 *
 * Module    : MVM
 * Copyright : (c) [2009..2010] Trevor L. McDonell
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#include "mvm.h"
#include "utils.h"
#include "texture.h"
#include "algorithms.h"

#include <stdint.h>

#if 0
__device__ static float
dotp_f4(const float4 &a, const float4 &b)
{
    return (a.w*b.w + a.x*b.x + a.y*b.y + a.z*b.z);
}
#endif


/*
 * Matrix-vector multiplication: y = A*x
 * The matrix is stored in row-major order.
 */
template <bool UseCache>
__global__ static void
mvm_core
(
    float               *d_y,
    const uint32_t      *d_A,
    const float         *d_x,
    const uint32_t      rows,
    const uint32_t      cols
)
{
    __shared__ float    s_data[BLOCKDIM_X * BLOCKDIM_Y];

    const uint32_t      tid    = threadIdx.x + threadIdx.y * blockDim.x;
    const uint32_t      rowIdx = threadIdx.y + blockIdx.y  * blockDim.y;
    const uint32_t      full   = cols & ~(blockDim.x * blockDim.y - 1);
    const uint32_t      *d_row = &d_A[cols * rowIdx];

    float               sum    = 0.0f;

    /*
     * Loop over sub-blocks of the matrix A across its width. At the end, we
     * will have a square block of partial sums for each row.
     */
    for (uint32_t j = 0; j < full;  j += BLOCKDIM_X * BLOCKDIM_Y)
    {
        s_data[tid] = fetch_x<UseCache>(tid + j, d_x);
        __syncthreads();

        /*
         * All threads in this block must remain active to share loading chunks
         * of the x-vector, but we still need to be careful that they don't read
         * out-of-bounds
         */
        if (rowIdx < rows)
        {
            /*
             * Every thread fetches data from the x-vector, so we can step over
             * several sub-blocks.
             */
#pragma unroll
            for (uint32_t i = threadIdx.x; i < BLOCKDIM_X * BLOCKDIM_Y; i += BLOCKDIM_X)
            {
                /*
                 * TODO: change BLOCKDIM_X to 8 and read d_A as float4 ~~> 128bytes/row
                 */
                sum += d_row[j+i] * s_data[i];
            }
        }

        /*
         * Make sure all threads in the block have finished processing this
         * chunk of the x-vector.
         */
        __syncthreads();
    }

    /*
     * Add remaining chunks, when cols is not evenly divisible by our chunk size
     */
    if (full != cols)
    {
        if (tid + full < cols)
            s_data[tid] = fetch_x<UseCache>(tid + full, d_x);

        __syncthreads();

        if (rowIdx < rows)
        {
            for (uint32_t i = threadIdx.x;i < (cols - full); i += BLOCKDIM_X)
            {
                sum+= d_row[full+i] * s_data[i];
            }
        }
        __syncthreads();
    }

    /*
     * Reduce the partial sums to the final result vector. No synchronisation is
     * required since we support blockDim.x <= 64 only, meaning warps process an
     * entire row(s) in lockstep.
     */
#ifndef __DEVICE_EMULATION__
    s_data[tid] = sum;
    if (BLOCKDIM_X >= 64) s_data[tid] = sum = sum + s_data[tid + 32];
    if (BLOCKDIM_X >= 32) s_data[tid] = sum = sum + s_data[tid + 16];
    if (BLOCKDIM_X >= 16) s_data[tid] = sum = sum + s_data[tid +  8];
    if (BLOCKDIM_X >=  8) s_data[tid] = sum = sum + s_data[tid +  4];
    if (BLOCKDIM_X >=  4) s_data[tid] = sum = sum + s_data[tid +  2];
    if (BLOCKDIM_X >=  2) s_data[tid] = sum = sum + s_data[tid +  1];
#else
    s_data[tid] = sum;                                                             __EMUSYNC;
    if (BLOCKDIM_X >= 64 && threadIdx.x < 32) { s_data[tid] += s_data[tid + 32]; } __EMUSYNC;
    if (BLOCKDIM_X >= 32 && threadIdx.x < 16) { s_data[tid] += s_data[tid + 16]; } __EMUSYNC;
    if (BLOCKDIM_X >= 16 && threadIdx.x <  8) { s_data[tid] += s_data[tid +  8]; } __EMUSYNC;
    if (BLOCKDIM_X >=  8 && threadIdx.x <  4) { s_data[tid] += s_data[tid +  4]; } __EMUSYNC;
    if (BLOCKDIM_X >=  4 && threadIdx.x <  2) { s_data[tid] += s_data[tid +  2]; } __EMUSYNC;
    if (BLOCKDIM_X >=  2 && threadIdx.x <  1) { s_data[tid] += s_data[tid +  1]; } __EMUSYNC;
#endif

    /*
     * If the number of rows is not evenly divisible by the block size, these
     * threads still participate in loading the (shared) x-vector values, but
     * don't write their final result.
     */
    if (threadIdx.x == 0 && rowIdx < rows)
        d_y[rowIdx] = s_data[tid];
}


static void
mvm_control(const uint32_t m, const uint32_t n, dim3 &blocks, dim3 &threads)
{
    (void) n;

    blocks  = dim3(1, ((m + BLOCKDIM_Y - 1) / BLOCKDIM_Y), 1);
    threads = dim3(BLOCKDIM_X, BLOCKDIM_Y, 1);
}


template <bool UseCache>
static void
mvm
(
    float               *d_y,
    const uint32_t      *d_A,
    const float         *d_x,
    uint32_t            rows,
    uint32_t            cols
)
{
    dim3 blocks;
    dim3 threads;

    if (UseCache)
        bind_x(d_x);

    mvm_control(rows, cols, blocks, threads);
    mvm_core<UseCache><<<blocks,threads>>>(d_y, d_A, d_x, rows, cols);

    if (UseCache)
        unbind_x(d_x);
}


/* -----------------------------------------------------------------------------
 * Instances
 * ---------------------------------------------------------------------------*/

void
mvm_if(float *d_y, const uint32_t *d_A, const float *d_x, const uint32_t m, const uint32_t n)
{
    mvm<false>(d_y, d_A, d_x, m, n);
}

