/* -----------------------------------------------------------------------------
 *
 * Module    : MVM
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#include "mvm.h"
#include "utils.h"
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

    float               sum    = 0;

    /*
     * Loop over sub-blocks of the matrix A across its width. At the end, we
     * will have a square block of partial sums for each row.
     */
    for (uint32_t j = 0; j < full;  j += BLOCKDIM_X * BLOCKDIM_Y)
    {
        s_data[tid] = d_x[tid + j];
        __syncthreads();

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
        s_data[tid] = d_x[tid + full];
        __syncthreads();

        for (uint32_t i = threadIdx.x; i < (cols - full); i += BLOCKDIM_X)
        {
            sum+= d_row[full+i] * s_data[i];
        }
        __syncthreads();
    }

    /*
     * Reduce the partial sums to the final result vector. No synchronisation is
     * required since we support blockDim.x <= 64 only, meaning warps process an
     * entire row(s) in lockstep.
     */
    s_data[tid] = sum;                                                  __EMUSYNC;
    if (64 <= BLOCKDIM_X) { s_data[tid] = sum = sum + s_data[tid + 32]; __EMUSYNC; }
    if (32 <= BLOCKDIM_X) { s_data[tid] = sum = sum + s_data[tid + 16]; __EMUSYNC; }
    if (16 <= BLOCKDIM_X) { s_data[tid] = sum = sum + s_data[tid +  8]; __EMUSYNC; }
    if ( 8 <= BLOCKDIM_X) { s_data[tid] = sum = sum + s_data[tid +  4]; __EMUSYNC; }
    if ( 4 <= BLOCKDIM_X) { s_data[tid] = sum = sum + s_data[tid +  2]; __EMUSYNC; }
    if ( 2 <= BLOCKDIM_X) { s_data[tid] = sum = sum + s_data[tid +  1]; __EMUSYNC; }

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

    mvm_control(rows, cols, blocks, threads);
    mvm_core<<<blocks,threads>>>(d_y, d_A, d_x, rows, cols);
}


/* -----------------------------------------------------------------------------
 * Instances
 * ---------------------------------------------------------------------------*/

void
mvm_if(float *d_y, const uint32_t *d_A, const float *d_x, const uint32_t m, const uint32_t n)
{
    mvm(d_y, d_A, d_x, m, n);
}

