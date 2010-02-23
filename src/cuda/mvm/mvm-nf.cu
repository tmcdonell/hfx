/* -----------------------------------------------------------------------------
 *
 * Module    : MVM
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 *
 * Noriyuki Fujimoto, "Faster Matrix-Vector Multiplication on GeForce 8800GTX"
 * Proceedings of the IEEE Parallel & Distributed Processing Symposium 2008
 *
 * ---------------------------------------------------------------------------*/

#include "utils.h"
#include "algorithms.h"

#include <stdint.h>
#include <cuda_runtime.h>

texture<uint4, 2, cudaReadModeElementType> texRef_A;

__global__
static void
mvm_core
(
    float               *d_y,
    const uint32_t      *d_A,
    const float         *d_x,
    const uint32_t      rows,
    const uint32_t      cols
)
{
    __shared__ float    xs[16][16];
    __shared__ float    Ps[16][16];

    uint4               a;
    float               *Psptr = (float*) Ps + (threadIdx.y << 4) + threadIdx.x;
    uint32_t            ay     = (blockIdx.x << 4) + threadIdx.y;
    float               *xptr  = (float*) d_x + (threadIdx.y << 4) + threadIdx.x;
    float               *xsptr = (float*) xs + (threadIdx.x << 2);

    *Psptr = 0.0f;
    uint32_t i;
    for (i = 0; i < (cols & ~255); i += 256, xptr += 256)
    {
        xs[threadIdx.y][threadIdx.x] = *xptr;
        __syncthreads();

        int ax = threadIdx.x + (i >> 2);
        a = tex2D(texRef_A, ax,    ay);    *Psptr += a.x*xsptr[  0] + a.y*xsptr[  1] + a.z*xsptr[  2] + a.w*xsptr[  3];
        a = tex2D(texRef_A, ax+16, ay);    *Psptr += a.x*xsptr[ 64] + a.y*xsptr[ 65] + a.z*xsptr[ 66] + a.w*xsptr[ 67];
        a = tex2D(texRef_A, ax+32, ay);    *Psptr += a.x*xsptr[128] + a.y*xsptr[129] + a.z*xsptr[130] + a.w*xsptr[131];
        a = tex2D(texRef_A, ax+48, ay);    *Psptr += a.x*xsptr[192] + a.y*xsptr[193] + a.z*xsptr[194] + a.w*xsptr[195];
        __syncthreads();
    }

    if (i + (threadIdx.y << 4) + threadIdx.x < cols)
    {
        xs[threadIdx.y][threadIdx.x] = *xptr;
    }
    __syncthreads();

    uint32_t j;
    for (j = 0; j < ((cols - i) >> 6); j++, xsptr += 62)
    {
        a = tex2D(texRef_A, threadIdx.x + (i >> 2) + (j << 4), ay);
        *Psptr += a.x* *xsptr++ + a.y* *xsptr++ + a.z* *xsptr++ + a.w* *xsptr++;
    }
    __syncthreads();

    uint32_t remain = (cols-i) & 63;
    if ((threadIdx.x << 2) < remain)
    {
        a = tex2D(texRef_A, threadIdx.x + (i >> 2) + (j << 4), ay);
        *Psptr += a.x* *xsptr++;
    }

    if ((threadIdx.x << 2) + 1 < remain) *Psptr += a.y* *xsptr++;
    if ((threadIdx.x << 2) + 2 < remain) *Psptr += a.z* *xsptr++;
    if ((threadIdx.x << 2) + 3 < remain) *Psptr += a.w* *xsptr;
    __syncthreads();

    if (threadIdx.x < 8) *Psptr += *(Psptr + 8);
    if (threadIdx.x < 4) *Psptr += *(Psptr + 4);
    if (threadIdx.x < 2) *Psptr += *(Psptr + 2);
    if (threadIdx.x < 1) *Psptr += *(Psptr + 1);
    __syncthreads();

    if (threadIdx.y == 0 && (blockIdx.x << 4) + threadIdx.x < rows)
        d_y[(blockIdx.x << 4) + threadIdx.x] = Ps[threadIdx.x][0];
}


void
mvm_if(float *d_y, const uint32_t *d_A, const float *d_x, const uint32_t m, const uint32_t n)
{
    uint32_t blkNum = (m >> 4) + ((m & 15) ? 1 : 0);
    uint32_t height = blkNum << 4;
    uint32_t width  = (n & 255) ? (256 * ((n >> 8) + 1)) : n;

    dim3 threads(16,16);
    dim3 blocks(blkNum,1);

    size_t                offset = size_t (-1);
    cudaChannelFormatDesc desci4 = cudaCreateChannelDesc<uint4>();

    cudaBindTexture2D(&offset, texRef_A, d_A, desci4, width >> 2, height, width*sizeof(uint32_t));
    assert(offset == 0 || !"memory is not aligned, cannot use texture cache");

    mvm_core<<<blocks,threads>>>(d_y, d_A, d_x, m, n);

    cudaUnbindTexture(texRef_A);
}

