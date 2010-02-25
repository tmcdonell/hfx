/* -----------------------------------------------------------------------------
 *
 * Module    : Ion Series
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#include "utils.h"
#include "ion_series.h"
#include "algorithms.h"

#include <stdint.h>


__inline__ __device__ static float
ionMZ(const float m, const float c)
{
    return __fdividef(m + MASS_H * c, c);
}

__inline__ __device__ static uint32_t
binMZ(const float mz)
{
    return rintf(__fdividef(mz, BIN_WIDTH_MONO));
}

__inline__ __device__ static void
addIon(uint32_t *d_spec, const uint32_t N, const int32_t x, const uint32_t y)
{
    if (0 <= x && x < N) atomicMax(&d_spec[x], y);
}


template <uint32_t charge>
__device__ void
addIonsAB(uint32_t *d_spec, const uint32_t N, const float mass)
{
    float   m;
    int32_t x;

    // A-ions
    addIon(d_spec, N, binMZ(ionMZ(mass - MASS_CO, charge)), 10);

    // B-ions
    m = ionMZ(mass, charge);
    x = binMZ(m);

    addIon(d_spec, N, x,   50);
    addIon(d_spec, N, x+1, 25); // technically, should be binMZ(m+1)
    addIon(d_spec, N, x-1, 25);

    addIon(d_spec, N, binMZ(m - __fdividef(MASS_H2O, charge)), 10);
    addIon(d_spec, N, binMZ(m - __fdividef(MASS_NH3, charge)), 10);
}


template <uint32_t charge>
__device__ void
addIonsY(uint32_t *d_spec, const uint32_t N, const float mass)
{
    float   m = ionMZ(mass + MASS_H2O, charge);
    int32_t x = binMZ(m);

    // Y-ions
    addIon(d_spec, N, x,   50);
    addIon(d_spec, N, x+1, 25);
    addIon(d_spec, N, x-1, 25);

    addIon(d_spec, N, binMZ(m - __fdividef(MASS_NH3, charge)), 10);
}


template <uint32_t charge>
__device__ void
addIons_k(uint32_t *d_spec, const uint32_t N, const float b_mass, const float y_mass)
{
    addIonsAB<charge>(d_spec, N, b_mass);
    addIonsY <charge>(d_spec, N, y_mass);
}


/*
 * Generate theoretical spectra for a collection of peptide fragments. The
 * yIonLadder array contains data for all fragments in the database, although we
 * are only interested in those beginning at the inRangeIdx positions.
 *
 * A warp of threads iterates over the fragment masses for a peptide, issuing a
 * (long) sequence of (slow) global atomic update requests. The input spectra
 * matrix should be initially zero, and on output will contain the theoretical
 * spectral peaks in a square (but mostly sparse) matrix.
 */
template <uint32_t BlockSize, uint32_t MaxCharge>
__global__ static void
addIons_core
(
    uint32_t            *d_spec,
    const float         *d_residual,
    const float         *d_yIonLadder,
    const uint32_t      *d_rowPtr,
    const uint32_t      *d_inRangeIdx,
    const uint32_t      num_inRange,
    const uint32_t      len_spec
)
{
    assert(BlockSize % WARP_SIZE == 0);

    const uint32_t vectorsPerBlock = BlockSize / WARP_SIZE;
    const uint32_t numVectors      = vectorsPerBlock * gridDim.x;
    const uint32_t thread_id       = BlockSize * blockIdx.x + threadIdx.x;
    const uint32_t vector_id       = thread_id / WARP_SIZE;
    const uint32_t thread_lane     = threadIdx.x & (WARP_SIZE-1);
    const uint32_t vector_lane     = threadIdx.x / WARP_SIZE;

    __shared__ volatile uint32_t s_ptrs[vectorsPerBlock][2];

    for (uint32_t row = vector_id; row < num_inRange; row += numVectors)
    {
        const uint32_t idx      = d_inRangeIdx[row];
        const float    residual = d_residual[idx];
        uint32_t       *spec    = &d_spec[row * len_spec];

        /*
         * Use two threads to fetch the indices of the start and end of this
         * segment. This is a single coalesced (unaligned) global read.
         */
        if (thread_lane < 2)
            s_ptrs[vector_lane][thread_lane] = d_rowPtr[idx + thread_lane];

        __EMUSYNC;
        const uint32_t row_start = s_ptrs[vector_lane][0];
        const uint32_t row_end   = s_ptrs[vector_lane][1];

        /*
         * Have all threads read in values for this segment, writing the
         * spectral peaks out to global memory (very, very slowly...)
         */
        for (uint32_t j = row_start + thread_lane; j < row_end; j += WARP_SIZE)
        {
            const float y_mass = d_yIonLadder[j];
            const float b_mass = residual - y_mass;

            if (1 <= MaxCharge) addIons_k<1>(spec, len_spec, b_mass, y_mass);
            if (2 <= MaxCharge) addIons_k<2>(spec, len_spec, b_mass, y_mass);
            if (3 <= MaxCharge) addIons_k<3>(spec, len_spec, b_mass, y_mass);
            if (4 <= MaxCharge) addIons_k<4>(spec, len_spec, b_mass, y_mass);
        }
    }
}

/*
 * Select a number of threads and blocks. Each block will have at least one full
 * warp, as required by the core kernel
 */
static void
addIons_control(uint32_t N, uint32_t &blocks, uint32_t &threads)
{
    threads = (N < MAX_THREADS) ? max(WARP_SIZE, ceilPow2(N)) : MAX_THREADS;
    blocks  = (N + threads - 1) / threads;
    blocks  = min(blocks, MAX_BLOCKS);
}


template <uint32_t MaxCharge>
static void
addIons_dispatch
(
    uint32_t            *d_spec,
    const float         *d_residual,
    const float         *d_ladder,
    const uint32_t      *d_rowPtr,
    const uint32_t      *d_inRangeIdx,
    const uint32_t      num_inRange,
    const uint32_t      len_spec
)
{
    uint32_t blocks;
    uint32_t threads;

    addIons_control(num_inRange, blocks, threads);
    switch (threads)
    {
    case 512: addIons_core<512,MaxCharge><<<blocks,threads>>>(d_spec, d_residual, d_ladder, d_rowPtr, d_inRangeIdx, num_inRange, len_spec); break;
    case 256: addIons_core<256,MaxCharge><<<blocks,threads>>>(d_spec, d_residual, d_ladder, d_rowPtr, d_inRangeIdx, num_inRange, len_spec); break;
    case 128: addIons_core<128,MaxCharge><<<blocks,threads>>>(d_spec, d_residual, d_ladder, d_rowPtr, d_inRangeIdx, num_inRange, len_spec); break;
    case  64: addIons_core< 64,MaxCharge><<<blocks,threads>>>(d_spec, d_residual, d_ladder, d_rowPtr, d_inRangeIdx, num_inRange, len_spec); break;
    case  32: addIons_core< 32,MaxCharge><<<blocks,threads>>>(d_spec, d_residual, d_ladder, d_rowPtr, d_inRangeIdx, num_inRange, len_spec); break;
    default:
        assert(!"Non-exhaustive patterns in match");
    }
}


void
addIons
(
    uint32_t            *d_spec,
    const float         *d_residual,
    const float         *d_ladder,
    const uint32_t      *d_rowPtr,
    const uint32_t      *d_inRangeIdx,
    const uint32_t      num_inRange,
    const uint32_t      max_charge,
    const uint32_t      len_spec
)
{
    switch (max_charge)
    {
    case 1: addIons_dispatch<1>(d_spec, d_residual, d_ladder, d_rowPtr, d_inRangeIdx, num_inRange, len_spec); break;
    case 2: addIons_dispatch<2>(d_spec, d_residual, d_ladder, d_rowPtr, d_inRangeIdx, num_inRange, len_spec); break;
    case 3: addIons_dispatch<3>(d_spec, d_residual, d_ladder, d_rowPtr, d_inRangeIdx, num_inRange, len_spec); break;
    case 4: addIons_dispatch<4>(d_spec, d_residual, d_ladder, d_rowPtr, d_inRangeIdx, num_inRange, len_spec); break;
    default:
        assert(!"Non-exhaustive patterns in match");
    }
}


