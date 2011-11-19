/* -----------------------------------------------------------------------------
 *
 * Module    : Ion Series
 * Copyright : (c) [2009..2010] Trevor L. McDonell
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#include "utils.h"
#include "device.h"
#include "texture.h"
#include "ion_series.h"
#include "algorithms.h"

#include <stdint.h>


/*
 * Scan a warp-sized chunk of data. Because warps execute instructions in SIMD
 * fashion, there is no need to synchronise in order to share data. The most
 * efficient algorithm is the step-efficient method of Hillis & Steele that
 * takes log(N) steps, rather than the work-efficient tree-based algorithm
 * described by Blelloch that takes 2 * log(N) steps.
 */
template <class T, bool inclusive>
static __device__ T
scan_warp(T val, volatile T* s_data)
{
    const uint32_t idx  = threadIdx.x;
    const uint32_t lane = threadIdx.x & (WARP_SIZE-1);

    /*
     * If we double the size of the s_data array and pad the bottom half with
     * zero, then we can avoid branching (although there is plenty already).
     *
     * In device emulation mode, the warp size is 1 and so sync-less operation
     * does not work.
     */
    s_data[idx] = val;                                                        __EMUSYNC;
#ifdef __DEVICE_EMULATION__
    val = (lane >=  1) ? s_data[idx -  1] : 0; __EMUSYNC; s_data[idx] += val; __EMUSYNC;
    val = (lane >=  2) ? s_data[idx -  2] : 0; __EMUSYNC; s_data[idx] += val; __EMUSYNC;
    val = (lane >=  4) ? s_data[idx -  4] : 0; __EMUSYNC; s_data[idx] += val; __EMUSYNC;
    val = (lane >=  8) ? s_data[idx -  8] : 0; __EMUSYNC; s_data[idx] += val; __EMUSYNC;
    val = (lane >= 16) ? s_data[idx - 16] : 0; __EMUSYNC; s_data[idx] += val; __EMUSYNC;
#else
    if (lane >=  1) s_data[idx] = val = val + s_data[idx -  1];
    if (lane >=  2) s_data[idx] = val = val + s_data[idx -  2];
    if (lane >=  4) s_data[idx] = val = val + s_data[idx -  4];
    if (lane >=  8) s_data[idx] = val = val + s_data[idx -  8];
    if (lane >= 16) s_data[idx] = val = val + s_data[idx - 16];
#endif

    if (inclusive) return s_data[idx];
    else           return (lane > 0) ? s_data[idx - 1] : 0;
}


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
 * Return the mass of an amino acid residue in atomic mass units, for the given
 * short abbreviation.
 */
template <bool UseCache>
__device__ float
getAAMass(const float *d_mass, const char aa)
{
    return fetch_x<UseCache>(aa - 'A', d_mass);
}


/*
 * Generate theoretical spectra for a collection of peptide fragments. The
 * 'ions' array contains the individual amino-acid masses for the database
 * entries. We are interested in the sequences generated between the terminal
 * indices (tc,tn) of the locations specified in the 'idx' array.
 *
 * A warp of threads iterates between the (tc,tn) indices, generating the b- and
 * y-ion mass ladders. A (long) sequence of (slow) global atomic update requests
 * is subsequently issued. The input d_spec should be initially zero, and on
 * output will contain the theoretical spectra peaks in a dense (although
 * mostly zero) matrix.
 */
template <uint32_t BlockSize, uint32_t MaxCharge, bool UseCache>
__global__ static void
addIons_core
(
    uint32_t            *d_spec,
    const float         *d_residual,    // peptide residual mass
    const float         *d_mass,        // lookup table for ion character codes ['A'..'Z']
    const uint8_t       *d_ions,        // individual ion character codes (the database)
    const uint32_t      *d_tc,          // c-terminal indices
    const uint32_t      *d_tn,          // n-terminal indices
    const uint32_t      *d_idx,         // The indices of the data for peptides under consideration
    const uint32_t      num_idx,
    const uint32_t      len_spec
)
{
    assert(BlockSize % WARP_SIZE == 0);

    const uint32_t vectorsPerBlock = BlockSize / WARP_SIZE;
    const uint32_t numVectors      = vectorsPerBlock * gridDim.x;
    const uint32_t thread_id       = BlockSize * blockIdx.x + threadIdx.x;
    const uint32_t vector_id       = thread_id / WARP_SIZE;
    const uint32_t thread_lane     = threadIdx.x & (WARP_SIZE-1);

    __shared__ volatile float s_data[BlockSize];

    for (uint32_t row = vector_id; row < num_idx; row += numVectors)
    {
        const uint32_t idx       = d_idx[row];
        const uint32_t row_start = d_tc[idx];
        const uint32_t row_end   = d_tn[idx];
        const float    residual  = d_residual[idx];

        uint32_t       *spec     = &d_spec[row * len_spec];
        float          b_mass;
        float          y_mass;

        s_data[threadIdx.x]      = 0;

        /*
         * Have all threads read in values for this segment, writing the
         * spectral peaks out to global memory (very, very slowly...)
         */
        for (uint32_t j = row_start + thread_lane; j < row_end; j += WARP_SIZE)
        {
            /*
             * Load the ion mass, and propagate the partial scan results
             */
            b_mass = getAAMass<UseCache>(d_mass, d_ions[j]);

            if (thread_lane == 0)
                b_mass += s_data[threadIdx.x + (WARP_SIZE-1)];

            /*
             * Generate fragment mass ladder
             */
            b_mass = scan_warp<float,true>(b_mass, s_data);
            y_mass = residual - b_mass;

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


template <uint32_t MaxCharge, bool UseCache>
static void
addIons_dispatch
(
    uint32_t            *d_spec,
    const float         *d_residual,
    const float         *d_mass,
    const uint8_t       *d_ions,
    const uint32_t      *d_tc,
    const uint32_t      *d_tn,
    const uint32_t      *d_idx,
    const uint32_t      num_idx,
    const uint32_t      len_spec
)
{
    uint32_t blocks;
    uint32_t threads;

    if (UseCache)
        bind_x(d_mass);

    addIons_control(num_idx, blocks, threads);
    switch (threads)
    {
//  case 512: addIons_core<512,MaxCharge,UseCache><<<blocks,threads>>>(d_spec, d_residual, d_mass, d_ions, d_tc, d_tn, d_idx, num_idx, len_spec); break;
//  case 256: addIons_core<256,MaxCharge,UseCache><<<blocks,threads>>>(d_spec, d_residual, d_mass, d_ions, d_tc, d_tn, d_idx, num_idx, len_spec); break;
    case 128: addIons_core<128,MaxCharge,UseCache><<<blocks,threads>>>(d_spec, d_residual, d_mass, d_ions, d_tc, d_tn, d_idx, num_idx, len_spec); break;
    case  64: addIons_core< 64,MaxCharge,UseCache><<<blocks,threads>>>(d_spec, d_residual, d_mass, d_ions, d_tc, d_tn, d_idx, num_idx, len_spec); break;
    case  32: addIons_core< 32,MaxCharge,UseCache><<<blocks,threads>>>(d_spec, d_residual, d_mass, d_ions, d_tc, d_tn, d_idx, num_idx, len_spec); break;
    default:
        assert(!"Non-exhaustive patterns in match");
    }

    if (UseCache)
      unbind_x(d_mass);
}


void
addIons
(
    uint32_t            *d_spec,
    const float         *d_residual,
    const float         *d_mass,
    const uint8_t       *d_ions,
    const uint32_t      *d_tc,
    const uint32_t      *d_tn,
    const uint32_t      *d_idx,
    const uint32_t      num_idx,
    const uint32_t      max_charge,
    const uint32_t      len_spec
)
{
    switch (max_charge)
    {
    case 1: addIons_dispatch<1,true>(d_spec, d_residual, d_mass, d_ions, d_tc, d_tn, d_idx, num_idx, len_spec); break;
    case 2: addIons_dispatch<2,true>(d_spec, d_residual, d_mass, d_ions, d_tc, d_tn, d_idx, num_idx, len_spec); break;
    case 3: addIons_dispatch<3,true>(d_spec, d_residual, d_mass, d_ions, d_tc, d_tn, d_idx, num_idx, len_spec); break;
    case 4: addIons_dispatch<4,true>(d_spec, d_residual, d_mass, d_ions, d_tc, d_tn, d_idx, num_idx, len_spec); break;
    default:
        assert(!"Non-exhaustive patterns in match");
    }
}

