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


/*
 * Reduce a warp-sized chunk of data in shared memory.
 */
template <class T>
static __device__ T
reduce_warp(T sum, volatile T* s_data)
{
    s_data[threadIdx.x] = sum;                                  __EMUSYNC;
    s_data[threadIdx.x] = sum = sum + s_data[threadIdx.x + 16]; __EMUSYNC;
    s_data[threadIdx.x] = sum = sum + s_data[threadIdx.x +  8]; __EMUSYNC;
    s_data[threadIdx.x] = sum = sum + s_data[threadIdx.x +  4]; __EMUSYNC;
    s_data[threadIdx.x] = sum = sum + s_data[threadIdx.x +  2]; __EMUSYNC;
    s_data[threadIdx.x] = sum = sum + s_data[threadIdx.x +  1]; __EMUSYNC;

    return sum;
}


/*
 * Histogram binning
 */
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

/*
 * Dissociation products
 */
template <bool UseCache>
__inline__ __device__ static void
addIon(float &sum, const float *d_spec, const uint32_t N, const int32_t x, const float y)
{
    if (0 <= x && x < N)
        sum += y * fetch_x<UseCache>(x, d_spec);
}


template <uint32_t charge, bool UseCache>
__device__ void
addIonsAB(float &sum, const float *d_spec, const uint32_t N, const float mass)
{
    float   m;
    int32_t x;

    // A-ions
    addIon<UseCache>(sum, d_spec, N, binMZ(ionMZ(mass - MASS_CO, charge)), 10.0f);

    // B-ions
    m = ionMZ(mass, charge);
    x = binMZ(m);

    addIon<UseCache>(sum, d_spec, N, x,   50.0f);
    addIon<UseCache>(sum, d_spec, N, x+1, 25.0f); // technically, should be binMZ(m+1)
    addIon<UseCache>(sum, d_spec, N, x-1, 25.0f);

    addIon<UseCache>(sum, d_spec, N, binMZ(m - __fdividef(MASS_H2O, charge)), 10.0f);
    addIon<UseCache>(sum, d_spec, N, binMZ(m - __fdividef(MASS_NH3, charge)), 10.0f);
}


template <uint32_t charge, bool UseCache>
__device__ void
addIonsY(float &sum, const float *d_spec, const uint32_t N, const float mass)
{
    float   m = ionMZ(mass + MASS_H2O, charge);
    int32_t x = binMZ(m);

    // Y-ions
    addIon<UseCache>(sum, d_spec, N, x,   50.0f);
    addIon<UseCache>(sum, d_spec, N, x+1, 25.0f);
    addIon<UseCache>(sum, d_spec, N, x-1, 25.0f);

    addIon<UseCache>(sum, d_spec, N, binMZ(m - __fdividef(MASS_NH3, charge)), 10.0f);
}


template <uint32_t charge, bool UseCache>
__device__ void
addIons_k(float &sum, const float *d_spec, const uint32_t N, const float b_mass, const float y_mass)
{
    addIonsAB<charge,UseCache>(sum, d_spec, N, b_mass);
    addIonsY <charge,UseCache>(sum, d_spec, N, y_mass);
}


/* -----------------------------------------------------------------------------
 * Sequest theoretical spectrum and cross-correlation analysis
 * -----------------------------------------------------------------------------
 *
 * Generate theoretical spectra for a collection of peptide fragments. The
 * 'ions' array contains the individual amino-acid masses for the database
 * entries. We are interested in the sequences generated between the terminal
 * indices (tc,tn) of the locations specified in the 'idx' array.
 *
 * A warp of threads iterates between the (tc,tn) indices, generating the b- and
 * y-ion mass ladders. These fragment locations, together with those
 * corresponding to neutral losses of H2O and NH3, are the spectral peaks that
 * will be combined with the pre-processed experimental spectrum to calculate
 * the sequest correlation score.
 *
 * Optionally, the texture cache may be used for accessing the d_spec
 * (experimental spectrum) vector. This generally shows good improvements.
 *
 */
template <uint32_t BlockSize, uint32_t MaxCharge, bool UseCache>
__global__ static void
addIons_core
(
    float               *d_score,       // output array of scores, length num_idx
    const float         *d_spec,        // experimental spectrum
    const float         *d_residual,    // peptide residual mass
    const float         *d_ions,        // individual ion masses
    const uint32_t      *d_tc,          // c-terminal indices
    const uint32_t      *d_tn,          // n-terminal indices
    const uint32_t      *d_idx,         // The indices of the sequences under consideration
    const uint32_t      num_idx,
    const uint32_t      len_spec
)
{
    /*
     * Require at least a full warp for each row. This could be relaxed by
     * modifying the cooperative reduction step
     */
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
        float          sum;
        float          b_mass;
        float          y_mass;

        s_data[threadIdx.x]      = 0.0f;

        /*
         * Have all threads read in mass values for this segment, calculating
         * dissociation products and partial dot-product sums to the
         * pre-processed experimental spectrum.
         */
        for (uint32_t j = row_start + thread_lane; j < row_end; j += WARP_SIZE)
        {
            /*
             * Load the ion mass, and propagate the partial scan results
             */
            b_mass = d_ions[j];

            if (thread_lane == 0)
                b_mass += s_data[threadIdx.x + (WARP_SIZE-1)];

            /*
             * Generate fragment mass ladder
             */
            b_mass = scan_warp<float,true>(b_mass, s_data);
            y_mass = residual - b_mass;
            sum    = 0.0f;

            if (1 <= MaxCharge) addIons_k<1,UseCache>(sum, d_spec, len_spec, b_mass, y_mass);
            if (2 <= MaxCharge) addIons_k<2,UseCache>(sum, d_spec, len_spec, b_mass, y_mass);
            if (3 <= MaxCharge) addIons_k<3,UseCache>(sum, d_spec, len_spec, b_mass, y_mass);
            if (4 <= MaxCharge) addIons_k<4,UseCache>(sum, d_spec, len_spec, b_mass, y_mass);
        }

        /*
         * Reduce the partial dot-product results from all threads, and write
         * the result for this sequence
         */
        reduce_warp<float>(sum, s_data);

        if (thread_lane == 0)
            d_score[row] = s_data[threadIdx.x];
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
    float               *d_score,
    const float         *d_spec,
    const float         *d_residual,
    const float         *d_ions,
    const uint32_t      *d_tc,
    const uint32_t      *d_tn,
    const uint32_t      *d_idx,
    const uint32_t      num_idx,
    const uint32_t      len_spec
)
{
    uint32_t blocks;
    uint32_t threads;

    addIons_control(num_idx, blocks, threads);
    switch (threads)
    {
//  case 512: addIons_core<512,MaxCharge,UseCache><<<blocks,threads>>>(d_score, d_spec, d_residual, d_ions, d_tc, d_tn, d_idx, num_idx, len_spec); break;
//  case 256: addIons_core<256,MaxCharge,UseCache><<<blocks,threads>>>(d_score, d_spec, d_residual, d_ions, d_tc, d_tn, d_idx, num_idx, len_spec); break;
    case 128: addIons_core<128,MaxCharge,UseCache><<<blocks,threads>>>(d_score, d_spec, d_residual, d_ions, d_tc, d_tn, d_idx, num_idx, len_spec); break;
    case  64: addIons_core< 64,MaxCharge,UseCache><<<blocks,threads>>>(d_score, d_spec, d_residual, d_ions, d_tc, d_tn, d_idx, num_idx, len_spec); break;
    case  32: addIons_core< 32,MaxCharge,UseCache><<<blocks,threads>>>(d_score, d_spec, d_residual, d_ions, d_tc, d_tn, d_idx, num_idx, len_spec); break;
    default:
        assert(!"Non-exhaustive patterns in match");
    }
}


void
addIons_inplace
(
    float               *d_score,
    const float         *d_spec,
    const float         *d_residual,
    const float         *d_ions,
    const uint32_t      *d_tc,
    const uint32_t      *d_tn,
    const uint32_t      *d_idx,
    const uint32_t      num_idx,
    const uint32_t      max_charge,
    const uint32_t      len_spec
)
{
    bind_x(d_spec);

    switch (max_charge)
    {
    case 1: addIons_dispatch<1,true>(d_score, d_spec, d_residual, d_ions, d_tc, d_tn, d_idx, num_idx, len_spec); break;
    case 2: addIons_dispatch<2,true>(d_score, d_spec, d_residual, d_ions, d_tc, d_tn, d_idx, num_idx, len_spec); break;
    case 3: addIons_dispatch<3,true>(d_score, d_spec, d_residual, d_ions, d_tc, d_tn, d_idx, num_idx, len_spec); break;
    case 4: addIons_dispatch<4,true>(d_score, d_spec, d_residual, d_ions, d_tc, d_tn, d_idx, num_idx, len_spec); break;
    default:
        assert(!"Non-exhaustive patterns in match");
    }

    unbind_x(d_spec);
}

