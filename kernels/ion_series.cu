/*
 * Module    : IonSeries
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 */

#include "mass.h"
#include "utils.h"
#include "kernels.h"


/*
 * Convert a given mass into a mass/charge ratio
 * Locate the appropriate spectrum bin for a peak.
 */
__device__ float ionMZ(float m, float c) { return (m + massH * c) / c; }
__device__ int   bin(float x) { return rintf(x / binWidthMono); }


/*
 * Add a spectral peak for each fragment ion location, as well as the peaks
 * corresponding to the neutral losses of H2O and NH3.
 */
__device__ void
addIonsAB(float mass, float charge, int *spec, unsigned int N)
{
    int   idx;
    float m;

    // A
    idx = bin(ionMZ(mass - massCO, charge));
    if (0 <= idx && idx < N) atomicMax(&spec[idx], 10);

    // B
    m   = ionMZ(mass, charge);
    idx = bin(m);

    if (1 <= idx && idx < N-1)
    {
        atomicMax(&spec[idx],   50);
        atomicMax(&spec[idx+1], 25);
        atomicMax(&spec[idx-1], 25);
    }

    idx = bin(m - massH2O/charge);
    if (0 <= idx && idx < N) atomicMax(&spec[idx], 10);

    idx = bin(m - massNH3/charge);
    if (0 <= idx && idx < N) atomicMax(&spec[idx], 10);
}


__device__ void
addIonsY(float mass, float charge, int *spec, unsigned int N)
{
    float m   = ionMZ(mass + massH2O, charge);
    int   idx = bin(m);

    if (1 <= idx && idx < N-1)
    {
        atomicMax(&spec[idx],   50);
        atomicMax(&spec[idx+1], 25);
        atomicMax(&spec[idx-1], 25);
    }

    idx = bin(m - massNH3/charge);
    if (0 <= idx && idx < N) atomicMax(&spec[idx], 10);
}


/*
 * Add a spectral peak for each fragment ion location. The output spectrum array
 * must exist and be initialised to zero.
 */
template <bool lengthIsPow2>
__global__ static void
addIons_core
(
    int          max_charge,
    float        *b_ions,
    float        *y_ions,
    int          *spec,
    unsigned int len_ions,
    unsigned int len_spec
)
{
    unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x;

    if (lengthIsPow2 || idx < len_ions)
    {
        int   charge = 1;
        float b_mass = b_ions[idx];
        float y_mass = y_ions[idx];

        do
        {
            addIonsAB(b_mass, (float) charge, spec, len_spec);
            addIonsY (y_mass, (float) charge, spec, len_spec);
        }
        while (++charge < max_charge);
    }
}


void
addIons
(
    int          max_charge,
    float        *b_ions,
    float        *y_ions,
    int          *spec,
    unsigned int len_ions,
    unsigned int len_spec
)
{
    unsigned int threads = min(ceilPow2(len_ions), 512);
    unsigned int blocks  = (len_ions + threads - 1) / threads;

    if (isPow2(len_ions))
        addIons_core<true><<<blocks,threads>>>(max_charge, b_ions, y_ions, spec, len_ions, len_spec);
    else
        addIons_core<false><<<blocks,threads>>>(max_charge, b_ions, y_ions, spec, len_ions, len_spec);
}

