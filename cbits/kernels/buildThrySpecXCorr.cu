/*
 * Module    : IonSeries
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 */

#include <assert.h>

#include <host_defines.h>
#include <device_functions.h>

#include "cbits/mass.h"
#include "cbits/kernels.h"


/*
 * Convert a given mass into a mass/charge ratio, and locate the appropriate
 * spectrum bin for the peak.
 */
__device__ int
binIonMZ(float mass, int charge)
{
    int bin = rintf((mass + massH * charge) / (charge * binWidthMono));
#ifdef __DEVICE_EMULATION__
    assert(bin >= 0 && bin < 2048 && "index out of bounds");
#endif
    return bin;
}


/*
 * Add a spectral peak for each fragment ion location, as well as the peaks
 * corresponding to the neutral losses of H2O and NH3.
 */
__device__ void
addIonsAB(float mass, int charge, int *spec)
{
    // A
    atomicMax(&spec[binIonMZ(mass - massCO, charge)], 10);

    // B
    int m = binIonMZ(mass, charge);

    atomicMax(&spec[m],   50);
    atomicMax(&spec[m+1], 25);
    atomicMax(&spec[m-1], 25);
    atomicMax(&spec[binIonMZ(mass - massH2O, charge)], 10);
    atomicMax(&spec[binIonMZ(mass - massNH3, charge)], 10);
}


__device__ void
addIonsY(float mass, int charge, int *spec)
{
    int m = binIonMZ(mass + massH2O, charge);

    atomicMax(&spec[m],   50);
    atomicMax(&spec[m+1], 25);
    atomicMax(&spec[m-1], 25);
    atomicMax(&spec[binIonMZ(mass - massNH3, charge)], 10);
}


/*
 * Add a spectral peak for each fragment ion location. The output spectrum array
 * must exist and be initialised to zero.
 */
__global__ void
buildThrySpecXCorr_kernel
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

    if (idx < len_ions)
    {
        int   charge = 1;
        float b_mass = b_ions[idx];
        float y_mass = y_ions[idx];

        do
        {
            addIonsAB(b_mass, charge, spec);
            addIonsY (y_mass, charge, spec);
        }
        while (++charge < max_charge);
    }
}


/*
 * Kernel wrapper to be called from C. The Ion and Spectrum arrays are in the
 * device memory space.
 */
void
buildThrySpecXCorr
(
    int          charge,
    float        *b_ions,
    float        *y_ions,
    int          *spec,
    unsigned int len_ions,
    unsigned int len_spec
)
{
    int threads = min(len_ions, 64);
    int blocks  = (len_ions + threads - 1) / threads;

    buildThrySpecXCorr_kernel<<<blocks,threads>>>(charge, b_ions, y_ions, spec, len_ions, len_spec);
}

