/*
 * Module    : Spectrum
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 *
 * Manipulate the results of a mass-spectroscopy experiment
 */

#include "utils.h"
#include "kernels.h"
#include "kernels-priv.h"


/*
 * Normalise each element of the input array according to the maximum value in
 * each of 10 equally sized windows.
 */
void
normaliseByRegion
(
    float       *spec,
    int         cutoff
)
{
    assert(!"Not implemented yet");
}


/*
 * Calculate the sequest cross-correlation function for the normalised
 * experimental spectrum
 */
void
calculateXCorr
(
    float       *spec,
    int         len
)
{
    assert(!"Not implemented yet");
}

