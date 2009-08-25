/*
 * Module    : IonSeries
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 */

#include "utils.h"
#include "kernels.h"
#include "kernels-priv.h"


/*
 * Generate the theoretical spectral representation of a peptide from its ion
 * fragment sequences. The ion and spectrum arrays exit in device memory, and
 * the latter is assumed initialised to zero.
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

