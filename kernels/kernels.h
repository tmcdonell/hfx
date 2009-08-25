/*
 * Module    : Kernels
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 *
 * Functions which are internally implemented using CUDA, public interface.
 */

#ifndef __KERNELS_H__
#define __KERNELS_H__

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Ion Series
 */
void
buildThrySpecXCorr(int, float*, float*, int*, unsigned int, unsigned int);

/*
 * Spectrum
 */
void
normaliseByRegion(float *, int);

void
calculateXCorr(float *, int);


#ifdef __cplusplus
}
#endif
#endif
/*
 * vim: syn=cuda
 */
