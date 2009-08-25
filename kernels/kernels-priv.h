/*
 * Module    : Kernels
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 *
 * Functions which are internally implemented using CUDA, internal interface
 */

#ifndef __KERNELS_PRIV_H__
#define __KERNELS_PRIV_H__

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Ion Series
 */
__global__ void
buildThrySpecXCorr_kernel(int, float*, float*, int*, unsigned int, unsigned int);


#ifdef __cplusplus
}
#endif
#endif
/*
 * vim: syn=cuda
 */
