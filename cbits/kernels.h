/*
 * Module    : Kernels
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 *
 * Functions which are internally implemented using CUDA
 *
 */

#ifndef __KERNELS_H__
#define __KERNELS_H__

#ifdef __cplusplus
extern "C" {
#endif

void
buildThrySpecXCorr(int charge, float *b_ions, float *y_ions, int *spec, unsigned int len_ions, unsigned int len_spec);


#ifdef __cplusplus
}
#endif
#endif
