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
 * test.cu
 */
int  reducePlusi(int *xs, int N);

void scanl1Plusi(int *in, int *out, int N);
void scanr1Plusi(int *in, int *out, int N);

void zipWithPlusi(int *xs, int *ys, int *zs, int N);
void zipWithMaxf(float *xs, float *ys, float *zs, int N);


#ifdef __cplusplus
}
#endif
#endif
/*
 * vim: syn=cuda
 */
