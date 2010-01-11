/* -----------------------------------------------------------------------------
 *
 * Module    : SMVM
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#ifndef __SMVM_PRIV_H__
#define __SMVM_PRIV_H__

/*
 * Optimised for Tesla C1060 (compute 1.3)
 * Maximum performance for your card may be achieved with different values.
 *
 * http://developer.download.nvidia.com/compute/cuda/CUDA_Occupancy_calculator.xls
 */
#define MAX_THREADS             128
#define MAX_BLOCKS_PER_SM       8
#define MAX_BLOCKS              (MAX_BLOCKS_PER_SM * 30)
#define WARP_SIZE               32

#endif
