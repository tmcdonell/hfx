/* -----------------------------------------------------------------------------
 *
 * Module    : Device
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#ifndef __DEVICE_H__
#define __DEVICE_H__

/*
 * Sensible default limits for thread and block sizes. Based on testing with a
 * Tesla C1060 (compute 1.3), but may differ for other processors / kernels.
 *
 * http://developer.download.nvidia.com/compute/cuda/CUDA_Occupancy_calculator.xls
 */
#define MAX_BLOCKS_PER_MP       8
#define NUM_MULTI_PROCESSORS    30

#define WARP_SIZE               32
#define MAX_THREADS             128
#define MAX_BLOCKS              (MAX_BLOCKS_PER_MP * NUM_MULTI_PROCESSORS)

#endif

