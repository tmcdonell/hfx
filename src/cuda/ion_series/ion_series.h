/* -----------------------------------------------------------------------------
 *
 * Module    : Ion Series
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/


#ifndef __ION_SERIES_H__
#define __ION_SERIES_H__

/*
 * Based on Tesla C1060 (compute 1.3)
 * http://developer.download.nvidia.com/compute/cuda/CUDA_Occupancy_calculator.xls
 */
#define WARP_SIZE               32
#define MAX_THREADS             128
#define MAX_BLOCKS_PER_SM       8
#define MAX_BLOCKS              (MAX_BLOCKS_PER_SM * 30)

/*
 * Spectrum bin width
 */
#define BIN_WIDTH_MONO  1.0005079f
#define BIN_WIDTH_AVG   1.0011413f

/*
 * The monoisotopic mass of several elements and molecules
 */
#define MASS_H2O        18.01056f
#define MASS_NH3        17.02655f
#define MASS_CO         27.9949f
#define MASS_O          16.0013f
#define MASS_H          1.0078246f

#endif
