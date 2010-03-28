/* -----------------------------------------------------------------------------
 *
 * Module    : MVM
 * Copyright : (c) [2009..2010] Trevor L. McDonell
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#ifndef __MVM_PRIV_H__
#define __MVM_PRIV_H__

/*
 * Prefer the y-dimension to be larger, and require the x-dimension <= 64
 */
#define BLOCKDIM_X      16
#define BLOCKDIM_Y      16

#define MAX_THREADS     undefined
#define MAX_BLOCKS      undefined
#define WARP_SIZE       32

#endif


