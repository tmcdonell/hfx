/* -----------------------------------------------------------------------------
 *
 * Module    : Utils
 * Copyright : (c) [2009..2010] Trevor L. McDonell
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#ifndef __CUDPP_UTILS__
#define __CUDPP_UTILS__

#include <cudpp.h>

template <typename T> CUDPPDatatype inline getType();
template <> CUDPPDatatype inline getType<float>()        { return CUDPP_FLOAT; }
template <> CUDPPDatatype inline getType<unsigned int>() { return CUDPP_UINT;  }

#endif

