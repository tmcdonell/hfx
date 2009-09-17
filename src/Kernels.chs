{-# LANGUAGE CPP, ForeignFunctionInterface #-}
--------------------------------------------------------------------------------
-- |
-- Module    : Kernels
-- Copyright : (c) 2009 Trevor L. McDonell
-- License   : BSD
--
-- Foreign function bindings to CUDA kernels
--
--------------------------------------------------------------------------------


module Kernels where

import C2HS
import Foreign.CUDA (DevicePtr, withDevicePtr)

#include "cuda/kernels.h"


--------------------------------------------------------------------------------
-- Ion Series
--------------------------------------------------------------------------------

{# fun unsafe addIons
    { cIntConv          `Int'              ,
      withDevicePtr*    `DevicePtr CFloat' ,
      withDevicePtr*    `DevicePtr CFloat' ,
      withDevicePtr*    `DevicePtr CInt'   ,
                        `Int'              ,
                        `Int'              } -> `()' #}

--------------------------------------------------------------------------------
-- Prelude
--------------------------------------------------------------------------------

{# fun unsafe zipWith_timesif
    { withDevicePtr* `DevicePtr CInt'   ,
      withDevicePtr* `DevicePtr CFloat' ,
      withDevicePtr* `DevicePtr CFloat' ,
                     `Int'              } -> `()' #}

{# fun unsafe fold_plusf
    { withDevicePtr* `DevicePtr CFloat' ,
                     `Int'              } -> `Float' #}

{# fun unsafe scanl1Seg_plusf
    { withDevicePtr* `DevicePtr CFloat' ,
      withDevicePtr* `DevicePtr CUInt'  ,
      withDevicePtr* `DevicePtr CFloat' ,
                     `Int'              } -> `()' #}

{# fun unsafe scanr1Seg_plusf
    { withDevicePtr* `DevicePtr CFloat' ,
      withDevicePtr* `DevicePtr CUInt'  ,
      withDevicePtr* `DevicePtr CFloat' ,
                     `Int'              } -> `()' #}

{# fun unsafe permute_f
    { withDevicePtr* `DevicePtr CFloat' ,
      withDevicePtr* `DevicePtr CFloat' ,
      withDevicePtr* `DevicePtr CInt'   ,
                     `Int'              } -> `()' #}

