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
-- ZipWith
--------------------------------------------------------------------------------

{# fun unsafe zipWith_timesif
    { withDevicePtr* `DevicePtr CInt'   ,
      withDevicePtr* `DevicePtr CFloat' ,
      withDevicePtr* `DevicePtr CFloat' ,
                     `Int'              } -> `()' #}


--------------------------------------------------------------------------------
-- Fold
--------------------------------------------------------------------------------

{# fun unsafe fold_plusf
    { withDevicePtr* `DevicePtr CFloat' ,
                     `Int'              } -> `Float' #}


--------------------------------------------------------------------------------
-- Scan
--------------------------------------------------------------------------------

{# fun unsafe scanl_plusui
    { withDevicePtr* `DevicePtr CUInt' ,
      withDevicePtr* `DevicePtr CUInt' ,
                     `Int'             } -> `()' #}

{# fun unsafe scanr_plusui
    { withDevicePtr* `DevicePtr CUInt' ,
      withDevicePtr* `DevicePtr CUInt' ,
                     `Int'             } -> `()' #}

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


--------------------------------------------------------------------------------
-- Permute
--------------------------------------------------------------------------------

{# fun unsafe permute_ui
    { withDevicePtr* `DevicePtr CUInt' ,
      withDevicePtr* `DevicePtr CUInt' ,
      withDevicePtr* `DevicePtr CUInt' ,
                     `Int'             } -> `()' #}

{# fun unsafe bpermute_f
    { withDevicePtr* `DevicePtr CFloat' ,
      withDevicePtr* `DevicePtr CFloat' ,
      withDevicePtr* `DevicePtr CUInt'  ,
                     `Int'              } -> `()' #}

{# fun unsafe compact_f
    { withDevicePtr* `DevicePtr CFloat' ,
      withDevicePtr* `DevicePtr CFloat' ,
      withDevicePtr* `DevicePtr CUInt'  ,
                     `Int'              } -> `Int' #}

--------------------------------------------------------------------------------
-- Map
--------------------------------------------------------------------------------

{# fun unsafe map_getAAMass
    { withDevicePtr* `DevicePtr CChar'  ,
      withDevicePtr* `DevicePtr CFloat' ,
                     `Int'              } -> `()' #}

--------------------------------------------------------------------------------
-- Replicate
--------------------------------------------------------------------------------

{# fun unsafe replicate_ui
    { withDevicePtr* `DevicePtr CUInt' ,
                     `Int'             ,
                     `Int'             } -> `()' #}


--------------------------------------------------------------------------------
-- Sort
--------------------------------------------------------------------------------

{# fun unsafe sort_f
    { withDevicePtr* `DevicePtr CFloat' ,
                     `Int'              } -> `()' #}

{# fun unsafe sortPairs_f
    { withDevicePtr* `DevicePtr CFloat' ,
      withDevicePtr* `DevicePtr CFloat' ,
                     `Int'              } -> `()' #}

