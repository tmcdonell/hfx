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
import Foreign.CUDA (DevicePtr)

#include "cuda/kernels.h"

--
-- XXX: Need to replace this with a rewrite rules to achieve the same effect
--
withDevicePtr :: DevicePtr a -> (Ptr b -> IO c) -> IO c
withDevicePtr =  withForeignPtr . castForeignPtr

--------------------------------------------------------------------------------
-- Ion Series
--------------------------------------------------------------------------------

{# fun unsafe addIons
    { cIntConv          `Int'             ,
                        `Float'           ,
      withDevicePtr*    `DevicePtr Float' ,
      withDevicePtr*    `DevicePtr Int'   ,
                        `Int'             ,
                        `Int'             ,
                        `Int'             } -> `()' #}

--------------------------------------------------------------------------------
-- ZipWith
--------------------------------------------------------------------------------

{# fun unsafe zipWith_timesif
    { withDevicePtr* `DevicePtr Int'   ,
      withDevicePtr* `DevicePtr Float' ,
      withDevicePtr* `DevicePtr Float' ,
                     `Int'             } -> `()' #}


--------------------------------------------------------------------------------
-- Fold
--------------------------------------------------------------------------------

{# fun unsafe fold_plusf
    { withDevicePtr* `DevicePtr Float' ,
                     `Int'             } -> `Float' #}


--------------------------------------------------------------------------------
-- Scan
--------------------------------------------------------------------------------

{# fun unsafe scanl_plusui
    { withDevicePtr* `DevicePtr Int' ,
      withDevicePtr* `DevicePtr Int' ,
                     `Int'           } -> `()' #}

{# fun unsafe scanr_plusui
    { withDevicePtr* `DevicePtr Int' ,
      withDevicePtr* `DevicePtr Int' ,
                     `Int'           } -> `()' #}

{# fun unsafe scanl1Seg_plusf
    { withDevicePtr* `DevicePtr Float' ,
      withDevicePtr* `DevicePtr Int'   ,
      withDevicePtr* `DevicePtr Float' ,
                     `Int'             } -> `()' #}

{# fun unsafe scanr1Seg_plusf
    { withDevicePtr* `DevicePtr Float' ,
      withDevicePtr* `DevicePtr Int'   ,
      withDevicePtr* `DevicePtr Float' ,
                     `Int'             } -> `()' #}


--------------------------------------------------------------------------------
-- Permute
--------------------------------------------------------------------------------

{# fun unsafe permute_ui
    { withDevicePtr* `DevicePtr Int' ,
      withDevicePtr* `DevicePtr Int' ,
      withDevicePtr* `DevicePtr Int' ,
                     `Int'           } -> `()' #}

{# fun unsafe bpermute_f
    { withDevicePtr* `DevicePtr Float' ,
      withDevicePtr* `DevicePtr Float' ,
      withDevicePtr* `DevicePtr Int'   ,
                     `Int'             } -> `()' #}

{# fun unsafe compact_f
    { withDevicePtr* `DevicePtr Float' ,
      withDevicePtr* `DevicePtr Float' ,
      withDevicePtr* `DevicePtr Int'   ,
                     `Int'             } -> `Int' #}

--------------------------------------------------------------------------------
-- Map
--------------------------------------------------------------------------------

{# fun unsafe map_getAAMass
    { withDevicePtr* `DevicePtr Char'  ,
      withDevicePtr* `DevicePtr Float' ,
                     `Int'             } -> `()' #}

--------------------------------------------------------------------------------
-- Replicate
--------------------------------------------------------------------------------

{# fun unsafe replicate_ui
    { withDevicePtr* `DevicePtr Int' ,
                     `Int'           ,
                     `Int'           } -> `()' #}


--------------------------------------------------------------------------------
-- Sort
--------------------------------------------------------------------------------

{# fun unsafe sort_f
    { withDevicePtr* `DevicePtr Float' ,
                     `Int'             } -> `()' #}

{# fun unsafe sortPairs_f
    { withDevicePtr* `DevicePtr Float' ,
      withDevicePtr* `DevicePtr Float' ,
                     `Int'             } -> `()' #}


--------------------------------------------------------------------------------
-- groupBy
--------------------------------------------------------------------------------

{# fun unsafe group_f
    { withDevicePtr* `DevicePtr Float' ,
      withDevicePtr* `DevicePtr Int'   ,
                     `Int'             } -> `()' #}

