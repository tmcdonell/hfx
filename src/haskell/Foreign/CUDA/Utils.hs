--------------------------------------------------------------------------------
-- |
-- Module    : Foreign.CUDA.Utils
-- Copyright : (c) 2010 Trevor L. McDonell
-- License   : BSD
--
--------------------------------------------------------------------------------

module Foreign.CUDA.Utils where

import Foreign
import Control.Exception

import Data.Vector.Generic       as G
import Data.Vector.Fusion.Stream as Stream
import Foreign.CUDA              as CUDA


--
-- Copy a vector to the device and perform a computation, returning the result.
-- This requires the data to be marshalled to a heap allocated array so that it
-- can be copied to the device.
--
withVector :: (G.Vector v a, Storable a) => v a -> (DevicePtr a -> IO b) -> IO b
withVector vec action = let l = G.length vec in
  CUDA.allocaArray l  $ \d_ptr -> do

    -- Transfer via a pinned heap array, which the device can read directly at
    -- much higher bandwidth than pageable memory.
    --
    bracket (CUDA.mallocHostArray [] l) CUDA.freeHost $
      flip withHostPtr $ \ptr ->
        Stream.foldM' (\p e -> poke p e >> return (p `advancePtr` 1)) ptr (G.stream vec) >>
        CUDA.pokeArray l ptr d_ptr

    -- Release host memory and perform computation on device array
    --
    action d_ptr

