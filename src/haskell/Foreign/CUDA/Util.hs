--------------------------------------------------------------------------------
-- |
-- Module    : Foreign.CUDA.Util
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
--------------------------------------------------------------------------------

module Foreign.CUDA.Util where

import Foreign
import Control.Exception

import Data.Vector.Storable      as S
import Data.Vector.Generic       as G
import Data.Vector.Fusion.Stream as Stream
import Foreign.CUDA              as CUDA


--
-- Return the size of each element in a device array
--
sizeOfPtr :: Storable a => DevicePtr a -> Int
sizeOfPtr =  sizeOf . (undefined :: DevicePtr a -> a)


--
-- Copy a vector to the device and perform a computation, returning the result.
-- This requires the data to be marshalled to a heap allocated array so that it
-- can be copied to the device.
--
{-# INLINE withVector #-}
withVector :: (G.Vector v a, Storable a) => v a -> (DevicePtr a -> IO b) -> IO b
withVector vec action = let l = G.length vec in
  CUDA.allocaArray l  $ \d_ptr -> do

    -- Transfer via a pinned heap array, which the device can read directly at
    -- much higher bandwidth than pageable memory.
    --
    bracket (CUDA.mallocHostArray [] l) CUDA.freeHost $ \h_ptr ->
      withHostPtr h_ptr $ \ptr ->
        Stream.foldM' (\p e -> poke p e >> return (p `advancePtr` 1)) ptr (G.stream vec) >>
        CUDA.pokeArrayAsync l h_ptr d_ptr Nothing

    -- Release host memory and perform computation on device array
    --
    action d_ptr


--
-- Copy the contents of a storable-based vector to the device and perform a
-- computation, returning the result. As storable vectors are stored as foreign
-- pointers which are visible to the C heap, it is not necessary to marshal the
-- data via a temporary array.
--
{-# INLINE withVectorS #-}
withVectorS :: Storable a => S.Vector a -> (DevicePtr a -> IO b) -> IO b
withVectorS vec action = let l = S.length vec in
  CUDA.allocaArray l $ \d_ptr -> do
    S.unsafeWith vec $ \p -> CUDA.pokeArray l p d_ptr
    action d_ptr

