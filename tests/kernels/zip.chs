{-# LANGUAGE ForeignFunctionInterface #-}

module Main where

import Foreign.CUDA (DevicePtr, withDevicePtr)
import qualified Foreign.CUDA as CUDA

import Foreign
import Foreign.C


#include "../kernels/kernels.h"

{# fun unsafe zipWithMaxf
    { withDevicePtr* `DevicePtr CFloat' ,
      withDevicePtr* `DevicePtr CFloat' ,
      withDevicePtr* `DevicePtr CFloat' ,
      fromIntegral   `Int'              } -> `()' #}


main :: IO ()
main =
    let nums1  = map sin [1..1024]
        nums2  = map cos [1..1024]
        cpu_r = zipWith max nums1 nums2
    in
        CUDA.withArrayLen nums1 $ \l xs -> do
        CUDA.withArray    nums2 $ \ys   -> do
        CUDA.allocaBytes (fromIntegral l*4) $ \zs -> do
        zipWithMaxf xs ys zs l >> do
        CUDA.forceEither `fmap` CUDA.peekArray l zs >>= \gpu_r ->
            if gpu_r == cpu_r
            then putStrLn  "Test: PASSED"
            else putStrLn ("Test: FAILED (" ++ show (take 10 gpu_r) ++ " ... )")

