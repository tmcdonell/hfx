{-# LANGUAGE ForeignFunctionInterface #-}

module Main where

import Foreign.CUDA (DevicePtr, withDevicePtr)
import qualified Foreign.CUDA as CUDA

import Foreign
import Foreign.C


#include "../kernels/kernels.h"

{# fun unsafe scanr1Plusi
    { withDevicePtr* `DevicePtr CInt' ,
      withDevicePtr* `DevicePtr CInt' ,
      fromIntegral   `Int'              } -> `()' #}


main :: IO ()
main =
    let nums  = [1..473]
        cpu_r = scanr1 (+) nums
    in
        CUDA.withArrayLen nums $ \l xs -> do
        CUDA.allocaBytes (fromIntegral l*4) $ \zs -> do
        scanr1Plusi xs zs l >> do
        CUDA.forceEither `fmap` CUDA.peekArray l zs >>= \gpu_r ->
            if gpu_r == cpu_r
            then putStrLn  "Test: PASSED"
            else putStrLn ("Test: FAILED (" ++ show (take 10 gpu_r) ++ " ... )")

