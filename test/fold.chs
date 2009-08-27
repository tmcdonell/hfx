{-# LANGUAGE ForeignFunctionInterface #-}

module Main where

import Foreign.CUDA (DevicePtr, withDevicePtr)
import qualified Foreign.CUDA as CUDA

import Foreign
import Foreign.C

#if 0
int
foldPlus(const int *xs, int n)
{
    return fold< Plus<int> >(xs, n);
}
#endif

#include "../kernels/kernels.h"

{# fun unsafe foldPlusi
    { withDevicePtr* `DevicePtr CInt' ,
      fromIntegral   `Int'            } -> `Int' fromIntegral #}


--
-- Result should be 524800
--
main :: IO ()
main =
    let nums  = [1..2048]
        cpu_r = fromIntegral (foldl (+) 0 nums)
    in do
        CUDA.withArrayLen nums $ \l xs -> do
        foldPlusi xs l >>= \gpu_r      ->
            if gpu_r == cpu_r
            then putStrLn "Test: PASSED"
            else putStrLn ("Test: FAILED (" ++ (show gpu_r) ++ ")")

