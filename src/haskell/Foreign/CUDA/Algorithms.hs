{-# LANGUAGE ForeignFunctionInterface #-}
--------------------------------------------------------------------------------
-- |
-- Module    : Foreign.CUDA.Algorithms
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- FFI Bindings to a collection of algorithms, implemented to run on the
-- graphics processor using CUDA.
--
--------------------------------------------------------------------------------

module Foreign.CUDA.Algorithms
  (
    findIndicesInRange,
    addIons, addIonsIP,
    rsort,
    mvm
  )
  where

import Util.C2HS
import Data.Word
import Foreign
import Foreign.CUDA                             (DevicePtr, withDevicePtr)


findIndicesInRange :: DevicePtr Float -> DevicePtr Word32 -> Int -> Float -> Float -> IO Int
findIndicesInRange a1 a2 a3 a4 a5 =
  withDevicePtr a1 $ \a1' ->
  withDevicePtr a2 $ \a2' ->
  cIntConv `fmap` findIndicesInRange'_ a1' a2' (cIntConv a3) a4 a5

foreign import ccall unsafe "algorithms.h findIndicesInRange_f"
  findIndicesInRange'_ :: Ptr Float -> Ptr Word32 -> Word32 -> Float -> Float -> IO Word32


addIons :: DevicePtr Word32 -> DevicePtr Float -> DevicePtr Float -> (DevicePtr Word32, DevicePtr Word32) -> DevicePtr Word32 -> Int -> Int -> Int -> IO ()
addIons a1 a2 a3 (a4,a5) a6 a7 a8 a9 =
  withDevicePtr a1 $ \a1' ->
  withDevicePtr a2 $ \a2' ->
  withDevicePtr a3 $ \a3' ->
  withDevicePtr a4 $ \a4' ->
  withDevicePtr a5 $ \a5' ->
  withDevicePtr a6 $ \a6' ->
  addIons'_ a1' a2' a3' a4' a5' a6' (cIntConv a7) (cIntConv a8) (cIntConv a9)

foreign import ccall unsafe "algorithms.h addIons"
  addIons'_ :: Ptr Word32 -> Ptr Float -> Ptr Float -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Word32 -> Word32 -> Word32 -> IO ()


addIonsIP
    :: DevicePtr Float                          -- [out] sequence scores
    -> DevicePtr Float                          -- experimental spectrum
    -> DevicePtr Float                          -- residual masses
    -> DevicePtr Float                          -- individual ion masses
    -> (DevicePtr Word32, DevicePtr Word32)     -- c- and n- terminal indices
    -> DevicePtr Word32                         -- indices of the sequences under consideration
    -> Int                                      -- number of sequences to consider
    -> Int                                      -- peptide charge state
    -> Int                                      -- length of input spectrum
    -> IO ()
addIonsIP d_score d_spec d_residual d_ions (d_tc, d_tn) d_idx num_idx max_charge len_spec =
  withDevicePtr d_score    $ \a1' ->
  withDevicePtr d_spec     $ \a2' ->
  withDevicePtr d_residual $ \a3' ->
  withDevicePtr d_ions     $ \a4' ->
  withDevicePtr d_tc       $ \a5' ->
  withDevicePtr d_tn       $ \a6' ->
  withDevicePtr d_idx      $ \a7' ->
  addIons_ip'_ a1' a2' a3' a4' a5' a6' a7' (cIntConv num_idx) (cIntConv max_charge) (cIntConv len_spec)

foreign import ccall unsafe "algorithms.h addIons_inplace"
  addIons_ip'_ :: Ptr Float -> Ptr Float -> Ptr Float -> Ptr Float -> Ptr Word32 -> Ptr Word32 -> Ptr Word32 -> Word32 -> Word32 -> Word32 -> IO ()



rsort :: DevicePtr Float -> DevicePtr Word32 -> Int -> IO ()
rsort a1 a2 a3 =
  withDevicePtr a1 $ \a1' ->
  withDevicePtr a2 $ \a2' ->
  rsort'_ a1' a2' (cIntConv a3)

foreign import ccall unsafe "algorithms.h sort_rf"
  rsort'_ :: Ptr Float -> Ptr Word32 -> Word32 -> IO ()


mvm :: DevicePtr Float -> DevicePtr Word32 -> DevicePtr Float -> Int -> Int -> IO ()
mvm a1 a2 a3 a4 a5 =
  withDevicePtr a1 $ \a1' ->
  withDevicePtr a2 $ \a2' ->
  withDevicePtr a3 $ \a3' ->
  mvm'_ a1' a2' a3' (fromIntegral a4) (fromIntegral a5)

foreign import ccall unsafe "algorithms.h mvm_if"
  mvm'_ :: Ptr Float -> Ptr Word32 -> Ptr Float -> Word32 -> Word32 -> IO ()


