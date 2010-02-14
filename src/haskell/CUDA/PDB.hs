{-# LANGUAGE ForeignFunctionInterface #-}
--------------------------------------------------------------------------------
-- |
-- Module    : CUDA.PDB
-- Copyright : (c) 2009 Trevor L. McDonell
-- License   : BSD
--
-- Sequest matching by constructing and sampling theoretical spectra in-situ
--
--------------------------------------------------------------------------------

module CUDA.PDB where

import Time
import Mass
import C2HS
--import Utils
import Config
import Spectrum
import Sequest                                  (MatchCollection, Match(..))
import CUDA.Database                            (ProteinDatabase(..))

import Data.Word
import Data.List                                (intersperse)
import Control.Monad                            (when)
import System.IO
import Foreign.Ptr
import Foreign.CUDA
import Foreign.Storable
import Foreign.ForeignPtr
import qualified Data.Vector                    as V
import qualified Data.Vector.Storable           as S



--------------------------------------------------------------------------------
-- Main
--------------------------------------------------------------------------------

searchForMatches :: ConfigParams Float -> ProteinDatabase Float -> Spectrum Float -> IO (MatchCollection Float)
searchForMatches cp db dta =
  withVector spec         $ \d_specExp ->
  allocaArray numpeptides $ \d_pepIdx  -> do

    -- Find peptides in the database with a residual mass close to the spectral
    -- precursor mass
    --
    (tf,nr) <- bracketTime $ findIndicesInRange (residuals db) d_pepIdx numpeptides (mass-delta) (mass+delta)
    whenVerbose cp ["> findIndices: " ++ showTime tf, shows nr " found"]

    -- Generate theoretical spectra for each of these candidates
    --
    allocaArray nr             $ \d_score    -> do
    allocaArray (nr * specLen) $ \d_specThry -> do
      memset d_specThry (fromIntegral $ specLen * nr * sizeOfPtr d_specThry) 0
      (ta,_) <- bracketTime $ addIons d_specThry (yIonLadders db) (rowOffsets db) d_pepIdx nr (round chrg) specLen
      whenVerbose cp ["> addIons: " ++ showTime ta]

      -- Score each peptide
      --
      (tm,_) <- bracketTime $ mvm d_score d_specThry d_specExp nr specLen
      let elm = 2 * (sizeOfPtr d_specThry) * (nr * specLen + specLen)
      whenVerbose cp ["> mvm: " ++ showTime tm, showRateSI elm tm "FLOPS"]

      -- Sort results and retrieve the best matches
      --
      (ts,_) <- bracketTime $ radixsort d_score d_pepIdx nr
      whenVerbose cp ["> sort: " ++ showTime ts, showRateSI nr ts "pairs"]

      sc <- peekListArray results (d_score  `advanceDevPtr` (nr-results))
      ix <- peekListArray results (d_pepIdx `advanceDevPtr` (nr-results))

      finish sc ix
  where
    numpeptides  = V.length (peptides db)
    results      = max (numMatches cp) (numMatchesDetail cp)

    spec         = buildExpSpecXCorr cp dta
    specLen      = S.length spec
    chrg         = max 1 (charge dta - 1)

    finish sc ix = return . reverse $ zipWith (\s i -> Match (peptides db V.! cIntConv i) (s/10000)) sc ix

    -- Residual mass plus the mass of the water molecule released when forming
    -- the peptide bond
    --
    mass         = (precursor dta * charge dta) - ((charge dta * massH) - 1) - (massH + massH2O)
    delta        = massTolerance cp


--------------------------------------------------------------------------------
-- Utils
--------------------------------------------------------------------------------

whenVerbose      :: ConfigParams a -> [String] -> IO ()
whenVerbose cp s =  when (verbose cp) (hPutStrLn stderr . concat . intersperse ", " $ s)

sizeOfPtr :: Storable a => DevicePtr a -> Int
sizeOfPtr =  sizeOf . (undefined :: DevicePtr a -> a)

withVector :: Storable a => S.Vector a -> (DevicePtr a -> IO b) -> IO b
withVector (S.Vector _ l p) f =
  allocaArray l      $ \dp -> do
    withForeignPtr p $ \p' -> pokeArray l p' dp
    f dp


--------------------------------------------------------------------------------
-- FFI Bindings (put these somewhere else...)
--------------------------------------------------------------------------------

findIndicesInRange :: DevicePtr Float -> DevicePtr Word32 -> Int -> Float -> Float -> IO Int
findIndicesInRange a1 a2 a3 a4 a5 =
  withDevicePtr a1 $ \a1' ->
  withDevicePtr a2 $ \a2' ->
  findIndicesInRange'_ a1' a2' (cIntConv a3) a4 a5 >>= \res ->
  return (cIntConv res)

foreign import ccall unsafe "algorithms.h findIndicesInRange_f"
  findIndicesInRange'_ :: Ptr Float -> Ptr Word32 -> Word32 -> Float -> Float -> IO Word32


addIons :: DevicePtr Word32 -> DevicePtr Float -> DevicePtr Word32 -> DevicePtr Word32 -> Int -> Int -> Int -> IO ()
addIons a1 a2 a3 a4 a5 a6 a7 =
  withDevicePtr a1 $ \a1' ->
  withDevicePtr a2 $ \a2' ->
  withDevicePtr a3 $ \a3' ->
  withDevicePtr a4 $ \a4' ->
  addIons'_ a1' a2' a3' a4' (cIntConv a5) (cIntConv a6) (cIntConv a7)

foreign import ccall unsafe "algorithms.h addIons"
  addIons'_ :: Ptr Word32 -> Ptr Float -> Ptr Word32 -> Ptr Word32 -> Word32 -> Word32 -> Word32 -> IO ()


radixsort :: Storable a => DevicePtr Float -> DevicePtr a -> Int -> IO ()
radixsort a1 a2 a3 =
  withDevicePtr a1 $ \a1' ->
  withDevicePtr a2 $ \a2' ->
  radixsort'_ a1' (castPtr a2') (cIntConv a3)

foreign import ccall unsafe "algorithms.h radixsort_f"
  radixsort'_ :: Ptr Float -> Ptr () -> Word32 -> IO ()


mvm :: DevicePtr Float -> DevicePtr Word32 -> DevicePtr Float -> Int -> Int -> IO ()
mvm a1 a2 a3 a4 a5 =
  withDevicePtr a1 $ \a1' ->
  withDevicePtr a2 $ \a2' ->
  withDevicePtr a3 $ \a3' ->
  mvm'_ a1' a2' a3' (fromIntegral a4) (fromIntegral a5)

foreign import ccall unsafe "algorithms.h mvm_if"
  mvm'_ :: Ptr Float -> Ptr Word32 -> Ptr Float -> Word32 -> Word32 -> IO ()

