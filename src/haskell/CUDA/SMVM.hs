{-# LANGUAGE ForeignFunctionInterface, BangPatterns #-}
--------------------------------------------------------------------------------
-- |
-- Module    : CUDA.SMVM
-- Copyright : (c) 2009 Trevor L. McDonell
-- License   : BSD
--
-- Sequest matching based on sparse-matrix vector multiplication
--
--------------------------------------------------------------------------------

module CUDA.SMVM (searchForMatches) where

import Time
import Config
import Protein
import Spectrum
import IonSeries
import Sequest                          (MatchCollection, Match(..), findCandidates)

import Control.Monad                    (when, foldM)
import Control.Exception.Extensible     (assert)
import Data.List
import System.IO
import Foreign
import Foreign.C.Types
import Foreign.CUDA                     (DevicePtr, withDevicePtr)

import qualified Foreign.CUDA           as CUDA
import qualified Data.Vector.Generic    as V
import qualified Data.Vector.Storable   as S


--------------------------------------------------------------------------------
-- Data marshalling
--------------------------------------------------------------------------------

data Vec a = Vec {-# UNPACK #-} !Int            -- length
                 {-# UNPACK #-} !(Ptr a)        -- payload

newV :: Storable a => Int -> IO (Vec a)
newV l = Vec l `fmap` mallocArray l

--
-- Write an element, expanding the vector if necessary
--
writeV :: Storable a => Vec a -> Int -> a -> IO (Vec a)
writeV !(Vec l p) n e
  | n < l     = do pokeElemOff p n e >> return (Vec l p)
  | otherwise = do let l' = ceiling (growthFactor * fromIntegral l)
                   p' <- reallocArray p l'
                   pokeElemOff p' n e
                   return (Vec l' p')
  where
    growthFactor = 1.5 :: Double

--
-- bracketed allocating/deallocating of the marshalling structures
--
withDB :: [XCorrSpecThry CUInt CFloat]
         -> (Int -> Int -> DevicePtr CUInt -> DevicePtr CUInt -> DevicePtr CFloat -> IO b)
         -> IO b
withDB ss f = do
  (numRows, numNZ, (Vec _ rowPtr), (Vec _ colIdx), (Vec _ vals)) <- copyDB ss
  CUDA.allocaArray (numRows+1) $ \d_rowPtr -> do
  CUDA.allocaArray numNZ       $ \d_colIdx -> do
  CUDA.allocaArray numNZ       $ \d_vals   -> do
    CUDA.pokeArray (numRows+1) rowPtr d_rowPtr >> free rowPtr
    CUDA.pokeArray numNZ       colIdx d_colIdx >> free colIdx
    CUDA.pokeArray numNZ       vals   d_vals   >> free vals

    f numRows numNZ d_rowPtr d_colIdx d_vals

--
-- Marshal the intensity spectra into a compressed sparse-row format suitable
-- for processing on the device, without requiring the lengths to be known
-- a-priori.
--
copyDB :: [XCorrSpecThry CUInt CFloat] -> IO (Int, Int, Vec CUInt, Vec CUInt, Vec CFloat)
copyDB ss = do
  rowPtr <- newV 10000          -- number of candidate peptides in a database
  colIdx <- newV 1400000        -- approx. 140 peaks per theoretical spectra
  vals   <- newV 1400000
  finish =<< foldM k (0, 0, rowPtr, colIdx, vals) ss

  where
    finish (!nr,!nnz,!rv,!cv,!vv) = do
      rv' <- writeV rv nr (cIntConv nnz)
      return (nr,nnz,rv',cv,vv)

    k (!nr,!nnz,!rv,!cv,!vv) s = do
      rv'            <- writeV rv nr (cIntConv nnz)
      (nnz',cv',vv') <- foldM k' (nnz,cv,vv) s
      return (nr+1,nnz',rv',cv',vv')

    k' (!nnz,!cv,!vv) (x,y) = do
      cv' <- writeV cv nnz x
      vv' <- writeV vv nnz y
      return (nnz+1,cv',vv')


--------------------------------------------------------------------------------
-- Database search, CUDA version
--------------------------------------------------------------------------------

searchForMatches :: ConfigParams CFloat -> [Protein CFloat] -> Spectrum CFloat -> IO (MatchCollection CFloat)
searchForMatches cp db spec = do
  tc1 <- getTime

  withDB specThry       $ \numRows numNZ d_rowPtr d_colIdx d_vals -> do
  withVector specExp    $ \d_x -> do
  CUDA.allocaArray numRows $ \d_i -> do
  CUDA.allocaArray numRows $ \d_y -> do

    -- Estimate the time and amount of data used by each kernel, abusing the
    -- fact that all types on the device must be 32-bit. This is only used in
    -- verbose mode.
    --
    tc2 <- getTime
    let tc  = elapsedTime tc1 tc2
        mbc = (2*numNZ + V.length specExp + numRows + 1) * 4
        mbm = (3*numNZ + 3*numRows) * 4
        gfm = 2*numNZ
        mps = numRows

    vprint cp ["Setup: " ++ showTime tc, showRateSI mbc tc "b", shows numRows " peptides", shows numNZ " non-zeros"]

    -- score each candidate peptide, and sort based on that ranking
    --
    (tm,_) <- bracketTime $ cu_smvm_f d_y d_x d_vals d_rowPtr d_colIdx numRows
    (te,_) <- bracketTime $ cu_enum_i d_i 0 (numRows-1)
    (ts,_) <- bracketTime $ cu_sort_f d_y d_i numRows

    vprint cp ["SMVM:  " ++ showTime tm, showRateSI gfm tm "FLOPS", showRateSI mbm tm "b"]
    vprint cp ["Enum:  " ++ showTime te]
    vprint cp ["Sort:  " ++ showTime ts, showRateSI mps ts "Pairs"]

    -- retrieve the highest ranked results (ascending sort)
    --
    sc <- CUDA.peekListArray n (d_y `CUDA.advanceDevPtr` (numRows-n))
    ix <- CUDA.peekListArray n (d_i `CUDA.advanceDevPtr` (numRows-n))
    finish sc ix
  where
    specExp      = buildExpSpecXCorr cp spec
    specThry     = map (buildThrySpecXCorr cp (charge spec) bnds) candidates
    candidates   = concatMap fragments . findCandidates cp spec . map (digestProtein cp) $ db

    finish sc ix = return . reverse $ zipWith (\s i -> Match (candidates !! cIntConv i) (s/10000)) sc ix
    n            = max (numMatches cp) (numMatchesDetail cp)
    bnds         = (0, cIntConv (V.length specExp -1))


--------------------------------------------------------------------------------
-- Misc/Info
--------------------------------------------------------------------------------

vprint :: ConfigParams a -> [String] -> IO ()
vprint cp s = when (verbose cp) (hPutStrLn stderr . concat . intersperse ", " $ s)

cIntConv :: (Integral a, Integral b) => a -> b
cIntConv =  fromIntegral

withVector :: Storable a => S.Vector a -> (DevicePtr a -> IO b) -> IO b
withVector vec action = let l = V.length vec in
  S.unsafeWith vec   $ \p  ->
  CUDA.allocaArray l $ \dp -> do
    CUDA.pokeArray l p dp
    action dp


--------------------------------------------------------------------------------
-- FFI Bindings
--------------------------------------------------------------------------------

sizeOfPtr :: Storable a => DevicePtr a -> Int
sizeOfPtr =  sizeOf . (undefined :: DevicePtr a -> a)

cu_sort_f :: Storable a => DevicePtr CFloat -> DevicePtr a -> Int -> IO ()
cu_sort_f d_keys d_vals l =
  assert (sizeOfPtr d_vals == 4) $
  withDevicePtr d_keys $ \keys ->
  withDevicePtr d_vals $ \vals -> do
    cu_sort_f' keys (castPtr vals) (cIntConv l)


cu_smvm_f :: DevicePtr CFloat -> DevicePtr CFloat -> DevicePtr CFloat -> DevicePtr CUInt -> DevicePtr CUInt -> Int -> IO ()
cu_smvm_f d_y d_x d_vals d_ptr d_idx num_rows =
  withDevicePtr d_y    $ \p_y    ->
  withDevicePtr d_x    $ \p_x    ->
  withDevicePtr d_vals $ \p_vals ->
  withDevicePtr d_ptr  $ \p_ptr  ->
  withDevicePtr d_idx  $ \p_idx  -> do
    cu_smvm_f' p_y p_x p_vals p_ptr p_idx (cIntConv num_rows)


cu_enum_i :: DevicePtr CUInt -> Int -> Int -> IO ()
cu_enum_i d_x l u = withDevicePtr d_x $ \p_x -> cu_enum_i' p_x (cIntConv l) (cIntConv u)


foreign import ccall unsafe "algorithms.h radixsort_f" cu_sort_f'
  :: Ptr CFloat -> Ptr () -> CUInt -> IO ()

foreign import ccall unsafe "algorithms.h smvm_f" cu_smvm_f'
  :: Ptr CFloat -> Ptr CFloat -> Ptr CFloat -> Ptr CUInt -> Ptr CUInt -> CUInt -> IO ()

foreign import ccall unsafe "algorithms.h enumFromTo_i" cu_enum_i'
  :: Ptr CUInt -> CUInt -> CUInt -> IO ()

