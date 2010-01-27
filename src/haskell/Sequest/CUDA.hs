{-# LANGUAGE ForeignFunctionInterface, BangPatterns #-}
--------------------------------------------------------------------------------
-- |
-- Module    : Sequest.CUDA
-- Copyright : (c) 2009 Trevor L. McDonell
-- License   : BSD
--
--------------------------------------------------------------------------------

module Sequest.CUDA (searchForMatches) where

import Time
import Config
import Protein
import Spectrum
import IonSeries
import Sequest.Base                     (MatchCollection, Match(..), findCandidates)

import Numeric
import Control.Monad                    (when, foldM)
import Control.Exception.Extensible     (assert)
import Data.List
import Data.Vector.Storable             (Vector(..))
import System.IO
import Foreign
import Foreign.C
import Foreign.CUDA (DevicePtr, withDevicePtr)

import qualified Foreign.CUDA           as C
import qualified Data.Vector.Storable   as V


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
  C.allocaArray (numRows+1) $ \d_rowPtr -> do
  C.allocaArray numNZ       $ \d_colIdx -> do
  C.allocaArray numNZ       $ \d_vals   -> do
    C.pokeArray (numRows+1) rowPtr d_rowPtr >> free rowPtr
    C.pokeArray numNZ       colIdx d_colIdx >> free colIdx
    C.pokeArray numNZ       vals   d_vals   >> free vals

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

searchForMatches :: ConfigParams CFloat -> ProteinDatabase CFloat -> Spectrum CFloat -> IO (MatchCollection CFloat)
searchForMatches cp db spec = do
  tc1 <- getTime

  withDB specThry       $ \numRows numNZ d_rowPtr d_colIdx d_vals -> do
  withVector specExp    $ \d_x -> do
  C.allocaArray numRows $ \d_i -> do
  C.allocaArray numRows $ \d_y -> do

    -- Estimate the time and amount of data used by each kernel, abusing the
    -- fact that all types on the device must be 32-bit. This is only used in
    -- verbose mode.
    --
    tc2 <- getTime
    let tc  = elapsedTime tc1 tc2
        mbc = fromIntegral (2*numNZ + V.length specExp + numRows + 1) * 4 / 1E6
        mbm = fromIntegral (3*numNZ + 3*numRows) * 4 / 1E6
        gfm = fromIntegral (2*numNZ) / 1E9
        mps = fromIntegral numRows / 1E6

    vprint cp ["Setup: " ++ showTime tc, showRate mbc tc " MB", shows numRows " peptides", shows numNZ " non-zeros"]

    -- score each candidate peptide, and sort based on that ranking
    --
    (tm,_) <- bracketTime $ cu_smvm_f d_y d_x d_vals d_rowPtr d_colIdx numRows
    (te,_) <- bracketTime $ cu_enum_i d_i 0 (numRows-1)
    (ts,_) <- bracketTime $ cu_sort_f d_y d_i numRows

    vprint cp ["SMVM:  " ++ showTime tm, showRate gfm tm " GFLOPS", showRate mbm tm " MB"]
    vprint cp ["Enum:  " ++ showTime te]
    vprint cp ["Sort:  " ++ showTime ts, showRate mps ts " MPairs"]

    -- retrieve the highest ranked results (ascending sort)
    --
    sc <- C.peekListArray n (d_y `C.advanceDevPtr` (numRows-n))
    ix <- C.peekListArray n (d_i `C.advanceDevPtr` (numRows-n))
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

showRate :: Double -> Time -> String -> String
showRate _ (Time 0) u = "--" ++ u ++ "/s"
showRate n t        u = showFFloat (Just 3) (1000 * n / (fromInteger (timeIn millisecond t))) (u++"/s")

cIntConv :: (Integral a, Integral b) => a -> b
cIntConv =  fromIntegral

withVector :: Storable a => Vector a -> (DevicePtr a -> IO b) -> IO b
withVector (Vector _ l p) f =
  C.allocaArray l $ \dp -> do
    withForeignPtr p $ \p' -> C.pokeArray l p' dp
    f dp


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

