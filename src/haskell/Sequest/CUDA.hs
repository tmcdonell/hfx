{-# LANGUAGE ForeignFunctionInterface #-}
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
import Control.Monad                    (when)
import Control.Exception.Extensible     (assert)
import Data.List
import System.IO
import Foreign
import Foreign.C
import Foreign.CUDA (DevicePtr, withDevicePtr)
import qualified Foreign.CUDA           as C
import qualified Data.Vector.Storable   as V


--------------------------------------------------------------------------------
-- Database search, CUDA version
--------------------------------------------------------------------------------

searchForMatches :: ConfigParams CFloat -> ProteinDatabase CFloat -> Spectrum CFloat -> IO (MatchCollection CFloat)
searchForMatches cp db spec = do
  -- setup
  --
  t1 <- getTime
  C.withListArray offsets            $ \d_rowPtr -> do
  C.withListArray col                $ \d_colIdx -> do
  C.withListArray val                $ \d_data   -> do
  C.withListArray ident              $ \d_i      -> do
  C.withListArray (V.toList specExp) $ \d_x      -> do
  C.allocaArray   num_peptides       $ \d_y      -> do
    t2 <- getTime
    let tc = elapsedTime t1 t2
    vprint cp ["Setup: " ++ showTime tc, showRate mb_copy tc " MB", shows num_peptides " peptides", shows num_nonzeros " non-zeros"]

    -- score
    --
    (tm,_) <- bracketTime $ cu_smvm_f d_y d_x d_data d_rowPtr d_colIdx num_peptides
    (ts,_) <- bracketTime $ cu_sort_f d_y d_i num_peptides

    vprint cp ["SMVM:  " ++ showTime tm, showRate gflops tm " GFLOPS", showRate mb_smvm tm " MB"]
    vprint cp ["Sort:  " ++ showTime ts, showRate mpairs ts " MPairs"]

    -- retrieve results
    --
    r <- C.peekListArray n (d_y `C.advanceDevPtr` (num_peptides-n))
    i <- C.peekListArray n (d_i `C.advanceDevPtr` (num_peptides-n))
    return $ reverse (zipWith k i r)
  where
    n            = max (numMatches cp) (numMatchesDetail cp)
    k i r        = Match (concatMap fragments candidates !! cIntConv i) (r/10000)
    bnds         = (0, fromIntegral (V.length specExp -1))
    offsets      = scanl (+) 0 . map (fromIntegral . length) $ specThry
    (col,val)    = unzip . concat $ specThry
    num_peptides = length offsets - 1
    ident        = enumFromTo 0 (cIntConv num_peptides - 1) :: [CUInt]

    candidates   = findCandidates cp spec . map (digestProtein cp) $ db
    specExp      = buildExpSpecXCorr cp spec
    specThry     = [ buildThrySpecXCorr cp (charge spec) bnds peptide
                      | protein <- candidates
                      , peptide <- fragments protein ]

    -- misc / info
    --
    mb_copy      = (\x -> fromIntegral x/1E6)
                 $ (length offsets + length col + num_peptides) * sizeOf (undefined::CUInt)
                 + (V.length specExp + length val)              * sizeOf (undefined::CFloat)

    num_nonzeros = fromIntegral $ last offsets
    mb_smvm      = (\x -> fromIntegral x/1E6)
                 $ 2 * num_peptides * sizeOf (undefined::CUInt)       -- row pointer
                 + 1 * num_nonzeros * sizeOf (undefined::CUInt)       -- column index
                 + 2 * num_nonzeros * sizeOf (undefined::CFloat)      -- values, A(i,j) and x(j)
                 + 1 * num_peptides * sizeOf (undefined::CFloat)      -- y(i) = ...

    mpairs       = fromIntegral num_peptides / 1E6
    gflops       = (2 * fromIntegral num_nonzeros) / 1E9


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
    cu_sort_f' keys (castPtr vals) (fromIntegral l)


cu_smvm_f :: DevicePtr CFloat -> DevicePtr CFloat -> DevicePtr CFloat -> DevicePtr CUInt -> DevicePtr CUInt -> Int -> IO ()
cu_smvm_f d_y d_x d_vals d_ptr d_idx num_rows =
  withDevicePtr d_y    $ \p_y    ->
  withDevicePtr d_x    $ \p_x    ->
  withDevicePtr d_vals $ \p_vals ->
  withDevicePtr d_ptr  $ \p_ptr  ->
  withDevicePtr d_idx  $ \p_idx  -> do
    cu_smvm_f' p_y p_x p_vals p_ptr p_idx (fromIntegral num_rows)


foreign import ccall unsafe "algorithms.h radixsort_f" cu_sort_f'
  :: Ptr CFloat -> Ptr () -> CUInt -> IO ()

foreign import ccall unsafe "algorithms.h smvm_f" cu_smvm_f'
  :: Ptr CFloat -> Ptr CFloat -> Ptr CFloat -> Ptr CUInt -> Ptr CUInt -> CUInt -> IO ()

