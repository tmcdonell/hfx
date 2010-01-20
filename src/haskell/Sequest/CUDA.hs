{-# LANGUAGE ForeignFunctionInterface #-}
--------------------------------------------------------------------------------
-- |
-- Module    : Sequest.CUDA
-- Copyright : (c) 2009 Trevor L. McDonell
-- License   : BSD
--
--------------------------------------------------------------------------------

module Sequest.CUDA where

import Time
import Config
import Protein
import Spectrum
import IonSeries
import Sequest.Base                     (MatchCollection, Match(..), findCandidates)

import Control.Monad                    (when)
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
  t1 <- getTime

  -- setup
  --
  C.withListArray offsets                   $ \d_rowPtr -> do
  C.withListArray ix                        $ \d_colIdx -> do
  C.withListArray val                       $ \d_data   -> do
  C.withListArray (V.toList specExp)        $ \d_x      -> do
  C.allocaArray num_peptides                $ \d_y      -> do
  C.withListArray [0 .. (num_peptides - 1)] $ \d_i      -> do

    t2 <- getTime
    when (verbose cp) (hPutStrLn stderr $ "Setup: " ++ showTime (elapsedTime t1 t2))

    -- score
    --
    (tm,_) <- bracketTime $ cu_smvm_f d_y d_x d_data d_rowPtr d_colIdx num_peptides
    (ts,_) <- bracketTime $ cu_sort_f d_y d_i num_peptides

    when (verbose cp) (hPutStrLn stderr $ "SMVM:  " ++ showTime tm)
    when (verbose cp) (hPutStrLn stderr $ "Sort:  " ++ showTime ts)

    -- retrieve results
    --
    r <- C.peekListArray n (d_y `C.plusDevPtr` ((num_peptides-n) * sizeOf (undefined::CFloat)))
    i <- C.peekListArray n (d_i `C.plusDevPtr` ((num_peptides-n) * sizeOf (undefined::CFloat)))
    return $ reverse (zipWith k i r)
  where
    n            = max (numMatches cp) (numMatchesDetail cp)
    k i r        = Match (concatMap fragments candidates !! i) (r/10000)
    bnds         = (0, fromIntegral (V.length specExp -1))
    offsets      = scanl (+) 0 . map (fromIntegral . length) $ specThry
    (ix,val)     = unzip . concat $ specThry
    num_peptides = length offsets - 1

    candidates   = findCandidates cp spec . map (digestProtein cp) $ db
    specExp      = buildExpSpecXCorr cp spec
    specThry     = [ buildThrySpecXCorr cp (charge spec) bnds peptide
                      | protein <- candidates
                      , peptide <- fragments protein ]


--------------------------------------------------------------------------------
-- FFI Bindings
--------------------------------------------------------------------------------

cu_sort_f :: DevicePtr CFloat -> DevicePtr a -> Int -> IO ()
cu_sort_f d_keys d_vals l =
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


foreign import ccall unsafe "algorithms.h sort_f" cu_sort_f'
  :: Ptr CFloat -> Ptr () -> CUInt -> IO ()

foreign import ccall unsafe "algorithms.h smvm_f" cu_smvm_f'
  :: Ptr CFloat -> Ptr CFloat -> Ptr CFloat -> Ptr CUInt -> Ptr CUInt -> CUInt -> IO ()

