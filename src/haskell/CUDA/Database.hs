{-# LANGUAGE ForeignFunctionInterface #-}
--------------------------------------------------------------------------------
-- |
-- Module    : CUDA.Database
-- Copyright : (c) 2009 Trevor L. McDonell
-- License   : BSD
--
-- Store a collection of data relating to protein fragments on the device
--
--------------------------------------------------------------------------------

module CUDA.Database (ProteinDatabase(..), withPDB) where

import Config
import Protein

import Foreign
import Foreign.C
import Control.Exception.Extensible             (assert, bracket)
import Foreign.CUDA                             (DevicePtr, withDevicePtr)
import Data.Vector                              (Vector)
import qualified Data.Vector                    as V
import qualified Foreign.CUDA                   as C
import qualified Data.ByteString.Lazy.Char8     as L

import Unsafe.Coerce


--------------------------------------------------------------------------------
-- Data Structures
--------------------------------------------------------------------------------

data ProteinDatabase a = ProteinDatabase
  {
    peptides     :: Vector (Peptide a), -- The peptide fragments
    rowOffsets   :: DevicePtr CUInt,    -- index of each of the head flags
    headFlags    :: DevicePtr CUInt,    -- start of a peptide block
    yIonLadders  :: DevicePtr a,        -- mass ladder, y-ion
    residuals    :: DevicePtr a         -- sum of residual masses for each peptide
  }


--------------------------------------------------------------------------------
-- Marshalling
--------------------------------------------------------------------------------

withPDB :: (Fractional a, Ord a, Storable a) => ConfigParams a -> [Protein a] -> (ProteinDatabase a -> IO b) -> IO b
withPDB cp ps = bracket (newPDB cp ps) freePDB

freePDB :: ProteinDatabase a -> IO ()
freePDB (ProteinDatabase _  rp hf yi rm) =
  C.free rp >> C.free hf >> C.free yi >> C.free rm


--
-- Generate a new protein database on the graphics device
--
newPDB :: (Fractional a, Ord a, Storable a) => ConfigParams a -> [Protein a] -> IO (ProteinDatabase a)
newPDB cp ps = do

  -- The offset index for the start of each peptide sequence, and corresponding
  -- head-flags array.
  --
  rowP  <- C.newListArray . scanl (+) 0 . map (cIntConv . length) $ seqs
  headF <- C.mallocArray seqLen

  C.allocaArray numPeps $ \ones -> do
    C.memset headF (fromIntegral (seqLen * sizeOfPtr headF)) 0
    cu_replicate ones 1 numPeps
    cu_permute ones headF rowP numPeps

  -- Copy the y-ion mass ladder and corresponding total residual mass. This
  -- could also be calculated from the ion masses with segscan / bpermute.
  --
  yIons <- C.newListArray . concatMap (scanr1 (+)) $ seqs
  resi  <- C.newListArray . map sum $ seqs

  return (ProteinDatabase peps rowP headF yIons resi)
  where
    -- Do some extra conversion work to return the peptides in an O(1) indexable form
    --
    peps    = V.fromList . concatMap fragments . map (digestProtein cp) $ ps
    numPeps = V.length peps
    seqs    = V.toList . V.map (map (getAAMass cp) . L.unpack . lyse) $ peps
    seqLen  = sum . map length $ seqs


--------------------------------------------------------------------------------
-- FFI Bindings
--------------------------------------------------------------------------------

sizeOfPtr :: Storable a => DevicePtr a -> Int
sizeOfPtr =  sizeOf . (undefined :: DevicePtr a -> a)

cIntConv :: (Integral a, Integral b) => a -> b
cIntConv =  fromIntegral

-- Permute
--
cu_permute :: Storable a => DevicePtr a -> DevicePtr a -> DevicePtr CUInt -> Int -> IO ()
cu_permute d_src d_dst d_idx n =
  assert (sizeOfPtr d_src == 4) $
  withDevicePtr d_src $ \src ->
  withDevicePtr d_dst $ \dst ->
  withDevicePtr d_idx $ \idx ->
    cu_permute' (castPtr src) (castPtr dst) idx (cIntConv n)

foreign import ccall unsafe "algorithms.h permute" cu_permute'
  :: Ptr () -> Ptr () -> Ptr CUInt -> CUInt -> IO ()

-- Replicate
--
cu_replicate :: Storable a => DevicePtr a -> a -> Int -> IO ()
cu_replicate d_xs x n =
  assert (sizeOf x == 4) $
  withDevicePtr d_xs     $ \xs ->
    cu_replicate' (castPtr xs) (unsafeCoerce x) (cIntConv n)

foreign import ccall unsafe "algorithms.h replicate" cu_replicate'
  :: Ptr () -> CUInt -> CUInt -> IO ()

