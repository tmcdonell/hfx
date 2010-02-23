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

import C2HS
import Time
import Config
import Protein

import Foreign.Storable
import System.IO
import Data.Word
import Data.List                                (intersperse)
import Control.Monad                            (when)
import Control.Exception.Extensible             (bracket)
import Foreign.CUDA                             (DevicePtr)
import Data.Vector                              (Vector, fromList)
import qualified Foreign.CUDA                   as CUDA
import qualified Data.ByteString.Lazy.Char8     as L


{-
import System.IO.Unsafe
import Control.Monad.State
import Foreign
import Foreign.CUDA                             (DevicePtr)
import qualified Foreign.CUDA                   as CUDA

data Builder a = B !Int !Int !(Ptr a)   -- number of elements, maximum array size, payload

newListArrayLen :: Storable a => [a] -> IO (DevicePtr a, Int)
newListArrayLen vals = finish =<< foldM loop initialState vals
  where
    initialState = B 0 0 nullPtr
    initialSize  = 4 * 1024 * 1024
    growthFactor = 1.25 :: Double

    finish (B n _ h_arr) = do
        d_arr <- CUDA.mallocArray n
        CUDA.pokeArray n h_arr d_arr >> free h_arr
        return (d_arr, n)

    loop (B n l arr) x = unsafeInterleaveIO $
        if n < l
          then pokeElemOff arr n x >> return (B (n+1) l arr)
          else let l' = max initialSize (ceiling (growthFactor * fromIntegral l)) in
               reallocArray arr l' >>= \arr' -> loop (B n l' arr') x
-}
{-
newListArrayLen :: Storable a => [a] -> IO (DevicePtr a, Int)
newListArrayLen vals = evalStateT (loop vals) initialState
  where
    initialState = B 0 0 nullPtr
    initialSize  = 4 * 1024 * 1024      -- 16 megabytes for 32-bit words
    growthFactor = 1.5 :: Double

    -- Finished consuming the list; copy to the device and free host array
    --
    loop []      = do
      (B n _ arr) <- get
      d_arr       <- liftIO $ CUDA.mallocArray n
      h_arr       <- liftIO $ evaluate arr
      liftIO $ CUDA.pokeArray n h_arr d_arr >> free h_arr
      return (d_arr, n)

    -- If the temporary array is large enough, add the next element.
    -- Otherwise resize the array and try again.
    --
    loop (x:xs)  = do
      (B n l arr) <- get
      if n < l then do liftIO . unsafeInterleaveIO $ pokeElemOff arr n x
                       put (B (n+1) l arr)
                       loop xs
               else do let l' = max initialSize (ceiling (growthFactor * fromIntegral l))
                       arr'  <- liftIO $ reallocArray arr l'
                       put (B n l' arr')
                       loop (x:xs)
-}


--------------------------------------------------------------------------------
-- Data Structures
--------------------------------------------------------------------------------

data ProteinDatabase a = ProteinDatabase
  {
    peptides     :: Vector (Peptide a), -- The peptide fragments
    rowOffsets   :: DevicePtr Word32,   -- index of each of the head flags
    yIonLadders  :: DevicePtr a,        -- mass ladder, y-ion
    residuals    :: DevicePtr a         -- sum of residual masses for each peptide
  }


--------------------------------------------------------------------------------
-- Marshalling
--------------------------------------------------------------------------------

withPDB :: (Fractional a, Ord a, Storable a) => ConfigParams a -> [Protein a] -> (ProteinDatabase a -> IO b) -> IO b
withPDB cp ps = flip bracket freePDB $ do
  (t,db) <- bracketTime (newPDB cp ps)
  whenVerbose cp [ "> setup: " ++ showTime t ]
  return db

freePDB :: ProteinDatabase a -> IO ()
freePDB (ProteinDatabase _  rp yi rm) =
  CUDA.free rp >> CUDA.free yi >> CUDA.free rm

--
-- Generate a new protein database on the graphics device.
--
-- This consists of a flattened description of the fragment mass ladders for
-- each peptide in the database. This is the y-ion ladder ladder (C->N terminus)
-- except that we also include the mass of the unbroken peptide in the sequence.
--
-- Also includes the list of peptides in an O(1) indexable form, for fast
-- reverse lookup.
--
newPDB :: (Fractional a, Ord a, Storable a) => ConfigParams a -> [Protein a] -> IO (ProteinDatabase a)
newPDB cp proteins = do
  rowP   <- CUDA.newListArray . scanl (+) 0 . map (cIntConv . length) $ seqs
  ladder <- CUDA.newListArray . concatMap (scanr1 (+))                $ seqs
  resi   <- CUDA.newListArray . map sum                               $ seqs

  return $ ProteinDatabase (fromList peps) rowP ladder resi
  where
    peps = concatMap fragments . map (digestProtein cp) $ proteins
    seqs = map (map (getAAMass cp) . L.unpack . lyse)   $ peps

--
-- Generate a new protein database on the graphics device
--


{-
newPDB :: (Fractional a, Ord a, Storable a) => ConfigParams a -> [Protein a] -> IO (ProteinDatabase a)
newPDB cp ps = do
  t1 <- getTime

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

  t2 <- getTime
  let t = elapsedTime t1 t2
      b = fromIntegral $ (2*numPeps + 2*seqLen + 1) * sizeOf (undefined::Word32)
  whenVerbose cp [ "> setup: " ++ showTime t , shows numPeps " peptides" , showFFloatSI (b::Double) "B" ]

  return (ProteinDatabase peps rowP headF yIons resi)
  where
    -- Do some extra conversion work to return the peptides in an O(1) indexable form
    --
    peps    = V.fromList . concatMap fragments . map (digestProtein cp) $ ps
    numPeps = V.length peps
    seqs    = V.toList . V.map (map (getAAMass cp) . L.unpack . lyse) $ peps
    seqLen  = sum . map length $ seqs
-}


--------------------------------------------------------------------------------
-- FFI Bindings
--------------------------------------------------------------------------------

whenVerbose      :: ConfigParams a -> [String] -> IO ()
whenVerbose cp s =  when (verbose cp) (hPutStrLn stderr . concat . intersperse ", " $ s)

{-
sizeOfPtr :: Storable a => DevicePtr a -> Int
sizeOfPtr =  sizeOf . (undefined :: DevicePtr a -> a)

-- Permute
--
cu_permute :: Storable a => DevicePtr a -> DevicePtr a -> DevicePtr Word32 -> Int -> IO ()
cu_permute d_src d_dst d_idx n =
  assert (sizeOfPtr d_src == 4) $
  withDevicePtr d_src $ \src ->
  withDevicePtr d_dst $ \dst ->
  withDevicePtr d_idx $ \idx ->
    cu_permute'_ (castPtr src) (castPtr dst) idx (cIntConv n)

foreign import ccall unsafe "algorithms.h permute"
  cu_permute'_ :: Ptr () -> Ptr () -> Ptr Word32 -> Word32 -> IO ()

-- Replicate
--
cu_replicate :: Storable a => DevicePtr a -> a -> Int -> IO ()
cu_replicate d_xs x n =
  assert (sizeOf x == 4) $
  withDevicePtr d_xs     $ \xs ->
    cu_replicate'_ (castPtr xs) (unsafeCoerce x) (cIntConv n)

foreign import ccall unsafe "algorithms.h replicate"
  cu_replicate'_ :: Ptr () -> Word32 -> Word32 -> IO ()

-- Segmented scan
--
cu_postsegscanr_plusf :: DevicePtr Float -> DevicePtr Word32 -> DevicePtr Float -> Int -> IO ()
cu_postsegscanr_plusf d_src d_flags d_dst n =
  withDevicePtr d_src   $ \d_src'   ->
  withDevicePtr d_dst   $ \d_dst'   ->
  withDevicePtr d_flags $ \d_flags' ->
    cu_postsegscanr_plusf'_ d_src' d_flags' d_dst' (cIntConv n)

foreign import ccall unsafe "algorithms.h postsegscanr_plusf"
  cu_postsegscanr_plusf'_ :: Ptr Float -> Ptr Word32 -> Ptr Float -> Word32 -> IO ()
-}
