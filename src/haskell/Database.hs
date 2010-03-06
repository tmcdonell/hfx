--------------------------------------------------------------------------------
-- |
-- Module    : Database
-- Copyright : (c) 2009 Trevor L. McDonell
-- License   : BSD
--
-- Store a collection of data relating to protein fragments on the device
--
--------------------------------------------------------------------------------

module Database (DevicePDB(..), withPDB) where

import Config
import Protein

import Foreign.Storable
import Data.Word
import Data.Vector                              (Vector)
import Control.Exception.Extensible             (bracket)
import Foreign.CUDA                             (DevicePtr)
import qualified Foreign.CUDA                   as CUDA

import qualified Data.Vector.Storable           as U
import qualified Data.Vector.Generic            as G
import qualified Data.Vector.Generic.Mutable    as M


--------------------------------------------------------------------------------
-- Data Structures
--------------------------------------------------------------------------------

data DevicePDB a = DevicePDB
  {
    peptides     :: Vector (Peptide a), -- the peptide fragments
    rowOffsets   :: DevicePtr Word32,   -- index of the start of each sequence
    yIonLadders  :: DevicePtr a,        -- concatenated y-ion mass ladders for each peptide
    residuals    :: DevicePtr a         -- sum of residual masses for each peptide
  }


--------------------------------------------------------------------------------
-- Marshalling
--------------------------------------------------------------------------------

withPDB :: (Fractional a, Ord a, Storable a)
        => ConfigParams a -> Vector (Protein a) -> (DevicePDB a -> IO b) -> IO b
withPDB cp ps = bracket (newPDB cp ps) freePDB

freePDB :: DevicePDB a -> IO ()
freePDB db = do
  CUDA.free (rowOffsets  db)
  CUDA.free (yIonLadders db)
  CUDA.free (residuals   db)

newFromVector :: Storable a => U.Vector a -> IO (DevicePtr a)
newFromVector v =
  let l = G.length v in
  do
    dp <- CUDA.mallocArray l
    U.unsafeWith v $ \p -> CUDA.pokeArray l p dp
    return dp


--
-- Generate a new protein database on the graphics device.
--
-- This consists of a flattened description of the fragment mass ladders for
-- each peptide in the database. Also includes the list of peptides in an O(1)
-- indexable form, for fast reverse lookup.
--
newPDB :: (Fractional a, Storable a) => ConfigParams a -> Vector (Protein a) -> IO (DevicePDB a)
newPDB cp proteins = do
  resi   <- newFromVector $ G.generate (G.length peps) (\i -> residual (peps G.! i))
  rowP   <- newFromVector offsets

  mv     <- M.unsafeNew (fromIntegral (G.last offsets))
  let fill i  =
        let o = fromIntegral (offsets G.! i)
            s = bIonLadder cp (peps G.! i)
        in mapM_ (\n -> M.unsafeWrite mv (o+n) (s G.! n)) [0 .. G.length s - 1]

  mapM_ fill [ 0 .. G.length peps - 1 ]
  ladder <- newFromVector =<< G.unsafeFreeze mv

  return $ DevicePDB peps rowP ladder resi
  where
    peps    = G.concatMap fragments proteins
    offsets = G.scanl' (+) 0 $ G.generate (G.length peps) seql
--    offsets = G.generate (G.length peps + 1) $ \i ->
--      case i of
--        0 -> 0
--        n -> (offsets G.! (n-1) + (seql (n-1)))

    seql i  = fromIntegral . (\(c,n) -> n-c) . terminals $ peps G.! i

