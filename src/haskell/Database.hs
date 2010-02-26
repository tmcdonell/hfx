--------------------------------------------------------------------------------
-- |
-- Module    : Database
-- Copyright : (c) 2009 Trevor L. McDonell
-- License   : BSD
--
-- Store a collection of data relating to protein fragments on the device
--
--------------------------------------------------------------------------------

module Database (ProteinDatabase(..), withPDB) where

import Time
import Config
import Protein

import Foreign.Storable
import System.IO
import Data.Word
import Data.List                                (intersperse)
import Data.Vector                              (Vector)
import Control.Monad                            (when)
import Control.Exception.Extensible             (bracket)
import Foreign.CUDA                             (DevicePtr)
import qualified Foreign.CUDA                   as CUDA

import qualified Data.Vector.Storable           as U
import qualified Data.Vector.Generic            as G
import qualified Data.Vector.Generic.Mutable    as M


--------------------------------------------------------------------------------
-- Data Structures
--------------------------------------------------------------------------------

data ProteinDatabase a = ProteinDatabase
  {
    peptides     :: Vector (Peptide a), -- the peptide fragments
    rowOffsets   :: DevicePtr Word32,   -- index of the start of each sequence
    yIonLadders  :: DevicePtr a,        -- concatenated y-ion mass ladders for each peptide
    residuals    :: DevicePtr a         -- sum of residual masses for each peptide
  }


--------------------------------------------------------------------------------
-- Marshalling
--------------------------------------------------------------------------------

withPDB :: (Fractional a, Ord a, Storable a) => ConfigParams a -> Vector (Protein a) -> (ProteinDatabase a -> IO b) -> IO b
withPDB cp ps = flip bracket freePDB $ do
  (t,db) <- bracketTime (newPDB cp ps)
  whenVerbose cp [ "> setup: " ++ showTime t, shows (G.length (peptides db)) " peptides" ]
  return db

freePDB :: ProteinDatabase a -> IO ()
freePDB db = do
  CUDA.free (rowOffsets  db)
  CUDA.free (yIonLadders db)
  CUDA.free (residuals   db)

whenVerbose      :: ConfigParams a -> [String] -> IO ()
whenVerbose cp s =  when (verbose cp) (hPutStrLn stderr . concat . intersperse ", " $ s)

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
newPDB :: (Fractional a, Storable a) => ConfigParams a -> Vector (Protein a) -> IO (ProteinDatabase a)
newPDB cp proteins = do
  resi   <- newFromVector $ G.generate (G.length peps) (\i -> residual (peps G.! i))
  rowP   <- newFromVector offsets

  mv     <- M.unsafeNew (fromIntegral $ G.last offsets)
  let fill i  =
        let o = fromIntegral (offsets G.! i)
            s = bIonLadder cp (peps G.! i)
        in mapM_ (\n -> M.unsafeWrite mv (o+n) (s G.! n)) [0 .. G.length s - 1]

  mapM_ fill [ 0 .. G.length peps - 1 ]
  ladder <- newFromVector =<< G.unsafeFreeze mv

  return $ ProteinDatabase peps rowP ladder resi
  where
    peps    = G.concatMap fragments proteins
    offsets = G.scanl' (+) 0 $ G.generate (G.length peps) seql
    seql i  = fromIntegral . (\(c,n) -> n-c) . terminals $ peps G.! i

