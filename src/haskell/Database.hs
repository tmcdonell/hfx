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

withPDB :: (Fractional a, Ord a, Storable a) => ConfigParams a -> [Protein a] -> (ProteinDatabase a -> IO b) -> IO b
withPDB cp ps = flip bracket freePDB $ do
  (t,db) <- bracketTime (newPDB cp ps)
  whenVerbose cp [ "> setup: " ++ showTime t ]
  return db

freePDB :: ProteinDatabase a -> IO ()
freePDB db = do
  CUDA.free (rowOffsets  db)
  CUDA.free (yIonLadders db)
  CUDA.free (residuals   db)

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

whenVerbose      :: ConfigParams a -> [String] -> IO ()
whenVerbose cp s =  when (verbose cp) (hPutStrLn stderr . concat . intersperse ", " $ s)

