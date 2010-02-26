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
import Config
import Spectrum
import Sequest                                  (MatchCollection, Match(..))
import Database                                 (ProteinDatabase(..))

import Data.Word
import Data.List                                (intersperse)
import Control.Monad                            (when)
import System.IO
import Foreign
import Foreign.CUDA                             (DevicePtr, withDevicePtr)
import qualified Foreign.CUDA                   as CUDA
import qualified Data.Vector.Generic            as V
import qualified Data.Vector.Storable           as S

import Foreign.CUDA.Algorithms


--------------------------------------------------------------------------------
-- Main
--------------------------------------------------------------------------------

searchForMatches :: ConfigParams Float -> ProteinDatabase Float -> Spectrum Float -> IO (MatchCollection Float)
searchForMatches cp db dta =
  withVector spec              $ \d_specExp ->
  CUDA.allocaArray numPeptides $ \d_pepIdx  -> do

    -- Find peptides in the database with a residual mass close to the spectral
    -- precursor mass
    --
    (tf,nr) <- bracketTime $ findIndicesInRange (residuals db) d_pepIdx numPeptides (mass-delta) (mass+delta)
    whenVerbose cp ["> findIndices: " ++ showTime tf, shows nr " found"]

    -- Generate theoretical spectra for each of these candidates
    --
    CUDA.allocaArray nr             $ \d_score    -> do
    CUDA.allocaArray (nr * specLen) $ \d_specThry -> do
      CUDA.memset d_specThry (fromIntegral $ specLen * nr * sizeOfPtr d_specThry) 0
      (ta,_) <- bracketTime $ addIons d_specThry (residuals db) (yIonLadders db) (rowOffsets db) d_pepIdx nr (round chrg) specLen
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

      sc <- CUDA.peekListArray numResults (d_score  `CUDA.advanceDevPtr` (nr-numResults))
      ix <- CUDA.peekListArray numResults (d_pepIdx `CUDA.advanceDevPtr` (nr-numResults))

      finish sc ix
  where
    numPeptides  = V.length (peptides db)
    numResults   = max (numMatches cp) (numMatchesDetail cp)

    spec         = buildExpSpecXCorr cp dta
    specLen      = V.length spec
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
withVector vec action = let l = V.length vec in
  S.unsafeWith vec   $ \p  ->
  CUDA.allocaArray l $ \dp -> do
    CUDA.pokeArray l p dp
    action dp

