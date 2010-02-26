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

import Mass
import C2HS
import Config
import Spectrum
import Sequest                                  (MatchCollection, Match(..))
import Database                                 (ProteinDatabase(..))
import qualified Foreign.CUDA.Algorithms        as CUDA

import Foreign
import Foreign.CUDA                             (DevicePtr)
import qualified Foreign.CUDA                   as CUDA
import qualified Data.Vector.Generic            as G
import qualified Data.Vector.Storable           as U


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
    nr <- CUDA.findIndicesInRange (residuals db) d_pepIdx numPeptides (mass-delta) (mass+delta)

    -- Generate theoretical spectra for each of these candidates
    --
    CUDA.allocaArray nr             $ \d_score    -> do
    CUDA.allocaArray (nr * specLen) $ \d_specThry -> do
      CUDA.memset  d_specThry (fromIntegral  (specLen * nr * sizeOfPtr d_specThry)) 0
      CUDA.addIons d_specThry (residuals db) (yIonLadders db) (rowOffsets db) d_pepIdx nr (round chrg) specLen

      -- Score each peptide
      -- Sort results and retrieve the best matches (reverse order)
      --
      CUDA.mvm d_score d_specThry d_specExp nr specLen
      CUDA.radixsort d_score d_pepIdx nr

      sc <- CUDA.peekListArray numResults (d_score  `CUDA.advanceDevPtr` (nr-numResults))
      ix <- CUDA.peekListArray numResults (d_pepIdx `CUDA.advanceDevPtr` (nr-numResults))

      return . reverse $ zipWith (\s i -> Match (peptides db G.! cIntConv i) (s/10000)) sc ix
  where
    numPeptides  = G.length (peptides db)
    numResults   = max (numMatches cp) (numMatchesDetail cp)

    spec         = buildExpSpecXCorr cp dta
    specLen      = G.length spec
    chrg         = max 1 (charge dta - 1)

    -- Residual mass plus the mass of the water molecule released when forming
    -- the peptide bond
    --
    mass         = (precursor dta * charge dta) - ((charge dta * massH) - 1) - (massH + massH2O)
    delta        = massTolerance cp


--------------------------------------------------------------------------------
-- Utils
--------------------------------------------------------------------------------

sizeOfPtr :: Storable a => DevicePtr a -> Int
sizeOfPtr =  sizeOf . (undefined :: DevicePtr a -> a)

withVector :: Storable a => U.Vector a -> (DevicePtr a -> IO b) -> IO b
withVector vec action = let l = G.length vec in
  U.unsafeWith vec   $ \p  ->
  CUDA.allocaArray l $ \dp -> do
    CUDA.pokeArray l p dp
    action dp

