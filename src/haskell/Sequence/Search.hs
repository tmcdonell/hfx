{-# LANGUAGE TupleSections #-}
--------------------------------------------------------------------------------
-- |
-- Module    : Sequence.Search
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- An implementation of the SEQUEST algorithm for fast cross-correlation based
-- identification of protein sequences.
--
--
-- References:
--
-- [1] J. K. Eng, B. Fischer, J. Grossmann, and M. J. MacCoss. "A fast sequest
--     cross correlation algorithm." Journal of Proteome Research,
--     7(10):4598-4602, 2008.
--
-- [2] J. K. Eng, A. L. McCormack, and I. John R. Yates. "An approach to
--     correlate tandem mass spectral data of peptides with amino acid sequences
--     in a protein database." Journal of the American Society for Mass
--     Spectrometry, 5(11):976-989, November 1994.
--
--------------------------------------------------------------------------------

module Sequence.Search (searchForMatches) where

import Mass
import Config
import Spectrum
import Sequence.Match
import Sequence.Fragment
import Sequence.Location

import Data.Word
import Data.Maybe
import Control.Monad
import System.IO
import Prelude                                  hiding (lookup)

import Foreign.CUDA (DevicePtr)
import qualified Data.Vector.Generic            as G
import qualified Foreign.CUDA                   as CUDA
import qualified Foreign.CUDA.Util              as CUDA
import qualified Foreign.CUDA.Algorithms        as CUDA


type Candidates = (DevicePtr Word32, Int)
type IonSeries  = DevicePtr Word32

--------------------------------------------------------------------------------
-- Search
--------------------------------------------------------------------------------

--
-- Search an amino acid sequence database to find the most likely matches to a
-- given experimental spectrum.
--
searchForMatches :: ConfigParams -> SequenceDB -> DeviceSeqDB -> MS2Data -> IO MatchCollection
searchForMatches cp sdb ddb ms2 =
  findCandidates cp ddb ms2                                  $ \candidates ->
  mkSpecXCorr ddb candidates (ms2charge ms2) (G.length specExp) $ \specThry   ->
  mapMaybe lookupF `fmap` sequestXC cp candidates specExp specThry
  where
    specExp       = sequestXCorr cp ms2
    lookupF (s,i) = liftM (flip Match s) (lookup sdb i)


--
-- Search the protein database for candidate peptides within the specified mass
-- tolerance that should be examined by spectral cross-correlation.
--
-- This generates a device array containing the indices of the relevant
-- sequences, together with the number of candidates found.
--
findCandidates :: ConfigParams -> DeviceSeqDB -> MS2Data -> (Candidates -> IO b) -> IO b
findCandidates cp db ms2 action =
  CUDA.allocaArray np $ \d_idx ->
  CUDA.findIndicesInRange (dbResidual db) d_idx np (mass-delta) (mass+delta) >>= action . (d_idx,)
  where
    np    = dbNumFrag db
    delta = massTolerance cp
    mass  = (ms2precursor ms2 * ms2charge ms2) - ((ms2charge ms2 * massH) - 1) - (massH + massH2O)


--
-- Generate a theoretical spectral representation for each of the specified
-- candidates. This generates spectral peaks for all of the A-, B- and Y-ions of
-- the given sequences, retaining only the most intense peak in each bin.
--
-- On the device, this is stored as a dense matrix, each row corresponding to a
-- single sequence.
--
mkSpecXCorr :: DeviceSeqDB -> Candidates -> Float -> Int -> (IonSeries -> IO b) -> IO b
mkSpecXCorr db (d_idx, nIdx) chrg len action =
  CUDA.allocaArray n $ \d_spec -> do
    let bytes = fromIntegral $ n * CUDA.sizeOfPtr d_spec

    CUDA.memset  d_spec bytes 0
    CUDA.addIons d_spec (dbResidual db) (dbIonMass db) (dbTerminal db) d_idx nIdx ch len

    action d_spec
  where
    n  = len * nIdx
    ch = round $ max 1 (chrg - 1)


--
-- Score each candidate sequence against the observed intensity spectra,
-- returning the most relevant results.
--
sequestXC :: ConfigParams -> Candidates -> Spectrum -> IonSeries -> IO [(Float,Int)]
sequestXC cp (d_idx,nIdx) expr d_thry = let n = max (numMatches cp) (numMatchesDetail cp) in
  CUDA.withVector  expr $ \d_expr  ->
  CUDA.allocaArray nIdx $ \d_score -> do
    when (verbose cp) $ hPutStrLn stderr ("Matched peptides: " ++ show nIdx)

    -- Score and rank each candidate sequence
    --
    CUDA.mvm       d_score d_thry d_expr nIdx (G.length expr)
    CUDA.radixsort d_score d_idx nIdx

    -- Retrieve the most relevant matches
    --
    sc <- CUDA.peekListArray n (d_score `CUDA.advanceDevPtr` (nIdx-n))
    ix <- CUDA.peekListArray n (d_idx   `CUDA.advanceDevPtr` (nIdx-n))

    return . reverse $ zipWith (\s i -> (s/10000,fromIntegral i)) sc ix

