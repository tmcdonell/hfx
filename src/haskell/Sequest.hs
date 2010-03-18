{-# LANGUAGE TupleSections #-}
--------------------------------------------------------------------------------
-- |
-- Module    : Sequest
-- Copyright : (c) 2009 Trevor L. McDonell
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

module Sequest where

import Mass
import Match
import Config
import Spectrum
import Sequence

import Data.Word
import Data.Maybe
import Control.Monad
import Prelude                                  hiding (lookup)

import Foreign.CUDA (DevicePtr)
import qualified Data.Vector.Generic            as G
import qualified Foreign.CUDA                   as CUDA
import qualified Foreign.CUDA.Utils             as CUDA
import qualified Foreign.CUDA.Algorithms        as CUDA


type Candidates      = (DevicePtr Word32, Int)
type XCorrSpecThry   = DevicePtr Word32

--------------------------------------------------------------------------------
-- Search
--------------------------------------------------------------------------------

--
-- Search an amino acid sequence database to find the most likely matches to a
-- given experimental spectrum.
--
searchForMatches :: ConfigParams -> SequenceDB -> DeviceSeqDB -> Spectrum -> IO MatchCollection
searchForMatches cp sdb ddb dta =
  findCandidates cp ddb dta                                  $ \candidates ->
  mkSpecXCorr ddb candidates (charge dta) (G.length specExp) $ \specThry   ->
  mapMaybe lookupF `fmap` sequestXC cp candidates specExp specThry
  where
    specExp       = buildExpSpecXCorr cp dta
    lookupF (s,i) = liftM (flip Match s) (lookup sdb i)


--
-- Search the protein database for candidate peptides within the specified mass
-- tolerance that should be examined by spectral cross-correlation.
--
-- This generates a device array containing the indices of the relevant
-- sequences, together with the number of candidates found.
--
findCandidates :: ConfigParams -> DeviceSeqDB -> Spectrum -> (Candidates -> IO b) -> IO b
findCandidates cp db dta action =
  CUDA.allocaArray np $ \d_idx ->
  CUDA.findIndicesInRange (dbResidual db) d_idx np (mass-delta) (mass+delta) >>= action . (d_idx,)
  where
    np    = dbNumFrag db
    delta = massTolerance cp
    mass  = (precursor dta * charge dta) - ((charge dta * massH) - 1) - (massH + massH2O)


--
-- Generate a theoretical spectral representation for each of the specified
-- candidates. This generates spectral peaks for all of the A-, B- and Y-ions of
-- the given sequences, retaining only the most intense peak in each bin.
--
-- On the device, this is stored as a dense matrix, each row corresponding to a
-- single sequence.
--
mkSpecXCorr :: DeviceSeqDB -> Candidates -> Float -> Int -> (XCorrSpecThry -> IO b) -> IO b
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
sequestXC :: ConfigParams -> Candidates -> XCorrSpecExp -> XCorrSpecThry -> IO [(Float,Int)]
sequestXC cp (d_idx,nIdx) expr d_thry = let n = max (numMatches cp) (numMatchesDetail cp) in
  CUDA.withVectorS expr $ \d_expr  ->
  CUDA.allocaArray nIdx $ \d_score -> do

    -- Score and rank each candidate sequence
    --
    CUDA.mvm       d_score d_thry d_expr nIdx (G.length expr)
    CUDA.radixsort d_score d_idx nIdx

    -- Retrieve the most relevant matches
    --
    sc <- CUDA.peekListArray n (d_score `CUDA.advanceDevPtr` (nIdx-n))
    ix <- CUDA.peekListArray n (d_idx   `CUDA.advanceDevPtr` (nIdx-n))

    return . reverse $ zipWith (\s i -> (s/10000,fromIntegral i)) sc ix


{-
module Sequest
  (
    Match(..),
    MatchCollection,
    findCandidates,
    searchForMatches
  )
  where

import Mass
import Match
import Config
import Protein
import Spectrum
import IonSeries

import Data.List
import Data.Maybe
import Data.Ord
import Data.Vector                      (Vector)
import Data.Vector.Storable             (Storable)
import qualified Data.Vector.Generic    as G


--------------------------------------------------------------------------------
-- Database search
--------------------------------------------------------------------------------

--
-- Search the database for amino acid sequences within a defined mass tolerance.
-- Only peptides which fall within this range will be considered.
--
searchForMatches :: (RealFrac a, Floating a, Enum a, Storable a)
                 => ConfigParams a -> Vector (Protein a) -> Spectrum a -> MatchCollection a
searchForMatches cp database spec
    = finish
    . foldl' record nomatch
    . G.toList
    . G.map score
    $ G.concatMap fragments (candidates database)
    where
        specExp    = buildExpSpecXCorr  cp spec
        specThry   = buildThrySpecXCorr cp (charge spec) (G.length specExp)
        candidates = findCandidates cp spec
        finish     = reverse . catMaybes

        record l x = tail $ insertBy (comparing (fmap scoreXC)) (Just x) l
        n          = max (numMatches cp) (numMatchesDetail cp)
        nomatch    = replicate n Nothing

        score p    = Match {
            scoreXC   = sequestXC specExp (specThry p),
            candidate = p
          }

--
-- Search the protein database for candidate peptides within the specified mass
-- tolerance that should be examined by spectral cross-correlation.
--
findCandidates :: (Fractional a, Ord a) => ConfigParams a -> Spectrum a -> Vector (Protein a) -> Vector (Protein a)
findCandidates cp spec
    = G.filter (not . G.null . fragments)
    . G.map (\p -> p {fragments = G.filter inrange (fragments p)})
    where
        inrange p = (mass - limit) <= pmass p && pmass p <= (mass + limit)
        mass      = (precursor spec * charge spec) - ((charge spec * massH) - 1)
        limit     = massTolerance cp


--------------------------------------------------------------------------------
-- Scoring
--------------------------------------------------------------------------------

--
-- Score a peptide against the observed intensity spectrum. The sequest cross
-- correlation is the dot product between the theoretical representation and the
-- preprocessed experimental spectra.
--
sequestXC :: (Fractional a, Storable a) => XCorrSpecExp a -> XCorrSpecThry a -> a
sequestXC x y = G.sum ( G.zipWith (*) x y ) / 10000
-}

