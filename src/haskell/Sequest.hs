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

module Sequest
  (
    Match(..),
    MatchCollection,
    findCandidates,
    searchForMatches
  )
  where

import Mass
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
-- Data Structures
--------------------------------------------------------------------------------

type MatchCollection a = [Match a]

--
-- A structure to store the result of a peptide/spectrum similarity test
--
data Match a = Match
    {
        candidate :: Peptide a,         -- The fragment that was examined
        scoreXC   :: a                  -- Sequest cross-correlation score
--        scoreSP   :: (Int, Int)         -- Matched ions / total ions
    }
    deriving (Eq, Show)


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

