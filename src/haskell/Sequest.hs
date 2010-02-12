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
import Data.Vector.Storable (Storable)
import qualified Data.Vector.Storable as V

{-
import GHC.Conc                                 (numCapabilities)
import Control.Parallel
import Control.Parallel.Strategies
-}

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
-- Parallel MapReduce
--------------------------------------------------------------------------------
{-
mapReduce
  :: Strategy b -> (a -> b)     -- evaluation strategy for mapping function
  -> Strategy c -> ([b] -> c)   -- evaluation strategy for reduction
  -> [a]                        -- list to map/reduce over
  -> c
mapReduce mapStrat mapFunc reduceStrat reduceFunc input =
  mapResult `pseq` reduceResult
  where
    mapResult    = parMap mapStrat mapFunc input
    reduceResult = reduceFunc mapResult `using` reduceStrat
-}

--------------------------------------------------------------------------------
-- Database search
--------------------------------------------------------------------------------

--
-- Search the database for amino acid sequences within a defined mass tolerance.
-- Only peptides which fall within this range will be considered.
--
searchForMatches :: (RealFrac a, Floating a, Enum a, Storable a)
                 => ConfigParams a -> [Protein a] -> Spectrum a -> MatchCollection a
searchForMatches cp database spec
    = finish
--    . mapReduce rwhnf score rwhnf (foldl' record nomatch)
--    . concatMap fragments
--    $ candidates database

    . foldl' record nomatch
    $ [ score peptide | protein <- candidates database, peptide <- fragments protein ]
    where
        specExp    = buildExpSpecXCorr  cp spec
        specThry   = buildThrySpecXCorr cp (charge spec) (0, V.length specExp - 1)
        candidates = findCandidates cp spec . map (digestProtein cp)
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
findCandidates :: (Fractional a, Ord a) => ConfigParams a -> Spectrum a -> [Protein a] -> [Protein a]
findCandidates cp spec
    = filter (not . null . fragments)
    . map (\p -> p {fragments = filter inrange (fragments p)})
--    . parBuffer numCapabilities rwhnf
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
sequestXC :: (Fractional a, Storable a) => XCorrSpecExp a -> XCorrSpecThry Int a -> a
sequestXC v sv = sum [ x * v V.! i | (i,x) <- sv ] / 10000
