{-
 - An implementation of the SEQUEST algorithm for fast cross-correlation based
 - identification of protein sequences.
 -
 -
 - References:
 -
 - [1] J. K. Eng, B. Fischer, J. Grossmann, and M. J. MacCoss. "A fast sequest
 -     cross correlation algorithm." Journal of Proteome Research,
 -     7(10):4598-4602, 2008.
 -
 - [2] J. K. Eng, A. L. McCormack, and I. John R. Yates. "An approach to
 -     correlate tandem mass spectral data of peptides with amino acid sequences
 -     in a protein database." Journal of the American Society for Mass
 -     Spectrometry, 5(11):976-989, November 1994.
 -}

module Sequest where

import Config
import Protein
import Spectrum
import IonSeries
import AminoAcid

import Data.List
import Data.Maybe
import Data.Array.Unboxed

--------------------------------------------------------------------------------
-- Data Structures
--------------------------------------------------------------------------------

type MatchCollection = [Match]

--
-- A structure to store the result of a peptide/spectrum similarity test
--
data Match = Match
    {
--        scoreSP   :: (Int, Int),        -- Matched ions / total ions
        scoreXC   :: Float,             -- Sequest cross-correlation score
        candidate :: Peptide            -- The fragment that was examined
    }
    deriving (Eq, Show)


--------------------------------------------------------------------------------
-- Database search
--------------------------------------------------------------------------------

--
-- Search the database for amino acid sequences within a defined mass tolerance.
-- Only peptides which fall within this range will be considered.
--
searchForMatches :: ConfigParams -> ProteinDatabase -> Spectrum -> MatchCollection
searchForMatches cp database spec = finish $
    foldl' record nomatch [ sequest cp experimental peptide cutoff cr |
                              protein <- candidates database,
                              peptide <- fragments protein
                          ]
    where
        experimental = buildExpSpecXCorr cp spec
        candidates   = findCandidates cp spec . map (digestProtein cp)
        finish       = reverse . catMaybes

        record l     = tail . flip (insertBy cmp) l . Just
        n            = max (numMatches cp) (numMatchesDetail cp)
        nomatch      = replicate n Nothing

        mass         = precursor spec
        cr           = charge spec
        cutoff       = 50 + mass * cr

        cmp (Just x) (Just y) = compare (scoreXC x) (scoreXC y)
        cmp _        _        = GT


--
-- Search the protein database for candidate peptides within the specified mass
-- tolerance that should be examined by spectral cross-correlation.
--
findCandidates :: ConfigParams -> Spectrum -> ProteinDatabase -> ProteinDatabase
findCandidates cp spec =
    filter (not.null.fragments) . map (\p -> p {fragments = filter inrange (fragments p)})
    where
        inrange p = (mass - limit) <= pmass p && pmass p <= (mass + limit)
        mass      = (precursor spec * charge spec) - ((charge spec * massH) - 1)
        limit     = massTolerance cp


--------------------------------------------------------------------------------
-- Scoring
--------------------------------------------------------------------------------

sequest :: ConfigParams -> XCorrSpecExp -> Peptide -> Float -> Float -> Match
sequest cp spec pep cutoff cr = Match
    {
        scoreXC   = (sequestXC cp spec (buildThrySpecXCorr cp pep cutoff cr)) / 10000,
        candidate = pep
    }


--
-- Explicitly de-forested array dot product [P. Walder, Deforestation, 1988]
--
dot     :: (Ix a, Num a, Num e) => Array a e -> Array a e -> e
dot v w =  loop m 0
    where
        (m,n)                  = bounds v
        loop i acc | i > n     = acc
                   | otherwise = loop (i+1) (v!i * w!i + acc)

--
-- Score a peptide against the observed intensity spectrum. The sequest cross
-- correlation is the dot product between the theoretical representation and the
-- preprocessed experimental spectra.
--
sequestXC :: ConfigParams -> XCorrSpecExp -> XCorrSpecThry -> Float
sequestXC cp v sv = dot v w
    where
        w      = accumArray max 0 (bounds v) [(bin i,e) | (i,e) <- sv, inRange (bounds v) (bin i)]
        bin mz = round (mz / width)
        width  = if aaMassTypeMono cp then 1.0005079 else 1.0011413

