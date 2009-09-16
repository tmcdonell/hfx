{-# LANGUAGE CPP, ForeignFunctionInterface #-}
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
    searchForMatches
  )
  where

#include <cublas.h>
#include "kernels/kernels.h"

import Mass
import Config
import Protein
import Spectrum
import IonSeries

import Data.List
import Data.Maybe

import C2HS
import Foreign.CUDA (DevicePtr, withDevicePtr)
import qualified Foreign.CUDA as CUDA


--------------------------------------------------------------------------------
-- Data Structures
--------------------------------------------------------------------------------

type MatchCollection = [Match]

--
-- A structure to store the result of a peptide/spectrum similarity test
--
data Match = Match
    {
        candidate :: Peptide,           -- The fragment that was examined
        scoreXC   :: Float              -- Sequest cross-correlation score
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
searchForMatches :: ConfigParams -> ProteinDatabase -> Spectrum -> IO MatchCollection
searchForMatches cp database spec = do
    specExp <- buildExpSpecXCorr cp spec
    scores  <- sequence [ score specExp peptide |
                            protein <- candidates database,
                            peptide <- fragments protein
                        ]

    return . finish $ foldl' record nomatch scores
    where
        specThry b = buildThrySpecXCorr cp b (round (charge spec))
        candidates = findCandidates cp spec . map (digestProtein cp)
        finish     = reverse . catMaybes

        record l   = tail . flip (insertBy cmp) l . Just
        n          = max (numMatches cp) (numMatchesDetail cp)
        nomatch    = replicate n Nothing

        score e p  = Match p `fmap` (specThry (bnds e) p >>= sequestXC cp e)

        bnds (XCorrSpecExp b _) = b
        cmp (Just x) (Just y)   = compare (scoreXC x) (scoreXC y)
        cmp _        _          = GT


#if 0
searchForMatches :: ConfigParams -> ProteinDatabase -> Spectrum -> MatchCollection
searchForMatches cp database spec = finish $
    foldl' record nomatch [ score peptide |
                              protein <- candidates database,
                              peptide <- fragments protein
                          ]
    where
        specExp    = buildExpSpecXCorr  cp spec
        specThry   = buildThrySpecXCorr cp (charge spec)
        candidates = findCandidates cp spec . map (digestProtein cp)
        finish     = reverse . catMaybes

        record l   = tail . flip (insertBy cmp) l . Just
        n          = max (numMatches cp) (numMatchesDetail cp)
        nomatch    = replicate n Nothing

        score p    = Match {
            scoreXC   = sequestXC cp spec specExp (specThry p),
            candidate = p
          }

        cmp (Just x) (Just y) = compare (scoreXC x) (scoreXC y)
        cmp _        _        = GT
#endif


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

--
-- Score a peptide against the observed intensity spectrum. The sequest cross
-- correlation is the dot product between the theoretical representation and the
-- preprocessed experimental spectra.
--
sequestXC :: ConfigParams -> XCorrSpecExp -> XCorrSpecThry -> IO Float
sequestXC _cp (XCorrSpecExp (m,n) d_exp) (XCorrSpecThry _ d_thry) = do
    CUDA.allocaBytes bytes $ \d_xs -> do
    map_fromIntegralf d_thry d_xs len
    cublasSdot len d_exp 1 d_xs 1 >>= \x -> return (x / 10000)
    where
      len   = (n-m)
      bytes = fromIntegral len * fromIntegral (sizeOf (undefined::CFloat))


{# fun unsafe zipWithTimesif
    { withDevicePtr* `DevicePtr CInt'   ,
      withDevicePtr* `DevicePtr CFloat' ,
      withDevicePtr* `DevicePtr CFloat' ,
                     `Int'              } -> `()' #}

{# fun unsafe reducePlusf
    { withDevicePtr* `DevicePtr CFloat' ,
                     `Int'              } -> `Float' #}

{# fun unsafe map_fromIntegralf
    { withDevicePtr* `DevicePtr CInt'   ,
      withDevicePtr* `DevicePtr CFloat' ,
                     `Int'              } -> `()' #}

{# fun unsafe cublasSdot
    {                `Int'              ,
      withDevicePtr* `DevicePtr CFloat' ,
                     `Int'              ,
      withDevicePtr* `DevicePtr CFloat' ,
                     `Int'              } -> `Float' #}


#if 0
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
sequestXC :: ConfigParams -> Spectrum -> XCorrSpecExp -> XCorrSpecThry -> Float
sequestXC cp spec v sv = (dot v w) / 10000
    where
        w      = accumArray max 0 bnds [(bin i,e) | (i,e) <- sv, inRange bnds (bin i)]
        bin mz = round (mz / width)
        width  = if aaMassTypeMono cp then 1.0005079 else 1.0011413

        cutoff = 50 + precursor spec * charge spec
        (m,n)  = bounds v
        bnds   = (m, min n (bin cutoff))
#endif

