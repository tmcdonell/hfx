--------------------------------------------------------------------------------
-- |
-- Module    : IonSeries
-- Copyright : (c) 2009 Trevor L. McDonell
-- License   : BSD
--
-- Functions to work with a series of ions, and generate theoretical spectra
-- suitable for matching against experimental results with the sequest
-- cross-correlation algorithm.
--
--------------------------------------------------------------------------------

module IonSeries
  (
    XCorrSpecThry,
    buildThrySpecXCorr
  )
  where

import Mass
import Config
import Protein
import Spectrum                         (specBinWidth)

import Data.Ix
import Data.Vector.Storable             (Vector, Storable)
import qualified Data.Vector.Generic    as G


--------------------------------------------------------------------------------
-- Data structures
--------------------------------------------------------------------------------

--
-- The mz/intensity spectrum array for the theoretical spectrum. Although the
-- vector is fairly sparse, store in dense format for simplicity.
--
type XCorrSpecThry a = Vector a


--------------------------------------------------------------------------------
-- Theoretical Spectrum
--------------------------------------------------------------------------------

--
-- Generate the theoretical spectral representation of a peptide from its
-- character code sequence.
--
-- This generates spectral peaks for all of the A-, B- and Y-ions for the given
-- peptide, retaining only the most intense peak in each bin.
--
buildThrySpecXCorr :: (RealFrac a, Enum a, Storable a)
                   => ConfigParams a -> a -> Int -> Peptide a -> XCorrSpecThry a
buildThrySpecXCorr cp charge len peptide =
  G.accum max (G.replicate len 0) ions
  where
    ions      = filter (inRange (0,len-1) . fst)
              . map (\(x,y) -> (bin x, y))
              . concatMap addIons $ [ 1 .. max 1 (charge-1) ]

    addIons c = concatMap (addIonsAB c) b_ions ++ concatMap (addIonsY c) y_ions
    b_ions    = G.toList $ bIonLadder cp peptide
    y_ions    = G.toList $ yIonLadder cp peptide

    bin mz    = round (mz / specBinWidth cp)


--
-- Convert mass to mass/charge ratio
--
ionMZ :: Fractional a => a -> a -> a
ionMZ m c = (m + massH*c) / c


--
-- Add a spectral peak for each fragment location, as well as the peaks
-- corresponding to neutral losses of H2O and NH3.
--
-- The factors that contribute to the collision induced dissociation process
-- used in tandem-MS experiments are not completely understood, so accurate
-- prediction of fragment ion abundances are not possible. Magnitude components
-- are assigned based on empirical knowledge.
--
addIonsAB, addIonsY :: Fractional a => a -> a -> [(a, a)]
addIonsAB charge mass = addIonsA : addIonsB
  where
    addIonsA = let m = ionMZ (mass - massCO) charge in (m, 10)
    addIonsB = let m = ionMZ mass charge in
      [
        (m,50), (m+1,25), (m-1,25),
        (m - massH2O/charge, 10),
        (m - massNH3/charge, 10)
      ]

addIonsY charge mass =
  let m = ionMZ (mass + massH2O) charge in
    [
      (m,50), (m+1,25), (m-1,25),
      (m - massNH3/charge, 10)
    ]

