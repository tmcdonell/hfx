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
import Spectrum                                 (specBinWidth)

import Data.List
import Data.Array
import Data.Function


--------------------------------------------------------------------------------
-- Data structures
--------------------------------------------------------------------------------

--
-- The mz/intensity spectrum array for the theoretical spectrum. Keep this as a
-- sparse array, as there are so few elements compared to the experimental data.
--
type XCorrSpecThry = [(Int, Float)]


--------------------------------------------------------------------------------
-- Theoretical Spectrum
--------------------------------------------------------------------------------

--
-- Generate the theoretical spectral representation of a peptide from its
-- character code sequence.
--
-- This generates spectral peaks for all of the A-, B- and Y-ions for the given
-- peptide, retaining only the most intense peak in each bin. Incidentally, the
-- output is also sorted by bin index.
--
buildThrySpecXCorr :: ConfigParams -> Float -> (Int,Int) -> Peptide -> XCorrSpecThry
buildThrySpecXCorr cp charge bnds peptide =
  finish [ (bin x,y) | ions  <- map addIons [1.. max 1 (charge-1)]
                     , (x,y) <- ions ]
  where
    addIons c = concatMap (addIonsAB c) b_ions ++ concatMap (addIonsY c) y_ions
    b_ions    = bIonLadder cp peptide
    y_ions    = yIonLadder cp peptide

    bin mz    = round (mz / specBinWidth cp)
    finish    = filter (inRange bnds . fst) . map (foldl1' max) . groupBy ((==) `on` fst) . sort

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

