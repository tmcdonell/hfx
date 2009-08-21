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

module IonSeries where

import Config
import Protein
import AminoAcid


--------------------------------------------------------------------------------
-- Data structures
--------------------------------------------------------------------------------

--
-- The mz/intensity spectrum array for the theoretical spectrum. Keep this as a
-- sparse array, as there are so few elements compared to the experimental data.
--
type XCorrSpecThry = [(Float, Float)]


--------------------------------------------------------------------------------
-- Theoretical Spectrum
--------------------------------------------------------------------------------

--
-- Generate the theoretical spectral representation of a peptide from its
-- character code sequence.
--
buildThrySpecXCorr :: ConfigParams -> Float -> Peptide -> XCorrSpecThry
buildThrySpecXCorr _cp charge peptide =
    concatMap addIons $ [1 .. (max 1 (charge-1))]
    where
        addIons c = concatMap (addIonsAB c) b_ions ++ concatMap (addIonsY c) y_ions
        b_ions    = bIonLadder peptide
        y_ions    = yIonLadder peptide


--
-- Convert mass to mass/charge ratio
--
ionMZ :: Float -> Float -> Float
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
addIonsAB, addIonsY :: Float -> Float -> [(Float, Float)]
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

