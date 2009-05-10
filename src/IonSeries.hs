{-
 - Functions to work with a series of ions, and generate theoretical spectra
 - suitable for matching against experimental results with the sequest
 - cross-correlation algorithm.
 -}

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
-- XXX: Apparently there are different types of ions other than the B- Y- and
-- A-series we consider. Also, assumes a charge state of one.
--
buildThrySpecXCorr :: ConfigParams -> Peptide -> XCorrSpecThry
buildThrySpecXCorr cp peptide =
    concatMap (addIonsAB 1.0) b_ions ++ concatMap (addIonsY 1.0) y_ions
    where
        b_ions = init $ ladder peptide
        y_ions = map (\x -> residual peptide - x) b_ions


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
