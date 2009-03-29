{-
 - Data structure to store results of a mass spectroscopy experiment
 -}

module Spectrum where

import Data.Array.Unboxed

--------------------------------------------------------------------------------
-- Data Structures
--------------------------------------------------------------------------------

--
-- The actual mass/charge ratio intensity measurements
--
type Spectrum = Array Int Double

--
-- Output a list of binned (mz,intensity) pairs to an immutable array
--
mkSpectrum      :: (Int, Int) -> [(Int, Double)] -> Spectrum
mkSpectrum bnds =  accumArray (+) 0 bnds . filter (inRange bnds . fst)


--
-- The mass spectroscopy data generated from a protein identification
-- experiment. It consists of:
--
data MS2Data = MS2
        Double                  -- The singly protonated peptide mass
        Double                  -- Peptide charge state
        [(Double,Double)]       -- The actual mass/charge ratio intensity measurements
    deriving (Eq, Show)

--
-- Extract readings from the experimental spectrum
--
getPrecursor, getParentMass :: MS2Data -> Double
getData                     :: MS2Data -> [(Double, Double)]

getPrecursor  (MS2 p c _) = (p + c - 1) / c
getParentMass (MS2 p _ _) = p
getData       (MS2 _ _ s) = s

