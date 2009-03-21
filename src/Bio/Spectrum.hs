{-
 - Data structure to store results of a mass spectroscopy experiment
 -}

module Bio.Spectrum where


--------------------------------------------------------------------------------
-- Data Structures
--------------------------------------------------------------------------------

--
-- The mass spectroscopy data generated from a protein identification
-- experiment. It consists of:
--
data Spectrum = Spec
        Double                  -- The singly protonated peptide mass
        Double                  -- Peptide charge state
        [(Double, Double)]      -- The actual mass/charge ratio intensity measurements
    deriving (Eq, Show)

--
-- Extract readings from the experimental spectrum
--
getPrecursor  (Spec p c _) = (p + c - 1) / c
getParentMass (Spec p _ _) = p
getSpectrum   (Spec _ _ s) = s

