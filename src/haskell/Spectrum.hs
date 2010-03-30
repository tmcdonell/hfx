--------------------------------------------------------------------------------
-- |
-- Module    : Spectrum
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- Data structures and functions to store and manipulate the results of a
-- mass-spectroscopy experiment.
--
--------------------------------------------------------------------------------

module Spectrum
  (
    MS2Data(..),
    Spectrum, SpectrumCollection(..),
  )
  where

import Data.ByteString.Lazy             (ByteString)
import qualified Data.Vector.Unboxed    as U


--
-- A group of multiple spectra, read from the same experimental results file.
--
data SpectrumCollection = SpectrumCollection
  {
    scSpectra :: [Spectrum],            -- Collection of measurements
    scHeader  :: ByteString             -- Description/header for this collection
  }
  deriving (Eq, Show)

--
-- The mass spectroscopy data generated from an experiment
--
-- A mass spectrum consists mainly of a list of detected spectra peaks,
-- consisting of a m/z location and intensity, along with some identifying
-- information. A single spectrum is generated from one or more "scans" of the
-- mass spectrometer.
--
-- In addition to scan information, a tandem fragmentation mass spectrum has
-- information about the m/z of the intact ion that generated the spectrum,
-- known as the "precursor".
--
type Peak     = (Float, Float)
type Spectrum = U.Vector Float

data MS2Data  = MS2Data
  {
    ms2info      :: ByteString,
    ms2precursor :: Float,              -- The precursor mass; (M+zH)/z
    ms2charge    :: Float,              -- Peptide charge state
    ms2data      :: U.Vector Peak       -- The actual mass/charge ratio intensity measurements
  }

