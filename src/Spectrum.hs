{-
 - Data structures and functions to store and manipulate the results of a
 - mass-spectroscopy experiment. This will generate an intensity array suitable
 - for sequest cross-correlation analysis.
 -}

module Spectrum where

import Config
import Utils

import Data.List
import Data.Array.Unboxed
import Data.ByteString.Lazy (ByteString)


--------------------------------------------------------------------------------
-- Data Structures
--------------------------------------------------------------------------------

--
-- A group of multiple spectra, read from the same experimental results file.
--
data SpectrumCollection = SpectrumCollection
    {
        spectra     :: [Spectrum],      -- Collection of measurements
        description :: ByteString       -- Description/header from the input file
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
type Peak = (Float, Float)

data Spectrum = Spectrum
    {
        precursor :: Float,        -- The singly protonated peptide mass
        charge    :: Float,        -- Peptide charge state
        peaks     :: [Peak]        -- The actual mass/charge ratio intensity measurements
    }
    deriving (Eq, Show)

--
-- Find the m/z range of the recorded peaks
--
mzRange      :: Spectrum -> (Float, Float)
mzRange spec =  minmax (peaks spec)
    where
        minmax []       = error "Spectrum.mzRange: empty list"
        minmax (x:xs)   = foldl' cmp x xs
        cmp (a,_) (b,_) = (min a b, max a b)

--
-- A mz/intensity spectrum array of experimental data suitable for sequest
-- cross-correlation ranking
--
type XCorrSpecExp = Array Int Float


--------------------------------------------------------------------------------
-- Experimental Intensity Spectrum
--------------------------------------------------------------------------------

--
-- Process the observed spectral peaks and generate an array suitable for
-- Sequest cross correlation analysis.
--
buildExpSpecXCorr    :: ConfigParams -> Spectrum -> XCorrSpecExp
buildExpSpecXCorr cp =  calculateXCorr . normaliseByRegion . observedIntensity cp


--
-- Generate an intensity array for the observed spectrum. The square root of the
-- input intensity is taken, and only the most intense sample in a given bin is
-- recorded.
--
-- The bin width index for a given mass/charge ratio is stolen from crux...(?)
--
observedIntensity :: ConfigParams -> Spectrum -> Array Int Float
observedIntensity cp spec =
    accumArray max 0 bnds [(bin x,sqrt y) | (x,y) <- filter limits (peaks spec)]
    where
        bin mz = round (mz / width)
        mass   = precursor spec
        cutoff = 50 + mass * (charge spec)
        bnds   = let (m,n) = mzRange spec in (floor m, ceiling n)

        width  = if aaMassTypeMono cp then 1.0005079 else 1.0011413

        limits (x,_) = if removePrecursorPeak cp
                       then x <= cutoff && (x < (mass-5) || (mass+5) < x)
                       else x <= cutoff


--
-- Normalise each element of the input array according to the maximum value in
-- each of 10 equally sized windows. Although the input array stores only the
-- minimum range of values, the algorithm requires that the normalisation
-- windows are from zero.
--
-- The bins are sized such that the range of the last may extend past the end of
-- the input array. This means that fewer values may be placed into the final
-- bin, but all values from the input array are guaranteed to be considered.
--
normaliseByRegion :: Array Int Float -> Array Int Float
normaliseByRegion a = array (bounds a) [(i,norm i) | i <- indices a]
    where
        rgn_max :: Array Int Float
        rgn_max = accumArray max 0 (0,9) [(rgn i,e) | (i,e) <- assocs a]
        norm i  = let m = rgn_max ! (rgn i) in
                  if  m > 1E-6 then 50 * ((a!i) / m) else 0

        rgn i   = i `div` sel
        sel     = (9 + snd (bounds a)) `div` 10


--
-- Calculate the sequest cross-correlation function for the given input array.
-- Each sequest matching score is then a dot product between a theoretical input
-- and this pre-processed spectrum.
--
calculateXCorr :: Array Int Float -> XCorrSpecExp
calculateXCorr a = array (bounds a) [(i,xcorr i e) | (i,e) <- assocs a]
    where
        (m,n)     = bounds a
        xcorr i e = e - (subFoldA1' (+) a (xrange i)) / 150
        xrange i  = [(max m (i-75))..(i-1)] ++ [(i+1)..(min n (i+75))]

