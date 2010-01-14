--------------------------------------------------------------------------------
-- |
-- Module    : Spectrum
-- Copyright : (c) 2009 Trevor L. McDonell
-- License   : BSD
--
-- Data structures and functions to store and manipulate the results of a
-- mass-spectroscopy experiment. This will generate an intensity array suitable
-- for sequest cross-correlation analysis.
--
--------------------------------------------------------------------------------

module Spectrum
  (
    Spectrum(..),
    SpectrumCollection(..),
    XCorrSpecExp,

    buildExpSpecXCorr,
    specBinWidth
  )
  where

import Config
import Utils

import Data.ByteString.Lazy (ByteString)
import Data.Vector.Storable (Vector, Storable)
import qualified Data.Vector.Storable as V


--------------------------------------------------------------------------------
-- Data Structures
--------------------------------------------------------------------------------

--
-- A group of multiple spectra, read from the same experimental results file.
--
data SpectrumCollection a = SpectrumCollection
    {
        spectra :: [Spectrum a],        -- Collection of measurements
        header  :: ByteString           -- Description/header from the input file
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
type Peak a = (a, a)

data Spectrum a = Spectrum
    {
        precursor :: a,                 -- The precursor mass; (M+zH)/z
        charge    :: a,                 -- Peptide charge state
        peaks     :: [Peak a]           -- The actual mass/charge ratio intensity measurements
    }
    deriving (Eq, Show)

--
-- Find the m/z range of the recorded peaks
--
mzRange      :: Ord a => Spectrum a -> (a, a)
mzRange spec =  minmax (peaks spec)
    where
        minmax []       = error "Spectrum.mzRange: empty list"
        minmax (x:xs)   = foldl cmp x xs
        cmp (a,_) (b,_) = (min a b, max a b)

--
-- A mz/intensity spectrum array of experimental data suitable for sequest
-- cross-correlation ranking
--
type XCorrSpecExp a = Vector a


--------------------------------------------------------------------------------
-- Experimental Intensity Spectrum
--------------------------------------------------------------------------------

--
-- Width of the mass/charge ratio bins used to discretize the spectra. These are
-- stolen from Crux.
--
specBinWidth    :: Fractional a => ConfigParams a -> a
specBinWidth cp =  if aaMassTypeMono cp then 1.0005079 else 1.0011413

--
-- Process the observed spectral peaks and generate an array suitable for
-- Sequest cross correlation analysis
--
buildExpSpecXCorr    :: (Floating a, RealFrac a, Storable a) => ConfigParams a -> Spectrum a -> XCorrSpecExp a
buildExpSpecXCorr cp =  calculateXCorr . normaliseByRegion . observedIntensity cp


--
-- Generate an intensity array for the observed spectrum. The square root of the
-- input intensity is taken, and only the most intense sample in a given bin is
-- recorded.
--
-- Some slight-of-hand going on here. The boundaries of the array are not
-- entirely clear, and determined by several limits:
--   1. An explicit cutoff at (50 + maximum recorded M/Z)
--   2. An implicit cutoff due to the way normaliseByRegion must calculate the
--      window ranges (a "feature" of crux)
--   3. Since this must mimic a FFT in real space, we need to include space for
--      the "wings" in the (-75,+75) region in calculateXCorr
--
observedIntensity :: (Floating a, RealFrac a, Storable a) => ConfigParams a -> Spectrum a -> XCorrSpecExp a
observedIntensity cp spec =
  V.accum max (V.replicate (75 + ceiling cutoff) 0) [(bin x,sqrt y) | (x,y) <- filter limits (peaks spec)]
  where
    bin mz = round (mz / specBinWidth cp)
    mass   = precursor spec
    (_,n)  = mzRange spec

    cutoff = let mzcut = 50 + mass * charge spec
                 mzlim = fromInteger (truncate (n/10) * 10)
             in  min mzcut mzlim

    limits (x,_) = if removePrecursorPeak cp
                   then x <= cutoff && (x <= (mass-15) || (mass+15) <= x)
                   else x <= cutoff

--
-- Normalise each element of the input array according to the maximum value in
-- each of 10 equally sized windows. Although the input array stores only the
-- minimum range of values, the algorithm requires that the normalisation
-- windows are from zero.
--
-- The bins are sized such that the set of bins fits inside the input array.
-- This means that some values from the input will not be considered, and be set
-- to zero.
--
--normaliseByRegion :: UArray Int Float -> UArray Int Float
normaliseByRegion :: (Fractional a, Ord a, Storable a) => XCorrSpecExp a -> XCorrSpecExp a
normaliseByRegion a = V.zipWith norm ix a
  where
    rgn_max  = V.accum max (V.replicate 11 0) [(rgn i, a V.! i) | i <- [0 .. V.length a - 1]]

    rgn i    = i `div` sel
    sel      = cutoff `div` 10
    cutoff   = V.length a - 75

    ix       = V.enumFromTo 0 (V.length a - 1)
    norm i e = let m = rgn_max V.! (rgn i)  in
               if  m > 1E-6 && i < sel*10 then 50 * (e / m) else 0


--
-- Calculate the sequest cross-correlation function for the given input array.
-- Each sequest matching score is then a dot product between a theoretical input
-- and this pre-processed spectrum.
--
calculateXCorr :: (Fractional a, Storable a) => XCorrSpecExp a -> XCorrSpecExp a
calculateXCorr a = V.zipWith xcorr ix a
  where
    ix        = V.enumFromTo 0 (ceilPow2 (V.length a) - 1)
    xcorr i e = e - (V.foldl' (+) 0 (subv i)) / 150
    subv  i   = V.take 151 . V.drop (i-75) $ a

