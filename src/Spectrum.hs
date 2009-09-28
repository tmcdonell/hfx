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
    XCorrSpecExp(..),

    buildExpSpecXCorr
  )
  where

import Config
import Utils

import Data.List
import Data.Array.Unboxed
import Data.ByteString.Lazy (ByteString)

import C2HS
import qualified Foreign.CUDA as G


--------------------------------------------------------------------------------
-- Data Structures
--------------------------------------------------------------------------------

--
-- A group of multiple spectra, read from the same experimental results file.
--
data SpectrumCollection = SpectrumCollection
    {
        spectra :: [Spectrum],          -- Collection of measurements
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
type Peak = (Float, Float)

data Spectrum = Spectrum
    {
        precursor :: Float,             -- The precursor mass; (M+zH)/z
        charge    :: Float,             -- Peptide charge state
        peaks     :: [Peak]             -- The actual mass/charge ratio intensity measurements
    }
    deriving (Eq, Show)

--
-- Find the m/z range of the recorded peaks
--
mzRange      :: Spectrum -> (Float, Float)
mzRange spec =  minmax (peaks spec)
    where
        minmax []       = error "Spectrum.mzRange: empty list"
        minmax (x:xs)   = foldl cmp x xs
        cmp (a,_) (b,_) = (min a b, max a b)

--
-- A mz/intensity spectrum array of experimental data suitable for sequest
-- cross-correlation ranking
--
data XCorrSpecExp = XCorrSpecExp
        (Int,Int)                       -- bounds of the array
        (G.DevicePtr Float)             -- array data, stored on the device


--------------------------------------------------------------------------------
-- Experimental Intensity Spectrum
--------------------------------------------------------------------------------

--
-- Process the observed spectral peaks and generate an array suitable for
-- Sequest cross correlation analysis
--
-- TODO: This goes via a temporary IArray, but it might have been possible to
--       use a StorableArray and copy directly from that memory. Otherwise,
--       investigate the use of an STUArray for fast construction.
--
buildExpSpecXCorr :: ConfigParams
                  -> Spectrum
                  -> IO XCorrSpecExp
buildExpSpecXCorr cp spec =
    let sp = calculateXCorr . normaliseByRegion . observedIntensity cp $ spec
    in do
      d_ptr <- G.newArray . elems $ sp
      return $ XCorrSpecExp (bounds sp) d_ptr


--
-- Generate an intensity array for the observed spectrum. The square root of the
-- input intensity is taken, and only the most intense sample in a given bin is
-- recorded.
--
-- The bin width index for a given mass/charge ratio is stolen from crux...(?)
--
-- Some slight-of-hand going on here. The boundaries of the array are not
-- entirely clear, and determined by several limits:
--   1. An explicit cutoff at (50 + maximum recorded M/Z)
--   2. An implicit cutoff due to the way normaliseByRegion must calculate the
--      window ranges (a "feature" of crux)
--   3. Since this must mimic a FFT in real space, we need to include space for
--      the "wings" in the (-75,+75) region in calculateXCorr
--
observedIntensity :: ConfigParams -> Spectrum -> Array Int Float
observedIntensity cp spec =
    accumArray max 0 bnds [(bin x,sqrt y) | (x,y) <- filter limits (peaks spec)]
    where
        bin mz = round (mz / width)
        mass   = precursor spec
        (m,n)  = mzRange spec
        bnds   = (max 0 ((floor m)-75), 75 + ceiling cutoff)

        cutoff = let mzcut = 50 + mass * charge spec
                     mzlim = fromInteger (truncate (n/10) * 10)
                 in  min mzcut mzlim

        width  = if aaMassTypeMono cp then 1.0005079 else 1.0011413

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
normaliseByRegion :: Array Int Float -> Array Int Float
normaliseByRegion a = array (bounds a) [ (i,norm i) | i <- indices a ]
    where
        rgn_max :: Array Int Float
        rgn_max =  accumArray max 0 (0,10) [(rgn i,e) | (i,e) <- assocs a]

        norm i  = let m = rgn_max ! (rgn i)  in
                  if  m > 1E-6 && i < sel*10 then 50 * ((a!i) / m) else 0

        rgn i   = i `div` sel
        sel     = cutoff `div` 10
        cutoff  = (snd (bounds a)) - 75

--
-- Calculate the sequest cross-correlation function for the given input array.
-- Each sequest matching score is then a dot product between a theoretical input
-- and this pre-processed spectrum.
--
calculateXCorr :: Array Int Float -> Array Int Float
calculateXCorr a  = listArray (0, ceilPow2 n - 1) (repeat 0) // [(i,xcorr i e) | (i,e) <- assocs a]
    where
        (m,n)     = bounds a
        xcorr i e = e - (subFoldA1' (+) a (xrange i)) / 150
        xrange i  = [(max m (i-75)) .. (min n (i+75))]

