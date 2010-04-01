--------------------------------------------------------------------------------
-- |
-- Module    : Spectrum.Correlation
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- SEQUEST cross-correlation analysis of mass-spectroscopy spectra
--
--------------------------------------------------------------------------------

module Spectrum.Correlation (sequestXCorr) where

import Config
import Spectrum.Data

import Control.Arrow
import qualified Data.Vector.Generic as G


--------------------------------------------------------------------------------
-- Experimental Intensity Spectrum
--------------------------------------------------------------------------------

--
-- Width of the mass/charge ratio bins used to discretize the spectra. These are
-- stolen from Crux.
--
binWidth :: ConfigParams -> Float
binWidth cp = if aaMassTypeMono cp then 1.0005079 else 1.0011413


--
-- Process the observed spectral peaks and generate an array suitable for
-- SEQUEST cross-correlation analysis
--
sequestXCorr :: ConfigParams -> MS2Data -> Spectrum
sequestXCorr cp = crossCorrolation . normaliseByRegion . observedIntensity cp


--
-- Generate an intensity array for the observed spectrum. The square root of the
-- input intensity is taken, and only the most intense sample in a given bin is
-- recorded.
--
-- Some slight-of-hand going on here. The boundaries of the array are not
-- entirely clear, and determined by several limits:
--
--   1. An explicit cutoff at (50 + maximum recorded M/Z)
--
--   2. An implicit cutoff due to the way normaliseByRegion must calculate the
--      window ranges (a "feature" of crux)
--
--   3. Since this must mimic a FFT in real space, we need to include space for
--      the "wings" in the (-75,+75) region in calculateXCorr
--
{-# INLINE observedIntensity #-}
observedIntensity :: ConfigParams -> MS2Data -> Spectrum
observedIntensity cp ms2 =
  G.accumulate max zeros . G.map (bin *** sqrt) $ G.filter limits (ms2data ms2)
  where
    zeros  = G.replicate (75 + ceiling cutoff) 0
    pcr    = ms2precursor ms2
    bin mz = round (mz / binWidth cp)
    cutoff = let mzcut = 50 + pcr * ms2charge ms2
                 mzlim = fromInteger (truncate (fst (G.last (ms2data ms2))/10) * 10)
             in  min mzcut mzlim

    limits (x,_) = if removePrecursorPeak cp
                     then x <= cutoff && (x <= (pcr-15) || (pcr+15) <= x)
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
{-# INLINE normaliseByRegion #-}
normaliseByRegion :: Spectrum -> Spectrum
normaliseByRegion s = G.zipWith norm ix s
  where
    zeros    = G.replicate 12 0
    rgn_max  = G.accumulate_ max zeros (G.map rgn ix) s

    rgn i    = i `div` sel
    sel      = cutoff `div` 10
    cutoff   = G.length s - 75

    ix       = G.enumFromN 0 (G.length s)
    norm i e = let m = rgn_max G.! rgn i  in
               if  m > 1E-6 && i < sel*10 then 50 * (e / m) else 0


--
-- Calculate the SEQUEST cross-correlation function for the given input array.
-- Each matching score is then a dot product between a theoretical input and
-- this pre-processed spectrum.
--
{-# INLINE crossCorrolation #-}
crossCorrolation :: Spectrum -> Spectrum
crossCorrolation a = G.imap xcorr a
  where
    xcorr i e = e - (G.sum (slice' i a)) / 150
    slice' i  = G.take 151 . G.drop (i-75)

