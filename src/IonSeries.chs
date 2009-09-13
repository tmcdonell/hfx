{-# LANGUAGE CPP, ForeignFunctionInterface #-}
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
    XCorrSpecThry(..),
    buildThrySpecXCorr
  )
  where

#include "kernels/kernels.h"

import Config
import Protein

import C2HS
import Foreign.CUDA (DevicePtr, withDevicePtr)
import qualified Foreign.CUDA as CUDA


--------------------------------------------------------------------------------
-- Data structures
--------------------------------------------------------------------------------

--
-- The mz/intensity spectrum array for the theoretical spectrum. Keep this as a
-- sparse array, as there are so few elements compared to the experimental data.
--
-- XXX: Changed to a dense array on the device, to facilitate eventual dot
-- product operation (and because my cuda-fu is weak...)
--
data XCorrSpecThry = XCorrSpecThry
        (Int,Int)
        (CUDA.DevicePtr CInt)


--------------------------------------------------------------------------------
-- Theoretical Spectrum
--------------------------------------------------------------------------------

--
-- Generate the theoretical spectral representation of a peptide from its
-- character code sequence, and do something useful with it. The device memory
-- is deallocated once the action completes.
--
buildThrySpecXCorr :: ConfigParams
                   -> (Int,Int)                 -- ^ bounds of the output array
                   -> Int                       -- ^ precursor charge state
                   -> Peptide                   -- ^ peptide to build spectrum for
                   -> (XCorrSpecThry -> IO b)   -- ^ action to perform
                   -> IO b
buildThrySpecXCorr _cp (m,n) chrg pep action =
    CUDA.allocaBytesMemset bytes 0     $ \spec     ->
    CUDA.withArrayLen (map cFloatConv . bIonLadder $ pep) $ \l b_ions ->
    CUDA.withArray    (map cFloatConv . yIonLadder $ pep) $ \y_ions   ->
    addIons chrg b_ions y_ions spec l len >>

    action (XCorrSpecThry (m,n) spec)

    where
      len   = n - m + 1
      bytes = fromIntegral len * fromIntegral (sizeOf (undefined::Int))


{# fun unsafe addIons
    { cIntConv          `Int'              ,
      withDevicePtr*    `DevicePtr CFloat' ,
      withDevicePtr*    `DevicePtr CFloat' ,
      withDevicePtr*    `DevicePtr CInt'   ,
                        `Int'              ,
                        `Int'              } -> `()' #}

#if 0
buildThrySpecXCorr :: ConfigParams -> Int -> Int -> Peptide -> IO XCorrSpecThry
buildThrySpecXCorr _cp len_spec charge peptide = do
    spec <- CUDA.forceEither `fmap` CUDA.malloc bytes
    rv   <- CUDA.memset spec bytes 0

    case rv of
      Just e  -> error e
      Nothing -> CUDA.withArrayLen (bIonLadder peptide) $ \len_ions b_ions -> do
                 CUDA.withArray    (yIonLadder peptide) $ \y_ions -> do
                   addIons charge b_ions y_ions spec len_ions len_spec

    return $ XCorrSpecThry (0,len_spec) spec

    where
      bytes = fromIntegral (len_spec * sizeOf (undefined::Int))
#endif
#if 0
--
-- Sequential version
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
#endif

