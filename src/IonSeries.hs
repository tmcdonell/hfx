{-# LANGUAGE CPP #-}
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
    XCorrSpecThry,
    buildThrySpecXCorr
  )
  where


import Config
import Kernels
import Protein

--import Mass
--import Data.Array.Unboxed

import C2HS
import qualified Foreign.CUDA as G


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
type XCorrSpecThry = G.DevicePtr CInt


--data XCorrSpecThry = XCorrSpecThry
--  {
--    smdata   :: G.DevicePtr CFloat,     -- ^ sparse matrix data
--    smidx    :: G.DevicePtr CUInt,      -- ^ corresponding index of each non-zero element
--    smlen    :: Int,                    -- ^ number of non-zero elements
--    smrowlen :: G.DevicePtr CUInt       -- ^ number of non-zero elements in each row
--  }
--  deriving (Show)


--------------------------------------------------------------------------------
-- Theoretical Spectrum
--------------------------------------------------------------------------------

{-

map (foldl1 (\(x,y) (_,z) -> (x, max y z)))
    . groupBy (\x y -> fst x == fst y)
    . sortBy (\x y -> compare (fst x) (fst y))
    . map (\(x,y) -> (bin x,y)) $ spec

 -}


--
-- Generate the theoretical spectral representation of a peptide from its
-- character code sequence
--
buildThrySpecXCorr :: ConfigParams
                   -> (Int,Int)                 -- ^ bounds of the output array
                   -> Int                       -- ^ precursor charge state
                   -> Peptide                   -- ^ peptide to build spectrum for
                   -> IO XCorrSpecThry
buildThrySpecXCorr _cp (p,q) cg pep = do
  spec <- G.forceEither `fmap` G.malloc bytes
  rv   <- G.memset spec bytes 0
  case rv of
    Just e  -> error e
    Nothing -> addIons charge (residual pep) (seqladder (parent pep)) spec len_i len_s (offset pep) >>
               return spec
  where
    charge = max 1 (cg-1)
    len_i  = let (c,n) = terminals pep in fromIntegral (n-c)
    len_s  = q - p + 1
    bytes  = fromIntegral len_s * fromIntegral (sizeOf (undefined::CInt))


#if 0
buildThrySpecXCorr _cp (m,n) chrg pep =
    G.allocaBytesMemset bytes 0     $ \spec     ->
    G.withArrayLen (map cFloatConv . bIonLadder $ pep) $ \l b_ions ->
    G.withArray    (map cFloatConv . yIonLadder $ pep) $ \y_ions   ->
    addIons chrg b_ions y_ions spec l len >>

    return (XCorrSpecThry (m,n) spec)

    where
      len   = n - m + 1
      bytes = fromIntegral len * fromIntegral (sizeOf (undefined::Int))
#endif
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
#endif
#if 0
buildThrySpecXCorr :: ConfigParams
                   -> (Int,Int)                 -- ^ bounds of the output array
                   -> Int                       -- ^ precursor charge state
                   -> Peptide                   -- ^ peptide to build spectrum for
                   -> (XCorrSpecThry -> IO b)   -- ^ action to perform
                   -> IO b
buildThrySpecXCorr cp bnds chrg pep action =
    CUDA.withArray (map cIntConv . elems $ spec) $ \s' -> action (XCorrSpecThry bnds s')
    where
        spec      :: Array Int Int
        spec      = accumArray max 0 bnds [(bin i,e) | (i,e) <- concatMap addIons [1 .. (max 1 (fromIntegral chrg-1))], inRange bnds (bin i)]
        bin mz    = round (mz / width)
        width     = if aaMassTypeMono cp then 1.0005079 else 1.0011413

        addIons c = concatMap (addIonsAB c) b_ions ++ concatMap (addIonsY c) y_ions
        b_ions    = bIonLadder pep
        y_ions    = yIonLadder pep
#endif
#if 0
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
addIonsAB, addIonsY :: Num a => Float -> Float -> [(Float, a)]
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

