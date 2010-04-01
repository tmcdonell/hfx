--------------------------------------------------------------------------------
-- |
-- Module    : Spectrum
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- Meta-module reexporting spectrum related stuff
--
--------------------------------------------------------------------------------

module Spectrum
  (
    -- Data Structures
    MS2Data(..),
    Spectrum,

    -- File Formats
    readMS2Data,
    readDTA,
    readMGF,

    -- Spectrum analysis
    module Spectrum.Correlation
  )
  where

import Spectrum.Data
import Spectrum.DTA
import Spectrum.MGF
import Spectrum.Correlation

import Data.Char (toLower)


--
-- Read MS/MS sample data, based on file extension
--
readMS2Data :: FilePath -> IO (Either String [MS2Data])
readMS2Data fp =
  case suffix fp of
    "dta" -> either Left (Right . unit) `fmap` readDTA fp
    "mgf" -> Right `fmap` readMGF fp
    _     -> return $ Left ("Unsupported file type: " ++ show fp)

  where
    unit x = [x]
    suffix = map toLower . reverse . takeWhile (/= '.') . reverse

