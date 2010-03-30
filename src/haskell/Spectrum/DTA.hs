--------------------------------------------------------------------------------
-- |
-- Module    : Spectrum.DTA
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- Parse a DTA LC-MS/MS results file.
--
-- The file format is very simple. The first line contains the singly protonated
-- peptide mass (MH+) and the peptide charge state as a pair of space separated
-- values. Subsequent lines contain space separated pairs of fragment ion m/z
-- ratio and intensity values. Note that the precursor peptide mass is
-- independent of the charge state.
--
-- The filename usually used to identify the dataset, and each file contains
-- only a single MS/MS sample set.
--
--------------------------------------------------------------------------------

module Spectrum.DTA (readDTA) where

import Mass
import Spectrum
import Spectrum.Parsec

import Text.ParserCombinators.Parsec


--------------------------------------------------------------------------------
-- DTA File Parser/Lexer
--------------------------------------------------------------------------------

--
-- The DTA file contains at least one line, each of which is terminated by an
-- end-of-line character (eol)
--
dtaFile :: RealFrac a => Parser [(a,a)]
dtaFile =  endBy line eol

--
-- Encase the values read from the DTA file into a data structure
--
mkSpec :: FilePath -> [(Float,Float)] -> Either String Spectrum
mkSpec name []      =  Left ("Error parsing file: " ++ show name ++ "\nempty spectrum")
mkSpec name ((m,c):ss)
    | trunc' c /= c =  Left ("Error parsing file: " ++ show name ++ "\ninvalid peptide charge state\nexpecting integer")
    | otherwise     =  Right (Spectrum pcr c ss)
    where
        pcr    = (m - 1) / c + massH
        trunc' = fromInteger . truncate


--------------------------------------------------------------------------------
-- File I/O
--------------------------------------------------------------------------------

--
-- Read the given file and return either an error or the MS/MS spectrum data.
--
readDTA :: FilePath -> IO (Either String Spectrum)
readDTA name =  do
    dta <- parseFromFile dtaFile name
    case dta of
        Left  e -> return (Left ("Error parsing file: " ++ show e))
        Right s -> return (mkSpec name s)

