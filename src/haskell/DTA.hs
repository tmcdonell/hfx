--------------------------------------------------------------------------------
-- |
-- Module    : DTA
-- Copyright : (c) 2009 Trevor L. McDonell
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

module DTA
    (
      Spectrum(..),     -- Data structure
      readDTA           -- File formats
    ) where

import Mass
import Spectrum

import Numeric
import Control.Monad (liftM2)
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
-- Each line contains exactly two data values, separated by white space. These
-- are returned as a pair of (mass/charge ratio, intensity) values. Detecting
-- signed values isn't really necessary, but is done for completeness.
--
line :: RealFrac a => Parser (a,a)
line =  liftM2 (,) fval fval
    where fval = (fst . head . readSigned readFloat) `fmap` value

--
-- Each value is a floating point number. Discard any leading white space
-- encountered.
--
value :: Parser String
value =  skipMany (oneOf " \t") >> getValue
    where getValue =  many1 (digit <|> char '.' <|> char '-')
                  <?> "floating-point number"

--
-- The end of line character. Different operating systems use different
-- characters to mark the end-of-line, so just look for all combinations
--
eol :: Parser String
eol =  try (string "\n\r")
   <|> try (string "\r\n")
   <|> string "\r"
   <|> string "\n"
   <?> "end of line"


--------------------------------------------------------------------------------
-- Spectrum
--------------------------------------------------------------------------------

--
-- Encase the values read from the DTA file into a data structure
--
mkSpec :: RealFrac a => [(a,a)] -> Either String (Spectrum a)
mkSpec []           =  Left "Error: empty spectrum"
mkSpec ((m,c):ss)
    | trunc' c /= c =  Left "Error: invalid peptide charge state\nexpecting integer"
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
readDTA :: RealFrac a => FilePath -> IO (Either String (Spectrum a))
readDTA name =  do
    dta <- parseFromFile dtaFile name
    case dta of
        Left  e -> return (Left ("Error parsing file: " ++ show e))
        Right s -> return (mkSpec s)
