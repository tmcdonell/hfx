--------------------------------------------------------------------------------
-- |
-- Module    : Spectrum.MGF
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- Parse a MGF MS/MS results file, which may contain many sample sets. Not all
-- features of this format are implemented. Specifically, only the parameters
-- TITLE, CHARGE and PEPMASS are registered, the latter ignoring the (optional)
-- intensity component. Default / global values are not allowed.
--
-- <http://www.matrixscience.com/help/data_file_help.html#GEN>
--
--------------------------------------------------------------------------------

module Spectrum.MGF (readMGF) where

import Util.Parsec
import Spectrum.Data

import Numeric
import Bio.Util.Parsex
import Control.Monad
import Control.Applicative                      hiding (many)
import Text.ParserCombinators.Parsec

import qualified Data.Vector.Unboxed            as U
import qualified Data.ByteString.Lazy.Char8     as L

instance Applicative (GenParser tok st) where
    (<*>) = ap
    pure  = return


--
-- Read the given file and return either an error or a list of the contained
-- MS/MS spectral samples.
--
-- The file is read lazily, so any errors in the file will cause an exception to
-- be thrown once they are encountered.
--
readMGF :: FilePath -> IO [MS2Data]
readMGF name = lazyMany sample name <$> readFile name


--
-- Each sample represents a complete MS/MS sample. Local parameters must be
-- specified before the ions list. We require that the peptide mass and
-- precursor charge state be provided, and optionally, a sample title. All other
-- parameters are currently ignored.
--
-- XXX: read global configuration options, optional charge state
--
sample :: Parser MS2Data
sample = do
  kv <-         string "BEGIN IONS" *> eol *> params
  pk <- ions <* string "END IONS"   <* eol

  let title     = maybe L.empty L.pack (lookup "TITLE" kv)
      charge    = maybe 0 (fst . head . readSigned readFloat) (lookup "CHARGE"  kv)
      pepmass   = maybe 0 (fst . head . readSigned readFloat) (lookup "PEPMASS" kv)

  if charge == 0 || pepmass == 0
    then unexpected "missing charge state or peptide mass"
    else return $ MS2Data title pepmass charge pk

--
-- Each sample block may contain a sequence of key/value pairs.
--
params :: Parser [(String, String)]
params = many (try keyval)

--
-- Each line contains a pair of whitespace separated floating point values,
-- representing the (mass/charge ratio, intensity) measurements.
--
ions :: Parser (U.Vector Peak)
ions = liftM U.fromList (endBy readF2 eol)

