--------------------------------------------------------------------------------
-- |
-- Module    : Spectrum.Parsec
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- Combinators for parsing spectrum files
--
--------------------------------------------------------------------------------

module Spectrum.Parsec
  where

import Numeric
import Control.Monad (liftM2)
import Text.ParserCombinators.Parsec


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


