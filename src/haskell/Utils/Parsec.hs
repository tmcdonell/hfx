--------------------------------------------------------------------------------
-- |
-- Module    : Utils.Parsec
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- Combinators for parsing files
--
--------------------------------------------------------------------------------

module Utils.Parsec
  where

import Numeric
import Data.Char
import Control.Monad (liftM2)
import Text.ParserCombinators.Parsec


--
-- Each line contains exactly two data values, separated by white space. These
-- are returned as a pair of (mass/charge ratio, intensity) values. Detecting
-- signed values isn't really necessary, but is done for completeness.
--
readF2 :: RealFrac a => Parser (a,a)
readF2 =  liftM2 (,) fval fval
    where fval = (fst . head . readSigned readFloat) `fmap` float

--
-- Each value is a floating point number. Discard any leading white space
-- encountered.
--
float :: Parser String
float =  skipMany (oneOf " \t") >> getValue
    where getValue =  many1 (digit <|> char '.' <|> char '-')
                  <?> "floating-point number"

--
-- The left hand side of an identifier is a C-like string, which must begin with
-- a letter or underscore, and may subsequently also contain numbers and dashes.
--
identifier :: Parser String
identifier =
    do  c  <- letter <|> char '_'
        cs <- many (letter <|> digit <|> char '_' <|> char '-')
        return (c:cs)
    <?> "identifier"

--
-- The right hand side is any sequence of characters until the end of line or
-- the newline symbol, stripping trailing whitespace
--
value :: Parser String
value =
    do  c  <- noneOf "#\n\r"
        cs <- manyTill anyChar (try comment <|> try eol <|> eof)
        return (c:cs)
    <?> "value"

--
-- Comments begin with the pound sign, and continue to the end of the line
--
comment :: Parser ()
comment =  char '#' >> skipMany (noneOf "\n\r")
       <?> "comment"

--
-- A configuration (identifier) key/value pair, separated by an equals sign with
-- any surrounding whitespace removed.
--
keyval :: Parser (String, String)
keyval = liftM2 (\k v -> (k, rstrip v)) key val
    where
        key    = skipMany space >> identifier
        val    = skipMany space >> char '=' >> skipMany space >> value
        rstrip = reverse . dropWhile isSpace . reverse

--
-- The end of line character. Different operating systems use different
-- characters to mark the end-of-line, so just look for all combinations
--
eol :: Parser ()
eol =  many1 (oneOf "\n\r") >> return ()
   <?> "end of line"

