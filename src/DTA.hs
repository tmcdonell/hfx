{-
 - Parse a DTA LC-MS/MS results file.
 -
 - The file format is very simple. The first line contains the singly protonated
 - peptide mass (MH+) and the peptide charge state as a pair of space separated
 - values. Subsequent lines contain space separated pairs of fragment ion m/z
 - ratio and intensity values. Note that the precursor peptide mass is
 - independent of the charge state.
 -
 - The filename usually used to identify the dataset, and each file contains
 - only a single MS/MS sample set.
 -}

module MS2.DTA where

import Numeric
import ApplicativeParsec


--------------------------------------------------------------------------------
-- Data Structures
--------------------------------------------------------------------------------

--
-- The mass spectroscopy data generated from a protein identification
-- experiment. It consists of:
--
data Spectrum = Spec
        Double                  -- The singly protonated peptide mass
        Int                     -- Peptide charge state
        [(Double, Double)]      -- The actual mass/charge ratio intensity measurements
    deriving (Eq, Show)


--------------------------------------------------------------------------------
-- DTA File Parser/Lexer
--------------------------------------------------------------------------------

--
-- The DTA file contains at least one line, each of which is terminated by an
-- end-of-line character (eol)
--
dtaFile = endBy line eol

-- 
-- Each line contains exactly two data values, separated by white space. These
-- are returned as a pair of (mass/charge ratio, intensity) values
--
line = liftA2 (,) fval fval
    where fval = fst . head . readSigned readFloat <$> value

--
-- Each value is a (possibly signed) floating point number. Discard any leading
-- white space encountered
--
value = skipMany (oneOf " \t") >> getValue
    where getValue =  many1 (oneOf (['0'..'9']++"-."))
                  <?> "floating-point number"

--
-- The end of line character. Different operating systems use different
-- characters to mark the end-of-line, so just look for all combinations
--
eol =  try (string "\n\r")
   <|> try (string "\r\n")
   <|> string "\r"
   <|> string "\n"
   <?> "end of line"

--
-- Parse the input file
--
parseDTA       :: [Char] -> Either ParseError [(Double, Double)]
parseDTA input =  parse dtaFile "(unknown)" input


--------------------------------------------------------------------------------
-- Spectrum
--------------------------------------------------------------------------------



--------------------------------------------------------------------------------
-- File I/O
--------------------------------------------------------------------------------

readDTA   :: FilePath -> Spectrum
readDTA _ =  error "not implemented yet"

{-
    doReadFile ::  FilePath -> (Handle -> IO c) -> IO c
    doReadFile filename readHandler=
        bracket (openFile filename ReadMode) hClose readHandler

    --
    -- Read the spectrum data from a MS experiment results file.
    --
    --readDTA :: FilePath -> IO Spectrum
    readDTA filename = do
        hdl <- openFile filename ReadMode
        top <- hGetLine hdl
        bot <- hGetContents hdl

        case parse dtaFile "" bot of
            Left  e -> error ("Error parsing input: " ++ show e)
    --        Right s -> liftM2 (,) (many1 val) (many1 val)
            Right s -> mapM_ print s --return Spec mass charge spectrum
-}
