{-
 - A simple test file for the peptide sequence matching algorithm
 -}

module Main where

--
-- Custom libraries
--
import Config
import DTA
import Sequest
import Utils

--
-- System libraries
--
import Bio.Sequence (readFasta)
import System.Environment (getArgs)

--------------------------------------------------------------------------------
-- Main
--------------------------------------------------------------------------------
help, version :: String
help    = "No help for you!"
version = "2.71"

main :: IO ()
main = do
    [dta]    <- getArgs
    spectrum <- readDTA dta
    config   <- readParams "sequest.params"

    let cp    = forceEither config
        ms    = forceEither spectrum
        fasta = databasePath cp

    database     <- readFasta fasta
    printResults $! findMatch cp ms database


--
-- Pulls a "Right" value out of an Either construct. If the either is a "Left",
-- raises an exception with that string.
--
forceEither :: Either String a -> a
forceEither (Left  e) = error e
forceEither (Right x) = x

