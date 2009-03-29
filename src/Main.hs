{-
 - A simple test file for the peptide sequence matching algorithm
 -}

module Main where

--
-- Custom libraries
--
import DTA
import Sequest

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
    [dta,fasta] <- getArgs
    spectrum    <- readDTA dta
    database    <- readFasta fasta

    case spectrum of
        Left  err -> putStrLn err
        Right ms2 -> printResults $ take 5 $
            findMatch (getParentMass ms2) (mkXCorrSpec (getData ms2)) database

