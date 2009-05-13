{-
 - A simple test file for the peptide sequence matching algorithm
 -}

module Main where

--
-- Custom libraries
--
import DTA
import Utils
import Config
import Protein
import Sequest
import PrettyPrint

--
-- System libraries
--
import System.Environment (getArgs)


--------------------------------------------------------------------------------
-- Program Defaults
--------------------------------------------------------------------------------

defaultConfigFile :: FilePath
defaultConfigFile =  "sequest.params"


--------------------------------------------------------------------------------
-- Main
--------------------------------------------------------------------------------

main :: IO ()
main = do
    argv           <- getArgs
    (cp, dtaFiles) <- sequestConfig defaultConfigFile argv

    proteins <- case databasePath cp of
        Nothing -> error "Protein database not specified"
        Just fp -> readFasta fp

    mapM_ (search cp proteins) dtaFiles


--
-- Search the protein database for a match to the experimental spectra
--
search :: ConfigParams -> ProteinDatabase -> FilePath -> IO ()
search cp proteins fp = do
    dta         <- readDTA fp
    let spec     = forceEitherStr dta
        matches  = searchForMatches cp proteins spec

    printConfig cp fp spec

    printResults        $! (take (numMatches cp)       matches)
    printResultsDetail  $! (take (numMatchesDetail cp) matches)

