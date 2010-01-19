--------------------------------------------------------------------------------
-- |
-- Module    : Main
-- Copyright : (c) 2009 Trevor L. McDonell
-- License   : BSD
--
-- A simple test file for the peptide sequence matching algorithm
--
--------------------------------------------------------------------------------

module Main where

--
-- Custom libraries
--
import DTA
import Utils
import Config
import Protein
import PrettyPrint

import Sequest.Base
import qualified Sequest.CUDA as C

--
-- System libraries
--
import Foreign.C.Types
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

    mapM_ (search cp) dtaFiles


--
-- Search the protein database for a match to the experimental spectra
--
search :: ConfigParams CFloat -> FilePath -> IO ()
search cp fp = do
    dta         <- readDTA fp
    proteins    <- case databasePath cp of
        Nothing -> error "Protein database not specified"
        Just db -> readFasta db

    let spec = forceEitherStr dta
        ref  = searchForMatches cp proteins spec

    printConfig cp fp spec
    matches <- C.searchForMatches cp proteins spec

    printResults       $! (take (numMatches cp)       matches)
    printResultsDetail $! (take (numMatchesDetail cp) matches)

