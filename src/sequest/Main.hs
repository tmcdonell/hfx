{-# LANGUAGE ScopedTypeVariables #-}
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
import Sequest
import PrettyPrint

--
-- System libraries
--
import Data.Vector.Storable (Storable)
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

    mapM_ (search (cp :: ConfigParams Float)) dtaFiles


--
-- Search the protein database for a match to the experimental spectra
--
search :: (RealFloat a, Storable a, Read a, Enum a) => ConfigParams a -> FilePath -> IO ()
search cp fp = do
    dta         <- readDTA fp
    proteins    <- case databasePath cp of
        Nothing -> error "Protein database not specified"
        Just db -> readFasta db

    let spec    = forceEitherStr dta
        matches = searchForMatches cp proteins spec

    printConfig cp fp spec

    printResults       $! (take (numMatches cp)       matches)
    printResultsDetail $! (take (numMatchesDetail cp) matches)

