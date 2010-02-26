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

import Database
import qualified CUDA.PDB       as PDB
--import qualified CUDA.SMVM      as SMVM

--
-- System libraries
--
import Control.Monad                                    (when)
import System.Environment                               (getArgs)
import System.IO
import qualified Foreign.CUDA   as CUDA


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
  (cp,files) <- sequestConfig defaultConfigFile =<< getArgs
  proteins   <- maybe (error "Protein database not specified") (digestFasta cp) (databasePath cp)

  when (verbose cp && not (useCPU cp)) $ do
    dev   <- CUDA.get
    props <- CUDA.props dev
    hPutStrLn stderr $ "> device " ++ shows dev ": " ++ CUDA.deviceName props

  withPDB cp proteins $ \pdb -> mapM_ (search cp pdb) files


--
-- Search the protein database for a match to the experimental spectra
--
search :: ConfigParams Float -> ProteinDatabase Float -> FilePath -> IO ()
search cp pdb fp = do
  dta     <- forceEitherStr `fmap` readDTA fp
  matches <- PDB.searchForMatches cp pdb dta

--  when (verbose cp) $
--    let ref = searchForMatches cp proteins dta in
--    hPutStrLn stderr $ "> verify: " ++ shows (zipWith (==) ref matches) "\n"

  printConfig cp fp dta
  printResults       $! (take (numMatches cp)       matches)
  printResultsDetail $! (take (numMatchesDetail cp) matches)

