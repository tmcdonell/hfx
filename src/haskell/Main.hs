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

import Sequest
import CUDA.Database
import qualified CUDA.PDB       as PDB
--import qualified CUDA.SMVM      as SMVM

--
-- System libraries
--
import Control.Monad                                    (when)
import System.Environment                               (getArgs)
import System.IO
import qualified Foreign.CUDA as C


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
  proteins   <- maybe (error "Protein database not specified") readFasta (databasePath cp)

  when (verbose cp && not (useCPU cp)) $ do
    dev   <- C.get
    props <- C.props dev
    hPutStrLn stderr $ "> device " ++ shows dev ": " ++ C.deviceName props

  withPDB cp proteins $ \pdb -> mapM_ (search cp proteins pdb) files


--
-- Search the protein database for a match to the experimental spectra
--
search :: ConfigParams Float -> [Protein Float] -> ProteinDatabase Float -> FilePath -> IO ()
search cp proteins pdb fp = do
  dta     <- forceEitherStr `fmap` readDTA fp
  matches <- PDB.searchForMatches cp pdb dta

  when (verbose cp) $
    let ref = searchForMatches cp proteins dta in
    hPutStrLn stderr $ "> verify: " ++ shows (zipWith verify ref matches) "\n"

  printConfig cp fp dta
  printResults       $! (take (numMatches cp)       matches)
  printResultsDetail $! (take (numMatchesDetail cp) matches)

  where
    verify (Match p s) (Match p' s') = p == p' && (s-s')/(s+s'+0.0005) < 0.0005

