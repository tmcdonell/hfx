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
import Sequest
import Sequence
import PrettyPrint

--
-- System libraries
--
import Control.Monad            (when)
import System.Environment       (getArgs)
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
  sequences  <- maybe (error "Protein database not specified") readFasta (databasePath cp)

  when (verbose cp && not (useCPU cp)) $ do
    dev   <- CUDA.get
    props <- CUDA.props dev
    hPutStrLn stderr $ "> device " ++ shows dev ": "
                                   ++ CUDA.deviceName props
                                   ++ " (compute " ++ shows (CUDA.computeCapability props) ")"

  withSeqs cp sequences $ \pdb -> mapM_ (search cp pdb) files


--
-- Search the protein database for a match to the experimental spectra
--
search :: ConfigParams -> ProteinDatabase -> FilePath -> IO ()
search cp pdb fp = do
  dta     <- forceEitherStr `fmap` readDTA fp
  matches <- searchForMatches cp pdb dta

  printConfig cp fp dta
  printResults       $! (take (numMatches cp)       matches)
  printResultsDetail $! (take (numMatchesDetail cp) matches)

