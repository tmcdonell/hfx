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
    argv           <- getArgs
    (cp, dtaFiles) <- sequestConfig defaultConfigFile argv

    when (verbose cp && not (useCPU cp)) $ do
      props <- C.get >>= C.props
      hPutStrLn stderr $ "Device      :: " ++ C.deviceName props

    mapM_ (search cp) dtaFiles


--
-- Search the protein database for a match to the experimental spectra
--
search :: ConfigParams CFloat -> FilePath -> IO ()
search cp fp = do
    dta <- forceEitherStr `fmap` readDTA fp
    db  <- maybe (error "Protein database not specified") readFasta (databasePath cp)

    printConfig cp fp dta
    let ref = searchForMatches cp db dta

    matches <- if useCPU cp
      then return ref
      else C.searchForMatches cp db dta

    when (verbose cp && not (useCPU cp)) $
      hPutStrLn stderr ("Valid: " ++ shows (zipWith verify ref matches) "\n")

    printResults       $! (take (numMatches cp)       matches)
    printResultsDetail $! (take (numMatchesDetail cp) matches)

    where
      verify (Match p s) (Match p' s') = p == p' && (s-s')/(s+s'+0.0005) < 0.0005

