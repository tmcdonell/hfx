{-# LANGUAGE TupleSections #-}
--------------------------------------------------------------------------------
-- |
-- Module    : Main
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- A simple test file for the peptide sequence matching algorithm
--
--------------------------------------------------------------------------------

module Main where

--
-- Custom libraries
--
import Config
import Sequence
import Spectrum
import Util.Time
import Util.Show
import Util.PrettyPrint

--
-- System libraries
--
import Data.Char
import Data.Maybe
import Control.Monad
import System.Environment
import System.IO
import Prelude                          hiding (lookup)

import qualified Data.Vector.Generic    as G
import qualified Foreign.CUDA           as CUDA


--------------------------------------------------------------------------------
-- Program Defaults
--------------------------------------------------------------------------------

defaultConfigFile :: FilePath
defaultConfigFile =  "hfx.params"


--------------------------------------------------------------------------------
-- Main
--------------------------------------------------------------------------------

main :: IO ()
main = do
  (cp,dta) <- sequestConfig defaultConfigFile =<< getArgs
  let fp   = fromMaybe (error "Protein database not specified") (databasePath cp)

  when (verbose cp && not (useCPU cp)) $ do
    dev   <- CUDA.get
    props <- CUDA.props dev
    hPutStrLn stderr $ "Device " ++ shows dev ": "
                                 ++ CUDA.deviceName props ++ ", "
                                 ++ "compute v" ++ shows (CUDA.computeCapability props) ", "
                                 ++ "global memory: " ++ showFFloatSIBase 1024 (fromIntegral $ CUDA.totalGlobalMem props :: Double) "B, "
                                 ++ "core clock: "    ++ showFFloatSI (fromIntegral $ 1000 * CUDA.clockRate props :: Double) "Hz"

  --
  -- Load the proteins from file, marshal to the device, and then get to work!
  --
  (cp',db) <- loadDatabase cp fp
  withDeviceDB cp' db $ forM_ dta . search cp' db


{-# INLINE loadDatabase #-}
loadDatabase :: ConfigParams -> FilePath -> IO (ConfigParams, SequenceDB)
loadDatabase cp fp = do
  (cp',db) <- case suffix fp of
    "fasta" -> (cp,) `fmap` makeSeqDB cp fp
    "index" -> readIndex cp fp
    _       -> error ("Unsupported database type: " ++ show fp)

  when (verbose cp) $ do
    hPutStrLn stderr $ "Database: " ++ fp
    hPutStrLn stderr $ " # amino acids: " ++ (show . G.length . dbIon    $ db)
    hPutStrLn stderr $ " # proteins:    " ++ (show . G.length . dbHeader $ db)
    hPutStrLn stderr $ " # peptides:    " ++ (show . G.length . dbFrag   $ db)

  return (cp',db)
  where
    suffix = map toLower . reverse . takeWhile (/= '.') . reverse


--
-- Search the protein database for a match to the experimental spectra
--
search :: ConfigParams -> SequenceDB -> DeviceSeqDB -> FilePath -> IO ()
search cp db dev fp =
  readMS2Data fp >>= \r -> case r of
    Left  s -> hPutStrLn stderr s
    Right d -> forM_ d $ \ms2 -> do
      (t,matches) <- bracketTime $ searchForMatches cp db dev ms2
      when (verbose cp) $ hPutStrLn stderr ("Elapsed time: " ++ showTime t)

      printConfig cp fp ms2
      printResults           $! take (numMatches cp)       matches
      printResultsDetail     $! take (numMatchesDetail cp) matches
      printIonMatchDetail cp $! take (numMatchesIon cp)    matches

