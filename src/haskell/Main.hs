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
defaultConfigFile =  "sequest.params"


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
  sdb <- makeSeqDB' cp fp
  withDeviceDB cp sdb $ forM_ dta . search cp sdb


makeSeqDB' :: ConfigParams -> FilePath -> IO SequenceDB
{-# INLINE makeSeqDB' #-}
makeSeqDB' cp fp = do
  when (verbose cp) $ hPutStr stderr "Loading database ... " >> hFlush stdout
  (t,sdb) <- bracketTime $ makeSeqDB cp fp

  when (verbose cp) $ do
    let lf = G.length (dbFrag   sdb)
        li = G.length (dbIon    sdb)
        lp = G.length (dbIonSeg sdb) - 1
    hPutStrLn stderr $ "done (" ++ showTime t ++ ")"
    hPutStrLn stderr $ "  " ++ shows lp " proteins"
    hPutStrLn stderr $ "  " ++ shows li " amino acids, " ++ showFFloatSIBase 1024 (fromIntegral (li * 4)  :: Double) "B"
    hPutStrLn stderr $ "  " ++ shows lf " peptides, "    ++ showFFloatSIBase 1024 (fromIntegral (lf * 12) :: Double) "B"

  return sdb


--
-- Search the protein database for a match to the experimental spectra
--
search :: ConfigParams -> SequenceDB -> DeviceSeqDB -> FilePath -> IO ()
search cp sdb ddb fp =
  readDTA fp >>= \r -> case r of
    Left  s   -> hPutStrLn stderr s
    Right dta -> do
      (t,matches) <- bracketTime $ searchForMatches cp sdb ddb dta
      when (verbose cp) $ hPutStrLn stderr ("Elapsed time: " ++ showTime t)

      printConfig cp fp dta
      printResults       $! take (numMatches cp)       matches
      printResultsDetail $! take (numMatchesDetail cp) matches

