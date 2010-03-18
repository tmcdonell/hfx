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
import Time
import Config
import Sequest
import Sequence
import PrettyPrint

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
    hPutStrLn stderr $ "device " ++ shows dev ": "
                                 ++ CUDA.deviceName props
                                 ++ ", compute " ++ show (CUDA.computeCapability props)

  --
  -- Load the proteins from file, marshal to the device, and then get to work!
  --
  sdb <- makeSeqDB' cp fp
  withDeviceDB cp sdb $ forM_ dta . search cp sdb


makeSeqDB' :: ConfigParams -> FilePath -> IO SequenceDB
{-# INLINE makeSeqDB' #-}
makeSeqDB' cp fp = do
  when (verbose cp) $ hPutStr stderr "loading database ... " >> hFlush stdout
  (t,sdb) <- bracketTime $ makeSeqDB cp fp

  when (verbose cp) $ do
    hPutStrLn stderr $ "done (" ++ showTime t ++ ")"
    hPutStrLn stderr $ "  " ++ shows (G.length (dbHeader sdb)) " proteins"
    hPutStrLn stderr $ "  " ++ shows (G.length (dbFrag   sdb)) " fragments"

  return sdb


--
-- Search the protein database for a match to the experimental spectra
--
search :: ConfigParams -> SequenceDB -> DeviceSeqDB -> FilePath -> IO ()
search cp sdb ddb fp = do
  readDTA fp >>= \r -> case r of
    Left  s   -> hPutStrLn stderr s
    Right dta -> do
      matches <- searchForMatches cp sdb ddb dta

      printConfig cp fp dta
      printResults       $! take (numMatches cp)       matches
      printResultsDetail $! take (numMatchesDetail cp) matches

