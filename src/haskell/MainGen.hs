--------------------------------------------------------------------------------
-- |
-- Module    : MainGen
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- Generate a list of peptides that meet the digestion criteria
--
--------------------------------------------------------------------------------

module Main where

import Config
import Sequence.Index
import Sequence.Fragment

import Data.Maybe
import System.IO
import System.Environment
import qualified Data.Vector.Generic as G


main :: IO ()
main = do
  args   <- getArgs
  (cp,f) <- sequestConfig "hfx.params" args
  let fp =  fromMaybe (error "Protein database not specified") (databasePath cp)

  db <- makeSeqDB cp fp
  if null f
     then writeIndex stdout cp db
     else withFile (head f) WriteMode (\h -> writeIndex h cp db)

  hPutStrLn stderr $ "Database: " ++ fp
  hPutStrLn stderr $ " # amino acids: " ++ (show . G.length . dbIon    $ db)
  hPutStrLn stderr $ " # proteins:    " ++ (show . G.length . dbHeader $ db)
  hPutStrLn stderr $ " # peptides:    " ++ (show . G.length . dbFrag   $ db)

