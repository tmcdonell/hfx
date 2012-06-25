{-# LANGUAGE TupleSections #-}
--------------------------------------------------------------------------------
-- |
-- Module    : Sequence.Index
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- Reading and writing digested sequence database
--
--------------------------------------------------------------------------------

module Sequence.Index where

import Mass
import Config
import Sequence.Fragment

import Data.Ix
import Data.Char
import Data.List
import Data.Maybe
import Data.Binary
import Data.Binary.Get

import qualified Codec.Compression.GZip     as GZip
import qualified Data.ByteString.Lazy.Char8 as L
import qualified Data.Vector.Generic        as G


--
-- Output the digested sequence database together with the relevant digestion
-- rules used to generate it, in a form suitable for the configuration parser.
--
-- Annoyingly, binary instances for generic vector match with String, so pack
-- into a ByteString first. Probably not too bad, as 'encode' would have done
-- this anyway.
--
writeIndex :: ConfigParams -> FilePath -> SequenceDB -> IO ()
writeIndex cp fp db
  = L.writeFile fp
  . compress
  $ L.concat [header, content]
  where
    compress    = if ".gz" `isSuffixOf` fp then GZip.compress else id
    content     = encode db
    header      = encode . L.pack . ('\n':) . unlines $
      [ "database         = " ++ (fromJust (databasePath cp))
      , "missed-cleavages = " ++ (show (missedCleavages cp))
      , "digestion-rule   = " ++ (takeWhile isDigit . dropWhile isSpace . snd $ digestionRule cp)
      , "min-peptide-mass = " ++ (show (minPeptideMass cp))
      , "max-peptide-mass = " ++ (show (maxPeptideMass cp))
      , "aa-mass-type     = " ++ (if aaMassTypeMono cp then "Mono" else "Average")
      ]
      ++ G.ifoldr' aaMod [] (aaMassTable cp)

    aaMod i m mods =
      let aa  = chr (ord 'A' + i)
          det = m - if aaMassTypeMono cp then getAAMassMono aa else getAAMassAvg aa
      in
      if det == 0 then mods
                  else ("add_" ++ [aa] ++ " = " ++ show det) : mods


--
-- Read an indexed sequence database. Update the configuration options with
-- those used to generate the database.
--
readIndex :: ConfigParams -> FilePath -> IO (ConfigParams, SequenceDB)
{-# INLINE readIndex #-}
readIndex cp fp = do
  let decompress = if ".gz" `isSuffixOf` fp then GZip.decompress else id

  f     <- decompress `fmap` L.readFile fp
  let (opt,f',_) = runGetState get f 0
      db         = decode f'
      table      = G.replicate (rangeSize ('A','Z')) 0

  cp'   <- readConfig (L.unpack opt) fp (cp {aaMassTable = table, databasePath = Just fp})
  cp' `seq` db `seq` return (cp', db)

