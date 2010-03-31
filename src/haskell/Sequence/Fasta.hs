--------------------------------------------------------------------------------
-- |
-- Module    : Sequence.Fasta
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- Reading sequence data in the FASTA format. Each sequence consists of a head
-- (prefixed by a ''>'') and a set of lines containing the sequence data.
--
--------------------------------------------------------------------------------

module Sequence.Fasta where

import Data.List                        (isSuffixOf)
import Control.Applicative              ((<$>))

import qualified Bio.Sequence                   as F
import qualified Bio.Sequence.Fasta             as F
import qualified Data.ByteString.Lazy.Char8     as L
import qualified Codec.Compression.GZip         as GZip


type Protein = F.Sequence F.Amino

--
-- Lazily read sequences from a FASTA-formatted file. This is identical to the
-- code of the bio package, except the sequence type is cast on output and GZip
-- compressed files are deflated as necessary.
--
readFasta :: FilePath -> IO [Protein]
{-# INLINE readFasta #-}
readFasta fp = map F.castToAmino . F.mkSeqs . L.lines . prepare <$> L.readFile fp
  where
    prepare = if ".gz" `isSuffixOf` fp then GZip.decompress
                                       else id

--
-- Count the number of sequences in a file. Each sequence consist of a header
-- and a set of lines containing the sequence data. GZip compressed files are
-- inflated if necessary, which is the only addition over the bio package
-- function of the same name.
--
countSeqs :: FilePath -> IO Int
{-# INLINE countSeqs #-}
countSeqs fp = length . describe . prepare <$> L.readFile fp
  where
    describe = filter (('>' ==) . L.head) . filter (not . L.null) . L.lines
    prepare  = if ".gz" `isSuffixOf` fp then GZip.decompress
                                        else id


