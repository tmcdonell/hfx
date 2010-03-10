--------------------------------------------------------------------------------
-- |
-- Module    : Location
-- Copyright : (c) 2010 Trevor L. McDonell
-- License   : BSD
--
-- Locate sequences fragments from index keys
--
--------------------------------------------------------------------------------

module Location (SeqMap(..), lookupSeq) where


import Mass
import Sequence

import Data.Int
import Data.Word
import Numeric.Search.Range

import qualified Bio.Sequence               as F
import qualified Data.ByteString.Lazy.Char8 as L
import qualified Data.Vector                as V
import qualified Data.Vector.Unboxed        as U
import qualified Data.Vector.Generic        as G


data SeqMap = SeqMap
  {
    smSeqs      :: V.Vector Protein,
    smZones     :: U.Vector Int,
    smResiduals :: U.Vector Float,
    smOffset    :: U.Vector Word32,
    smTerminals :: U.Vector (Word32, Word32)
  }


--
-- Find the last index in the ordered array whose value is less than or equal to
-- the given search element.
--
searchVector :: (Ord a, G.Vector v a) => v a -> a -> Maybe Int
searchVector vec x =
  searchFromTo (\i -> vec G.! i > x) 0 (G.length vec - 1)


--
-- Extract an amino acid sequence between the given terminals, including the
-- flanking residuals (if present)
--
slice :: L.ByteString -> (Int64,Int64) -> L.ByteString
slice s (c,n) = L.pack [ca,'.'] `L.append` aa `L.append` L.pack ['.',na]
  where
    aa = L.take (n-c+1) . L.drop c $ s
    l  = L.length s - 2
    ca = if c > 0 then L.index s (c-1) else '-'
    na = if n < l then L.index s (n+1) else '-'


--
-- Locate a fragment from its index key, returning the parent description text
-- and sequence representation, including flanking residuals.
--
lookupSeq :: SeqMap -> Int -> Maybe Fragment
lookupSeq sm key = do
  idx <- searchVector (smZones sm) key
  pro <- smSeqs sm `G.indexM` idx
  let o     = smOffset sm    G.! idx
      (c,n) = smTerminals sm G.! key
      res   = smResiduals sm G.! key + massH2O + massH
      aa    = slice (F.seqdata pro) (fromIntegral (c-o), fromIntegral (n-o))
      hdr   = F.seqheader pro

  return $ Fragment res hdr aa

