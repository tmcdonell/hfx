--------------------------------------------------------------------------------
-- |
-- Module    : Sequence.Location
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- Locate sequences fragments from index keys
--
--------------------------------------------------------------------------------

module Sequence.Location (SKey, lookup) where

import Mass
import Sequence.Fragment
import Util.Misc

import Prelude                          hiding (lookup)
import Numeric.Search.Range

import qualified Data.ByteString.Lazy   as L
import qualified Data.Vector.Generic    as G


--------------------------------------------------------------------------------
-- Database Search
--------------------------------------------------------------------------------

type SKey = Int

--
-- Find the last index in the ordered array whose value is less than or equal to
-- the given search element. Binary search, O(log n).
--
searchVector :: (Ord a, G.Vector v a) => v a -> a -> Maybe Int
searchVector vec x =
  searchFromTo (\i -> vec G.! i > x) 0 (G.length vec - 1)


--
-- Locate a particular sequence in the database
--
lookup :: SequenceDB -> SKey -> Maybe Fragment
lookup db k = do
  -- Index of the sequence this fragment derives from
  --
  seqIdx <- searchVector (G.tail (dbFragSeg db)) (fromIntegral k)

  -- Extract the supporting information for the fragment
  --
  let (res,c,n) = dbFrag db G.! k
      [a,b]     = G.toList $ G.slice seqIdx 2 (dbIonSeg db)
      hdr       = dbHeader db G.! seqIdx
      aa        = G.toList $ G.slice (fromIntegral c) (fromIntegral (n-c+1)) (dbIon db)
      ca        = if c > a   then dbIon db G.! (fromIntegral c-1) else c2w '-'
      na        = if n < b-1 then dbIon db G.! (fromIntegral n+1) else c2w '-'

  return $ Fragment (res + massH + massH2O) hdr (L.pack $ [ca,c2w '.'] ++ aa ++ [c2w '.',na])

