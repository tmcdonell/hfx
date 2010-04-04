--------------------------------------------------------------------------------
-- |
-- Module    : Sequence.Match
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
--------------------------------------------------------------------------------

module Sequence.Match where

import Sequence.Fragment
import qualified Data.ByteString.Lazy as L

--
-- A structure to store the result of a peptide/spectrum similarity test
--
type MatchCollection = [Match]

data Match = Match
  {
    fragment :: Fragment,         -- The sequence fragment (header and amino acid chain)
    scoreXC  :: Float,            -- Sequest cross-correlation score
    matchSP  :: ([Bool], [Bool])  -- Matched ion sequences
--  scoreSP  :: (Int, Int)        -- Matched ions / total ions
  }
  deriving (Show)

instance Eq Match where
  Match f s p == Match f' s' p' =  f == f' && p == p' && (s-s')/(s+s'+0.0005) < 0.0005


scoreSP :: Match -> (Int, Int)
scoreSP (Match f _ (b,y)) = (matched, total)
  where
    boolToInt True  = 1
    boolToInt False = 0

    total   = fromIntegral (L.length (fragdata f) - 5) * 2
    matched = sum . map boolToInt $ b ++ y

