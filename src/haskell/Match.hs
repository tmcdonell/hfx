--------------------------------------------------------------------------------
-- |
-- Module    : Match
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
--------------------------------------------------------------------------------

module Match where

import Sequence

--
-- A structure to store the result of a peptide/spectrum similarity test
--
type MatchCollection = [Match]

data Match = Match
  {
    fragment :: Fragment,       -- The sequence fragment (header and amino acid chain)
    scoreXC  :: Float           -- Sequest cross-correlation score
--  scoreSP  :: (Int, Int)      -- Matched ions / total ions
  }
  deriving (Show)

instance Eq Match where
  Match f s == Match f' s' = f == f' && (s-s')/(s+s'+0.0005) < 0.0005

