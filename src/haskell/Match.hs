--------------------------------------------------------------------------------
-- |
-- Module    : Match
-- Copyright : (c) 2009 Trevor L. McDonell
-- License   : BSD
--
--------------------------------------------------------------------------------

module Match where

import Protein

type MatchCollection a = [Match a]

--
-- A structure to store the result of a peptide/spectrum similarity test
--
data Match a = Match
    {
	parent    :: Protein a,		-- The parent protein
        candidate :: Peptide a,         -- The fragment that was examined
        scoreXC   :: a                  -- Sequest cross-correlation score
--        scoreSP   :: (Int, Int)         -- Matched ions / total ions
    }
    deriving (Show)

instance (Fractional a, Ord a) => Eq (Match a) where
  Match p f s == Match p' f' s' = p == p' && f == f' && (s-s')/(s+s'+0.0005) < 0.0005



