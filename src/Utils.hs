{-
 - General utilities that don't really belong anywhere
 -}

module Utils where

import Text.Printf
import Bio.Sequence

import Data.Int
import Data.Array
import qualified Data.ByteString.Lazy.Char8 as L

--------------------------------------------------------------------------------
-- Constants
--------------------------------------------------------------------------------

--
-- The monoisotopic mass of several elements and molecules
--
massH2O, massNH3, massCO, massO, massH :: Float
massH2O = 18.01056
massNH3 = 17.02655
massCO  = 27.9949
massO   = 16.0013
massH   = 1.0078246

--------------------------------------------------------------------------------
-- Arrays
--------------------------------------------------------------------------------

--
-- Strict left fold over a given subset of a one-dimensional array. The given
-- indices must lie within the bounds of the array.
--
subFoldA' :: Ix k => (a -> b -> a) -> a -> Array k b -> [k] -> a
subFoldA' f q a i = go q i
    where go s (j:js) = let s' = f s (a ! j)
                        in  s' `seq` go s' js
          go s _      = s

--
-- Strict left fold over a given subset of an array, using the first index as
-- the starting value.
--
subFoldA1' :: Ix k => (a -> a -> a) -> Array k a -> [k] -> a
subFoldA1' f a i = subFoldA' f (a ! head i) a (tail i)


--------------------------------------------------------------------------------
-- Lazy ByteString
--------------------------------------------------------------------------------

--
-- A more generic byte string scan. The scanning function is not limited to
-- conversions between the Char type, but as a consequence lifts the result to a
-- standard list
--
scanlBS                 :: (a -> Char -> a) -> a -> (Int64, Int64) -> L.ByteString -> [a]
scanlBS f q (c,n) bs
    | L.null bs         = []
    | otherwise         =
    q : (case c <= n && c < L.length bs of
            False -> []
            True  -> scanlBS f (f q (L.index bs c)) (c+1,n) bs)


--------------------------------------------------------------------------------
-- Results
--------------------------------------------------------------------------------

-- 
-- Dumb print wrapper (i.e. doesn't calculate optimum column widths, etc)
--
printResults :: [(Double, Sequence)] -> IO ()
printResults =  mapM_ printMatch 
    where
        printMatch (score, match) = printf "%7.3f  %20s  %s\n"
            score
            (toStr (seqdata match))
            (toStr (seqheader match))

