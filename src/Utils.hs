{-
 - General utilities that don't really belong anywhere
 -}

module Utils where

import Data.Int
import Data.Array
import qualified Data.ByteString.Lazy.Char8 as L


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
-- Strict left scan over a given subset of a lazy byte string. The scanning
-- function is not limited to conversions between the Char type, but as a
-- consequence lifts the result to a standard list. The limiting indices are
-- inclusive.
--
-- The limiting indices are inclusive, and the seed element is discarded.
--
subScanBS' :: (a -> Char -> a) -> a -> L.ByteString -> (Int64, Int64) -> [a]
subScanBS' f q bs (c,n) = go q [c..n]
    where go s (i:is) = let s' = f s (L.index bs i)
                        in  s' `seq` s' : go s' is
          go _ []     = []

--
-- Strict left fold over a given subset of a lazy byte string. The limiting
-- indices are inclusive.
--
-- This is (somewhat surprisingly) faster than using an index-based method
-- similar to that above.
--
subFoldBS' :: (a -> Char -> a) -> a -> L.ByteString -> (Int64, Int64) -> a
subFoldBS' f q bs (c,n) = (L.foldl' f q . L.take (n-c) . L.drop c) bs


--------------------------------------------------------------------------------
-- Either
--------------------------------------------------------------------------------

--
-- Pulls a "Right" value out of an Either construct. If the either is a "Left",
-- raises an exception with that string.
--
forceEither :: Either String a -> a
forceEither (Left  e) = error e
forceEither (Right x) = x

