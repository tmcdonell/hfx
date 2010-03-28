{-# LANGUAGE CPP #-}
--------------------------------------------------------------------------------
-- |
-- Module    : Utils
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- General utilities that don't really belong anywhere
--
--------------------------------------------------------------------------------

module Utils where

import Data.Char
import Data.Word

#if defined(__GLASGOW_HASKELL__)
import GHC.Base                         (unsafeChr)
#endif

{-
import Data.Int
import Data.Bits
import Data.List
import Data.Array.IArray
import Data.ByteString.Lazy.Char8 (ByteString)
import qualified Data.ByteString.Lazy.Char8 as L


--------------------------------------------------------------------------------
-- Arrays
--------------------------------------------------------------------------------

--
-- Strict left fold over a given subset of a one-dimensional array. The given
-- indices must lie within the bounds of the array.
--
subFoldA' :: (IArray a e, Ix i) => (b -> e -> b) -> b -> a i e -> [i] -> b
subFoldA' f q a i = go q i
    where go s (j:js) = let s' = f s (a ! j)
                        in  s' `seq` go s' js
          go s _      = s

--
-- Strict left fold over a given subset of an array, using the first index as
-- the starting value.
--
subFoldA1' :: (IArray a e, Ix i) => (e -> e -> e) -> a i e -> [i] -> e
subFoldA1' f a i = subFoldA' f (a ! head i) a (tail i)


--------------------------------------------------------------------------------
-- Bits
--------------------------------------------------------------------------------

isPow2 :: Bits a => a -> Bool
isPow2 x = x .&. (x-1) == 0

ceilPow2 :: (Bits a, Integral a) => a -> a
ceilPow2 x | isPow2 x   = x
           | otherwise  = 1 `shiftL` ceiling (logBase 2 (fromIntegral x)::Double)


--------------------------------------------------------------------------------
-- Lazy ByteString
--------------------------------------------------------------------------------

scanlBS' :: (a -> Char -> a) -> a -> ByteString -> [a]
scanlBS' f z = snd . L.foldl' k (z,[z])
    where
        k (c,acc) a = let n = f c a in (n, acc ++ [n])

--
-- Strict left scan over a given subset of a lazy byte string. The scanning
-- function is not limited to conversions between the Char type, but as a
-- consequence lifts the result to a standard list.
--
-- The limiting indices are inclusive, and the seed element is discarded.
--
subScanBS' :: (a -> Char -> a) -> a -> ByteString -> (Int64, Int64) -> [a]
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
subFoldBS' :: (a -> Char -> a) -> a -> ByteString -> (Int64, Int64) -> a
subFoldBS' f q bs (c,n) = (L.foldl' f q . L.take (n-c+1) . L.drop c) bs


--------------------------------------------------------------------------------
-- Statistics
--------------------------------------------------------------------------------

--
-- Statistics about the elements of a list. Uses a single pass algorithm for
-- variance, which may be imprecise if the standard deviation is small relative
-- to the mean.
--
-- Calculates: (minimum, maximum, average, variance, standard deviation, length)
--
stats :: (Floating a, Ord a) => [a] -> (a,a,a,a,a,a)
stats []     = error "Utils.stats: empty list"
stats (x:xs) = finish . foldl' stats' (x,x,x,x*x,1) $ xs
  where
    stats' (mn,mx,s,ss,n) v = (min v mn, max v mx, s+v, ss+v*v, n+1)
    finish (mn,mx,s,ss,n)   = let av    = s/n
                                  var   = (1/(n-1))*ss - (n/(n-1))*av*av
                                  stdev = sqrt var
                              in (mn, mx, av, var, stdev, n)
-}


--
-- Extract components from a three-tuple
--
fst3 :: (a,b,c) -> a
fst3 (a,_,_) = a
{-# INLINE fst3 #-}

snd3 :: (a,b,c) -> b
snd3 (_,b,_) = b
{-# INLINE snd3 #-}

thd3 :: (a,b,c) -> c
thd3 (_,_,c) = c
{-# INLINE thd3 #-}

--
-- Conversion between 'Word8' and 'Char'. Should compile to a no-op.
--
w2c :: Word8 -> Char
#if !defined(__GLASGOW_HASKELL__)
w2c = chr . fromIntegral
#else
w2c = unsafeChr . fromIntegral
#endif
{-# INLINE w2c #-}

c2w :: Char -> Word8
c2w = fromIntegral . ord
{-# INLINE c2w #-}


--------------------------------------------------------------------------------
-- Either
--------------------------------------------------------------------------------

--
-- Pulls a "Right" value out of an Either construct. If the either is a "Left",
-- raises an exception with that string.
--
-- Stolen from the Missing-H package.
--
forceEither :: (Show e) => Either e a -> a
forceEither (Left  e) = error (show e)
forceEither (Right a) = a

forceEitherStr :: Either String a -> a
forceEitherStr (Left  e) = error e
forceEitherStr (Right a) = a

