{-
 - Useful data-mining functions
 -}

module Util.Stats (clamp, clamp1, stats) where

import Data.List (foldl')


--
-- Clamp a value between a low and high point, and map to the range [0,1]
--
clamp1 :: Floating a => a -> a -> a -> a
clamp1 lo hi val = (val - lo) / (hi - lo)

clamp  :: Floating a => a -> a -> [a] -> [a]
clamp lo hi = map (clamp1 lo hi)


--
-- Compute basic statistics: min, max, mean, variance, standard deviation, and
-- the number of elements in the list
--
stats :: (Floating a, Ord a) => [a] -> (a, a, a, a, a, a)
stats (x:xs) = finish . foldl' stats' (x, x, x, x*x, 1) $ xs

--
-- The incremental step which adds an element to the accumulator
--
stats' (mn, mx, s, ss, n) x = ( min x mn
                              , max x mx
                              , s  + x
                              , ss + x*x
                              , n  + 1 )

--
-- Calculate mean, variance and standard deviation from the values stored in the
-- accumulator
--
finish (mn, mx, s, ss, n)   = (mn, mx, av, va, stdev, n)
    where
        av    = s / n
        va    = (1 / (n-1)) * ss - (n / (n-1)) * av*av
        stdev = sqrt va

