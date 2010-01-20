--------------------------------------------------------------------------------
--
-- Module    : Time
-- Copyright : (c) 2009 Trevor L. McDonell
-- License   : BSD
--
-- Simple timing benchmarks
--
--------------------------------------------------------------------------------

module Time where

import Numeric
import System.CPUTime
import Control.Monad


-- Timing
--
data Time = Time { cpu_time :: Integer }

type TimeUnit = Integer -> Integer

picosecond, nanosecond, microsecond, millisecond, second :: TimeUnit
picosecond  n = n
nanosecond  n = n `div` 1000
microsecond n = n `div` 1000000
millisecond n = n `div` 1000000000
second      n = n `div` 1000000000000

getTime :: IO Time
getTime = Time `fmap` getCPUTime

timeIn :: TimeUnit -> Time -> Integer
timeIn u (Time t) = u t

elapsedTime :: Time -> Time -> Time
elapsedTime (Time t1) (Time t2) = Time (t2 - t1)

showTime :: Time -> String
showTime t
  | timeIn picosecond  t <= 1000 = shows (timeIn picosecond  t) " ps"
  | timeIn nanosecond  t <= 1000 = showT (timeIn picosecond  t) " ns"
  | timeIn microsecond t <= 1000 = showT (timeIn nanosecond  t) " us"
  | timeIn millisecond t <= 1000 = showT (timeIn microsecond t) " ms"
  | otherwise                    = showT (timeIn millisecond t) " s"
  where
    showT n = showFFloat (Just 3) (fromInteger n / 1000::Float)

-- Simple timing/benchmarking
--
bracketTime :: IO a -> IO (Time, a)
bracketTime f = do
  t1 <- getTime
  r  <- f
  t2 <- getTime
  return (elapsedTime t1 t2, r)

{-# NOINLINE benchmark #-}
benchmark
  :: Int                -- Number of times to repeat test
  -> IO a               -- Test to run
  -> IO b               -- Finaliser to before measuring elapsed time
  -> IO (Time,a)
benchmark n testee finaliser = do
  t1    <- getTime
  (r:_) <- replicateM n testee
  _     <- finaliser
  t2    <- getTime
  return (elapsedTime t1 t2, r)

