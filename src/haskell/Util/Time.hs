--------------------------------------------------------------------------------
--
-- Module    : Util.Time
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- Simple timing benchmarks
--
--------------------------------------------------------------------------------

module Util.Time where

import System.CPUTime
import Control.Monad


-- Timing
--
data Time = Time { cpuTime :: Integer }

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


-- Simple timing/benchmarking
--
bracketTime :: IO a -> IO (Time, a)
bracketTime f = do
  t1 <- getTime
  r  <- f
  t2 <- r `seq` getTime
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
  return (elapsedTime t1 t2 `divT` n, r)
  where
    divT (Time t) x = Time (t `div` toInteger x)

