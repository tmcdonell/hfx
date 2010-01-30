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
import Data.List
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


-- Show time in nearest SI unit
--
showTime :: Time -> String
showTime (Time t) = showsT . nubBy (\a b -> a == ' ' && b == ' ') $ ' ':si_unit:"s"
  where
    showsT  = showFFloat (Just 3) (fromInteger t / (1000 ^^ pow :: Double))
    pow     = min 4 . floor $ logBase 1000 (fromInteger t :: Double)
    si_unit = "pnum " !! pow


-- Shows with SI prefix
--
showFFloatSI :: RealFloat a => a -> ShowS
showFFloatSI n = showString . nubBy (\a b -> a == ' ' && b == ' ') $ showFFloat (Just 3) n' (' ':si_unit:[])
  where
    n'      = n / (1000 ^^ (pow-4))
    pow     = max 0 . min 8 . (+) 4 . floor $ logBase 1000 n
    si_unit = "pnum kMGT" !! pow


-- Show the rate of "things / second", with SI unit prefix
--
showRateSI :: Integral a => a -> Time -> String -> String
showRateSI _ (Time 0) unit = "-- " ++ unit ++ "/s"
showRateSI n (Time t) unit = showFFloatSI (fromInteger t * fromIntegral n / 1E12::Double) (unit ++ "/s")


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

