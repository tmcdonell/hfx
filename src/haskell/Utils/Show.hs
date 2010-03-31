--------------------------------------------------------------------------------
--
-- Module    : Utils.Show
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
--------------------------------------------------------------------------------

module Utils.Show where

import Numeric
import Data.List
import Utils.Time



--
-- Show time in seconds to nearest SI unit
--
showTime :: Time -> String
showTime (Time t) = showsT . nubBy (\a b -> a == ' ' && b == ' ') $ ' ':si_unit:"s"
  where
    showsT  = showFFloat (Just 3) (fromInteger t / (1000 ^^ pow :: Double))
    pow     = min 4 . floor $ logBase 1000 (fromInteger t :: Double)
    si_unit = "pnµm " !! pow


--
-- Shows with SI prefix
--
showFFloatSI :: RealFloat a => a -> ShowS
showFFloatSI =  showFFloatSIBase 1000

showFFloatSIBase :: RealFloat a => a -> a -> ShowS
showFFloatSIBase b n = showString . nubBy (\x y -> x == ' ' && y == ' ')
                     $ showFFloat (Just 3) n' [ ' ', si_unit ]
  where
    n'      = n / (b ^^ (pow-4))
    pow     = max 0 . min 8 . (+) 4 . floor $ logBase b n
    si_unit = "pnµm kMGT" !! pow


--
-- Show the rate of "things / second", with SI unit prefix
--
showRateSI :: Integral a => a -> Time -> String -> String
showRateSI _ (Time 0) unit = "-- " ++ unit ++ "/s"
showRateSI n (Time t) unit = showFFloatSI (fromIntegral n / fromInteger t * 1E12::Double) (unit ++ "/s")

