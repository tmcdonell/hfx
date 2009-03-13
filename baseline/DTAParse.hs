-- DTAParse.hs - Jason WH Wong (9/12/08)
-- Parser to read experimental ms/ms spectrum in DTA format
--
module DTAParse (parseSpec, getPrecur, getParent, getSpec) where

-- GHC libraries
import System.Environment (getArgs)
import Data.List (sortBy, maximum, unfoldr)
import Data.Ord(comparing)

-- Hackage libraries
import Text.ParserCombinators.Parsec
import qualified Text.ParserCombinators.Parsec.Token as P
import Text.ParserCombinators.Parsec.Language (emptyDef)

lexer = P.makeTokenParser emptyDef
float = P.float lexer
naturalOrFloat = P.naturalOrFloat lexer

-- Holds a spectrum of type (a b c) where a = parent mass
--                                        b = charge
--                                        c = actual spectrum
--
data Spectrum = Spectrum Double Integer [(Double,Double)]

instance Show Spectrum where
         show (Spectrum parent charge vals) = (show parent) ++ " " ++ (show charge) ++ "\n" ++ (show vals)

-- External methods for retrieving input experimental spectrum data
--
getPrecur :: Spectrum -> Double
getPrecur (Spectrum p c s) = (p+(fromInteger (c-1)))/(fromInteger c)

getParent :: Spectrum -> Double
getParent (Spectrum p c s) = p

getSpec :: Spectrum -> [(Double,Double)]
getSpec (Spectrum p c s) = s

-- Methods to process spectrum according to Jimmy Eng's paper - sqrtNorm -> Square roots all intensity values
--                                                              makeStdSpec -> bins mz values and convert to range (0..2000)
--                                                              normSpec -> Normalise spectrum to relative intensities in 200 Da windows
--                                                              xcorrSpec -> Convert to fast xcorr spectrum y'
--
processSpec :: [(Double,Either Integer Double)] -> Spectrum
processSpec i = Spectrum (fst (head i)) (snd(remLeft (head i))) (xcorrSpec (normSpec (makeStdSpec (map sqrtNorm (remLow (map remRight (tail i)) 200))) 200))


remRight :: (Double, Either Integer Double) -> (Double,Double)
remRight x = (fst x, case snd x of Left xs -> fromInteger xs
                                   Right xs -> xs)

remLeft :: (Double, Either Integer Double) -> (Double,Integer)
remLeft x = (fst x, case snd x of Left xs -> xs
                                  Right xs -> round xs)

makeStdSpec :: [(Double,Double)] -> [(Double,Double)]
makeStdSpec s = map (\x -> calcIntn x s) (take 2000 (map (\x -> (x,0)) (iterate (+1) 0)))

calcIntn :: (Double,Double) -> [(Double,Double)] -> (Double,Double)
calcIntn (x,y) s = (x,(maximum(map (\a -> bin (x,y) a) s)))

bin :: (Double,Double) -> (Double,Double) -> Double
bin (x,y) (mz,intn) = if truncate x == truncate mz then intn else y

remLow :: [(Double,Double)] -> Int -> [(Double,Double)]
remLow s l = (sortBy (comparing fst) (take l (reverse (sortBy (comparing snd) s))))

sqrtNorm :: (Double,Double) -> (Double,Double)
sqrtNorm s = (fst s, sqrt (snd s))

normSpec :: [(Double,Double)] -> Int -> [(Double,Double)]
normSpec c win = case (length c) of
                       0 -> []
                       _ -> (map (\x -> norm x (maximum (snd (unzip (take win c))))) (take win c))++(normSpec (drop win c) win)

norm :: (Double,Double) -> Double -> (Double,Double)
norm s m = case m of
                0 -> ((fst s), 0)
                _ -> ((fst s), ((snd s)/m))

xcorrSpec :: [(Double,Double)] -> [(Double,Double)]
xcorrSpec s = zip (fst(unzip s)) (zipWith (-) (snd(unzip s)) (map (/150) (xcorrSpec' (snd(unzip s)))))

xcorrSpec' :: [Double] -> [Double]
xcorrSpec' s = foldl (\x y -> makey' x s y) (take 2000 (repeat 0)) [-75..75]

makey' :: [Double] -> [Double] -> Int -> [Double]
makey' s orig shift= zipWith (+) s (shiftSpec orig shift)

shiftSpec :: [Double] -> Int -> [Double]
shiftSpec s shift = case shift of
                              x| x <0 -> (drop shift s)++(take (2000+shift) (repeat 0))
                               | x >0 -> (take shift (repeat 0))++(take (2000-shift) s)
                               | x == 0 -> take 2000 (repeat 0)

-- Parser methods
--
vals  :: Parser (Double,Either Integer Double)
vals  = do{ ds <- float
            ; ds2 <- naturalOrFloat
            ; return (ds, ds2)
            }

parseSpec    :: Parser Spectrum
parseSpec   = do{ v <- sepBy1 vals separator
                ; return (processSpec v)
                }

separator   :: Parser ()
separator   = skipMany newline
