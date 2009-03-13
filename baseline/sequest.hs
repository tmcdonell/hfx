-- sequest.hs - Jason WH Wong (9/12/08)
-- Processes database entries and runs Sequest
--

-- custom libaray
import DTAParse

-- GHC libraries
import Data.Maybe
import Data.List
import Data.Ord(comparing)
import System.Environment
import System.Environment (getArgs)

-- Hackage libraries
import Text.ParserCombinators.Parsec
import Bio.Sequence
import Bio.Sequence.SeqData

-- Holds spectral data (a, b) where a = mz
--                                  b = intensity
--
data Spec = Spec Double Double

instance Show Spec where
         show (Spec x y) = (show x) ++ " " ++ (show y) ++ "\n"

-- Holds peptide data (a, b, c) where a = peptide ID
--                                    b = Sequence
--                                    c = Xcorr score
--
data Peptide = Peptide String String Double

instance Show Peptide where
         show (Peptide id seq xcorr) = (show id) ++ "\t" ++ (show seq) ++ "\t"++ (show xcorr)++"\n"

-- Digests protein from database to make a list of peptides.
--
digestDB :: [Sequence] -> [Peptide]
digestDB d = concat (map digestProtein d)

-- peptides have to have a mass of at least 400 Da
--
digestProtein :: Sequence -> [Peptide]
digestProtein s = zipWith (\x y->Peptide x y (-1)) (repeat (toStr (seqheader s))) (filter (\x-> (calcMassParent x) > 400) (addMissedCleavage (doSimpleDigest ['K','R'] (toStr(seqdata s)))))

-- First argument are a list of amino acids where cleavage occurs. Does not currently support special digestion rules
--
doSimpleDigest :: [Char] -> String -> [String]
doSimpleDigest = unfoldr . doSimpleDigest'

doSimpleDigest' :: [Char] -> String -> Maybe (String, String)
doSimpleDigest' c l
  | null l = Nothing
  | otherwise = Just (if null t then h else insertBy back (head t) h, drop 1 t)
  where (h, t) = break (`elem` c)  l
        back a b  = GT

-- Currently only supports addition of 1 missed cleavage
--
addMissedCleavage :: [String] -> [String]
addMissedCleavage s = s++(map (\x -> concat [(s!!x),(s!!(x+1))]) [0..((length s)-2)])

-- Methods to make theoretical spectrum based on Jimmy Eng's paper
--
makeSpectrum :: String -> [Spec]
makeSpectrum y = addIons (map calcMassY (tails y))

makeSpectrum' :: Peptide -> [Spec]
makeSpectrum' (Peptide x y z) = addIons (map calcMassY (tails y))

addIons :: [Spec] -> [Spec]
addIons ss = ss++(map makeYm1Ions ss)++(map makeYp1Ions ss)++(map makeYmNHIons ss)++(map (\x -> makeBIons (head ss) x) ss)++(map (\x -> makeBm1Ions (head ss) x) ss)++(map (\x -> makeBp1Ions (head ss) x) ss)++(map (\x -> makeBmNHIons (head ss) x) ss)++(map (\x -> makeBmH2OIons (head ss) x) ss) -- ++(map (\x -> makeAIons (head ss) x) ss)

makeYm1Ions :: Spec -> Spec
makeYm1Ions (Spec a2 b2) = (Spec (a2-1) 0.25)
makeYp1Ions (Spec a2 b2) = (Spec (a2+1) 0.25)
makeYmNHIons (Spec a2 b2) = (Spec (a2-17.0278) 0.10)
makeBIons :: Spec -> Spec -> Spec
makeBIons (Spec a1 b1) (Spec a2 b2) = (Spec (a1-a2+1.0074) b1)
makeBm1Ions (Spec a1 b1) (Spec a2 b2) = (Spec (a1-a2+1.0074-1) 0.25)
makeBp1Ions (Spec a1 b1) (Spec a2 b2) = (Spec (a1-a2+1.0074+1) 0.25)
makeBmH2OIons (Spec a1 b1) (Spec a2 b2) = (Spec (a1-a2+1.0074-18.0105) 0.10)
makeBmNHIons (Spec a1 b1) (Spec a2 b2) = (Spec (a1-a2+1.0074-17.0278) 0.10)
makeAIons (Spec a1 b1) (Spec a2 b2) = (Spec (a1-a2+1.0074-27.9949) 0.10)

calcMassY :: String -> Spec
calcMassY s = Spec (sum(map getAAMass s)+19.017804) 0.5

getAAMass :: Char -> Double
getAAMass aa = case aa of
                   'A' -> 71.037114000
                   'R' -> 156.101111000
                   'N' -> 114.042927000
                   'D' -> 115.026943000
                   'C' -> 103.009185000+57   -- alkylation
                   'E' -> 129.042593000
                   'Q' -> 128.058578000
                   'G' -> 57.021464000
                   'H' -> 137.058912000
                   'I' -> 113.084064000
                   'L' -> 113.084064000
                   'K' -> 128.094963000
                   'M' -> 131.040485000
                   'F' -> 147.068414000
                   'P' -> 97.052764000
                   'S' -> 87.032028000
                   'T' -> 101.047679000
                   'W' -> 186.079313000
                   'Y' -> 163.063320000
                   'V' -> 99.068414000
                   'U' -> 168.053000
                   _ -> 0

-- Standardise spectrum to range (0..2000)
--
makeStdSpec :: [Spec] -> [Spec]
makeStdSpec s = map (\x -> calcIntn x s) (take 2000 (map (\x -> Spec x 0) (iterate (+1) 0)))

calcIntn :: Spec -> [Spec] -> Spec
calcIntn (Spec x y) s = Spec x (maximum(map (\a -> bin (Spec x y) a) s))

bin :: Spec -> Spec -> Double
bin (Spec x y) (Spec mz intn) = if truncate x == truncate mz then intn else y

calcMassParent :: String -> Double
calcMassParent s = sum(map getAAMass s)+19.017804

filterHits :: [Peptide] -> Double -> [Peptide]
filterHits orig h = filter (\(Peptide x y z) -> h > (calcMassParent y)-3 && h < (calcMassParent y)+3) orig

-- Runs the sequest algorithm for comparing
--
runSequest :: [Peptide] -> Double -> [(Double,Double)] -> [Peptide]
runSequest p h s = map (\x -> runSequest' x h s) p

runSequest' :: Peptide -> Double -> [(Double,Double)] -> Peptide
runSequest' (Peptide x y z) h s = (Peptide x y (calcXCorr (snd (unzip s)) (vectorise (makeStdSpec (makeSpectrum y)))))

-- Dot product calculation for Xcorr score... seems to be very slow!!
calcXCorr :: [Double] -> [Double] -> Double
calcXCorr e t = sum $ zipWith (*) e t

vectorise :: [Spec] -> [Double]
vectorise s = map (\(Spec x y) -> y) s

-- Sort final results
--
sortFinal :: [Peptide] -> [Peptide]
sortFinal s = reverse (sortBy (comparing getXcorr) s)

getXcorr :: Peptide -> Double
getXcorr (Peptide a b c) = c

main = do
    [input,fasta] <- getArgs
    spec <- parseFromFile parseSpec input
    db <- readFasta fasta

    let parent :: Double
        parent = case spec of Left x -> 0
                              Right x -> getParent x

    let pepList :: [Peptide]
        pepList = filterHits (digestDB db) parent

    print "Peptides to be searched:"
    print (length pepList)

    print "Searching..."
    case spec of Left x -> print "Error parsing"
                 Right x -> print (take 5 (sortFinal((runSequest pepList (getParent x) (getSpec x)))))

