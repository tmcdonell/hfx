{-
 - An implementation of the SEQUEST algorithm for fast cross-correlation based
 - identification of protein sequences.
 -
 -
 - References:
 -
 - [1] J. K. Eng, B. Fischer, J. Grossmann, and M. J. MacCoss. "A fast sequest
 -     cross correlation algorithm." Journal of Proteome Research,
 -     7(10):4598-4602, 2008.
 -
 - [2] J. K. Eng, A. L. McCormack, and I. John R. Yates. "An approach to
 -     correlate tandem mass spectral data of peptides with amino acid sequences
 -     in a protein database." Journal of the American Society for Mass
 -     Spectrometry, 5(11):976-989, November 1994.
 -}

module Sequest where

import Spectrum
import AminoAcid

import Bio.Sequence hiding ((!))
import Text.Printf
import Data.List
import Data.Array.Unboxed
import Data.ByteString.Internal (w2c)
import Data.ByteString.Lazy.Internal (ByteString(..))
import qualified Data.ByteString.Lazy as L


--------------------------------------------------------------------------------
-- Amino acid sequences
--------------------------------------------------------------------------------

--
-- Return the mass of the peptide, which is the sum of the amino acid residue
-- masses plus the mass of the water molecule released in forming the peptide
-- bond (plus one; from Eq. 1 of Eng.[1])
--
getPeptideMass :: SeqData -> Double
getPeptideMass =  L.foldl (flip $ (+) . getAAMass . w2c) (18.017804 + 1.0)

--
-- Scan a sequence from the database for linear combinations of amino acids,
-- proceeding from the N to the C terminus.
--
-- The first argument is the list of amino acids where cleavage occurs. This
-- currently does not support special digestion rules.
--
-- Additionally, only consider peptides must have a mass of at least 400 da
-- (atomic mass units)
--
digestProtein     :: String -> Sequence -> [Sequence]
digestProtein c s =  map (\x -> Seq name x Nothing) $ seqs s
    where
        name = seqheader s
        seqs = filter (\x -> getPeptideMass x >= 400) . addCleave . simpleDigest (`elem` c) . seqdata

--
-- Split at the given amino acids
--
simpleDigest         :: (Char -> Bool) -> SeqData -> [SeqData]
simpleDigest _ Empty =  []
simpleDigest p cs    =
    case L.findIndex (p . w2c) cs of
        Nothing                     -> [cs]
        Just n | n+1 == L.length cs -> [cs]
               | otherwise          -> let (a,b) = L.splitAt (n+1) cs in a : simpleDigest p b

--
-- Include the possibility of one missed cleavage
--
addCleave   :: [SeqData] -> [SeqData]
addCleave l =  l ++ zipper l (tail l)
    where
        zipper (a:as) (b:bs) = L.append a b : zipper as bs
        zipper _      _      = []


--------------------------------------------------------------------------------
-- Theoretical Spectrum
--------------------------------------------------------------------------------

--
-- Reconstruct a theoretical "spectrum" form a character-based amino acid
-- sequence. Contains values for the predicted mass-to-charge ratio of fragment
-- ions and a magnitude component.
--
mkAASpec         :: (Int, Int) -> (Double -> Int) -> Sequence -> Spectrum
mkAASpec bnds fn =  mkSpectrum bnds . bin . buildSeq . seqdata
    where
        bin      = map (\(x,y) -> (fn x,y))
        buildSeq = addIons . (map (\s -> (getPeptideMass s, 1.0))) . L.tails

--
-- The factors that contributed to the collision induced dissociation (CID)
-- process of peptide fragmentation in a tandem-MS experiment are not completely
-- understood, so accurate prediction of fragment ion abundances for a given
-- amino acid sequence are not yet possible.
--
-- As such, a magnitude component is assigned to the predicted mass-to-charge
-- ratio values of the fragment ions based on empirical knowledge.
--
-- Note that a factor of 50 has been removed, and the generated spectrum is
-- normalised to one.
--
addIons   :: [(Double, Double)] -> [(Double, Double)]
addIons s =  s ++
    concatMap (\fn -> map fn s)            [mkIonYm1, mkIonYp1, mkIonYmNH] ++
    concatMap (\fn -> map (fn (head s)) s) [mkIonB, mkIonBm1, mkIonBp1, mkIonBmH2O, mkIonBmNH, mkIonA]

    where
        mkIonYm1   (mz, _) = (mz - 1.0,     0.5)
        mkIonYp1   (mz, _) = (mz + 1.0,     0.5)
        mkIonYmNH  (mz, _) = (mz - 17.0278, 0.2)

        mkIonB     (mz, i) (mz', _) = (mz - mz' + 1.0074,           i)
        mkIonBm1   (mz, _) (mz', _) = (mz - mz' + 1.0074 - 1.0,     0.5)
        mkIonBp1   (mz, _) (mz', _) = (mz - mz' + 1.0074 + 1.0,     0.5)
        mkIonBmH2O (mz, _) (mz', _) = (mz - mz' + 1.0074 - 18.0105, 0.2)
        mkIonBmNH  (mz, _) (mz', _) = (mz - mz' + 1.0074 - 17.0278, 0.2)

        mkIonA     (mz, _) (mz', _) = (mz - mz' + 1.0074 - 27.9949, 0.2)


--------------------------------------------------------------------------------
-- Database search
--------------------------------------------------------------------------------

--
-- Search the database for amino acid sequences within a defined mass tolerance.
-- Only peptides which fall within this range will be considered.
--
findCandidates      :: Double -> [Sequence] -> [Sequence]
findCandidates mass =  filterDB . concatMap (digestProtein "KR")
    where
        filterDB = filter (limit . getPeptideMass . seqdata)
        limit x  = (mass - det) <= x && x <= (mass + det)
        det      = max 3 (0.05 / 100 * mass)


findMatch           :: Double -> Spectrum -> [Sequence] -> [(Double, Sequence)]
findMatch mass spec =  finish . map sequest . findCandidates mass
    where
        bnds      = bounds spec
        finish    = sortBy (\(x,_) (y,_) -> compare y x)
        sequest p = (score p, p)
        score     = dot spec . mkAASpec bnds truncate


--------------------------------------------------------------------------------
-- Cross Correlation
--------------------------------------------------------------------------------

--
-- Explicitly de-forested array dot product [P. Walder, Deforestation, 1988]
--
dot     :: (Ix a, Num a, Num e) => Array a e -> Array a e -> e
dot v w =  loop 0 0
    where
        n                      = snd (bounds v)
        loop i acc | i > n     = acc
                   | otherwise = loop (i+1) (v!i * w!i + acc)

--
-- Transform the data from a tandem-MS spectrum into the cross-correlation
-- spectrum that a theoretical peptide sequence will be matched against.
--
-- A spectrum histogram is created, where the mass/charge values are assigned
-- into unit bins with the square root of the input intensity. Measurements that
-- fall into the same bin are accumulated.
--
mkXCorrSpec        :: [(Double, Double)] -> Spectrum
mkXCorrSpec msdata =  mkSpectrum bnds . zip mz $ zipWith (-) itn rt
    where
        -- Could be extracted from the input data (minmax, as in normItn)
        bnds     = (0,2000)

        -- Bin and normalise the input data
        (mz,itn) = unzip . normItn bnds $ mkHist truncate msdata

        -- The processed spectrum:
        --   y' = y_0 - 1/150 ( \sum_{t=-75, t/=0}^{t=75} y_t )
        --      = y_0 - rt
        --
        shiftl _ _ []     = []
        shiftl n f (x:xs) = foldl f x (take n xs) : shiftl n f xs
        shiftr n f        = reverse . shiftl n f . reverse

        rt                = zipWith (\x y -> (x+y)/150) (shiftl 75 (+) itn) (shiftr 75 (+) itn)

--
-- Generate a histogram-like representation of the input spectrum. Values that
-- fall into the same bin have not yet been combined.
--
mkHist    :: (Double -> Int) -> [(Double, Double)] -> [(Int, Double)]
mkHist fn =  sort' . bin'
    where
        bin'    = map (\(x,y) -> (fn x, sqrt y))
        sort'   = sortBy (\(x,_) (y,_) -> compare x y)

--
-- Normalise the intensity measurements in ten equal windows across the
-- specified m/z bounding range. The intensity measurements are normalised to be
-- uniform the range [0,1] in each window.
--
normItn                :: (Int, Int) -> [(Int, Double)] -> [(Int, Double)]
normItn (lo,hi) msdata =  concat $ zipWith clamp limits windows
    where
        windows               = let n = (hi-lo+9) `div` 10 in mkWin ((<) . fst) [n,(n+n)..] msdata
        limits                = map minmax windows

        minmax  []            = let infinity = 1/0 in (infinity,-infinity)
        minmax  ((_,x):xs)    = foldl' minmax' (x, x) xs
        minmax' (mn,mx) (_,x) = (min x mn, max x mx)

        clamp lim             = map (\(x,y) -> (x, clamp' lim y))
        clamp' (l,u) v        = (v-l) / (u-l)

--
-- Split a list by repeatedly applying 'span' for each of the given splitting
-- locations.
--
mkWin                 ::  (a -> t -> Bool) -> [t] -> [a] -> [[a]]
mkWin _ _  []         =  []
mkWin _ [] ls         =  [ls]
mkWin p (n:ns) ls     =  let (xs,ys) = span (`p` n) ls in xs : mkWin p ns ys


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

