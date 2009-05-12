{-
 - Pretty printing with ugly code
 -}

module PrettyPrint where

import Sequest
import Protein

import Text.Printf


--------------------------------------------------------------------------------
-- Ugly Code
--------------------------------------------------------------------------------

printResults   :: MatchCollection -> IO ()
printResults m =  do
    putStrLn "  #     (M+H)+   deltCn  XCorr   Reference           Peptide"
    putStrLn " ---  ---------  ------  ------  ---------           -------"
    go 1 m
    where
        s0 = scoreXC (head m)

        go :: Int -> MatchCollection -> IO ()
        go _ []                         = putStrLn ""
        go _ ((Match _ NullPeptide):_)  = putStrLn ""
        go n ((Match score peptide):ms) = do
            printf " %2d.  %9.4f  %6.4f  %6.4f  %-18s  %s\n" n (pmass peptide) ((s0 - score)/s0) score (name (parent peptide)) (slice peptide)
            go (n+1) ms


printResultsDetail :: MatchCollection -> IO ()
printResultsDetail = go 1
    where
        go :: Int -> MatchCollection -> IO ()
        go _ [] = putStrLn ""
        go _ ((Match _ NullPeptide):_) = putStrLn ""
        go n ((Match _ peptide):ms) = do
            printf " %2d.  %s\n" n (Protein.description (parent peptide))
            go (n+1) ms


