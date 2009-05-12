{-
 - Pretty printing with ugly code
 -}

module PrettyPrint where

import Sequest
import Protein

import Numeric
import Data.List
import Text.PrettyPrint


--------------------------------------------------------------------------------
-- Results -> Render
--------------------------------------------------------------------------------

title :: [[Doc]]
title = map (map text) [[" # ", "  (M+H)+ ", "deltCn", "XCorr", "Reference", "Peptide"],
                        ["---", "---------", "------", "-----", "---------", "-------"]]

display :: [[Doc]] -> IO ()
display =  putStrLn . flip (++) "\n" . render . ppAsRows 1

toDoc :: Int -> Float -> Match -> [Doc]
toDoc _ _  (Match _ Null) = []
toDoc n s0 (Match sc pep) =
    [ space <> int n <> char '.'
    , float' (pmass pep)
    , float' ((s0 - sc)/s0)
    , float' sc
    , text  (name (parent pep))
    , text  (slice pep)
    ]
    where float' = text . flip (showFFloat (Just 4)) ""

toDocDetail :: Int -> Match -> [Doc]
toDocDetail _ (Match _ Null) = []
toDocDetail n (Match _ pep)  =
    [ space <> int n <> char '.'
    , text (description (parent pep))
    ]

printResults   :: MatchCollection -> IO ()
printResults m =  display . (++) title . snd . foldr k (length m,[]) $ m
    where
        s0           = scoreXC (head m)
        k z (n, acc) = (n-1, toDoc n s0 z : acc)

printResultsDetail   :: MatchCollection -> IO ()
printResultsDetail m =  display . snd . foldr k (length m,[]) $ m
    where
        k z (n, acc) = (n-1, toDocDetail n z : acc)


--------------------------------------------------------------------------------
-- Pretty Print
--------------------------------------------------------------------------------

--
-- Display the given grid of renderable data, given as either a list of rows or
-- columns, using the minimum size required for each column. An additional
-- parameter specifies extra space to be inserted between each column.
--
ppAsRows      :: Int -> [[Doc]] -> Doc
ppAsRows q    =  ppAsColumns q . transpose

ppAsColumns   :: Int -> [[Doc]] -> Doc
ppAsColumns q =  vcat . map hsep . transpose . map (\col -> pad (width col) col)
    where
        len   = length . render
        width = maximum . map len
        pad w = map (\x -> x <> (hcat $ replicate (w - (len x) + q) space))

