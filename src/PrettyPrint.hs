--------------------------------------------------------------------------------
-- |
-- Module    : PrettyPrint
-- Copyright : (c) 2009 Trevor L. McDonell
-- License   : BSD
--
-- Pretty printing with not-so-ugly code
--
--------------------------------------------------------------------------------

module PrettyPrint where

import Config
import Protein
import Sequest
import Spectrum
import AminoAcid

import Numeric
import Data.List
import Data.Maybe
import Text.PrettyPrint


--------------------------------------------------------------------------------
-- Doc -> IO
--------------------------------------------------------------------------------

displayIO :: Doc -> IO ()
displayIO =  putStrLn . flip (++) "\n" . render


--------------------------------------------------------------------------------
-- Configuration -> Render
--------------------------------------------------------------------------------

printConfig :: ConfigParams -> FilePath -> Spectrum -> IO ()
printConfig cp fp spec = displayIO . ppAsRows 0 . map (intersperse (text "::")) $
    [ [text "Spectrum"   , text fp]
    , [text "Database"   , text (fromJust (databasePath cp))]
    , [text "Enzyme"     , hsep $ map text (tail enzyme)]
    , [text "(M+H)+ Mass", float mass <+> text "~" <+> float (massTolerance cp)]
    ]
    where
        mass   = (precursor spec * charge spec) - ((charge spec * massH) - 1)
        enzyme = words . snd . digestionRule $ cp

--------------------------------------------------------------------------------
-- Results -> Render
--------------------------------------------------------------------------------

title :: [[Doc]]
title = map (map text) [[" # ", " (M+H)+  ", "deltCn", "XCorr", "Reference", "Peptide"],
                        ["---", "---------", "------", "-----", "---------", "-------"]]

toDoc :: Int -> Float -> Match -> [Doc]
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
toDocDetail n (Match _ pep)  =
    [ space <> int n <> char '.'
    , text (description (parent pep))
    ]

printResults   :: MatchCollection -> IO ()
printResults m =  displayIO . ppAsRows 1 . (++) title . snd . mapAccumL k 1 $ m
    where
        s0    = scoreXC (head m)
        k n z = (n+1, toDoc n s0 z)

printResultsDetail   :: MatchCollection -> IO ()
printResultsDetail m =  displayIO . ppAsRows 1 . snd . mapAccumL k 1 $ m
    where
        k n z = (n+1, toDocDetail n z)


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

