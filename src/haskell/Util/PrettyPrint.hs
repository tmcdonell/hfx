--------------------------------------------------------------------------------
-- |
-- Module    : Util.PrettyPrint
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- Pretty printing with not-so-ugly code
--
--------------------------------------------------------------------------------

module Util.PrettyPrint
  (
    Pretty(..),

    printConfig,
    printResults,
    printResultsDetail
  )
  where

import Mass
import Config
import Sequence.Match
import Sequence.Fragment
import Spectrum.Data

import Numeric
import Data.List
import Data.Maybe
import Text.PrettyPrint

import qualified Text.PrettyPrint.Boxes     as B
import qualified Data.ByteString.Lazy.Char8 as L

class Pretty a where ppr :: a -> Doc

instance Pretty Bool            where ppr = text . show
instance Pretty Char            where ppr = char
instance Pretty L.ByteString    where ppr = ppr . L.unpack
instance Pretty a => Pretty [a] where ppr = hcat . map ppr
--instance Pretty (Peptide a)     where ppr = text . slice

--------------------------------------------------------------------------------
-- Doc -> IO
--------------------------------------------------------------------------------

displayIO :: Doc -> IO ()
displayIO =  putStrLn . (++ "\n") . render

displayBoxIO :: B.Box -> IO ()
displayBoxIO =  putStrLn . B.render

{-
--
-- stolen from $fptools/ghc/compiler/utils/Pretty.lhs
--
-- This code has a BSD-style license
--
printDoc :: Mode -> Handle -> Doc -> IO ()
printDoc m hdl doc = do
  fullRender m cols 1.5 put done doc
  hFlush hdl
  where
    put (Chr c) next  = hPutChar hdl c >> next
    put (Str s) next  = hPutStr  hdl s >> next
    put (PStr s) next = hPutStr  hdl s >> next

    done = hPutChar hdl '\n'
    cols = 80
-}

--------------------------------------------------------------------------------
-- Configuration -> Render
--------------------------------------------------------------------------------

printConfig :: ConfigParams -> FilePath -> MS2Data -> IO ()
printConfig cp fp ms2 = displayIO . ppAsRows 0 . map (intersperse (text "::")) $
    [ [text "Title"      , ppr (ms2info ms2)]
    , [text "Spectrum"   , text fp]
    , [text "Database"   , text (fromJust (databasePath cp))]
    , [text "Enzyme"     , hsep $ map text (tail enzyme)]
    , [text "(M+H)+ Mass", float mass <+> text "~" <+> float (realToFrac $ massTolerance cp)]
    ]
    where
        mass   = realToFrac $ (ms2precursor ms2 * ms2charge ms2) - ((ms2charge ms2 * massH) - 1)
        enzyme = words . snd . digestionRule $ cp

--------------------------------------------------------------------------------
-- Results -> Render
--------------------------------------------------------------------------------

title :: [[Doc]]
title = map (map text) [[" # ", " (M+H)+  ", "deltCn", "XCorr", "Ions", "Reference", "Peptide"],
                        ["---", "---------", "------", "-----", "----", "---------", "-------"]]

toDoc :: Int -> Float -> Match -> [Doc]
toDoc n s0 m@(Match frag sc _) =
    [ space <> int n <> char '.'
    , float' (fragmass frag)
    , float' (realToFrac ((s0 - sc)/s0))
    , float' (realToFrac sc)
    , int    (fst sp) <> char '/' <> int (snd sp)
    , ppr    (fraglabel frag)
    , ppr    (fragdata  frag)
    ]
    where float' = text . flip (showFFloat (Just 4)) ""
          sp     = scoreSP m

toDocDetail :: Int -> Match -> B.Box
toDocDetail n (Match frag _ _) = B.hsep 2 B.top
    [ B.alignHoriz B.right 3 $ B.text (shows n ".")
    , B.para B.left cols     $ L.unpack (fragheader frag)
    ]
    where
        cols = 95       -- for a total width of 100 columns

printResults   :: MatchCollection -> IO ()
printResults m =  displayIO . ppAsRows 1 . (++) title . snd . mapAccumL k 1 $ m
    where
        s0    = scoreXC (head m)
        k n z = (n+1, toDoc n s0 z)

printResultsDetail :: MatchCollection -> IO ()
printResultsDetail =  displayBoxIO . B.vcat B.left . snd . mapAccumL k 1
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
        pad w = map (\x -> x <> hcat (replicate (w - len x + q) space))

