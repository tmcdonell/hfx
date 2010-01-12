--------------------------------------------------------------------------------
-- |
-- Module    : PrettyPrint
-- Copyright : (c) 2009 Trevor L. McDonell
-- License   : BSD
--
-- Pretty printing with not-so-ugly code
--
--------------------------------------------------------------------------------

module PrettyPrint
  (
    Pretty(..),

    printConfig,
    printResults,
    printResultsDetail
  )
  where

import Mass
import Config
import Protein
import Sequest
import Spectrum

import Numeric
import Data.List
import Data.Maybe
import Text.PrettyPrint
import System.IO

import qualified Data.ByteString.Lazy.Char8 as B

class Pretty a where ppr :: a -> Doc

instance Pretty Bool            where ppr = text . show
instance Pretty Char            where ppr = char
instance Pretty B.ByteString    where ppr = ppr . B.unpack
instance Pretty a => Pretty [a] where ppr = hcat . map ppr
instance Pretty Peptide         where ppr = text . slice

--------------------------------------------------------------------------------
-- Doc -> IO
--------------------------------------------------------------------------------

displayIO :: Doc -> IO ()
displayIO =  putStrLn . (++ "\n") . render

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

printConfig :: ConfigParams -> FilePath -> Spectrum -> IO ()
printConfig cp fp spec = displayIO . ppAsRows 0 . map (intersperse (text "::")) $
    [ [text "Spectrum"   , text fp]
    , [text "Database"   , text (fromJust (databasePath cp))]
    , [text "Enzyme"     , hsep $ map text (tail enzyme)]
    , [text "(M+H)+ Mass", float mass <+> text "~" <+> float (realToFrac $ massTolerance cp)]
    ]
    where
        mass   = realToFrac $ (precursor spec * charge spec) - ((charge spec * massH) - 1)
        enzyme = words . snd . digestionRule $ cp

--------------------------------------------------------------------------------
-- Results -> Render
--------------------------------------------------------------------------------

title :: [[Doc]]
title = map (map text) [[" # ", " (M+H)+  ", "deltCn", "XCorr", "Reference", "Peptide"],
                        ["---", "---------", "------", "-----", "---------", "-------"]]

toDoc :: Int -> Float -> Match -> [Doc]
toDoc n s0 (Match pep sc) =
    [ space <> int n <> char '.'
    , float' (pmass pep)
    , float' (realToFrac ((s0 - sc)/s0))
    , float' (realToFrac sc)
    , text   (name (parent pep))
    , text   (slice pep)
    ]
    where float' = text . flip (showFFloat (Just 4)) ""

toDocDetail :: Int -> Match -> [Doc]
toDocDetail n (Match pep _)  =
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

