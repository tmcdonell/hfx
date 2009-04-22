{-
 - General utilities that don't really belong anywhere
 -}

module Utils where

import Text.Printf
import Bio.Sequence

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

