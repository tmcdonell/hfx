--------------------------------------------------------------------------------
-- |
-- Module    : Sequence
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- Meta-module reexporting sequence related stuff
--
--------------------------------------------------------------------------------

module Sequence
  (
    -- Reading sequence files
    module Sequence.Fasta,

    -- Digesting protein sequences
    module Sequence.Index,
    module Sequence.Fragment,

    -- Locating sequences from index keys
    module Sequence.Location,

    -- Identifying and matching sequences
    module Sequence.Search,
    module Sequence.Match
  )
  where

import Sequence.Fasta
import Sequence.Index
import Sequence.Fragment
import Sequence.Location
import Sequence.Search
import Sequence.Match

