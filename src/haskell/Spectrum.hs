--------------------------------------------------------------------------------
-- |
-- Module    : Spectrum
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- Meta-module reexporting spectrum related spectrum related stuff
--
--------------------------------------------------------------------------------

module Spectrum
  (
    -- Data Structures
    MS2Data,
    Spectrum,

    -- File Formats
    readDTA,
    readMGF
  )
  where

import Spectrum.Data
import Spectrum.DTA
import Spectrum.MGF

