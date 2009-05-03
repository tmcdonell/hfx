{-
 - Functions to load user configuration parameters
 -}

module Config where

import AminoAcid

import Data.ConfigFile
import Data.Array.Unboxed
import Control.Monad.Error

--------------------------------------------------------------------------------
-- Data structures
--------------------------------------------------------------------------------

data ConfigParams = Params
    {
        databasePath        :: FilePath,        -- Path to the protein database file

        --
        -- Enzyme search parameters
        --
        massTolerence       :: Float,           -- Search peptides within ± this value of the spectrum mass
        removePrecursorPeak :: Bool,            -- Remove a ±5 da window surrounding the precursor mass
        missedCleavages     :: Int,             -- Number of missed cleavage sites to consider
        digestRule          :: (Char -> Bool),  -- Protein fragmentation rule
        minPeptideMass      :: Float,           -- Minimum mass of peptides to be considered
        maxPeptideMass      :: Float,           -- Maximum peptide mass

        --
        -- Allow static modifications to an amino acid mass, which affects every
        -- occurrence of that residue/terminus
        --
        aaMassTable         :: Array Char Float,
        aaMassTypeMono      :: Bool,

        --
        -- Output configuration
        --
        numMatches          :: Int,             -- Number of matches to show (summary statistics)
        numMatchesDetail    :: Int              -- Number of full protein descriptions to show
    }


--
-- The default configuration set
--
defaultParams :: ConfigParams
defaultParams =  Params
    {
        databasePath        = "data/uniprot_sprot_human+trypsin.fasta",

        massTolerence       = 3.0,
        removePrecursorPeak = True,
        missedCleavages     = 1,
        digestRule          = (`elem` "KR"),
        minPeptideMass      = 400,
        maxPeptideMass      = 7200,

        aaMassTable         = listArray ('A','Z') (repeat 0) // [('C',57.0)],
        aaMassTypeMono      = True,

        numMatches          = 5,
        numMatchesDetail    = 3
    }


--------------------------------------------------------------------------------
-- Read/Build
--------------------------------------------------------------------------------

--
-- Read a configuration file in a format that resembles that of an old-style
-- windows .INI file. This will return either a configuration object, or an
-- error string details what went wrong.
--
readParams          :: FilePath -> IO (Either String ConfigParams)
readParams filename =  do
    rv <- runErrorT $ (join . liftIO . readfile emptyCP) filename

    --
    -- Unwrap the result to expose the error string
    --
    return $ case rv of
        Left (ParseError s, _)           -> Left s
        Left (SectionAlreadyExists s, _) -> Left s
        Left (NoSection s, _)            -> Left s
        Left (NoOption s, _)             -> Left s
        Left (OtherProblem s, _)         -> Left s
        Left (InterpolationError s, _)   -> Left s
        Right config                     -> parseConfig config


--
-- Read values from the configuration file and populate the parameters data
-- structure, or return an appropriate error string
--
parseConfig   :: ConfigParser -> Either String ConfigParams
parseConfig _ =  Right (finish defaultParams)
    where
        finish                = initializeAAMasses

        --
        -- Add the base residue mass to the mass table, which currently only
        -- holds user-defined modifications (if any).
        --
        initializeAAMasses cp = cp { aaMassTable = accum (+) (aaMassTable cp) (zip alphabet (isolatedAAMass cp)) }
        isolatedAAMass cp     = map (aaMasses cp) alphabet
        aaMasses cp           = if aaMassTypeMono cp then getAAMassMono else getAAMassAvg
        alphabet              = ['A'..'Z']


--------------------------------------------------------------------------------
-- Extracting
--------------------------------------------------------------------------------

--
-- Chemical modifications to amino acids can be considered in the search by
-- changing the amino acid masses used to calculate the masses of the peptides.
-- The modified values are stored in the configuration parameters data
-- structure, and returned by the function below
--
-- XXX: Shouldn't live here...
--
getAAMass       :: ConfigParams -> Char -> Float
getAAMass cp aa =  aaMassTable cp ! aa

