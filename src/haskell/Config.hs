--------------------------------------------------------------------------------
-- |
-- Module    : Config
-- Copyright : (c) 2009 Trevor L. McDonell
-- License   : BSD
--
-- Functions to process user configuration parameters
--
--------------------------------------------------------------------------------

module Config
  (
    ConfigParams(..),
    sequestConfig,
    getAAMass
  )
  where

import Mass
import Utils

import Control.Monad
import Data.Char
import Data.List
import Data.Maybe
import System
import System.Console.GetOpt
import System.Directory
import System.IO
import Text.Show.Functions ()
import Text.ParserCombinators.Parsec

import Data.Ix
import Data.Vector.Unboxed (Vector)
import qualified Data.Vector.Unboxed as U


--------------------------------------------------------------------------------
-- Data structures
--------------------------------------------------------------------------------

--
-- A ginormous data structure to hold all of the configurable parameters
-- threaded throughout the program
--
data ConfigParams = ConfigParams
    {
        databasePath        :: Maybe FilePath,  -- Path to the protein database file

        --
        -- Enzyme search parameters
        --
        massTolerance       :: Float,                    -- Search peptides within ± this value of the spectrum mass
        removePrecursorPeak :: Bool,                     -- Remove a ±5 da window surrounding the precursor mass
        missedCleavages     :: Int,                      -- Number of missed cleavage sites to consider
        digestionRule       :: (Char -> Bool, String),   -- Protein fragmentation rule and description text
        minPeptideMass      :: Float,                    -- Minimum mass of peptides to be considered
        maxPeptideMass      :: Float,                    -- Maximum peptide mass

        --
        -- Allow static modifications to an amino acid mass, which affects every
        -- occurrence of that residue/terminus
        --
        aaMassTable         :: Vector Float,
        aaMassTypeMono      :: Bool,

        --
        -- Process configuration
        --
        useCPU              :: Bool,

        --
        -- Output configuration
        --
        verbose             :: Bool,
        numMatches          :: Int,             -- Number of matches to show (summary statistics)
        numMatchesDetail    :: Int              -- Number of full protein descriptions to show
    }
    deriving (Show)


--------------------------------------------------------------------------------
-- Helper Functions
--------------------------------------------------------------------------------

--
-- Chemical modifications to amino acids can be considered in the search by
-- changing the amino acid masses used to calculate the masses of the peptides.
-- The modified values are stored in the configuration parameters data
-- structure, and returned by the function below
--
-- XXX: Shouldn't really live here...
--
--getAAMass       :: (Fractional a, Storable a) => ConfigParams a -> Char -> a
getAAMass       :: ConfigParams -> Char -> Float
getAAMass cp aa =  aaMassTable cp U.! index ('A','Z') aa

{-# INLINE getAAMass #-}


--------------------------------------------------------------------------------
-- Options Processing
--------------------------------------------------------------------------------

--
-- The main option processing function, will process commands from file together
-- with the command line arguments. Returns a complete configuration object and
-- the list of non-option arguments, which is typically the list of input
-- spectrum files to analyse.
--
--sequestConfig :: (Fractional a, Storable a, Read a) => FilePath -> [String] -> IO (ConfigParams a, [String])
sequestConfig :: FilePath -> [String] -> IO (ConfigParams, [String])
sequestConfig fp argv = do
    p  <- doesFileExist fp
    cp <- if p then readConfigFile fp baseParams
               else return baseParams

    --
    -- Parse the command line options
    --
    case getOpt RequireOrder options argv of
        (a,n,[]) -> do cp' <- foldl' (>>=) (return cp) a
                       return (initializeAAMasses cp', n)
        (_,_,e)  -> error (unlines e)


--
-- Read a configuration file, returning a modified configuration set as well as
-- list of unprocessed options.
--
--readConfigFile :: (Fractional a, Storable a, Read a) => FilePath -> ConfigParams a -> IO (ConfigParams a)
readConfigFile :: FilePath -> ConfigParams -> IO ConfigParams
readConfigFile fp cp = do
    params   <- parseFromFile configFile fp

    --
    -- Apply amino-acid modifications. Theres are a special case for parameter
    -- files only.
    --
    let (cp', args) = processMods cp (forceEither params)

    --
    -- Parse file, getting a list of actions which we then apply to the supplied 
    -- configuration record. Any unrecognised options cause an exception
    --
    case getOpt RequireOrder (tail options) (map toOpt args) of
        (a,[],[]) -> foldl' (>>=) (return cp') a
        (_,n,e)   -> error . concat $ if (not.null) n then n else e

    where
        --
        -- Convert a (name,value) pair into the name=value form. This is
        -- required to handle parameters with optional arguments.
        --
        toOpt (k,v) = ("--"++k) ++ "=" ++ v


--
-- Traverse the list of parameters searching for instances of the 'add_x'
-- command to apply static modifications to amino acid masses. These are applied
-- immediately to the given configuration set, removing the command from the
-- argument key/value list
--
--processMods :: (Fractional a, Storable a, Read a) => ConfigParams a -> [(String, String)] -> (ConfigParams a, [(String, String)])
processMods :: ConfigParams -> [(String, String)] -> (ConfigParams, [(String, String)])
processMods config = foldr fn (config,[])
    where
        fn (k,v) (cp,acc) = case stripPrefix "add_" k of
                              Nothing     -> (cp, (k,v) : acc)
                              Just (c:[]) -> (apply cp c v, acc)
                              _           -> error ("Unrecognised modification: " ++ k)

        apply cp c v = cp { aaMassTable = aaMassTable cp U.// [(index ('A','Z') c, read v)] }


--------------------------------------------------------------------------------
-- Defaults
--------------------------------------------------------------------------------

--
-- The basic (almost empty) configuration set
--
--baseParams :: (Fractional a, Storable a) => ConfigParams a
baseParams :: ConfigParams
baseParams =  ConfigParams
    {
        databasePath        = Nothing,

        massTolerance       = 3.0,
        removePrecursorPeak = True,
        missedCleavages     = 2,
        digestionRule       = getDigestionRule 1, -- Trypsin
        minPeptideMass      = 400,
        maxPeptideMass      = 7200,

        aaMassTable         = U.replicate (rangeSize ('A','Z')) 0 U.// [(index ('A','Z') 'C',57.0)],
        aaMassTypeMono      = True,

        useCPU              = False,

        verbose             = False,
        numMatches          = 5,
        numMatchesDetail    = 3
    }

--
-- Add the base residue masses for all amino acid groups to the mass table,
-- which currently only holds user-defined modifications (if any).
--
--initializeAAMasses :: (Fractional a, Storable a) => ConfigParams a -> ConfigParams a
initializeAAMasses :: ConfigParams -> ConfigParams
initializeAAMasses cp = cp { aaMassTable = U.accum (+) (aaMassTable cp) (zipWith f alphabet isolatedAAMass) }
    where
        f c m          = (index ('A','Z') c, m)
        isolatedAAMass = map aaMasses alphabet
        aaMasses       = if aaMassTypeMono cp then getAAMassMono else getAAMassAvg
        alphabet       = ['A'..'Z']


--------------------------------------------------------------------------------
-- Actions
--------------------------------------------------------------------------------

--
-- The enormous processing options database, containing the parameter names,
-- option description, and state transformation action to apply when the given
-- command is encountered.
--
-- This is given to the get-opt command, together with the command line
-- arguments following the parameters read from file.
--
-- NOTE: the parameters option must remain at the head of the list, as this will
-- be dropped when reading parameters from file to prevent silly behaviour...
--
--options :: (Fractional a, Storable a, Read a) => [ OptDescr (ConfigParams a -> IO (ConfigParams a)) ]
options :: [ OptDescr (ConfigParams -> IO ConfigParams) ]
options =
    [ Option "p" ["parameters"]
        (ReqArg readConfigFile "FILE")
        "Read configuration options from file"

    , Option "d" ["database"]
        (ReqArg (\fp cp -> return cp { databasePath = Just fp }) "FILE")
        "Protein database to search"

    , Option "" ["add_[A..Z]"]
        (ReqArg (\_ cp -> return cp) "FLOAT")
        "Static modification to all occurrences of that residue"

    , Option "" ["mass-tolerance"]
        (ReqArg (\v cp -> return cp { massTolerance = read v }) "FLOAT")
        "Search for peptides within this value of the spectrum mass"

    , Option "" ["remove-precursor-peak"]
        (OptArg (\v cp -> return $ case v of
                            Nothing    -> cp { removePrecursorPeak = True }
                            Just []    -> cp { removePrecursorPeak = True }
                            Just (x:_) -> cp { removePrecursorPeak = toLower x `elem` "t1" }) "True|False")
        "Remove a 5 da window surrounding the precursor mass"

    , Option "" ["missed-cleavages"]
        (ReqArg (\v cp -> return cp { missedCleavages = read v }) "INT")
        "Number of missed cleavage sites to consider"

    , Option "" ["digestion-rule"]
        (ReqArg (\v cp -> return cp { digestionRule = getDigestionRule (read v) }) "INT")
        "Digestion rule number to use"

    , Option "" ["min-peptide-mass"]
        (ReqArg (\v cp -> return cp { minPeptideMass = read v }) "FLOAT")
        "Minimum mass of peptide to be considered"

    , Option "" ["max-peptide-mass"]
        (ReqArg (\v cp -> return cp { maxPeptideMass = read v}) "FLOAT")
        "Maximum peptide mass to be considered"

    , Option "" ["aa-mass-type"]
        (ReqArg (\(v:_) cp -> return cp { aaMassTypeMono = toLower v == 'm' }) "Mono|Average")
        "Use monoisotopic or average molecular mass of amino acids"

    , Option "n" ["num-matches"]
        (ReqArg (\v cp -> return cp { numMatches = read v }) "INT")
        "Number of peptide results to show"

    , Option "N" ["num-matches-detail"]
        (ReqArg (\v cp -> return cp { numMatchesDetail = read v}) "INT")
        "Number of full protein descriptions to show"

    , Option "" ["cpu"]
        (NoArg (\cp -> return cp { useCPU = True }))
        "Use CPU backend"

    , Option "v" ["verbose"]
        (NoArg (\cp -> return cp { verbose = True }))
        "Extra output on stderr"

    , Option "V" ["version"]
        (NoArg (\_ -> do hPutStrLn stderr "sequest version 2.71"
                         exitWith ExitSuccess))
        "Print version and exit"

    , Option "h" ["help"]
        (NoArg (\_ -> do prg <- getProgName
                         hPutStrLn stderr (usageInfo ("Usage: " ++ prg ++ " [OPTIONS...] dta_files...") options)
                         hPutStrLn stderr digestionRuleHelp
                         exitWith ExitSuccess))
        "This help text"
    ]


--
-- The list of support digestion rules, paired with a description to be output
-- as part of the help text.
--
digestionRuleHelp :: String
digestionRuleHelp = unlines $ "Digestion Rules:" : map (snd . getDigestionRule) [0..13]

getDigestionRule :: Int -> (Char -> Bool, String)
getDigestionRule 0  = (const False          , "  0:  No enzyme              0      -           -")
getDigestionRule 1  = ((`elem` "KR")        , "  1.  Trypsin                1      KR          P")
getDigestionRule 2  = ((`elem` "FWY")       , "  2.  Chymotrypsin           1      FWY         P")
getDigestionRule 3  = ((== 'R')             , "  3.  Clostripain            1      R           -")
getDigestionRule 4  = ((== 'M')             , "  4.  Cyanogen_Bromide       1      M           -")
getDigestionRule 5  = ((== 'W')             , "  5.  IodosoBenzoate         1      W           -")
getDigestionRule 6  = ((== 'P')             , "  6.  Proline_Endopept       1      P           -")
getDigestionRule 7  = ((== 'E')             , "  7.  Staph_Protease         1      E           -")
getDigestionRule 8  = ((== 'K')             , "  8.  Trypsin_K              1      K           P")
getDigestionRule 9  = ((== 'R')             , "  9.  Trypsin_R              1      R           P")
getDigestionRule 10 = ((== 'D')             , "  10. AspN                   0      D           -")
getDigestionRule 11 = ((`elem` "FWYL")      , "  11. Cymotryp/Modified      1      FWYL        P")
getDigestionRule 12 = ((`elem` "ALIV")      , "  12. Elastase               1      ALIV        P")
getDigestionRule 13 = ((`elem` "ALIVKRWFY") , "  13. Elastase/Tryp/Chymo    1      ALIVKRWFY   P")

getDigestionRule _  = error "getDigestionRule: unknown enzyme digestion rule"


--------------------------------------------------------------------------------
-- Parse a configuration file
--------------------------------------------------------------------------------

--
-- The left hand side of an identifier is a C-like string, which must begin with
-- a letter or underscore, and may subsequently also contain numbers and dashes.
--
identifier :: Parser String
identifier =
    do  c  <- letter <|> char '_'
        cs <- many (letter <|> digit <|> char '_' <|> char '-')
        return (c:cs)
    <?> "identifier"

--
-- The right hand side is any sequence of characters until the end of line or
-- the newline symbol, stripping trailing whitespace
--
value :: Parser String
value =
    do  c  <- noneOf "#\n\r"
        cs <- manyTill anyChar (try comment <|> try eol <|> eof)
        return (c:cs)
    <?> "value"

--
-- Comments begin with the pound sign, and continue to the end of the line
--
comment :: Parser ()
comment =  char '#' >> skipMany (noneOf "\n\r")
       <?> "comment"

--
-- The end of line character. Different operating systems use different
-- characters to mark the end-of-line, so just look for all combinations
--
eol :: Parser ()
eol =  many1 (oneOf "\n\r") >> return ()
   <?> "end of line"


--
-- A configuration key/value pair
--
item :: Parser (String, String)
item = liftM2 (\k v -> (k, rstrip v)) key val
    where
        key    = skipMany space >> identifier
        val    = skipMany space >> char '=' >> skipMany space >> value
        rstrip = reverse . dropWhile isSpace . reverse

--
-- A line can be a empty, contain a comment, or a configuration item
--
line :: Parser (Maybe (String, String))
line =  skipMany space >> content
    where content =  try (comment >> return Nothing)
                 <|> Just `fmap` item

--
-- Get just the name/value pairs from the file
--
configFile :: Parser [(String, String)]
configFile =  liftM catMaybes (many line)

