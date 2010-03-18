{-# LANGUAGE CPP          #-}
{-# LANGUAGE BangPatterns #-}
--------------------------------------------------------------------------------
-- |
-- Module    : Sequence
-- Copyright : (c) 2010 Trevor L. McDonell
-- License   : BSD
--
-- Reading and processing sequence related stuff
--
--------------------------------------------------------------------------------

module Sequence
  (
    SequenceDB(..), Fragment(..), SKey,
    makeSeqDB, lookup, fraglabel
  )
  where

import Mass
import Config

import Prelude                                  hiding (lookup)
import Data.Int
import Data.Char
import Data.List                        (unfoldr, isSuffixOf)
import Data.Word
import Control.Monad                    (foldM)
import Control.Applicative              ((<$>))

import qualified Bio.Sequence                   as F
import qualified Bio.Sequence.Fasta             as F
import qualified Data.ByteString.Lazy           as L
import qualified Data.ByteString.Lazy.Char8     as LC
import qualified Codec.Compression.GZip         as GZip

import qualified Data.Vector                    as V
import qualified Data.Vector.Unboxed            as U
import qualified Data.Vector.Generic            as G
import qualified Data.Vector.Generic.Mutable    as GM
import qualified Data.Vector.Fusion.Stream      as S
import qualified Data.Vector.Fusion.Stream.Size as S

#if defined(__GLASGOW_HASKELL__)
import GHC.Base                         (unsafeChr)
#endif

#define PHASE_STREAM [1]
#define PHASE_INNER  [0]

#define INLINE_STREAM INLINE PHASE_STREAM
#define INLINE_INNER  INLINE PHASE_INNER


--
-- Extract components from a three-tuple
--
fst3 :: (a,b,c) -> a
fst3 (a,_,_) = a
{-# INLINE fst3 #-}

snd3 :: (a,b,c) -> b
snd3 (_,b,_) = b
{-# INLINE snd3 #-}

thd3 :: (a,b,c) -> c
thd3 (_,_,c) = c
{-# INLINE thd3 #-}

--
-- Conversion between 'Word8' and 'Char'. Should compile to a no-op.
--
w2c :: Word8 -> Char
#if !defined(__GLASGOW_HASKELL__)
w2c = chr . fromIntegral
#else
w2c = unsafeChr . fromIntegral
#endif
{-# INLINE w2c #-}

c2w :: Char -> Word8
c2w = fromIntegral . ord
{-# INLINE c2w #-}



--------------------------------------------------------------------------------
-- Reading FASTA Files
--------------------------------------------------------------------------------

type Protein = F.Sequence F.Amino

--
-- Lazily read sequences from a FASTA-formatted file. This is identical to the
-- code of the bio package, except the sequence type is cast on output and GZip
-- compressed files are deflated as necessary.
--
readFasta :: FilePath -> IO [Protein]
{-# INLINE readFasta #-}
readFasta fp = map F.castToAmino . F.mkSeqs . LC.lines . prepare <$> L.readFile fp
  where
    prepare = if ".gz" `isSuffixOf` fp then GZip.decompress
                                       else id

--
-- Count the number of sequences in a file. Each sequence consist of a header
-- and a set of lines containing the sequence data. GZip compressed files are
-- inflated if necessary, which is the only addition over the bio package
-- function of the same name.
--
countSeqs :: FilePath -> IO Int
{-# INLINE countSeqs #-}
countSeqs fp = length . headers . prepare <$> L.readFile fp
  where
    headers = filter (('>' ==) . LC.head) . filter (not . L.null) . LC.lines
    prepare = if ".gz" `isSuffixOf` fp then GZip.decompress
                                       else id


--------------------------------------------------------------------------------
-- Sequences
--------------------------------------------------------------------------------

type SKey = Int32

--
-- A collection of protein sequences
--
data SequenceDB = SeqDB
  {
    dbHeader  :: V.Vector L.ByteString,          -- sequence ion description headers
    dbIon     :: U.Vector Word8,                 -- flattened array of amino character codes
    dbIonSeg  :: U.Vector Int32,                 -- segmenting information for ions
    dbFrag    :: U.Vector (Float, Int32, Int32), -- (residual mass, c-idx, n-idx)
    dbFragSeg :: U.Vector Int32                  -- fragment segmenting information for deriving parent
  }
  deriving Show

--
-- Locate a particular sequence in the database
--
lookup :: SKey -> SequenceDB -> Maybe Fragment
lookup =  undefined


--
-- Generate a sequence database from the given file. This is a flattened
-- representation of ion character codes and fragment information, together with
-- segmenting information describing the boundaries of the original Proteins
-- they derive from.
--
makeSeqDB :: ConfigParams -> FilePath -> IO SequenceDB
{-# INLINE makeSeqDB #-}
makeSeqDB cp fp = do
  ns <- countSeqs fp
  db <- readFasta fp

  let iseg = ionSeg  ns db
      hdr  = headers ns db
      ions = ionSeq (fromIntegral (U.last iseg)) db
      nf   = countFrags cp ions

  f  <- GM.new nf               -- maximum number as estimated above, trim later
  fs <- GM.new (ns+1)           -- include the zero index

  --
  -- OPT: make versions that operate directly from (boxed) streams?
  --
  let put !i !v = do GM.unsafeWrite f i v
                     return (i+1)

  let fill (!i,!n) (!x,!y) = do GM.unsafeWrite fs i (fromIntegral n)
                                n' <- foldM put n (digest cp ions (x,y))
                                return (i+1, n')

  --
  -- Now iterate over all of the protein sequences, keeping track of the segment
  -- number and number of fragments generated so far, so we can also fill in the
  -- fragment segmenting information.
  --
  nf' <- snd <$> foldM fill (0,0) (G.toList (G.zip iseg (G.tail iseg)))
  GM.unsafeWrite fs ns (fromIntegral nf')

  f'  <- G.unsafeFreeze (GM.take nf' f)
  fs' <- G.unsafeFreeze fs

  return $ SeqDB hdr ions iseg f' fs'


--------------------------------------------------------------------------------
-- Fragments
--------------------------------------------------------------------------------

--
-- Sequence fragments
--
data Fragment = Frag
  {
    fragmass   :: Float,                -- Total mass of this fragment
    fragheader :: L.ByteString,         -- Full header describing this fragment
    fragdata   :: L.ByteString          -- Fragment sequence data, including flanking residuals
  }
  deriving (Eq, Show)

--
-- Fragment label (first word of full header)
--
fraglabel :: Fragment -> L.ByteString
fraglabel = head . LC.words . fragheader


--------------------------------------------------------------------------------
-- Digestion
--------------------------------------------------------------------------------

--
-- Estimate the number of fragments that will be produced from a collection of
-- peptide sequences. This is an upper bound, as it does not take into account:
--   * fragments outside of the designated mass bounds
--   * split points at the end of a sequence, producing one fragment and not two
--
-- The latter accounts for a rather small discrepancy (<1%), while the former
-- may result in significant overestimation (>10%).
--
countFrags :: ConfigParams -> U.Vector Word8 -> Int
{-# INLINE countFrags #-}
countFrags cp = ((missedCleavages cp + 1) * ) . count
  where
    rule  = fst (digestionRule cp) . w2c
    count = U.length . U.filter rule


--
-- Extract sequence headers
--
headers :: Int -> [Protein] -> V.Vector L.ByteString
{-# INLINE headers #-}
headers n db
  = G.unstream
  $ S.fromList (map F.seqheader db) `S.sized` S.Exact n


--
-- Extract the amino acid character codes and segmenting information from the
-- given list of proteins
--
ionSeq :: Int -> [Protein] -> U.Vector Word8
{-# INLINE ionSeq #-}
ionSeq n db
  = G.unstream
  $ S.fromList (concatMap (L.unpack . F.seqdata) db) `S.sized` S.Exact n


ionSeg :: Int -> [Protein] -> U.Vector Int32
{-# INLINE ionSeg #-}
ionSeg n db
  = G.unstream
  $ S.fromList (scanl (+) 0 $ map (fromIntegral . L.length . F.seqdata) db) `S.sized` S.Exact (n+1)


--
-- Digest a single protein sequence with the given enzyme and cleavage rules,
-- removing any fragments outside the broad mass range of interest.
--
-- The mass of the peptide is the sum of the amino acid residue masses plus the
-- mass of the water molecule released in forming the peptide bond (plus one;
-- from Eq. 1 of Eng.[1])
--
digest :: ConfigParams -> U.Vector Word8 -> (Int32,Int32) -> [(Float,Int32,Int32)]
{-# INLINE digest #-}
digest cp ions (!x,!y) = filter (inrange . fst3) . splice cp . fragment cp x $ segment
  where
    segment   = G.unsafeSlice (fromIntegral x) (fromIntegral (y-x)) ions
    inrange m = let m' = m + massH2O + massH
                in  minPeptideMass cp <= m' && m' <= maxPeptideMass cp

--
-- Split a protein sequence with the given digestion rule. This outputs the
-- total residual mass of the fragment, together with the indices of the C- and
-- N- terminals in global index coordinates.
--
fragment :: ConfigParams -> Int32 -> U.Vector Word8 -> [(Float,Int32,Int32)]
{-# INLINE fragment #-}
fragment cp gid ions = unfoldr step (gid,ions)
  where
    rule = fst (digestionRule cp) . w2c
    mass = G.foldl' (+) 0 . G.map (getAAMass cp . w2c)

    step (!c,!v)
      | G.null v  = Nothing
      | otherwise = case G.findIndex rule v of
                      Nothing -> Just ((mass v, c, c + fromIntegral (G.length v)), (0, G.empty))
                      Just !i -> let n = fromIntegral i
                                     h = U.take (i+1) v
                                     t = U.drop (i+1) v
                                 in Just ((mass h, c, c+n), (c+n+1,t))


--
-- Generate additional sequences from missed cleavages of sequential, adjacent
-- peptide fragments.
--
splice :: ConfigParams -> [(Float,Int32,Int32)] -> [(Float,Int32,Int32)]
{-# INLINE splice #-}
splice cp = loop
  where
    loop []     = []
    loop (p:ps) = scanl join p (take n ps) ++ loop ps

    n        = missedCleavages cp
    join a b = (fst3 a + fst3 b, snd3 a, thd3 b)

