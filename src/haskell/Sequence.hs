{-# LANGUAGE BangPatterns, CPP #-}
--------------------------------------------------------------------------------
-- |
-- Module    : Sequence
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- Reading and processing sequence related stuff
--
--------------------------------------------------------------------------------

module Sequence
  (
    SequenceDB(..), DeviceSeqDB(..), Fragment(..),
    makeSeqDB, withDeviceDB, fraglabel
  )
  where

import Mass
import Config
import Util.Misc
import Sequence.Fasta

import Prelude                          hiding (lookup)
import Data.List                        (unfoldr)
import Data.Word
import Control.Monad                    (foldM)
import Control.Applicative              ((<$>))

import qualified Bio.Sequence                   as F
import qualified Data.ByteString.Lazy           as L
import qualified Data.ByteString.Lazy.Char8     as LC

import qualified Data.Vector                    as V
import qualified Data.Vector.Unboxed            as U
import qualified Data.Vector.Generic            as G
import qualified Data.Vector.Generic.Mutable    as GM
import qualified Data.Vector.Fusion.Stream      as S
import qualified Data.Vector.Fusion.Stream.Size as S

import qualified Foreign.CUDA                   as CUDA
import qualified Foreign.CUDA.Util              as CUDA

#define PHASE_STREAM [1]
#define PHASE_INNER  [0]

#define INLINE_STREAM INLINE PHASE_STREAM
#define INLINE_INNER  INLINE PHASE_INNER


--------------------------------------------------------------------------------
-- Sequences
--------------------------------------------------------------------------------

--
-- A collection of protein sequences
--
data SequenceDB = SeqDB
  {
    dbHeader  :: V.Vector L.ByteString,            -- sequence ion description headers
    dbIon     :: U.Vector Word8,                   -- flattened array of amino character codes
    dbIonSeg  :: U.Vector Word32,                  -- segmenting information for ions
    dbFrag    :: U.Vector (Float, Word32, Word32), -- (residual mass, c-idx, n-idx)
    dbFragSeg :: U.Vector Word32                   -- fragment segmenting information for deriving parent
  }
  deriving Show

--
-- An extension of the above referencing memory actually stored on the graphics
-- device. Meant to be used in tandem.
--
data DeviceSeqDB = DevDB
  {
    dbNumIon   :: Int,
    dbNumFrag  :: Int,
    dbIonMass  :: CUDA.DevicePtr Float,
    dbResidual :: CUDA.DevicePtr Float,
    dbTerminal :: (CUDA.DevicePtr Word32, CUDA.DevicePtr Word32)
  }
  deriving Show


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
  -- OPT: make versions that operate directly from (unboxed) streams?
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


--
-- Transfer sequence data to the graphics device and execute some computation
--
withDeviceDB :: ConfigParams -> SequenceDB -> (DeviceSeqDB -> IO a) -> IO a
{-# INLINE withDeviceDB #-}
withDeviceDB cp sdb action =
  let (r,c,n) = G.unzip3 (dbFrag sdb)
      ions    = G.unstream . S.map (getAAMass cp . w2c) . G.stream $ dbIon sdb :: U.Vector Float
      numIon  = G.length (dbIon  sdb)
      numFrag = G.length (dbFrag sdb)
  in
    CUDA.withVector r    $ \d_r    ->
    CUDA.withVector c    $ \d_c    ->
    CUDA.withVector n    $ \d_n    ->
    CUDA.withVector ions $ \d_ions ->
      action (DevDB numIon numFrag d_ions d_r (d_c, d_n))


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


ionSeg :: Int -> [Protein] -> U.Vector Word32
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
digest :: ConfigParams -> U.Vector Word8 -> (Word32,Word32) -> [(Float,Word32,Word32)]
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
fragment :: ConfigParams -> Word32 -> U.Vector Word8 -> [(Float,Word32,Word32)]
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
splice :: ConfigParams -> [(Float,Word32,Word32)] -> [(Float,Word32,Word32)]
{-# INLINE splice #-}
splice cp = loop
  where
    loop []     = []
    loop (p:ps) = scanl join p (take n ps) ++ loop ps

    n        = missedCleavages cp
    join a b = (fst3 a + fst3 b, snd3 a, thd3 b)

