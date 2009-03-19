{-
 - An implementation of the SEQUEST algorithm for fast cross-correlation based
 - identification of protein sequences.
 -
 -
 - References:
 -
 - [1] J. K. Eng, B. Fischer, J. Grossmann, and M. J. MacCoss. "A fast sequest
 -     cross correlation algorithm." Journal of Proteome Research,
 -     7(10):4598-4602, 2008.
 -
 - [2] J. K. Eng, A. L. McCormack, and I. John R. Yates. "An approach to
 -     correlate tandem mass spectral data of peptides with amino acid sequences
 -     in a protein database." Journal of the American Society for Mass
 -     Spectrometry, 5(11):976-989, November 1994.
 -}

module Bio.Sequest where

import Bio.AminoAcid


--
-- The factors that contributed to the collision induced dissociation (CID)
-- process of peptide fragmentation in a tandem-MS experiment are not completely
-- understood, so accurate prediction of fragment ion abundances for a given
-- amino acid sequence are not yet possible.
--
-- As such, a magnitude component is assigned to the predicted mass-to-charge
-- ratio values of the fragment ions based on empirical knowledge.
--
mkIonYm1   (mz, _) = (mz - 1.0,     0.25)
mkIonYp1   (mz, _) = (mz + 1.0,     0.25)
mkIonYmNH  (mz, _) = (mz - 17.0278, 0.10)

mkIonB     (mz, i) (mz', _) = (mz - mz' + 1.0074,           i)
mkIonBm1   (mz, _) (mz', _) = (mz - mz' + 1.0074 - 1.0,     0.25)
mkIonBp1   (mz, _) (mz', _) = (mz - mz' + 1.0074 + 1.0,     0.25)
mkIonBmH2O (mz, _) (mz', _) = (mz - mz' + 1.0074 - 18.0105, 0.10)
mkIonBmNH  (mz, _) (mz', _) = (mz - mz' + 1.0074 - 17.0278, 0.10)

mkIonA     (mz, _) (mz', _) = (mz - mz' + 1.0074 - 27.9949, 0.10)

