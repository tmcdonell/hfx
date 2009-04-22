{-
 - Utilities for proteins sequences and their amino acids
 -}

module AminoAcid where

--
-- Return the monoisotopic (ground state) mass of a peptide in atomic mass
-- units, for the given short abbreviation.
--
-- This is corresponds to the mass of the residue, where the mass of the
-- peptide/protein is the mass of the residue plus the mass of water.
--
-- Chemical modifications to amino acids can be considered in the search by
-- changing the amino acid masses used to calculate the masses of the peptides
-- (specifically, reading alterations from file, and defaulting to these
-- values).
--
getAAMass     :: Char -> Double
getAAMass 'A' =   71.037114         -- Alanine              Alg     C3H5NO
getAAMass 'R' =  156.101111         -- Arginine             Arg     C6H12N4O
getAAMass 'N' =  114.042927         -- Asparagine           Asn     C4H6N2O2
getAAMass 'D' =  115.026943         -- Aspartic acid        Asp     C4H5NO3
getAAMass 'C' =  103.009185 + 57.0  -- Cysteine             Cys     C3H5NOS (plus alkylation??)
getAAMass 'E' =  129.042593         -- Glutamic acid        Glu     C5H7NO3
getAAMass 'Q' =  128.058578         -- Glutamine            Gln     C5H8N2O2
getAAMass 'G' =   57.021464         -- Glycine              Gly     C2H3NO
getAAMass 'H' =  137.058912         -- Histidine            His     C6H11NO
getAAMass 'I' =  113.084064         -- Isoleucine           Ile     C6H11NO
getAAMass 'L' =  113.084064         -- Leucine              Leu     C6H11NO
getAAMass 'K' =  128.094963         -- Lysine               Lys     C6H12N2O
getAAMass 'M' =  131.040485         -- Methionine           Met     C5H9NOS
getAAMass 'F' =  147.068414         -- Phenylalanine        Phe     C9H9NO
getAAMass 'P' =   97.052764         -- Proline              Pro     C5H7NO
getAAMass 'O' =  255.15829          -- Pyrrolysine          Pyl     C12H21N3O3
getAAMass 'U' =  168.053            -- Selenocysteine       Sec     C3H5NOSe
getAAMass 'S' =   87.032028         -- Serine               Ser     C3H5NO2
getAAMass 'T' =  101.047679         -- Threonine            Thr     C4H7NO2
getAAMass 'W' =  186.079313         -- Tryptophan           Trp     C11H10N2O
getAAMass 'Y' =  163.06332          -- Tyrosine             Tyr     C9H9NO2
getAAMass 'V' =   99.068414         -- Valine               Val     C5H9NO

--
-- Ambiguous amino acids
--
getAAMass 'B' =  114.10272          -- Aspargine            Asx     C4H8N2O3
getAAMass 'J' =  113.16472          -- Leucine              Xle     C6H13NO2
getAAMass 'Z' =  128.12472          -- Glutamine            Glx     C5H10N2O3
getAAMass 'X' = 0                   -- Unknown              Xaa

getAAMass  x  = error $ "Unknown peptide abbreviation: " ++ [x]


{-
-- Explicitly building an unboxed-array for the lookup-table leads to
-- complicated core, which isn't actually any faster...
--
getAAMass   :: Char -> Double
getAAMass x =  amino ! (fromEnum x)
    where
        amino :: Array Int Double
        amino = listArray (fromEnum 'A', fromEnum 'Z') 
            [71.037114, 114.10272, 103.009185+57.0, 115.026943, 129.042593, 147.068414, 57.021464, 137.058912, 113.084064, 113.16472, 128.094963, 113.084064, 131.040485, 114.042927, 255.15829, 97.052764, 128.058578, 156.101111, 87.032028, 101.047679, 168.053, 99.068414, 186.079313, 0, 163.06332, 128.12472]
 -}
