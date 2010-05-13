--------------------------------------------------------------------------------
-- |
-- Module    : Mass
-- Copyright : (c) [2009..2010] Trevor L. McDonell
-- License   : BSD
--
-- Mass constants for elements, molecules and amino acid groups
--
-- TODO: Find a way to derive both the C and HS versions from a single source
--
--------------------------------------------------------------------------------

module Mass where


--
-- The monoisotopic mass of several elements and molecules
--
massH2O, massNH3, massCO, massO, massH :: Float
massH2O = 18.010565
massNH3 = 17.026549
massCO  = 27.994915
massO   = 15.994915
massH   = 1.007276


--
-- Return the monoisotopic (ground state) mass of an amino acid residue in
-- atomic mass units, for the given short abbreviation.
--
-- This corresponds to the mass of the residue, where the mass of the
-- peptide/protein is the mass of the residue plus the mass of water.
--
getAAMassMono     :: Char -> Float
getAAMassMono 'A' =   71.037114         -- Alanine              Alg     C3H5NO
getAAMassMono 'R' =  156.101111         -- Arginine             Arg     C6H12N4O
getAAMassMono 'N' =  114.042927         -- Asparagine           Asn     C4H6N2O2
getAAMassMono 'D' =  115.026943         -- Aspartic acid        Asp     C4H5NO3
getAAMassMono 'C' =  103.009185         -- Cysteine             Cys     C3H5NOS
getAAMassMono 'E' =  129.042593         -- Glutamic acid        Glu     C5H7NO3
getAAMassMono 'Q' =  128.058578         -- Glutamine            Gln     C5H8N2O2
getAAMassMono 'G' =   57.021464         -- Glycine              Gly     C2H3NO
getAAMassMono 'H' =  137.058912         -- Histidine            His     C6H11NO
getAAMassMono 'I' =  113.084064         -- Isoleucine           Ile     C6H11NO
getAAMassMono 'L' =  113.084064         -- Leucine              Leu     C6H11NO
getAAMassMono 'K' =  128.094963         -- Lysine               Lys     C6H12N2O
getAAMassMono 'M' =  131.040485         -- Methionine           Met     C5H9NOS
getAAMassMono 'F' =  147.068414         -- Phenylalanine        Phe     C9H9NO
getAAMassMono 'P' =   97.052764         -- Proline              Pro     C5H7NO
getAAMassMono 'O' =  114.07931          -- Pyrrolysine          Pyl     C12H21N3O3
getAAMassMono 'U' =  150.04344          -- Selenocysteine       Sec     C3H5NOSe
getAAMassMono 'S' =   87.032028         -- Serine               Ser     C3H5NO2
getAAMassMono 'T' =  101.047679         -- Threonine            Thr     C4H7NO2
getAAMassMono 'W' =  186.079313         -- Tryptophan           Trp     C11H10N2O
getAAMassMono 'Y' =  163.06332          -- Tyrosine             Tyr     C9H9NO2
getAAMassMono 'V' =   99.068414         -- Valine               Val     C5H9NO

--
-- Ambiguous amino acids
--
getAAMassMono 'B' =  114.53494          -- Aspargine            Asx     C4H8N2O3
getAAMassMono 'J' =  113.16472          -- Leucine              Xle     C6H13NO2
getAAMassMono 'Z' =  128.55059          -- Glutamine            Glx     C5H10N2O3
getAAMassMono 'X' =  113.08406          -- Unknown              Xaa

getAAMassMono  x  = error $ "Unknown amino acid abbreviation: " ++ [x]


--
-- Return the average (molecular) mass of an amino acid residue in atomic mass
-- units, for the given short abbreviation.
--
getAAMassAvg     :: Char -> Float
getAAMassAvg 'A' =   71.0788            -- Alanine              Alg
getAAMassAvg 'R' =  156.1875            -- Arginine             Arg
getAAMassAvg 'N' =  114.1038            -- Asparagine           Asn
getAAMassAvg 'D' =  115.0886            -- Aspartic acid        Asp
getAAMassAvg 'C' =  103.1388            -- Cysteine             Cys
getAAMassAvg 'E' =  129.1155            -- Glutamic acid        Glu
getAAMassAvg 'Q' =  128.1307            -- Glutamine            Gln
getAAMassAvg 'G' =   57.0519            -- Glycine              Gly
getAAMassAvg 'H' =  137.1411            -- Histidine            His
getAAMassAvg 'I' =  113.1594            -- Isoleucine           Ile
getAAMassAvg 'L' =  113.1594            -- Leucine              Leu
getAAMassAvg 'K' =  128.1741            -- Lysine               Lys
getAAMassAvg 'M' =  131.1926            -- Methionine           Met
getAAMassAvg 'F' =  147.1766            -- Phenylalanine        Phe
getAAMassAvg 'P' =   97.1167            -- Proline              Pro
getAAMassAvg 'O' =  114.1472            -- Pyrrolysine          Pyl
getAAMassAvg 'U' =  150.0388            -- Selenocysteine       Sec
getAAMassAvg 'S' =   87.0782            -- Serine               Ser
getAAMassAvg 'T' =  101.1051            -- Threonine            Thr
getAAMassAvg 'W' =  186.2132            -- Tryptophan           Trp
getAAMassAvg 'Y' =  163.1760            -- Tyrosine             Tyr
getAAMassAvg 'V' =   99.1326            -- Valine               Val

getAAMassAvg 'B' =  114.5962            -- Aspargine            Asx
getAAMassAvg 'J' =  113.240178          -- Leucine              Xle
getAAMassAvg 'Z' =  128.6231            -- Glutamine            Glx
getAAMassAvg 'X' =  113.1594            -- Unknown              Xaa

getAAMassAvg  x  = error $ "Unknown amino acid abbreviation: " ++ [x]

