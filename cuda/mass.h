/*
 * Module    : Mass
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 *
 * Mass constants for elements, molecules and amino acid groups
 *
 * TODO: Find a way to derive both the C and HS versions from a single source
 */

#ifndef __MASS_H__
#define __MASS_H__

#include "operator.h"

/*
 * Spectrum bin width
 */
#define binWidthMono 1.0005079f
#define binWidthAvg  1.0011413f


/*
 * The monoisotopic mass of several elements and molecules
 */
#define massH2O 18.01056f
#define massNH3 17.02655f
#define massCO  27.9949f
#define massO   16.0013f
#define massH   1.0078246f

extern __device__ float aa_table[26] =
{
    71.037114f,         //  A: Alanine              Alg     C3H5NO
    114.53494f,         //  B: Aspargine            Asx     C4H8N2O3
    103.009185f,        //  C: Cysteine             Cys     C3H5NOS
    115.026943f,        //  D: Aspartic acid        Asp     C4H5NO3
    129.042593f,        //  E: Glutamic acid        Glu     C5H7NO3
    147.068414f,        //  F: Phenylalanine        Phe     C9H9NO
    57.021464f,         //  G: Glycine              Gly     C2H3NO
    137.058912f,        //  H: Histidine            His     C6H11NO
    113.084064f,        //  I: Isoleucine           Ile     C6H11NO
    113.16472f,         //  J: Leucine              Xle     C6H13NO2
    128.094963f,        //  K: Lysine               Lys     C6H12N2O
    113.084064f,        //  L: Leucine              Leu     C6H11NO
    131.040485f,        //  M: Methionine           Met     C5H9NOS
    114.042927f,        //  N: Asparagine           Asn     C4H6N2O2
    114.07931f,         //  O: Pyrrolysine          Pyl     C12H21N3O3
    97.052764f,         //  P: Proline              Pro     C5H7NO
    128.058578f,        //  Q: Glutamine            Gln     C5H8N2O2
    156.101111f,        //  R: Arginine             Arg     C6H12N4O
    87.032028f,         //  S: Serine               Ser     C3H5NO2
    101.047679f,        //  T: Threonine            Thr     C4H7NO2
    150.04344f,         //  U: Selenocysteine       Sec     C3H5NOSe
    99.068414f,         //  V: Valine               Val     C5H9NO
    186.079313f,        //  W: Tryptophan           Trp     C11H10N2O
    113.08406f,         //  X: Unknown              Xaa
    163.06332f,         //  Y: Tyrosine             Tyr     C9H9NO2
    128.55059f          //  Z: Glutamine            Glx     C5H10N2O3
};


template <typename Ta, typename Tb>
class getAAMass : Functor<Ta, Tb>
{
public:
    static __device__ Tb apply(const Ta &x);
};


template <>
class getAAMass<int, float> : Functor<char, float>
{
public:
    static __device__ float apply(const char &x) { return aa_table[x - 'A']; }
};


#endif

