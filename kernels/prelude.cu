/*
 * Module    : Prelude
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 */


#include "scan.cuh"
#include "reduce.cuh"
#include "zipWith.cuh"

#include "utils.h"
#include "kernels.h"

// -----------------------------------------------------------------------------
// Map
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// Reduce
// -----------------------------------------------------------------------------
float reducePlusf(float *xs, int N)
{
    float result = reduce< Plus<float> >(xs, N);
    return result;
}


// -----------------------------------------------------------------------------
// Scan
// -----------------------------------------------------------------------------
void scanl1Plusi(int *in, int *out, int N)
{
    scan< Plus<int>, int, false, false >(in, out, N);
}

void scanr1Plusi(int *in, int *out, int N)
{
    scan< Plus<int>, int, true, false >(in, out, N);
}


// -----------------------------------------------------------------------------
// Segmented Scan
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// zipWith
// -----------------------------------------------------------------------------
void zipWithPlusif(int *xs, float *ys, float *zs, int N)
{
    zipWith< Plus<int, float, float> >(xs, ys, zs, N);
}

void zipWithTimesif(int *xs, float *ys, float *zs, int N)
{
    zipWith< Times<int, float, float> >(xs, ys, zs, N);
}

