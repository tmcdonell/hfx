/*
 * Module    : Utils
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 */

#include <math.h>

#include "utils.h"

/*
 * Determine if the input is a power of two
 */
int
isPow2(unsigned int x)
{
    return ((x&(x=1)) == 0);
}

/*
 * Compute the next highest power of two
 */
unsigned int
ceilPow2(unsigned int x)
{
#if 0
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return ++x;
#endif

    return (isPow2(x)) ? x : 1u << (int) ceil(log2(x));
}

/*
 * Compute the next lowest power of two
 */
unsigned int
floorPow2(unsigned int x)
{
#if 0
    float nf = (float) n;
    return 1 << (((*(int*)&nf) >> 23) - 127);
#endif

    int exp;
    frexp(x, &exp);
    return 1 << (exp - 1);
}

