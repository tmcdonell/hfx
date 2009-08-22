/*
 * Module    : Utils
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 */


#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdio.h>

#define assert(e)  \
    ((void) ((e) ? (void(0)) : __assert (#e, __FILE__, __LINE__)))
#define __assert(e, file, line) \
    ((void) printf ("%s:%u: failed assertion `%s'\n", file, line, e), abort())


#ifdef __cplusplus
extern "C" {
#endif

unsigned int ceilPow2(unsigned int x);

#ifdef __cplusplus
}
#endif
#endif

