/* -----------------------------------------------------------------------------
 *
 * Module    : Utils
 * Copyright : (c) [2009..2010] Trevor L. McDonell
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/


#ifndef __UTILS_H__
#define __UTILS_H__

#include <math.h>
#include <stdio.h>


/*
 * Core assert function. Don't let this escape...
 */
#if defined(__CUDACC__) || !defined(__DEVICE_EMULATION__)
#define __assert(e, file, line) ((void)0)
#else
#define __assert(e, file, line) \
    ((void) fprintf (stderr, "%s:%u: failed assertion `%s'\n", file, line, e), abort())
#endif

/*
 * Test the given expression, and abort the program if it evaluates to false.
 * Only available in debug mode.
 */
#ifndef _DEBUG
#define assert(e)               ((void)0)
#else
#define assert(e)  \
    ((void) ((e) ? (void(0)) : __assert (#e, __FILE__, __LINE__)))
#endif


/*
 * Macro to insert __syncthreads() in device emulation mode
 *
 * In emulation mode, this is required inside warp-synchronous code, to recreate
 * the behaviour of the warp threads executing in lock-step on the GPU.
 */
#ifdef __DEVICE_EMULATION__
#define __EMUSYNC               __syncthreads()
#else
#define __EMUSYNC
#endif


/*
 * Check the return status of CUDA API calls, and abort with an appropriate
 * error string on failure.
 */
#define CUDA_SAFE_CALL_NO_SYNC(call)                                           \
    do {                                                                       \
        cudaError err = call;                                                  \
        if(cudaSuccess != err) {                                               \
            const char *str = cudaGetErrorString(err);                         \
            __assert(str, __FILE__, __LINE__);                                 \
        }                                                                      \
    } while (0)

#define CUDA_SAFE_CALL(call)                                                   \
    do {                                                                       \
        CUDA_SAFE_CALL_NO_SYNC(call);                                          \
        CUDA_SAFE_CALL_NO_SYNC(cudaThreadSynchronize());                       \
    } while (0)


#ifdef __cplusplus
extern "C" {
#endif

/*
 * Determine if the input is a power of two
 */
inline bool
isPow2(unsigned int x)
{
    return ((x&(x-1)) == 0);
}


/*
 * Compute the next highest power of two
 */
inline unsigned int
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

    return (isPow2(x)) ? x : 1u << (int) ceil(log2((double)x));
}


/*
 * Compute the next lowest power of two
 */
inline unsigned int
floorPow2(unsigned int x)
{
#if 0
    float nf = (float) n;
    return 1 << (((*(int*)&nf) >> 23) - 127);
#endif

    int exp;
    frexp((double)x, &exp);
    return 1 << (exp - 1);
}

#undef __asert

#ifdef __cplusplus
}
#endif
#endif

