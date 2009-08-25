/*
 * Module    : Sequest
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 */

#include "utils.h"
#include "kernels.h"

#include "sequestXC_kernel.h"


/*
 * Score the peptide against the observed intensity spectrum, which corresponds
 * to a dot product between the theoretical representation and pre-processed
 * experimental spectra.
 */
float
sequestXC
(
    float       *spec_exp,
    int         *spec_thry,
    int         length
)
{
    int         maxThreads = 128;
    int         maxBlocks  = 64;

    int         threads = (length < maxBlocks*2) ? nextPow2((n+1)/2) : maxThreads;
    int         blocks  = (length + (threads * 2 - 1)) / (threads * 2);
    int         smem    = threads * sizeof(float);

    assert(isPow2(length));

    /*
     * Calculate partial sums for each block on the GPU
     */

    /*
     * Really want to call this recursively, with our partial sums
     */
    switch (threads)
    {
    case 512: xcorr_dotP<512, true><<< blocks, threads, smem >>>(spec_exp, spec_thry, o_data, length); break;
    case 256: xcorr_dotP<256, true><<< blocks, threads, smem >>>(spec_exp, spec_thry, o_data, length); break;
    case 128: xcorr_dotP<128, true><<< blocks, threads, smem >>>(spec_exp, spec_thry, o_data, length); break;
    case  64: xcorr_dotP< 64, true><<< blocks, threads, smem >>>(spec_exp, spec_thry, o_data, length); break;
    case  32: xcorr_dotP< 32, true><<< blocks, threads, smem >>>(spec_exp, spec_thry, o_data, length); break;
    case  16: xcorr_dotP< 16, true><<< blocks, threads, smem >>>(spec_exp, spec_thry, o_data, length); break;
    case   8: xcorr_dotP<  8, true><<< blocks, threads, smem >>>(spec_exp, spec_thry, o_data, length); break;
    case   4: xcorr_dotP<  4, true><<< blocks, threads, smem >>>(spec_exp, spec_thry, o_data, length); break;
    case   2: xcorr_dotP<  2, true><<< blocks, threads, smem >>>(spec_exp, spec_thry, o_data, length); break;
    case   1: xcorr_dotP<  1, true><<< blocks, threads, smem >>>(spec_exp, spec_thry, o_data, length); break;
    }

}

