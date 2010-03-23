/* -----------------------------------------------------------------------------
 *
 * Module    : Scan
 * Copyright : (c) 2010 Trevor L. McDonell
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#include "algorithms.h"
#include "cudpp/utils.h"

#include <cudpp.h>


/*
 * Apply a binary operator to an array similar to `fold', but return a
 * successive list of values reduced from the end. The reduction will take place
 * in parallel, so the operator must be associative.
 */
template <CUDPPOperator op, typename T, bool backward, bool exclusive>
void
scan
(
    const T             *d_in,
    T                   *d_out,
    const unsigned int  length
)
{

    CUDPPHandle         plan;
    CUDPPConfiguration  cp;

    cp.algorithm = CUDPP_SCAN;
    cp.datatype  = getType<T>();
    cp.op        = op;
    cp.options   = (backward  ? CUDPP_OPTION_BACKWARD  : CUDPP_OPTION_FORWARD)
                 | (exclusive ? CUDPP_OPTION_EXCLUSIVE : CUDPP_OPTION_INCLUSIVE);

    cudppPlan(&plan, cp, length, 1, 0);
    cudppScan(plan, d_out, d_in, length);

    cudppDestroyPlan(plan);
}


// -----------------------------------------------------------------------------
// Instances
// -----------------------------------------------------------------------------

void prescanl_plusui(const unsigned int *d_in, unsigned int *d_out, const unsigned int N)
{
    scan< CUDPP_ADD, unsigned int, false, true >(d_in, d_out, N);
}

void prescanr_plusui(const unsigned int *d_in, unsigned int *d_out, const unsigned int N)
{
    scan< CUDPP_ADD, unsigned int, true, true >(d_in, d_out, N);
}

