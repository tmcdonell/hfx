/* -----------------------------------------------------------------------------
 *
 * Module    : Filter
 * Copyright : (c) [2009..2011] Trevor L. McDonell
 * License   : BSD
 *
 * ---------------------------------------------------------------------------*/

#include <thrust/iterator/counting_iterator.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <stdint.h>

#include "algorithms.h"


template <typename T>
struct interval : public thrust::unary_function<T,bool>
{
    T min_val;
    T max_val;

    __host__ __device__
    interval(T _m, T _n) : min_val(_m), max_val(_n) {}

    __host__ __device__ bool operator() (T x)
    {
        return (min_val <= x && x <= max_val);
    }
};

uint32_t
findIndicesInRange_f
(
    const float         *d_in_raw,
    uint32_t            *d_indices_raw,
    const uint32_t      N,
    const float         min_val,
    const float         max_val
)
{
    thrust::device_ptr<const float>     d_in(d_in_raw);
    thrust::device_ptr<uint32_t>        d_indices(d_indices_raw);

    // define the sequence [0, N)
    thrust::counting_iterator<uint32_t> first(0);
    thrust::counting_iterator<uint32_t> last = first + N;

    // compute indices of elements in range
    thrust::device_ptr<uint32_t> indices_end =
        thrust::copy_if(first, last, d_in, d_indices, interval<const float>(min_val, max_val));

    return indices_end - d_indices;
}

