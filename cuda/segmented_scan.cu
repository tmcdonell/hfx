/*
 * Module    : Prelude
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 */

#include "utils.h"
#include "kernels.h"
#include "operator.h"

#include "cudpp/cudpp_globals.h"
#include "cudpp/segmented_scan_kernel.cu"
#include "cudpp/vector_kernel.cu"

template <typename T>
struct segscan_plan
{
    T                   **sums;
    unsigned int        **flags;
    unsigned int        **indices;
    unsigned int        num_levels;
};

static inline unsigned int
calc_num_blocks(unsigned int N)
{
    return max(1u, (unsigned int)ceil((double)N / (SEGSCAN_ELTS_PER_THREAD * CTA_SIZE)));
}

/*
 * Scans of large arrays must be split (possibly recursively) into a hierarchy
 * of block scans, where each block is processed by a single thread block. On
 * returning from a recursive call, the total sum of each block from the level
 * below is added to all elements of the first segment of the corresponding
 * block. This is the CPU-side workhorse that achieves this.
 */
template <class op, typename T, bool backward, bool exclusive, bool shift_flags>
static void
segscan_recursive
(
    const T             *in,
    const unsigned int  *flags,
    T                   *out,
    segscan_plan<T>     *plan,
    unsigned int        N,
    int                 level
)
{
    unsigned int        num_blocks = calc_num_blocks(N);
    unsigned int        per_block  = CTA_SIZE * 2;
    bool                is_full    = N == num_blocks * SEGSCAN_ELTS_PER_THREAD * CTA_SIZE;

    /*
     * Space to store flags in the shared memory. Two sets are required, one
     * gets modified and the other does not.
     */
    unsigned int        flag_space = per_block * sizeof(unsigned int);
    unsigned int        idx_space  = per_block * sizeof(unsigned int);

    dim3                grid(num_blocks, 1, 1);
    dim3                threads(CTA_SIZE, 1, 1);
    unsigned int        smem = sizeof(T) * per_block + flag_space + idx_space;

    /*
     * Check the hardware
     */
    int dev;
    cudaDeviceProp props;
    cudaGetDevice(&dev);
    cudaGetDeviceProperties(&props, dev);

    /*
     * Set up execution parameters, and execute the scan
     */
#define MULTIBLOCK      0x01
#define FULLBLOCK       0x02
#define SM12_HW         0x04
#define BACKWARD        0x08
    int traits = 0;
    if (num_blocks > 1)   traits |= MULTIBLOCK;
    if (is_full)          traits |= FULLBLOCK;
    if (props.minor >= 2) traits |= SM12_HW;

    switch (traits)
    {
    case 0:
        segmentedScan4
            < T, SegmentedScanTraits<T, op, backward, exclusive, shift_flags, false, false, false> >
            <<<grid, threads, smem>>>(out, in, flags, N);
        break;

    case MULTIBLOCK:
        segmentedScan4
            < T, SegmentedScanTraits<T, op, backward, exclusive, shift_flags, false, true, false> >
            <<<grid, threads, smem>>>(out, in, flags, N, plan->sums[level], plan->flags[level], plan->indices[level]);
        break;

    case FULLBLOCK:
        segmentedScan4
            < T, SegmentedScanTraits<T, op, backward, exclusive, shift_flags, true, false, false> >
            <<<grid, threads, smem>>>(out, in, flags, N);
        break;

    case SM12_HW:
        segmentedScan4
            < T, SegmentedScanTraits<T, op, backward, exclusive, shift_flags, false, false, true> >
            <<<grid, threads, smem>>>(out, in, flags, N);
        break;

    case MULTIBLOCK | FULLBLOCK:
        segmentedScan4
            < T, SegmentedScanTraits<T, op, backward, exclusive, shift_flags, true, true, false> >
            <<<grid, threads, smem>>>(out, in, flags, N, plan->sums[level], plan->flags[level], plan->indices[level]);
        break;

    case MULTIBLOCK | SM12_HW:
        segmentedScan4
            < T, SegmentedScanTraits<T, op, backward, exclusive, shift_flags, false, true, true> >
            <<<grid, threads, smem>>>(out, in, flags, N, plan->sums[level], plan->flags[level], plan->indices[level]);
        break;

    case FULLBLOCK | SM12_HW:
        segmentedScan4
            < T, SegmentedScanTraits<T, op, backward, exclusive, shift_flags, true, false, true> >
            <<<grid, threads, smem>>>(out, in, flags, N);
        break;

    case MULTIBLOCK | FULLBLOCK | SM12_HW:
        segmentedScan4
            < T, SegmentedScanTraits<T, op, backward, exclusive, shift_flags, true, true, true> >
            <<<grid, threads, smem>>>(out, in, flags, N, plan->sums[level], plan->flags[level], plan->indices[level]);
        break;

    default:
        assert(!"Non-exhaustive patterns in match");
    }

    /*
     * After scanning the sub-blocks, take all of the last values and
     * segment-scan those. This will give the new value which must be added to
     * the first segment of each block to get the final result.
     */
    if (num_blocks > 1)
    {
        T            *sums    = plan->sums[level];
        unsigned int *indices = plan->indices[level];

        segscan_recursive
            <op, T, backward, false, false>
            (sums, plan->flags[level], sums, plan, num_blocks, level + 1);

        traits = 0;
        if (is_full)  traits |= FULLBLOCK;
        if (backward) traits |= BACKWARD;

        switch (traits)
        {
        case 0:
            vectorSegmentedAddUniform4<T, op, false><<<grid, threads>>>(out, sums, indices, N, 0, 0);
            break;

        case FULLBLOCK:
            vectorSegmentedAddUniform4<T, op, true><<<grid, threads>>>(out, sums, indices, N, 0, 0);
            break;

        case BACKWARD:
            vectorSegmentedAddUniformToRight4<T, op, false><<<grid, threads>>>(out, sums, indices, N, 0, 0);
            break;

        case FULLBLOCK | BACKWARD:
            vectorSegmentedAddUniformToRight4<T, op, true><<<grid, threads>>>(out, sums, indices, N, 0, 0);
            break;

        default:
            assert(!"Non-exhaustive patterns in match");
        }
    }

#undef MULTIBLOCK
#undef FULLBLOCK
#undef SM12_HW
#undef BACKWARD
}


/*
 * Allocate temporary memory used by the segmented scan
 */
template <typename T>
static void
segscan_init(int N, segscan_plan<T> *plan)
{
    unsigned int level        = 0;
    unsigned int elements     = N;
    unsigned int num_blocks;

    /*
     * Determine how many intermediate block-level summations will be required
     */
    for (elements = N; elements > 1; elements = num_blocks)
    {
        num_blocks = calc_num_blocks(elements);

        if (num_blocks > 1)
            ++level;
    }

    plan->num_levels = level;
    plan->sums       = (T**) malloc(level * sizeof(T*));
    plan->flags      = (unsigned int**) malloc(level * sizeof(unsigned int *));
    plan->indices    = (unsigned int**) malloc(level * sizeof(unsigned int *));

    /*
     * Now allocate the necessary storage at each level
     */
    for (elements = N, level = 0; elements > 1; elements = num_blocks, level++)
    {
        num_blocks = calc_num_blocks(elements);

        if (num_blocks > 1)
        {
            cudaMalloc((void**) &plan->sums[level],    num_blocks * sizeof(T));
            cudaMalloc((void**) &plan->flags[level],   num_blocks * sizeof(unsigned int));
            cudaMalloc((void**) &plan->indices[level], num_blocks * sizeof(unsigned int));
        }
    }
}


/*
 * Clean up temporary memory
 */
template <typename T>
static void
segscan_finalise(segscan_plan<T> *p)
{
    for (unsigned int l = 0; l < p->num_levels; ++l)
    {
        cudaFree(p->sums[l]);
        cudaFree(p->flags[l]);
        cudaFree(p->indices[l]);
    }

    free(p->sums);
    free(p->flags);
    free(p->indices);
}


/*
 * Perform a segmented scan operation on the input array of data, much like
 * `scan', but with an additional input array of non-zero `flags' that demarcate
 * the first element of a segment.
 */
template <class op, typename T, bool backward, bool exclusive>
void
segmented_scan
(
    const T             *in,
    const unsigned int  *flags,
    T                   *out,
    unsigned int        length
)
{
    segscan_plan<T> plan;
    segscan_init<T>(length, &plan);

    segscan_recursive<op, T, backward, exclusive, backward>(in, flags, out, &plan, length, 0);

    segscan_finalise<T>(&plan);
}


// -----------------------------------------------------------------------------
// Instances
// -----------------------------------------------------------------------------

void scanl1Seg_plusf(float *in, unsigned int *flags, float *out, int N)
{
    segmented_scan< Plus<float>, float, false, false >(in, flags, out, N);
}

void scanr1Seg_plusf(float *in, unsigned int *flags, float *out, int N)
{
    segmented_scan< Plus<float>, float, true, false >(in, flags, out, N);
}

void scanl1Seg_plusui(unsigned int *in, unsigned int *flags, unsigned int *out, int N)
{
    segmented_scan< Plus<unsigned int>, unsigned int, false, false >(in, flags, out, N);
}

