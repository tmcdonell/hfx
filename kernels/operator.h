/*
 * Module    : Prelude
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 */

#ifndef __OPERATOR_H__
#define __OPERATOR_H__

/*
 * Template class for binary operators. Certain algorithms may require the
 * operator to be associative (that is, Ta == Tb), such as parallel scan and
 * reduction.
 *
 * As this is template code, it should compile down to something efficient...
 */
template <typename Ta, typename Tb=Ta, typename Tc=Ta>
class BinaryOp
{
public:
    /*
     * Apply the operation to the given operands.
     */
    static __device__ Tc apply(const Ta &a, const Tb &b);

    /*
     * Return an identity element for the type Tc.
     *
     * This may have special meaning for a given implementation, for example a
     * `max' operation over integers may want to return INT_MIN.
     */
    static __device__ Tc identity();
};



/*
 * Basic binary arithmetic operations. We take advantage of automatic type
 * promotion to keep the parameters general.
 */
#define BASIC_OP(name,expr,id)                                                 \
    template <typename Ta, typename Tb=Ta, typename Tc=Ta>                     \
    class name : BinaryOp<Ta, Tb, Tc>                                          \
    {                                                                          \
    public:                                                                    \
        static __device__ Tc apply(const Ta &a, const Tb &b) { return expr; }  \
        static __device__ Tc identity() { return id; }                         \
    };

BASIC_OP(Plus,  a+b,      0)
BASIC_OP(Minus, a-b,      0)
BASIC_OP(Max,   max(a,b), INT_MIN)
BASIC_OP(Min,   min(a,b), INT_MAX)

#endif

