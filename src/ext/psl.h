/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __psl__
#define __psl__

#include <limits.h>


#define SLC_E           2.71828182845904523536028747135

//pi, sqrt(pi), ln(pi), log(2*pi), log(sqrt(pi)), and ln(sqrt(2*pi))
#define SLC_PI          3.14159265358979323846264338328
#define SLC_SQRTPI      1.77245385090551602729816748334
#define SLC_LNPI        1.14472988584940017414342735135
#define SLC_LN2PI       1.8378770664093454835606594728111235279723
#define SLC_LNSQRTPI    0.57236494292470008706
#define SLC_LNSQRT2PI   0.9189385332046727418

// Euler constant
#define SLC_EULER       0.57721566490153286060651209008

//sqrt(2), and ln(2)
#define SLC_SQRT2       1.41421356237309504880168872421
#define SLC_LN2         0.69314718055994530941723212146
#define SLC_LN1K        6.90775527898213681510242167860

// machine dependent constants
#define SLC_DBL_EPSILON         ( 2.2204460492503131e-16 )
#define SLC_SQRT_DBL_EPSILON    ( 1.4901161193847656e-08 )
#define SLC_ROOT3_DBL_EPSILON   ( 6.0554544523933429e-06 )
#define SLC_ROOT4_DBL_EPSILON   ( 1.2207031250000000e-04 )
#define SLC_ROOT5_DBL_EPSILON   ( 7.4009597974140505e-04 )
#define SLC_ROOT6_DBL_EPSILON   ( 2.4607833005759251e-03 )
#define SLC_LOG_DBL_EPSILON     ( -3.6043653389117154e+01 )

#define SLC_DBL_MIN             ( 2.2250738585072014e-308 )
#define SLC_SQRT_DBL_MIN        ( 1.4916681462400413e-154 )
#define SLC_LOG_DBL_MIN         (-7.0839641853226408e+02)

#define SLC_DBL_MAX             ( 1.7976931348623157e+308 )
#define SLC_SQRT_DBL_MAX        ( 1.3407807929942596e+154 )
#define SLC_LOG_DBL_MAX         ( 7.0978271289338397e+02 )
//

// functions
#define SLC_POW2( X )   (( X )*( X ))
#define SLC_POW3( X )   (( X )*( X )*( X ))
#define SLC_POW4( X )   ( SLC_POW2( X ) * SLC_POW2( X ))
#define SLC_POW5( X )   ( SLC_POW3( X ) * SLC_POW2( X ))
#define SLC_POW6( X )   ( SLC_POW3( X ) * SLC_POW3( X ))

#define SLC_ODD( X )    (( X ) & 1 )
#define SLC_EVEN( X )   (!SLC_ODD( X ))
#define SLC_SIGN( X )   (( 0 <= ( X ))? 1: -1 )

#define SLC_MIN( X, Y ) ((( X ) < ( Y ))? ( X ): ( Y ))
#define SLC_MAX( X, Y ) ((( X ) < ( Y ))? ( Y ): ( X ))
#define SLC_MAX3( X1, X2, X3 )  SLC_MAX( X1, SLC_MAX( X2, X3 ))
#define SLC_MAX4( X1, X2, X3, X4 )  SLC_MAX( X1, SLC_MAX3( X2, X3, X4 ))
#define SLC_MAX5( X1, X2, X3, X4, X5 )  SLC_MAX( X1, SLC_MAX4( X2, X3, X4, X5 ))

//
#endif//__psl__
