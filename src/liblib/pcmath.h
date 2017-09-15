/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __pcmath_h__
#define __pcmath_h__


#define LN2 0.6931471805599453094

#ifndef SIZE_MAX
    #define SIZE_MAX    4294967295UL
#endif

#ifndef SHRT_MIN
    #define SHRT_MIN    -32768
#endif

#define SCORE_MIN ( SHRT_MIN )
#define LOG_PROB_MIN ( SHRT_MIN )

// avoid redefinition
#define PCMIN( a, b )   (( a ) < ( b ) ? ( a ) : ( b ))
#define TIMES2( arg )   (( arg )+( arg ))
#define TIMES3( arg )   ( TIMES2( arg )+( arg ))
#define TIMES4( arg )   ( TIMES3( arg )+( arg ))
#define TIMES5( arg )   ( TIMES4( arg )+( arg ))
#define SQUARE( arg )   (( arg )*( arg ))
#define CUBE( arg )     (( arg )*( arg )*( arg ))

// typedefs
typedef void ( *PCMathTFdF )( double, double*, double*, void* );

// interface
int r2int( register double );           //rounding to the nearest integer
int PC_Gcd( int, int );                 //the greatest common divisor

double LecuyerRand( long* );            //random number generator of L'Ecuyer
const char* root_by_NR_and_bisection(   //root finding by Newton-Raphson and bisection methods
        PCMathTFdF fdfunction,
        double x1,
        double x2,
        double xacc,
        int maxit,
        void* params,
        double* result );

#endif//__pcmath_h__
