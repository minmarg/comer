/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include "ext/psl.h"
#include "ext/pslcodes.h"
#include "ext/rng.h"
#include "rvexp.h"


// =========================================================================
// RVExp
//
RVExp::RVExp( Rng& rng )
:   RVar( rng ),
    scale_( 1.0 )
{
}

RVExp::~RVExp()
{
}

// -------------------------------------------------------------------------
// exponentially distributed random variable
//  pdf: 1/scale * exp(-x/scale)
//
int RVExp::rvgexp( double* result )
{
    if( !result )
        return PSL_ERR_ADDRESS;

    double scale = GetScale();

    if( scale < 0.0 )
        return PSL_ERR_DOMAIN;

    if( scale == 0.0 ) {
        *result = 0.0;
        return 0;
    }

    int err = rvgexprand( result );
    if( err )
        return err;
    *result *= scale;
    return 0;
}

// -------------------------------------------------------------------------
// rvgexprand: random variable from standard exponential distribution.
// Code by R.Ihaka(C)1998 under GNU General Public License
// Ref. Ahrens, Dieter. (1972) Computer methods for sampling from the 
// exponential and normal distributions. Comm. ACM 15, 873-82.
// (see also Devroye. Non Uniform Random Variate Generation. Springer, 1986)
//
int RVExp::rvgexprand( double* result )
{
    // q[k-1] = sum(log(2)^k / k!)  k=1,..,n,
    // The highest n (here 16) is determined by q[n-1] = 1.0
    // within standard precision
    const static double q[] = {
        0.6931471805599453,
        0.9333736875190459,
        0.9888777961838675,
        0.9984959252914960,
        0.9998292811061389,
        0.9999833164100727,
        0.9999985691438767,
        0.9999998906925558,
        0.9999999924734159,
        0.9999999995283275,
        0.9999999999728814,
        0.9999999999985598,
        0.9999999999999289,
        0.9999999999999968,
        0.9999999999999999,
        1.0000000000000000
    };

    const int   imax = 100;
    int     i;
    double  a = 0.;
    double  u = GetRng().GetDouble();    //precaution if u = 0 is ever returned

    if( !result )
        return PSL_ERR_ADDRESS;
    for( i = 0; i < imax &&( u <= 0. || 1 <= u ); ) u = GetRng().GetDouble();
    if( u <= 0. || 1 <= u )
        return PSL_MAXITERATS;

    for( ;; ) {
        u += u;
        if( 1. < u )
            break;
        a += q[0];
    }
    u -= 1.;

    if( u <= q[0]) {
        *result = a + u;
        return 0;
    }

    i = 0;
    double ustar = GetRng().GetDouble(), umin = ustar;

    do{ ustar = GetRng().GetDouble0();
        if( umin > ustar )
            umin = ustar;
        i++;
    } while( q[i] < u );

    *result = a + umin * q[0];
    return 0;
}
