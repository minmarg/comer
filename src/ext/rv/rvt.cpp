/***************************************************************************
 *   Copyright (C) 2010 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ext/psl.h"
#include "rvnorm.h"
#include "rvexp.h"
#include "rvchi2.h"
#include "rvt.h"

// =========================================================================
// RVt: Student's t random variable
//
RVt::RVt( Rng& rng, TRVt met )
:   RVar( rng ),
    met_( met ),
    df_( 1.0 ),
    mean_( 0.0 ),
    scl_( 1.0 )
{
}

RVt::~RVt()
{
}

// -------------------------------------------------------------------------
// t-distributed random variable;
//
//  pdf: Gamma((df+1)/2) / (Gamma(df/2)sqrt(pi df)scale) 
//          (1+((x-mean)/scale)^2/df)^(-(df+1)/2),  (df,scale>0)
//
int RVt::GenBasic( double* rv )
{
    if( !rv )
        return PSL_ERR_ADDRESS;

    double  df = GetDF();
    double  nn, c2;
    int err;

    if( df <= 0.0 )
        return PSL_ERR_INVALID;

    RVNorm  rnv( GetRng(), RVNorm::TRVN_Ziggurat_M00 );
    RVChi2  rc2( GetRng());

    rc2.SetDF( df );

    if(( err = rnv.Gen( &nn )) != PSL_OK )
        return err;
    if(( err = rc2.Gen( &c2 )) != PSL_OK )
        return err;

    if( c2 < 0.0 )
        return PSL_ERR_ILLEGAL;
    if( c2 == 0.0 )
        *rv = SLC_DBL_MAX;
    else
        *rv = nn * sqrt( df / c2 );
    return PSL_OK;
}

// -------------------------------------------------------------------------
// GenM80: generate Student's t random variate by the algorithm of
//  Marsaglia. (1980) Math. Comp. 34(149), 235-6.
//
int RVt::GenM80( double* rv )
{
    if( !rv )
        return PSL_ERR_ADDRESS;

    const int   maxn = 1000;
    double  df = GetDF();
    double  a, b, c;
    int n, err;

    if( df <= 2.0 )
        return PSL_ERR_INVALID;

    RVNorm  rnv( GetRng(), RVNorm::TRVN_Ziggurat_M00 );
    RVExp   rex( GetRng());

    rex.SetScale( 1.0 /( 0.5 * df - 1.0 ));

    for( n = 0; n < maxn; n++ ) {
        if(( err = rnv.Gen( &a )) != PSL_OK )
            return err;
        if(( err = rex.Gen( &c )) != PSL_OK )
            return err;
        b = a * a /( df - 2.0 );
        if( b < 1.0 && exp( -b - c ) < 1.0 - b )
            break;
    }
    if( maxn <= n )
        return PSL_MAXITERATS;
    *rv = a / sqrt(( 1.0 - 2.0 / df )*( 1.0 - b ));
    return PSL_OK;
}
