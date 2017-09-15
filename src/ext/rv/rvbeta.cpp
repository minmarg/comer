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
#include "rvbeta.h"


// =========================================================================
// RVBeta
//
RVBeta::RVBeta( Rng& rng, TRVB met, RVGamma::TRVG gmet )
:   RVar( rng ),
    met_( met ),
    rvg_( rng, gmet ),
    a_( 1.0 ),
    b_( 1.0 ),
    //Ch78 state vars
    Cap_( 0.0 ),
    Cbp_( 0.0 ),
    alpha_( 0.0 ), beta_( 0.0 ), gamma_( 0.0 ),
    delta_( 0.0 ), kap1_( 0.0 ), kap2_( 0.0 )
{
    rvg_.SetScale( 1.0 );
}

RVBeta::~RVBeta()
{
}

// -------------------------------------------------------------------------
// Beta distributed random variable;
//
//  pdf: Gamma(a+b) / (Gamma(a)Gamma(b)) x^(a-1) (1-x)^(b-1),
//      (0<=x<=1; a,b>0)
//
int RVBeta::GenBasic( double* result )
{
    if( !result )
        return PSL_ERR_ADDRESS;

    double  a = GetA();
    double  b = GetB();
    double  x, y;
    int err;

    if( a <= 0.0 || b <= 0.0 )
        return PSL_ERR_INVALID;

    rvg_.SetShape( a );
    if(( err = rvg_.Gen( &x )) != PSL_OK )
        return err;
    rvg_.SetShape( b );
    if(( err = rvg_.Gen( &y )) != PSL_OK )
        return err;
    *result = x /( x + y );
    return PSL_OK;
}

// -------------------------------------------------------------------------
// GenCh78: generate beta random variate by the algorithm of
//  Cheng. (1978) Communications of the ACM 21, 317-22.
//
int RVBeta::GenCh78( double* rv )
{
    if( !rv )
        return PSL_ERR_ADDRESS;

    const int   maxn = 100;
    double  a = GetA();
    double  b = GetB();
    double  aa, bb;
    double  r, s, t, u1, u2, v, w, y, z;
    int n, abchanged;

    if( a <= 0.0 || b <= 0.0 )
        return PSL_ERR_INVALID;

    abchanged = Cap_ != a || Cbp_ != b;
    if( abchanged ) {
        Cap_ = a;
        Cbp_ = b;
    }

    aa = SLC_MIN( a, b );
    bb = SLC_MAX( a, b );

    if( aa <= 1.0 ) {
        //BC algorithm for aa<=1.0
        //
        //change vars to correspond paper notation
        r = aa; aa = bb; bb = r;//aa=MAX, bb=MIN
        if( abchanged ) {
            //INITIALIZE
            alpha_ = aa + bb;
            beta_ = 1.0 / bb;
            delta_ = 1.0 + aa - bb;
            kap1_ = delta_ *( 0.0138889 + 0.0416667 * bb )/( aa * beta_ - 0.777778 );
            kap2_ = 0.25 + ( 0.5 + 0.25 / delta_ ) * bb;
        }
        for( n = 0; n < maxn; n++ ) {
            //step 1
            u1 = GetRng().GetDouble01();
            u2 = GetRng().GetDouble0();
            if( u1 < 0.5 ) {
                //step 2    
                y = u1 * u2;
                z = u1 * y;
                if( kap1_ <= 0.25 * u2 + z - y )
                    continue;
            }
            else {
                //step 3
                z = u1 * u1 * u2;
                if( z <= 0.25 ) {
                    v = beta_ * log( u1 /( 1.0 - u1 ));
                    w = log( aa ) + v;
                    if( w <= SLC_LOG_DBL_MIN )
                        w = 0.0;
                    else if( SLC_LOG_DBL_MAX <= w )
                        w = SLC_DBL_MAX;
                    else
                        w = exp( w );
                    break;
                }
                //step 4
                if( kap2_ <= z )
                    continue;
            }
            //step 5
            v = beta_ * log( u1 /( 1.0 - u1 ));
            w = log( aa ) + v;
            if( w <= SLC_LOG_DBL_MIN )
                w = 0.0;
            else if( SLC_LOG_DBL_MAX <= w )
                w = SLC_DBL_MAX;
            else
                w = exp( w );

            if( log( z ) <= alpha_ *( log( alpha_ /( bb + w )) + v ) - 1.3862944 )
                break;
        }
        //step 6
        if( aa == a ) 
            *rv = w /( bb + w );
        else
            *rv = bb /( bb + w );
        if( maxn <= n )
            return PSL_MAXITERATS;
        return PSL_OK;
    }

    //BB algorithm for aa>1.0
    //
    if( abchanged ) {
        //INITIALIZE
        alpha_ = aa + bb;
        beta_ = sqrt(( alpha_ - 2.0 )/( 2.0 * aa * bb - alpha_ ));
        gamma_ = aa + 1.0 / beta_;
    }
    for( n = 0; n < maxn; n++ ) {
        //step 1
        u1 = GetRng().GetDouble01();
        u2 = GetRng().GetDouble0();

        v = beta_ * log( u1 /( 1.0 - u1 ));
        w = log( aa ) + v;
        if( w <= SLC_LOG_DBL_MIN )
            w = 0.0;
        else if( SLC_LOG_DBL_MAX <= w )
            w = SLC_DBL_MAX;
        else
            w = exp( w );

        z = u1 * u1 * u2;
        r = gamma_ * v - 1.3862944;//const=log(4)
        s = aa + r - w;
        //step 2
        if( 5.0 * z <= s + 2.609438 )//const=1+log(5)
            break;
        //step 3
        t = log( z );
        if( t <= s )
            break;
        //step 4
        if( t <= r + alpha_ * log( alpha_ /( bb + w )))
            break;
    }
    //step 5
    if( aa == a ) 
        *rv = w /( bb + w );
    else
        *rv = bb /( bb + w );
    if( maxn <= n )
        return PSL_MAXITERATS;
    return PSL_OK;
}
