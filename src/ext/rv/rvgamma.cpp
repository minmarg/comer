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
#include "ext/log.h"
#include "ext/exp.h"
#include "rvexp.h"
#include "rvnorm.h"
#include "rvgamma.h"

// =========================================================================
// RVGamma
//
RVGamma::RVGamma( Rng& rng, TRVG met )
:   RVar( rng ),
    met_( met ),
    shape_( 1.0 ),
    scale_( 1.0 ),
    //initialization of the AD82 state vars.
    Aap_( 0.0 ),
    Aapp_( 0.0 ),
    As_( 0.0 ), As_2_( 0.0 ), Ad_( 0.0),
    Aq0_( 0.0 ), Ab_( 0.0 ), Asgm_( 0.0 ), Ac_( 0.0 ),
    //initialization of the MT00 state vars.
    Map_( 0.0 ),
    Md_( 0.0 ), Mc_( 0.0 )
{
}

RVGamma::~RVGamma()
{
}

// -------------------------------------------------------------------------
// GenAD82: generate random gamma variates using the 
//  Ahrens & Dieter algorithms;
//  (shape par. a>=1) Ahrens, Dieter. (1982) Comm. ACM 25, 47-54.
//  (shape par. a<1) Ahrens, Dieter. (1974) Computing 12, 223-46 (Alg. GS),
//      and Devroye (1986). Non Uniform Random Variate Generation, p.425.
//
//  pdf: 1 / (Gamma(shape) scale^shape) x^(shape-1) e^(-x/scale),
//      (x,shape,scale>0)
//
int RVGamma::GenAD82( double* rv, double shape )
{
    if( rv == NULL )
        return PSL_ERR_ADDRESS;

    static const bool   slcb_useMTalt1 = true;
    static const int    maxn = 100;
    static const double sqrt32 = 5.6568542494923805819;//sqrt(32)
    static const double exptom1 = 0.36787944117144232159;//1/e

    //q coef. to approximate q0 = sum(q[k]*a^(-k)), (|err|<3.2e-8)
    static const double q1 = 0.04166669;
    static const double q2 = 0.02083148;
    static const double q3 = 0.00801191;
    static const double q4 = 0.00144121;
    static const double q5 = -7.388e-5;
    static const double q6 = 2.4511e-4;
    static const double q7 = 2.424e-4;

    //a coef. to approximate q = q0+(t*t/2)*sum(a[k]*v^k), (|err|<1.4e-7)
    static const double a1 = 0.3333333;
    static const double a2 = -0.250003;
    static const double a3 = 0.2000062;
    static const double a4 = -0.1662921;
    static const double a5 = 0.1423657;
    static const double a6 = -0.1367177;
    static const double a7 = 0.1233795;

    //e coef. to approximate exp(q)-1 = sum(e[k]*q^k), (|err|<1.1e-7)
    static const double e1 = 1.0;
    static const double e2 = 0.4999897;
    static const double e3 = 0.166829;
    static const double e4 = 0.0407753;
    static const double e5 = 0.0102930;

    double  bb, e, p, q, r, t, u, v, w, wt, x, l1px;
    RVNorm  rnv( GetRng(), RVNorm::TRVN_Ziggurat_M00 );
    RVExp   rev( GetRng());
    int     n, err;

    double  a = shape;

    if( a <= 0.0 )
        return PSL_ERR_INVALID;

    //shape par. a < 1
    if( a < 1.0 ) {
        //Marsaglia & Tsang - 2000 (see below) 
        //recursive relation for shape par. a < 1
        if( slcb_useMTalt1 ) {
            u = GetRng().GetDouble0();
            //NOTE: RECURSION
            err = GenAD82( rv, 1.0 + a );
            if( err == PSL_OK )
                *rv *= exp( log(u)/ a );
            return err;
        }

        //GS algorithm for shape par. a < 1
        bb = 1.0 + exptom1 * a;
        for( n = 0; n < maxn; n++ ) {
            p = bb * GetRng().GetDouble0();
            //generate random exponential variate
            if(( err = rev.Gen( &e )) != PSL_OK )
                return err;
            if( p <= 1.0 ) {
                x = exp( log( p )/ a );
                //e=-log(U) is distributed exponentially;
                //thus it's equivalent testing U<=exp(-x) and e>=x
                if( x <= e )
                    break;
            }
            else {
                x = -log(( bb - p )/ a );
                //same as above: accept if U<=x^(a-1) or e>=(1-a)log(x)
                if(( 1.0 - a )* log( x ) <= e )
                    break;
            }
        }
        if( maxn <= n )
            return PSL_MAXITERATS;
        return PSL_OK;
    }

    //GD algorithm for shape par. a >= 1
    //
    //step 1: recalculation of s2, s, d if a has changed
    if( a != Aap_ ) {
        Aap_ = a;
        As_2_ = a - 0.5;
        As_ = sqrt( As_2_ );
        Ad_ = sqrt32 - 12.0 * As_;
    }

    //step 2: t - std. normal, x - (s,1/2) normal deviates;
    if(( err = rnv.Gen( &t )) != PSL_OK )
        return err;
    x = As_ + 0.5 * t;
    *rv = x * x;
    if( 0.0 <= t )
        //immediate acceptance (i)
        return PSL_OK;

    //step 3: u - uniform sample; squeeze acceptance (s)
    u = GetRng().GetDouble1();
    if( Ad_ * u <= t * t * t )
        return PSL_OK;

    //step 4: recalculations of q0, b, sgm, c if necessary
    if( a != Aapp_ ) {
        Aapp_ = a;
        r = 1.0 / a;
        Aq0_ = (((((( q7*r + q6 )*r + q5 )*r + q4 )*r + q3 )*r + q2 )*r + q1 )*r;

        //approximation depending on size of parameter a;
        //constants in expressions were established in AD82
        if( a <= 3.686 ) {
            Ab_ = 0.463 + As_ + 0.178 * As_2_;
            Asgm_ = 1.235;
            Ac_ = 0.195 / As_ - 0.079 + 0.16 * As_;
        } else if( a <= 13.022 ) {
            Ab_ = 1.654 + 0.0076 * As_2_;
            Asgm_ = 1.68 / As_ + 0.275;
            Ac_ = 0.062 / As_ + 0.024;
        } else {
            Ab_ = 1.77;
            Asgm_ = 0.75;
            Ac_ = 0.1515 / As_;
        }
    }

    //step 5: no quotient q test if x not positive
    if( 0.0 < x ) {
        //step 6: calculation of v and quotient q
        v = t /( As_ + As_ );
        if( fabs( v ) <= 0.25 )
            q = Aq0_ + 0.5*t*t * (((((( a7*v + a6 )*v + a5 )*v + a4 )*v + a3 )*v + a2 )*v + a1 )*v;
        else {
            if(( err = psl_log_1plusx_e( v, &l1px, NULL/*err*/ )) != PSL_OK )
                return err;
            q = Aq0_ - As_*t + 0.25*t*t + ( As_2_+As_2_ )* l1px;//l1px=log(1+v)
        }

        //step 7: quotient acceptance (q)
        if(( err = psl_log_1plusx_e( -u, &l1px, NULL/*err*/ )) != PSL_OK )
            return err;
        if( l1px <= q )//:log(1-u)<=q
            return PSL_OK;
    }

    for( n = 0; n < maxn; n++ ) {
        //step 8: e - std. exponential, u - uniform, 
        // t - (b,si) double exponential (laplace) deviates
        if(( err = rev.Gen( &e )) != PSL_OK )
            return err;
        u = GetRng().GetDouble();
        u = u + u - 1.0;
        if( u < 0.0 )
            t = Ab_ - e * Asgm_;
        else
            t = Ab_ + e * Asgm_;

        //step 9:  rejection if t <= tau_1 = -0.71874483771719
        if( -0.71874483771719 < t ) {
            //step 10: calculation of v and quotient q
            v = t /( As_ + As_ );
            if( fabs( v ) <= 0.25 )
                q = Aq0_ + 0.5*t*t * (((((( a7*v + a6 )*v + a5 )*v + a4 )*v + a3 )*v + a2 )*v + a1 )*v;
            else {
                if(( err = psl_log_1plusx_e( v, &l1px, NULL/*err*/ )) != PSL_OK )
                    return err;
                q = Aq0_ - As_*t + 0.25*t*t + ( As_2_+As_2_ )* l1px;//l1px=log(1+v)
            }

            //step 11: hat acceptance (h);
            //go to step 8 if q not positive
            if( 0.0 < q ) {
                wt = e - 0.5*t*t ;
                if( q <= 0.5 )
                    w = (((( e5*q + e4 )*q + e3 )*q + e2 )*q + e1 )*q;
                else if( SLC_LOG_DBL_MAX <= q || SLC_LOG_DBL_MAX <= q + wt ) {
                    w = SLC_DBL_MAX;
                    break;
                }
                else
                    if(( err = psl_expm1_e( q, &w, NULL/*err*/ )) != PSL_OK )
                        return err;
                if( wt <= SLC_LOG_DBL_MIN ) {
                    wt = 0.0;
                    continue;
                }
                wt = exp( wt );
                //sample again at step 8 if t is rejected
                if( Ac_ * fabs( u ) <= w * wt )
                    break;
            }
        }
    }//until t accept
    if( maxn <= n )
        return PSL_MAXITERATS;

    x = As_ + 0.5 * t;
    *rv = x * x;
    return PSL_OK;
}

// -------------------------------------------------------------------------
// GenMT00: generate random gamma variates using algorithms (a<1 & a>=1) in
//  Marsaglia & Tsang. (2000) 
//      ACM Transactions on Mathematical Software 26(3), 363-72.
//
int RVGamma::GenMT00( double* rv, double shape )
{
    if( rv == NULL )
        return PSL_ERR_ADDRESS;

    static const int maxn = 100;
    int     n, err;

    static const double v1o3 = 0.333333333333;
    double  x, v, u;
    double  a = shape;

    if( a <= 0.0 )
        return PSL_ERR_INVALID;

    if( a < 1 ) {
        u = GetRng().GetDouble0();
        //NOTE: RECURSION
        err = GenMT00( rv, 1.0 + a );
        if( err == PSL_OK )
            *rv *= exp( log(u)/ a );
        return err;
    }

    RVNorm  rnv( GetRng(), RVNorm::TRVN_Ziggurat_M00 );

    if( a != Map_ ) {
        Map_ = a;
        Md_ = a - v1o3;
        Mc_ = v1o3 / sqrt( Md_ );
    }

    for( n = 0; n < maxn; n++ ) {
        do {//generate std. normal x
            if(( err = rnv.Gen( &x )) != PSL_OK )
                return err;
            v = 1.0 + Mc_ * x;
        } while( v <= 0.0 );

        v = v * v * v;
        u = GetRng().GetDouble0();

        if( u < 1.0 - 0.0331 * x * x * x * x )
            break;

        if( log(u) < 0.5 * x * x + Md_ * ( 1.0 - v + log(v)))
            break;
    }
    if( maxn <= n )
        return PSL_MAXITERATS;

    *rv = Md_ * v;
    return PSL_OK;
}
