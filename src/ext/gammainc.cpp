/* Copyright (C) 2007 Brian Gough
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */
// *** the code adopted from GSL ***

#include <math.h>
#include "psl.h"
#include "pslcodes.h"
#include "log.h"
#include "exp.h"
#include "erf.h"
#include "expint.h"
#include "gamma.h"
#include "gammainc.h"


// -------------------------------------------------------------------------
// The dominant part:
// D(a,x) := x^a e^(-x) / Gamma(a+1)
//
static int gammainc_D( const double a, const double x, double* result, double* err )
{
    if( a < 10.0 ) {
        double  lnr;
        double  lg, lgerr;
        psl_lngamma_e( a + 1.0, &lg, &lgerr );
        lnr = a * log( x ) - x - lg;
        if( result ) {
            *result = exp( lnr );
            if( err )
                *err = 2.0 * SLC_DBL_EPSILON * ( fabs( lnr ) + 1.0 ) * fabs( *result );
        }
        return PSL_SUCCESS;
    }
    else {
        double  gstar, gstarerr;
        double  ln_term, ln_termerr;
        double  term1;
        if( x < 0.5 * a ) {
            double  u = x / a;   
            double  ln_u = log( u );
            ln_term = ln_u - u + 1.0;
            ln_termerr = ( fabs( ln_u ) + fabs( u ) + 1.0 ) * SLC_DBL_EPSILON;
        } else {
            double  mu = ( x - a ) / a;
            psl_log_1plusx_mx_e( mu, &ln_term, &ln_termerr );  // log(1+mu) - mu
            // * Propagate cancellation error from x-a, since the absolute
            // * error of mu=x-a is DBL_EPSILON
            ln_termerr += SLC_DBL_EPSILON * fabs( mu );
        };
        psl_gammastar_e( a, &gstar, &gstarerr );
        term1 = exp( a * ln_term ) / sqrt( 2.0 * SLC_PI * a );
        if( result ) {
            *result = term1 / gstar;
            if( err ) {
                *err = 2.0 * SLC_DBL_EPSILON * ( fabs( a * ln_term ) + 1.0 ) * fabs( *result );
                // * Include propagated error from log term
                *err += fabs( a ) * ln_termerr * fabs( *result );
                *err += gstarerr / fabs( gstar ) * fabs( *result );
            }
        }
        return PSL_SUCCESS;
    }
}


// -------------------------------------------------------------------------
// P, series representation.
//
static int gammainc_P_series( const double a, const double x, double* result, double* err )
{
    const int nmax = 10000;

    double  D, Derr;
    int     Dstatus = gammainc_D( a, x, &D, &Derr );

    // * Approximating the terms of the series using Stirling's
    // * approximation gives t_n = (x/a)^n * exp(-n(n+1)/(2a)), so the
    // * convergence condition is n^2 / (2a) + (1-(x/a) + (1/2a)) n >>
    // * -log(GSL_DBL_EPS) if we want t_n < O(1e-16) t_0. The condition
    // * below detects cases where the minimum value of n is > 5000

    if( 0.995 * a < x && 1e5 < a ) { // Difficult case: try continued fraction
        double  cf_res, cf_err;
        int     status = psl_exprel_n_CF_e( a, x, &cf_res, &cf_err );
        if( result ) {
            *result = D * cf_res;
            if( err )
                *err = fabs( D * cf_err ) + fabs( Derr * cf_res );
        }
        return status;
    }

    // Series would require excessive number of terms
    //
    if(( a + nmax ) < x )
        return PSL_ERR_GI_XRANGE;

    // Normal case: sum the series
    {
        double  sum  = 1.0;
        double  term = 1.0;
        double  remainder;
        int     n;

        // Handle lower part of the series where t_n is increasing, |x| > a+n
        int     nlow = ( a < x )? ( x - a ): 0;

        for( n = 1; n < nlow; n++ ) {
            term *= x / ( a + n );
            sum += term;
        }
        // Handle upper part of the series where t_n is decreasing, |x| < a+n
        for (/* n = previous n */ ; n < nmax; n++ )  {
            term *= x / ( a + n );
            sum += term;
            if( fabs( term / sum ) < SLC_DBL_EPSILON )
                break;
        }

        // Estimate remainder of series ~ t_(n+1)/(1-x/(a+n+1))
        {
            double    tnp1 = ( x / ( a + n )) * term;
            remainder = tnp1 / ( 1.0 - x / ( a + n + 1.0 ));
        }

        if( result ) {
            *result = D * sum;
            if( err ) {
                *err  = Derr * fabs( sum ) + fabs( D * remainder );
                *err += ( 1.0 + n ) * SLC_DBL_EPSILON * fabs( *result );
            }
        }
        if( n == nmax && SLC_SQRT_DBL_EPSILON < fabs( remainder / sum ))
            return PSL_ERR_GI_SMAXIT;
    }
    return Dstatus;
}


// -------------------------------------------------------------------------
// Q, large x, asymptotic
//
static int gammainc_Q_large_x( const double a, const double x, double* result, double* err )
{
    const int nmax = 5000;

    double  D, Derr;
    const int Dstatus = gammainc_D( a, x, &D, &Derr );

    double  sum  = 1.0;
    double  term = 1.0;
    double  last = 1.0;
    int     n;

    for( n = 1; n < nmax; n++ ) {
        term *= ( a - n ) / x;
        if( 1.0 < fabs( term / last ))
            break;
        if( fabs( term / sum ) < SLC_DBL_EPSILON )
            break;
        sum += term;
        last = term;
    }

    if( result ) {
        *result = D * ( a / x ) * sum;
        if( err ) {
            *err = Derr * fabs(( a / x ) * sum );
            *err += 2.0 * SLC_DBL_EPSILON * fabs( *result );
        }
    }
    if( n == nmax )
        return PSL_ERR_GI_AMAXIT;
    return Dstatus;
}


// -------------------------------------------------------------------------
// Uniform asymptotic for x near a, a and x large.
// See [Temme, p. 285]
//
static int gammainc_Q_asymp_unif( const double a, const double x, double* result, double* err )
{
    const double rta = sqrt( a );
    const double eps = ( x - a ) / a;

    double      ln_term, ln_term_err;
    const int   ln_status = psl_log_1plusx_mx_e( eps, &ln_term, &ln_term_err ); // log(1+eps) - eps
    const double eta  = SLC_SIGN( eps ) * sqrt( -2.0 * ln_term );

    double  erfc, erfc_err;

    double R;
    double c0, c1;

    psl_erfc_e( eta * rta / SLC_SQRT2, &erfc, &erfc_err );

    if( fabs( eps ) < SLC_ROOT5_DBL_EPSILON ) {
        c0 = -1.0 / 3.0 + eps * ( 1.0 / 12.0 - eps * ( 23.0 / 540.0 - eps * ( 353.0 / 12960.0 - eps * 589.0 / 30240.0 )));
        c1 = -1.0 / 540.0 - eps / 288.0;
    }
    else {  
        const double rt_term = sqrt( -2.0 * ln_term / ( eps * eps ));
        const double lam = x / a;
        c0 = ( 1.0 - 1.0 / rt_term ) / eps;
        c1 = -( SLC_POW3( eta ) * ( lam * lam + 10.0 * lam + 1.0 ) - 12.0 * SLC_POW3( eps )) /
              ( 12.0 * SLC_POW3( eta ) * SLC_POW3( eps ));
    }

    R = exp( -0.5 * a * eta * eta ) / ( SLC_SQRT2 * SLC_SQRTPI * rta ) * ( c0 + c1 / a );

    if( result ) {
        *result = 0.5 * erfc + R;
        if( err ) {
            *err = SLC_DBL_EPSILON * fabs( R * 0.5 * a * eta * eta ) + 0.5 * erfc_err;
            *err += 2.0 * SLC_DBL_EPSILON * fabs( *result );
        }
    }
    return ln_status;
}


// -------------------------------------------------------------------------
// * Continued fraction which occurs in evaluation
// * of Q(a,x) or Gamma(a,x).
// *
// *              1   (1-a)/x  1/x  (2-a)/x   2/x  (3-a)/x
// *   F(a,x) =  ---- ------- ----- -------- ----- -------- ...
// *             1 +   1 +     1 +   1 +      1 +   1 +
// *
// * Hans E. Plesser, 2002-01-22 (hans dot plesser at itf dot nlh dot no).
// *
// * Split out from gamma_inc_Q_CF() by GJ [Tue Apr  1 13:16:41 MST 2003].
// * See gamma_inc_Q_CF() below.
// *
static int gammainc_F_CF( const double a, const double x, double* result, double* err )
{
    const int    nmax  =  5000;
    const double small =  SLC_POW3( SLC_DBL_EPSILON );

    double  hn = 1.0; // convergent
    double  Cn = 1.0 / small;
    double  Dn = 1.0;
    int     n;

    // n == 1 has a_1, b_1, b_0 independent of a, x,
    // so that has been done by hand
    for( n = 2 ; n < nmax ; n++ )
    {
        double an;
        double delta;

        if( SLC_ODD( n ))
            an = 0.5 * ( n - 1 ) / x;
        else
            an = ( 0.5 * n - a ) / x;

        Dn = 1.0 + an * Dn;
        if( fabs( Dn ) < small )
            Dn = small;
        Cn = 1.0 + an / Cn;
        if ( fabs( Cn ) < small )
            Cn = small;
        Dn = 1.0 / Dn;
        delta = Cn * Dn;
        hn *= delta;
        if( fabs( delta - 1.0 ) < SLC_DBL_EPSILON )
            break;
    }

    if( result ) {
        *result = hn;
        if( err ) {
            *err = 2.0 * SLC_DBL_EPSILON * fabs( hn );
            *err += SLC_DBL_EPSILON * ( 2.0 + 0.5 * n ) * fabs( *result );
        }
    }

    if( n == nmax )
        return PSL_ERR_GI_FMAXIT;
    return PSL_SUCCESS;
}


// -------------------------------------------------------------------------
// * Continued fraction for Q.
// *
// * Q(a,x) = D(a,x) a/x F(a,x)
// *
// * Hans E. Plesser, 2002-01-22 (hans dot plesser at itf dot nlh dot no):
// *
// * Since the Gautschi equivalent series method for CF evaluation may lead
// * to singularities, I have replaced it with the modified Lentz algorithm
// * given in
// *
// * I J Thompson and A R Barnett
// * Coulomb and Bessel Functions of Complex Arguments and Order
// * J Computational Physics 64:490-509 (1986)
// *
// * In consequence, gamma_inc_Q_CF_protected() is now obsolete and has been
// * removed.
// *
// * Identification of terms between the above equation for F(a, x) and
// * the first equation in the appendix of Thompson&Barnett is as follows:
// *
// *    b_0 = 0, b_n = 1 for all n > 0
// *
// *    a_1 = 1
// *    a_n = (n/2-a)/x    for n even
// *    a_n = (n-1)/(2x)   for n odd
// *
static int gammainc_Q_CF( const double a, const double x, double* result, double* err )
{
    double  D, Derr;
    double  F, Ferr;
    const int   Dstatus = gammainc_D( a, x, &D, &Derr );
    const int   Fstatus = gammainc_F_CF( a, x, &F, &Ferr );

    if( result ) {
        *result = D * ( a / x ) * F;
        if( err )
            *err = Derr * fabs(( a / x ) * F ) + fabs( D * a / x * Ferr );
    }
    if( Fstatus != PSL_SUCCESS )
        return Fstatus;
    if( Dstatus != PSL_SUCCESS )
        return Dstatus;
    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// Q, series representation.
// Useful for small a and x. Handles the subtraction analytically.
//
static int gammainc_Q_series( const double a, const double x, double* result, double* err )
{
    double  term1;  // 1 - x^a/Gamma(a+1)
    double  sum;    // 1 + (a+1)/(a+2)(-x)/2! + (a+1)/(a+3)(-x)^2/3! + ...
    int     sum_status;
    double  term2;  // a temporary variable used at the end

    {
      // Evaluate series for 1 - x^a/Gamma(a+1), small a
      //
      const double pg21 = -2.404113806319188570799476;  // PolyGamma[2,1]
      const double lnx  = log( x );
      const double el   = SLC_EULER + lnx;
      const double c1 = -el;
      const double c2 = SLC_PI * SLC_PI / 12.0 - 0.5 * el * el;
      const double c3 = el * ( SLC_PI * SLC_PI / 12.0 - el * el / 6.0 ) + pg21 / 6.0;
      const double c4 = -0.04166666666666666667
                        * ( -1.758243446661483480 + lnx )
                        * ( -0.764428657272716373 + lnx )
                        * (  0.723980571623507657 + lnx )
                        * (  4.107554191916823640 + lnx );
      const double c5 = -0.0083333333333333333
                        * ( -2.06563396085715900 + lnx )
                        * ( -1.28459889470864700 + lnx )
                        * ( -0.27583535756454143 + lnx )
                        * (  1.33677371336239618 + lnx )
                        * (  5.17537282427561550 + lnx );
      const double c6 = -0.0013888888888888889
                        * ( -2.30814336454783200 + lnx )
                        * ( -1.65846557706987300 + lnx )
                        * ( -0.88768082560020400 + lnx )
                        * (  0.17043847751371778 + lnx )
                        * (  1.92135970115863890 + lnx )
                        * (  6.22578557795474900 + lnx );
      const double c7 = -0.00019841269841269841
                        * ( -2.5078657901291800 + lnx )
                        * ( -1.9478900888958200 + lnx )
                        * ( -1.3194837322612730 + lnx )
                        * ( -0.5281322700249279 + lnx )
                        * (  0.5913834939078759 + lnx )
                        * (  2.4876819633378140 + lnx )
                        * (  7.2648160783762400 + lnx );
      const double c8 = -0.00002480158730158730
                        * ( -2.677341544966400 + lnx )
                        * ( -2.182810448271700 + lnx )
                        * ( -1.649350342277400 + lnx )
                        * ( -1.014099048290790 + lnx )
                        * ( -0.191366955370652 + lnx )
                        * (  0.995403817918724 + lnx )
                        * (  3.041323283529310 + lnx )
                        * (  8.295966556941250 + lnx );
      const double c9 = -2.75573192239859e-6
                        * ( -2.8243487670469080 + lnx )
                        * ( -2.3798494322701120 + lnx )
                        * ( -1.9143674728689960 + lnx )
                        * ( -1.3814529102920370 + lnx )
                        * ( -0.7294312810261694 + lnx )
                        * (  0.1299079285269565 + lnx )
                        * (  1.3873333251885240 + lnx )
                        * (  3.5857258865210760 + lnx )
                        * (  9.3214237073814600 + lnx );
      const double c10 = -2.75573192239859e-7
                        * ( -2.9540329644556910 + lnx )
                        * ( -2.5491366926991850 + lnx )
                        * ( -2.1348279229279880 + lnx )
                        * ( -1.6741881076349450 + lnx )
                        * ( -1.1325949616098420 + lnx )
                        * ( -0.4590034650618494 + lnx )
                        * (  0.4399352987435699 + lnx )
                        * (  1.7702236517651670 + lnx )
                        * (  4.1231539047474080 + lnx )
                        * (  10.342627908148680 + lnx );

      term1 = 
      a * ( c1 + a *( c2 + a *( c3 + a *( c4 + a *( c5 + a *( c6 + a *( c7 + a *( c8 + a *( c9 + a * c10 )))))))));
    }
    {
      // Evaluate the sum
      //
      const int nmax = 5000;
      double    t = 1.0;
      int       n;

      sum = 1.0;

      for( n=1; n < nmax; n++ ) {
          t *= -x / ( n + 1.0 );
          sum += ( a + 1.0 ) / ( a + n + 1.0 ) * t;
          if( fabs( t / sum ) < SLC_DBL_EPSILON )
              break;
      }

      if( n == nmax )
          sum_status = PSL_ERR_GI_SMAXIT;
      else
          sum_status = PSL_SUCCESS;
    }

    term2 = ( 1.0 - term1 ) * a / ( a + 1.0 ) * x * sum;
    if( result ) {
        *result = term1 + term2;
        if( err ) {
            *err = SLC_DBL_EPSILON * ( fabs( term1 ) + 2.0 * fabs( term2 ));
            *err += 2.0 * SLC_DBL_EPSILON * fabs( *result );
        }
    }
    return sum_status;
}


// -------------------------------------------------------------------------
// unnormalized incomplete Gamma Function, series representation.
// series for small a and x, but not defined for a == 0
//
static int gammainc_series( double a, double x, double* result, double* err )
{
    double  Q, Qerr;
    double  G, Gerr;

    const int Qstatus = gammainc_Q_series( a, x, &Q, &Qerr );
    const int Gstatus = psl_gamma_e(a, &G, &Gerr );

    if( result ) {
        *result = Q * G;
        if( err ) {
            *err = fabs( Q * Gerr) + fabs( Qerr * G );
            *err += 2.0 * SLC_DBL_EPSILON * fabs( *result );
        }
    }
    if( Qstatus != PSL_SUCCESS )
        return Qstatus;
    if( Gstatus != PSL_SUCCESS )
        return Gstatus;
    return PSL_SUCCESS;
}


// -------------------------------------------------------------------------
// unnormalized incomplete Gamma Function; a > 0
//
static int gammainc_a_gt_0( double a, double x, double* result, double* err )
{
    // x > 0 and a > 0; use result for Q
    double  Q, Qerr;
    double  G, Gerr;
    const int Qstatus = psl_gammainc_Q_e( a, x, &Q, &Qerr );
    const int Gstatus = psl_gamma_e( a, &G, &Gerr );

    if( result ) {
        *result = G * Q;
        if( err ) {
            *err = fabs( G * Qerr ) + fabs( Gerr * Q );
            *err += 2.0 * SLC_DBL_EPSILON * fabs( *result );
        }
    }
    if( Gstatus != PSL_SUCCESS )
        return Gstatus;
    if( Qstatus != PSL_SUCCESS )
        return Qstatus;
    return PSL_SUCCESS;
}


// -------------------------------------------------------------------------
// unnormalized incomplete Gamma Function; Continued fraction representation
//
static int gammainc_CF( double a, double x, double* result, double* err )
{
    double  F, Ferr;
    double  pre, prerr;
    const double am1lgx = ( a - 1.0 ) * log( x );
    const int Fstatus = gammainc_F_CF( a, x, &F, &Ferr );
    const int Estatus = psl_exp_err_e( am1lgx - x, SLC_DBL_EPSILON * fabs( am1lgx ), &pre, &prerr );

    if( result ) {
        *result = F * pre;
        if( err ) {
            *err = fabs( Ferr * pre ) + fabs( F * prerr );
            *err += 2.0 * SLC_DBL_EPSILON * fabs( *result );
        }
    }
    if( Fstatus != PSL_SUCCESS )
        return Fstatus;
    if( Estatus != PSL_SUCCESS )
        return Estatus;
    return PSL_SUCCESS;
}


// -------------------------------------------------------------------------
// Gamma(0,x), x>0
#define GAMMA_INC_A_0( x, result, err ) psl_expint_E1_e( x, result, err )


// --------------------- Functions with Error Codes ------------------------
// Normalized incomplete Gamma Function Q (normalized upper inc. Gamma function)
//
int psl_gammainc_Q_e( const double a, const double x, double* result, double* err )
{
    if( a < 0.0 || x < 0.0 ) {
        return PSL_ERR_DOMAIN;
    }
    else if( x == 0.0 ) {
        if( result ) {
            *result = 1.0;
            if( err )
                *err = 0.0;
        }
        return PSL_SUCCESS;
    }
    else if( a == 0.0 ) {
        if( result ) {
            *result = 0.0;
            if( err )
                *err = 0.0;
        }
        return PSL_SUCCESS;
    }
    else if( x <= 0.5 * a ) {
        // * If the series is quick, do that. It is robust and simple
        //
        double  P, Perr;
        int     Pstatus = gammainc_P_series( a, x, &P, &Perr );
        if( result ) {
            *result = 1.0 - P;
            if( err ) {
                *err = Perr;
                *err += 2.0 * SLC_DBL_EPSILON * fabs( *result );
            }
        }
        return Pstatus;
    }
    else if( 1.0e+06 <= a && ( x - a )*( x - a ) < a ) {
        // * Then try the difficult asymptotic regime.
        // * This is the only way to do this region.
        return gammainc_Q_asymp_unif( a, x, result, err );
    }
    else if( a < 0.2 && x < 5.0 ) {
        // * Cancellations at small a must be handled analytically; 
        // * x should not be too big either since the 
        // * series terms grow with x and log(x).
        return gammainc_Q_series( a, x, result, err );
    }
    else if( a <= x ) {
        if( x <= 1.0e+06 ) {
            // * Continued fraction is excellent for x >~ a.
            // * We do not let x be too large when x > a since
            // * it is somewhat pointless to try this there;
            // * the function is rapidly decreasing for
            // * x large and x > a, and it will just
            // * underflow in that region anyway. We
            // * catch that case in the standard large-x method.
            return gammainc_Q_CF( a, x, result, err );
        }
        else {
            return gammainc_Q_large_x( a, x, result, err );
        }
    }
    else {
        if( a - sqrt( a ) < x ) {
            // * Continued fraction again. The convergence
            // * is a little slower here, but that is fine.
            // * We have to trade that off against the slow
            // * convergence of the series, which is the
            // * only other option.
            return gammainc_Q_CF( a, x, result, err );
        }
        else {
            double  P, Perr;
            int     Pstatus = gammainc_P_series( a, x, &P, &Perr );
            if( result ) {
                *result = 1.0 - P;
                if( err ) {
                    *err = Perr;
                    *err += 2.0 * SLC_DBL_EPSILON * fabs( *result );
                }
            }
            return Pstatus;
        }
    }
}


// -------------------------------------------------------------------------
// Normalized incomplete Gamma Function P 
// (normalized lower inc. Gamma function, complementary to Q: 1-Q)
//
int psl_gammainc_P_e( const double a, const double x, double* result, double* err )
{
    if( a <= 0.0 || x < 0.0 ) {
        return PSL_ERR_DOMAIN;
    }
    else if( x == 0.0 ) {
        if( result ) {
            *result = 0.0;
            if( err )
                *err = 0.0;
        }
        return PSL_SUCCESS;
    }
    else if( x < 20.0 || x < 0.5 * a ) {
        // * Do the easy series cases. Robust and quick
        return gammainc_P_series( a, x, result, err );
    }
    else if( 1.0e+06 < a && ( x - a )*( x - a ) < a ) {
        // * Crossover region. Note that Q and P are
        // * roughly the same order of magnitude here,
        // * so the subtraction is stable.
        double  Q, Qerr;
        int     Qstatus = gammainc_Q_asymp_unif( a, x, &Q, &Qerr );
        if( result ) {
            *result = 1.0 - Q;
            if( err ) {
                *err = Qerr;
                *err += 2.0 * SLC_DBL_EPSILON * fabs( *result );
            }
        }
        return Qstatus;
    }
    else if( a <= x ) {
        // * Q <~ P in this area, so the subtractions are stable
        double  Q, Qerr;
        int     Qstatus;
        if( 0.2 * x < a ) {
            Qstatus = gammainc_Q_CF( a, x, &Q, &Qerr );
        }
        else {
            Qstatus = gammainc_Q_large_x( a, x, &Q, &Qerr );
        }
        if( result ) {
            *result = 1.0 - Q;
            if( err ) {
                *err = Qerr;
                *err += 2.0 * SLC_DBL_EPSILON * fabs( *result );
            }
        }
        return Qstatus;
    }
    else {
        if(( x - a )*( x - a ) < a ) {
            // * This condition is meant to insure
            // * that Q is not very close to 1,
            // * so the subtraction is stable.
            double  Q, Qerr;
            int     Qstatus = gammainc_Q_CF( a, x, &Q, &Qerr );
            if( result ) {
                *result = 1.0 - Q;
                if( err ) {
                    *err = Qerr;
                    *err += 2.0 * SLC_DBL_EPSILON * fabs( *result );
                }
            }
            return Qstatus;
        }
        else {
            return gammainc_P_series( a, x, result, err );
        }
    }
}


// -------------------------------------------------------------------------
// Unnormalized incomplete Gamma Function
// (unnormalized upper inc. Gamma function)
//
int psl_gammainc_e( const double a, const double x, double* result, double* err )
{
    if( x < 0.0 ) {
        return PSL_ERR_DOMAIN;
    }
    else if( x == 0.0 ) {
        return psl_gamma_e( a, result, err );
    }
    else if( a == 0.0 ) {
        return GAMMA_INC_A_0( x, result, err );
    }
    else if( 0.0 < a ) {
        return gammainc_a_gt_0( a, x, result, err );
    }
    else if( 0.25 < x ) {
        // * continued fraction seems to fail for x too small; otherwise
        // * it is ok, independent of the value of |x/a|, because of the
        // * non-oscillation in the expansion, i.e. the CF is
        // * un-conditionally convergent for a < 0 and x > 0
        return gammainc_CF( a, x, result, err );
    }
    else if( fabs( a ) < 0.5 ) {
        return gammainc_series( a, x, result, err );
    }
    else {
        // * a = fa + da; da >= 0
        const double fa = floor( a );
        const double da = a - fa;

        double  g_da, g_da_err;
        const int g_da_status = (( 0.0 < da )? gammainc_a_gt_0( da, x, &g_da, &g_da_err ):
                                               GAMMA_INC_A_0( x, &g_da, &g_da_err ));
        double  alpha = da;
        double  gax = g_da;

        // * Gamma(alpha-1,x) = 1/(alpha-1) (Gamma(a,x) - x^(alpha-1) e^-x)
        do {
            const double shift = exp( -x + ( alpha - 1.0 ) * log( x ));
            gax = ( gax - shift ) / ( alpha - 1.0 );
            alpha -= 1.0;
        } while( a < alpha );

        if( result ) {
            *result = gax;
            if( err )
                *err = 2.0 * ( 1.0 + fabs( a )) * SLC_DBL_EPSILON * fabs( gax );
        }
        return g_da_status;
    }
}


// ----------------- Functions w/ Natural Prototypes -----------------------
//

// double psl_gammainc_P( const double a, const double x )
// {
//     double  result, err;
//     if( psl_gammainc_P_e( a, x, &result, &err ) != PSL_SUCCESS )
//         ;
//     return result;
// }
// 
// double psl_gammainc_Q( const double a, const double x )
// {
//     double  result, err;
//     if( psl_gammainc_Q_e( a, x, &result, &err ) != PSL_SUCCESS )
//         ;
//     return result;
// }
// 
// double psl_gammainc( const double a, const double x )
// {
//     double  result, err;
//     if( psl_gammainc_e( a, x, &result, &err ) != PSL_SUCCESS )
//         ;
//     return result;
// }

