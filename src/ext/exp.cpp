/* Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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
#include "gamma.h"
#include "exp.h"


// -------------------------------------------------------------------------
// * Evaluate the continued fraction for exprel.
// * [Abramowitz+Stegun, 4.2.41]
//
static int exprel_n_CF( const double N, const double x, double* result, double* err )
{
    const double RECUR_BIG = SLC_SQRT_DBL_MAX;
    const int maxiter = 5000;
    int     n = 1;
    double  Anm2 = 1.0;
    double  Bnm2 = 0.0;
    double  Anm1 = 0.0;
    double  Bnm1 = 1.0;
    double  a1 = 1.0;
    double  b1 = 1.0;
    double  a2 = -x;
    double  b2 = N + 1.0;
    double  an, bn;

    double fn;

    double An = b1 * Anm1 + a1 * Anm2;   // A1
    double Bn = b1 * Bnm1 + a1 * Bnm2;   // B1

    // * One explicit step, before we get to the main pattern
    n++;
    Anm2 = Anm1;
    Bnm2 = Bnm1;
    Anm1 = An;
    Bnm1 = Bn;
    An = b2 * Anm1 + a2 * Anm2;   // A2
    Bn = b2 * Bnm1 + a2 * Bnm2;   // B2

    fn = An / Bn;

    while( n < maxiter ) {
        double old_fn;
        double del;
        n++;
        Anm2 = Anm1;
        Bnm2 = Bnm1;
        Anm1 = An;
        Bnm1 = Bn;
        an = SLC_ODD( n )? (( n - 1 )/ 2 ) * x: -( N + ( n / 2 ) - 1 ) * x;
        bn = N + n - 1;
        An = bn * Anm1 + an * Anm2;
        Bn = bn * Bnm1 + an * Bnm2;

        if( RECUR_BIG < fabs( An ) || RECUR_BIG < fabs( Bn )) {
            An /= RECUR_BIG;
            Bn /= RECUR_BIG;
            Anm1 /= RECUR_BIG;
            Bnm1 /= RECUR_BIG;
            Anm2 /= RECUR_BIG;
            Bnm2 /= RECUR_BIG;
        }

        old_fn = fn;
        fn = An / Bn;
        del = old_fn / fn;

        if( fabs( del - 1.0 ) < 2.0 * SLC_DBL_EPSILON )
            break;
    }

    if( result ) {
        *result = fn;
        if( err )
            *err = 4.0 *( n + 1.0 )* SLC_DBL_EPSILON * fabs( fn );
    }
    if( n == maxiter )
        return PSL_MAXITERATS;
    return PSL_SUCCESS;
}


// --------------------- Functions with Error Codes ------------------------
// exp(x)
//
int psl_exp_e( const double x, double* result, double* err )
{
    if( SLC_LOG_DBL_MAX < x ) {
        return PSL_ERR_OVERFLOW;
    }
    else if( x < SLC_LOG_DBL_MIN ) {
        return PSL_ERR_UNDERFLOW;
    }
    else {
        if( result ) {
            *result = exp( x );
            if( err )
                *err = 2.0 * SLC_DBL_EPSILON * fabs( *result );
        }
        return PSL_SUCCESS;
    }
}

// -------------------------------------------------------------------------
// y exp(x)
//
int psl_exp_mult_e( const double x, const double y, double* result, double* err )
{
    const double ay = fabs( y );

    if( y == 0.0 ) {
        if( result ) {
            *result = 0.0;
            if( err )
                *err = 0.0;
        }
        return PSL_SUCCESS;
    }
    else if(( x < 0.5 * SLC_LOG_DBL_MAX   &&  0.5 * SLC_LOG_DBL_MIN < x ) &&
           ( ay < 0.8 * SLC_SQRT_DBL_MAX  &&  1.2 * SLC_SQRT_DBL_MIN < ay )) {
        const double ex = exp( x );
        if( result ) {
            *result = y * ex;
            if( err )
                *err = ( 2.0 + fabs( x )) * SLC_DBL_EPSILON * fabs( *result );
        }
        return PSL_SUCCESS;
    }
    else {
        const double ly  = log( ay );
        const double lnr = x + ly;

        if( SLC_LOG_DBL_MAX - 0.01 < lnr ) {
            return PSL_ERR_OVERFLOW;
        }
        else if( lnr < SLC_LOG_DBL_MIN + 0.01 ) {
            return PSL_ERR_UNDERFLOW;
        }
        else {
            const double sy   = SLC_SIGN( y );
            const double M    = floor( x );
            const double N    = floor( ly );
            const double a    = x  - M;
            const double b    = ly - N;
            const double berr = 2.0 * SLC_DBL_EPSILON * ( fabs( ly ) + fabs( N ));
            if( result ) {
                *result = sy * exp( M + N ) * exp( a + b );
                if( err ) {
                    *err = berr * fabs( *result );
                    *err += 2.0 * SLC_DBL_EPSILON * ( M + N + 1.0 ) * fabs( *result );
                }
            }
            return PSL_SUCCESS;
        }
    }
}


// -------------------------------------------------------------------------
// y exp(x) with add. error estimate
//
int psl_exp_mult_err_e( const double x, const double dx, const double y, const double dy, double* result, double* err )
{
    const double ay  = fabs(y);

    if( y == 0.0 ) {
        if( result ) {
            *result = 0.0;
            if( err )
                *err = fabs( dy * exp( x ));
        }
        return PSL_SUCCESS;
    }
    else if(( x < 0.5 * SLC_LOG_DBL_MAX   &&  0.5 * SLC_LOG_DBL_MIN < x ) &&
           ( ay < 0.8 * SLC_SQRT_DBL_MAX  &&  1.2 * SLC_SQRT_DBL_MIN < ay )) {
        double  ex = exp( x );
        if( result ) {
            *result = y * ex;
            if( err ) {
                *err = ex * ( fabs( dy ) + fabs( y * dx ));
                *err += 2.0 * SLC_DBL_EPSILON * fabs( *result );
            }
        }
        return PSL_SUCCESS;
    }
    else {
        const double ly  = log( ay );
        const double lnr = x + ly;

        if( SLC_LOG_DBL_MAX - 0.01 < lnr ) {
            return PSL_ERR_OVERFLOW;
        }
        else if( lnr < SLC_LOG_DBL_MIN + 0.01 ) {
            return PSL_ERR_UNDERFLOW;
        }
        else {
            const double sy  = SLC_SIGN( y );
            const double M   = floor( x );
            const double N   = floor( ly );
            const double a   = x  - M;
            const double b   = ly - N;
            const double eMN = exp( M + N );
            const double eab = exp( a + b );
            if( result ) {
                *result = sy * eMN * eab;
                if( err ) {
                    *err = eMN * eab * 2.0 * SLC_DBL_EPSILON;
                    *err += eMN * eab * fabs( dy / y );
                    *err += eMN * eab * fabs( dx );
                }
            }
            return PSL_SUCCESS;
        }
    }
}


// -------------------------------------------------------------------------
// exp(x)-1; accurate for small x
//
int psl_expm1_e( const double x, double* result, double* err )
{
    const double cut = 0.002;

    if( x < SLC_LOG_DBL_MIN ) {
        if( result ) {
            *result = -1.0;
            if( err )
                *err = SLC_DBL_EPSILON;
        }
        return PSL_SUCCESS;
    }
    else if( x < -cut ) {
        if( result ) {
            *result = exp( x ) - 1.0;
            if( err )
                *err = 2.0 * SLC_DBL_EPSILON * fabs( *result );
        }
        return PSL_SUCCESS;
    }
    else if( x < cut ) {
        if( result ) {
            *result = x *( 1.0 + 0.5 * x *( 1.0 + x / 3.0 *( 1.0 + 0.25 * x *( 1.0 + 0.2 * x ))));
            if( err )
                *err = 2.0 * SLC_DBL_EPSILON * fabs( *result );
        }
        return PSL_SUCCESS;
    } 
    else if( x < SLC_LOG_DBL_MAX ) {
        if( result ) {
            *result = exp( x ) - 1.0;
            if( err )
                *err = 2.0 * SLC_DBL_EPSILON * fabs( *result );
        }
        return PSL_SUCCESS;
    }
    else {
        return PSL_ERR_OVERFLOW;
    }
}


// -------------------------------------------------------------------------
// (exp(x)-1)/x; accurate for small x
//
int psl_exprel_e( const double x, double* result, double* err )
{
    const double cut = 0.002;

    if( x < SLC_LOG_DBL_MIN ) {
        if( result ) {
            *result = -1.0 / x;
            if( err )
                *err = SLC_DBL_EPSILON * fabs( *result );
        }
        return PSL_SUCCESS;
    }
    else if( x < -cut ) {
        if( result ) {
            *result = ( exp( x ) - 1.0 )/ x;
            if( err )
                *err = 2.0 * SLC_DBL_EPSILON * fabs( *result );
        }
        return PSL_SUCCESS;
    }
    else if( x < cut ) {
        if( result ) {
            *result = ( 1.0 + 0.5 * x *( 1.0 + x / 3.0 *( 1.0 + 0.25 * x *( 1.0 + 0.2 * x ))));
            if( err )
                *err = 2.0 * SLC_DBL_EPSILON * fabs( *result );
        }
        return PSL_SUCCESS;
    } 
    else if( x < SLC_LOG_DBL_MAX ) {
        if( result ) {
            *result = ( exp( x ) - 1.0 )/ x;
            if( err )
                *err = 2.0 * SLC_DBL_EPSILON * fabs( *result );
        }
        return PSL_SUCCESS;
    }
    else {
        return PSL_ERR_OVERFLOW;
    }
}


// -------------------------------------------------------------------------
// 2(exp(x)-1-x)/x^2; accurate for small x
//
int psl_exprel_2_e( double x, double* result, double* err )
{
    const double cut = 0.002;

    if( x < SLC_LOG_DBL_MIN ) {
        if( result ) {
            *result = -2.0 / x*( 1.0 + 1.0 / x );
            if( err )
                *err = 2.0 * SLC_DBL_EPSILON * fabs( *result );
        }
        return PSL_SUCCESS;
    }
    else if( x < -cut ) {
        if( result ) {
            *result = 2.0 *( exp( x ) - 1.0 - x )/( x * x );
            if( err )
                *err = 2.0 * SLC_DBL_EPSILON * fabs( *result );
        }
        return PSL_SUCCESS;
    }
    else if( x < cut ) {
        if( result ) {
            *result = ( 1.0 + 1.0 / 3.0 * x *( 1.0 + 0.25 * x *( 1.0 + 0.2 * x *( 1.0 + 1.0 / 6.0 * x ))));
            if( err )
                *err = 2.0 * SLC_DBL_EPSILON * fabs( *result );
        }
        return PSL_SUCCESS;
    } 
    else if( x < SLC_LOG_DBL_MAX ) {
        if( result ) {
            *result = 2.0 *( exp( x ) - 1.0 - x )/( x * x );
            if( err )
                *err = 2.0 * SLC_DBL_EPSILON * fabs( *result );
        }
        return PSL_SUCCESS;
    }
    else {
        return PSL_ERR_OVERFLOW;
    }
}


// -------------------------------------------------------------------------
// (exp(x)-1)/x; using continued fraction representation
//
int psl_exprel_n_CF_e( const double N, const double x, double* result, double* err )
{
    return exprel_n_CF( N, x, result, err );
}


// -------------------------------------------------------------------------
// N-relative exponential, n-th generalization of  exprel and exprel_2,
// N!/x^N (exp(x)-SUM(x^k/k!))
//
int psl_exprel_n_e( const int N, const double x, double* result, double* err )
{
    if( N < 0 ) {
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
    else if( fabs( x ) < SLC_ROOT3_DBL_EPSILON * N ) {
        if( result ) {
            *result = 1.0 + x /( N + 1 ) * ( 1.0 + x /( N + 2 ));
            if( err )
                *err = 2.0 * SLC_DBL_EPSILON;
        }
        return PSL_SUCCESS;
    }
    else if( N == 0 ) {
        return psl_exp_e( x, result, err );
    }
    else if( N == 1 ) {
        return psl_exprel_e( x, result, err );
    }
    else if( N == 2 ) {
        return psl_exprel_2_e( x, result, err );
    }
    else {
        if( N < x && ( -x + N *( 1.0 + log( x / N )) < SLC_LOG_DBL_EPSILON )) {
            // * x is much larger than N; ignore polynomial part, so
            // * exprel_N(x) ~= e^x N!/x^N
            double  lnf_N, lnf_Nerr;
            double  lnr_val;
            double  lnr_err;
            double  lnterm;
            psl_lnfact_e( N, &lnf_N, &lnf_Nerr );
            lnterm =  N * log( x );
            lnr_val = x + lnf_N - lnterm;
            lnr_err = SLC_DBL_EPSILON * ( fabs( x ) + fabs( lnf_N ) + fabs( lnterm ));
            lnr_err += lnf_Nerr;
            return psl_exp_err_e( lnr_val, lnr_err, result, err );
        }
        else if( N < x ) {
            // * Write the identity
            // *   exprel_n(x) = e^x n! / x^n (1 - Gamma[n,x]/Gamma[n])
            // * then use the asymptotic expansion
            // * Gamma[n,x] ~ x^(n-1) e^(-x) (1 + (n-1)/x + (n-1)(n-2)/x^2 + ...)
            double  ln_x = log( x );
            double  lnf_N, lnf_Nerr;
            double  lg_N;
            double  lnpre_val;
            double  lnpre_err;
            psl_lnfact_e( N, &lnf_N, &lnf_Nerr ); // log(N!)
            lg_N = lnf_N - log( N ); // log(Gamma(N))
            lnpre_val  = x + lnf_N - N * ln_x;
            lnpre_err  = SLC_DBL_EPSILON * ( fabs( x ) + fabs( lnf_N ) + fabs( N * ln_x ));
            lnpre_err += lnf_Nerr;
            if( lnpre_val < SLC_LOG_DBL_MAX - 5.0 ) {
                int     eGstatus;
                double  bigG_ratio, bigG_ratio_err;
                double  pre, prerr;
                int     exstatus = psl_exp_err_e( lnpre_val, lnpre_err, &pre, &prerr );
                double  ln_bigG_ratio_pre = -x + ( N - 1 )* ln_x - lg_N;
                double  bigGsum = 1.0;
                double  term = 1.0;
                int     k;
                for( k = 1; k < N; k++ ) {
                    term *= ( N - k )/ x;
                    bigGsum += term;
                }
                eGstatus = psl_exp_mult_e( ln_bigG_ratio_pre, bigGsum, &bigG_ratio, &bigG_ratio_err );
                if( eGstatus == PSL_SUCCESS ) {
                    if( result ) {
                        *result = pre * ( 1.0 - bigG_ratio );
                        if( err ) {
                            *err = pre * ( 2.0 * SLC_DBL_EPSILON + bigG_ratio_err );
                            *err += prerr * fabs( 1.0 - bigG_ratio );
                            *err += 2.0 * SLC_DBL_EPSILON * fabs( *result );
                        }
                    }
                    return exstatus;
                }
                else {
                    if( result ) {
                        *result = 0.0;
                        if( err )
                            *err = 0.0;
                    }
                    return eGstatus;
                }
            }
            else {
                return PSL_ERR_OVERFLOW;
            }
        }
        else if( -10.0 * N < x ) {
            return exprel_n_CF( N, x, result, err );
        }
        else {
            // * x -> -Inf asymptotic:
            // * exprel_n(x) ~ e^x n!/x^n - n/x (1 + (n-1)/x + (n-1)(n-2)/x + ...)
            // *             ~ - n/x (1 + (n-1)/x + (n-1)(n-2)/x + ...)
            double  sum  = 1.0;
            double  term = 1.0;
            int     k;
            for( k = 1; k < N; k++ ) {
                term *= ( N - k )/ x;
                sum += term;
            }
            if( result ) {
                *result = -N / x * sum;
                if( err )
                    *err = 2.0 * SLC_DBL_EPSILON * fabs( *result );
            }
            return PSL_SUCCESS;
        }
    }
}


// -------------------------------------------------------------------------
// exp(x) with add. error estimate
//
int psl_exp_err_e( const double x, const double dx, double* result, double* err )
{
    const double adx = fabs( dx );

    if( SLC_LOG_DBL_MAX < x + adx ) {
        return PSL_ERR_OVERFLOW;
    }
    else if( x - adx < SLC_LOG_DBL_MIN ) {
        return PSL_ERR_UNDERFLOW;
    }
    else {
        const double ex  = exp( x );
        const double edx = exp( adx );
        if( result ) {
            *result = ex;
            if( err ) {
                *err = ex * SLC_MAX( SLC_DBL_EPSILON, edx - 1.0 / edx );
                *err += 2.0 * SLC_DBL_EPSILON * fabs( *result );
            }
        }
        return PSL_SUCCESS;
    }
}


// -------------------------------------------------------------------------
