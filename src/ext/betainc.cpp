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
#include "beta.h"
#include "betainc.h"


// -------------------------------------------------------------------------
//
static int psl_loc_isnegint( const double x ) 
{
    return ( x < 0 ) && ( x == floor(x));
}

// -------------------------------------------------------------------------
//
static int psl_loc_beta_cont_frac( const double a, const double b, const double x, 
                                   double* result, double* err )
{
    const unsigned int max_iter = 512;          // * control iterations
    const double cutoff = 2.0 * SLC_DBL_MIN;    // * control the zero cutoff
    unsigned int iter_count = 0;
    double cf;

    // * standard initialization for continued fraction
    double num_term = 1.0;
    double den_term = 1.0 - (a+b)*x/(a+1.0);
    if( fabs(den_term) < cutoff ) 
        den_term = cutoff;
    den_term = 1.0 / den_term;
    cf = den_term;

    while( iter_count < max_iter ) {
        const int k  = iter_count + 1;
        double coeff = k*(b-k)*x/(((a-1.0)+2*k)*(a+2*k));
        double delta_frac;

        // * first step
        den_term = 1.0 + coeff*den_term;
        num_term = 1.0 + coeff/num_term;
        if( fabs(den_term) < cutoff ) den_term = cutoff;
        if( fabs(num_term) < cutoff ) num_term = cutoff;
        den_term = 1.0/den_term;

        delta_frac = den_term * num_term;
        cf *= delta_frac;

        coeff = -(a+k)*(a+b+k)*x/((a+2*k)*(a+2*k+1.0));

        // * second step
        den_term = 1.0 + coeff*den_term;
        num_term = 1.0 + coeff/num_term;
        if( fabs(den_term) < cutoff ) den_term = cutoff;
        if( fabs(num_term) < cutoff ) num_term = cutoff;
        den_term = 1.0/den_term;

        delta_frac = den_term * num_term;
        cf *= delta_frac;

        if( fabs(delta_frac-1.0) < 2.0*SLC_DBL_EPSILON )
            break;

        iter_count++;
    }

    if( result ) {
        *result = cf;
        if( err )
            *err = iter_count * 4.0 * SLC_DBL_EPSILON * fabs(cf);
    }

    if( max_iter <= iter_count )
        return PSL_MAXITERATS;
    return PSL_SUCCESS;
}



// --------------------- Functions with Error Codes ------------------------
//

int psl_betainc_e( const double a, const double b, const double x, 
                   double* result, double* err )
{
    if( x < 0.0 || 1.0 < x ) {
        return PSL_ERR_DOMAIN;
    } else if( psl_loc_isnegint(a) || psl_loc_isnegint(b)) {
        return PSL_ERR_DOMAIN;
    } else if( psl_loc_isnegint(a+b)) {
        return PSL_ERR_DOMAIN;
    } else if( x == 0.0 ) {
        if( result ) {
            *result = 0.0;
            if( err )
                *err = 0.0;
        }
        return PSL_SUCCESS;
    }
    else if( x == 1.0 ) {
        if( result ) {
            *result = 1.0;
            if( err )
                *err = 0.0;
        }
        return PSL_SUCCESS;
    }
    else if( a <= 0 || b <= 0 ) {
        double f, beta;
        double ferr, betaerr;
        int stat;
        const int fstatus = PSL_ERR_ILLEGAL;//TODO://gsl_sf_hyperg_2F1_e( a, 1-b, a+1, x, &f, &ferr );
        const int betastatus = psl_beta_e( a, b, &beta, &betaerr );
        double prefactor = pow(x,a) / a;
        if( result ) {
            *result = prefactor * f / beta;
            if( err )
                *err = fabs(prefactor) * ferr/ fabs(beta) + fabs(*result/beta) * betaerr;
        }
        if( fstatus != PSL_SUCCESS )
            return fstatus;
        if( betastatus != PSL_SUCCESS )
            return betastatus;
        if( result ) {
            if( fabs(*result) < SLC_DBL_MIN )
                return PSL_ERR_UNDERFLOW;
        }
        return PSL_SUCCESS;
    } 
    else {
        double ln_beta, ln_betaerr;
        double ln_x, ln_xerr;
        double ln_1mx, ln_1mxerr;
        double prefactor, prefactorerr;
        const int ln_betastatus = psl_lnbeta_e(a, b, &ln_beta, &ln_betaerr);
        const int ln_1mxstatus = psl_log_1plusx_e( -x, &ln_1mx, &ln_1mxerr );
        const int ln_xstatus = psl_log_e( x, &ln_x, &ln_xerr );

        if( result ) {
            *result = 0.0;
            if( err )
                *err = 0.0;
        }
        if( ln_betastatus != PSL_SUCCESS )
            return ln_betastatus;
        if( ln_1mxstatus != PSL_SUCCESS )
            return ln_1mxstatus;
        if( ln_xstatus != PSL_SUCCESS )
            return ln_xstatus;

        const double ln_pre_val = -ln_beta + a * ln_x + b * ln_1mx;
        const double ln_pre_err =  ln_betaerr + fabs(a*ln_xerr) + fabs(b*ln_1mxerr);
        const int expstatus = psl_exp_err_e( ln_pre_val, ln_pre_err, &prefactor, &prefactorerr );

        if( x < (a+1.0)/(a+b+2.0)) {
            // * apply continued fraction directly
            double cf, cferr;
            const int cfstatus = psl_loc_beta_cont_frac( a, b, x, &cf, &cferr );

            if( result ) {
                *result = prefactor * cf / a;
                if( err )
                    *err = ( fabs(prefactorerr*cf) + fabs(prefactor*cferr)) / a;
            }

            if( expstatus != PSL_SUCCESS )
                return expstatus;
            if( cfstatus != PSL_SUCCESS )
                return cfstatus;
            if( result ) {
                if( fabs(*result) < SLC_DBL_MIN )
                    return PSL_ERR_UNDERFLOW;
            }
            return PSL_SUCCESS;
        }
        else {
            // * apply continued fraction after hypergeometric transformation
            double cf, cferr;
            const int cfstatus = psl_loc_beta_cont_frac( b, a, 1.0-x, &cf, &cferr );

            const double term = prefactor * cf / b;

            if( result ) {
                *result = 1.0 - term;
                if( err ) {
                    *err  = fabs(prefactorerr*cf) / b;
                    *err += fabs(prefactor*cferr) / b;
                    *err += 2.0 * SLC_DBL_EPSILON * ( 1.0+fabs(term));
                }
            }
            // * since the prefactor term is subtracted from 1 we need to
            // * ignore underflow
            if( expstatus != PSL_ERR_UNDERFLOW ) {
                if( expstatus != PSL_SUCCESS )
                    return expstatus;
            }
            if( cfstatus != PSL_SUCCESS )
                return cfstatus;
            if( result ) {
                if( fabs(*result) < SLC_DBL_MIN )
                    return PSL_ERR_UNDERFLOW;
            }
            return PSL_SUCCESS;
        }
    }
}


// ----------------- Functions w/ Natural Prototypes -----------------------
//

// double psl_betainc( const double a, const double b, const double x )
// {
//     double  result, err;
//     if( psl_betainc_e( a, b, x, &result, &err ) != PSL_SUCCESS )
//         ;
//     return result;
// }
