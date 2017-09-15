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
#include "exp.h"
#include "log.h"
#include "digamma.h"
#include "gamma.h"
#include "beta.h"


// -------------------------------------------------------------------------
//
static double psl_loc_bt_isnegint( const double x )
{
    return ( x < 0 ) && ( x == floor(x));
}

// -------------------------------------------------------------------------
//
int psl_lnbeta_e( const double x, const double y, double* result, double* err )
{
    double sgn;
    int status = psl_lnbeta_sgn_e( x, y, result, err, &sgn );
    if( sgn == -1 )
        return PSL_ERR_DOMAIN;
    return status;
}

// -------------------------------------------------------------------------
//
int psl_lnbeta_sgn_e( const double x, const double y, 
                      double* result, double* err, double* sgn )
{
    if( x == 0.0 || y == 0.0 ) {
        *sgn = 0.0;
        return PSL_ERR_DOMAIN;
    }
    else if( psl_loc_bt_isnegint(x) || psl_loc_bt_isnegint(y)) {
        *sgn = 0.0;
        return PSL_ERR_DOMAIN; // * not defined for negative integers
    }

    // * see if we can handle the positive case with min/max < 0.2
    if( 0 < x && 0 < y )
    {
        const double max = SLC_MAX( x, y );
        const double min = SLC_MIN( x, y);
        const double rat = min/max;

        if(rat < 0.2) {
            // * min << max, so be careful with the subtraction
            double lnpre_val;
            double lnpre_err;
            double lnpow_val;
            double lnpow_err;
            double t1, t2, t3;
            double lnopr, lnoprerr;
            double gsx, gsy, gsxy;
            double gsxerr, gsyerr, gsxyerr;
            int status;
            if(( status = psl_gammastar_e( x, &gsx, &gsxerr )) != PSL_SUCCESS )
                return status;
            if(( status = psl_gammastar_e( y, &gsy, &gsyerr )) != PSL_SUCCESS )
                return status;
            if(( status = psl_gammastar_e( x+y, &gsxy, &gsxyerr )) != PSL_SUCCESS )
                return status;
            if(( status = psl_log_1plusx_e( rat, &lnopr, &lnoprerr )) != PSL_SUCCESS )
                return status;
            lnpre_val = log( gsx*gsy/gsxy * SLC_SQRT2*SLC_SQRTPI );
            lnpre_err = gsxerr/gsx + gsyerr/gsy + gsxyerr/gsxy;
            t1 = min*log(rat);
            t2 = 0.5*log(min);
            t3 = (x+y-0.5)*lnopr;
            lnpow_val  = t1 - t2 - t3;
            lnpow_err  = SLC_DBL_EPSILON * (fabs(t1) + fabs(t2) + fabs(t3));
            lnpow_err += fabs(x+y-0.5) * lnoprerr;
            if( result ) {
                *result = lnpre_val + lnpow_val;
                if( err ) {
                    *err  = lnpre_err + lnpow_err;
                    *err += 2.0 * SLC_DBL_EPSILON * fabs(*result);
                }
            }
            if( sgn )
                *sgn = 1.0;
            return PSL_SUCCESS;
        }
    }

    // * general case - fallback
    double lgxerr, lgyerr, lgxyerr;
    double lgx, lgy, lgxy;
    double sgx, sgy, sgxy, xy = x+y;
    int gxstatus  = psl_lngamma_sgn_e( x, &lgx, &lgxerr, &sgx );
    int gystatus  = psl_lngamma_sgn_e( y, &lgy, &lgyerr, &sgy );
    int gxystatus = psl_lngamma_sgn_e( xy, &lgxy, &lgxyerr, &sgxy );
    if( gxstatus != PSL_SUCCESS )
        return gxstatus;
    if( gystatus != PSL_SUCCESS )
        return gystatus;
    if( gxystatus != PSL_SUCCESS )
        return gxystatus;
    if( sgn )
        *sgn = sgx * sgy * sgxy;
    if( result ) {
        *result = lgx + lgy - lgxy;
        if( err ) {
            *err  = lgxerr + lgyerr + lgxyerr;
            *err += 2.0 * SLC_DBL_EPSILON * (fabs(lgx) + fabs(lgy) + fabs(lgxy));
            *err += 2.0 * SLC_DBL_EPSILON * fabs(*result);
        }
    }
    return PSL_SUCCESS;
}


// -------------------------------------------------------------------------
//
int psl_beta_e( const double x, const double y, double* result, double* err )
{
    if( result ) {
        *result = 0.0;
        if( err )
            *err = 0.0;
    }
    if( 0 < x && 0 < y && x < 50.0 && y < 50.0 ) {
        // * handle the easy case
        double gx, gy, gxy;
        double gxerr, gyerr, gxyerr;
        int status;
        if(( status = psl_gamma_e( x, &gx, &gxerr )) != PSL_SUCCESS )
            return status;
        if(( status = psl_gamma_e( y, &gy, &gyerr )) != PSL_SUCCESS )
            return status;
        if(( status = psl_gamma_e( x+y, &gxy, &gxyerr )) != PSL_SUCCESS )
            return status;
        if( result ) {
            *result = ( gx*gy )/ gxy;
            if( err ) {
                *err  = gxerr * fabs(gy/gxy);
                *err += gyerr * fabs(gx/gxy);
                *err += fabs((gx*gy)/(gxy*gxy)) * gxyerr;
                *err += 2.0 * SLC_DBL_EPSILON * fabs(*result);
            }
        }
        return PSL_SUCCESS;
    }
    else if( psl_loc_bt_isnegint(x) || psl_loc_bt_isnegint(y)) {
        return PSL_ERR_DOMAIN;
    }
    else if( psl_loc_bt_isnegint(x+y)) { // * infinity in the denominator
        return PSL_SUCCESS;
    }
    else {
        double lb, lberr;
        double sgn;
        int lbstatus = psl_lnbeta_sgn_e( x, y, &lb, &lberr, &sgn );
        if( lbstatus != PSL_SUCCESS ) 
            return lbstatus;
        int status = psl_exp_err_e( lb, lberr, result, err );
        if( result )
            *result *= sgn;
        return status;
    }
}


// ----------------- Functions w/ Natural Prototypes -----------------------
//

// double psl_lnbeta( const double x, const double y )
// {
//     double  result, err;
//     if( psl_lnbeta_e( x, y, &result, &err ) != PSL_SUCCESS )
//         ;
//     return result;
// }

// double gsl_sf_beta(const double x, const double y)
// {
//     double  result, err;
//     if( psl_beta_e( x, y, &result, &err ) != PSL_SUCCESS )
//         ;
//     return result;
// }
