/* Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
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

/* Author: G. Jungman */
// *** the code adopted from GSL ***

#include <math.h>
#include "psl.h"
#include "pslcodes.h"
#include "cheb.h"


// -------------------------------------------------------------------------
//
int cheb_eval_e( const Tcheb_series* cs, const double x, double* result, double* err )
{
    int     j;
    double  d  = 0.0;
    double  dd = 0.0;
    double  y  = ( 2.0 * x - cs->a - cs->b )/( cs->b - cs->a );
    double  y2 = 2.0 * y;
    double  e = 0.0;
    double  temp;

    for( j = cs->order; j >= 1; j-- ) {
        temp = d;
        d = y2 * d - dd + cs->c[j];
        e += fabs( y2 * temp ) + fabs( dd ) + fabs( cs->c[j] );
        dd = temp;
    }

        temp = d;
        d = y * d - dd + 0.5 * cs->c[0];
        e += fabs( y * temp ) + fabs( dd ) + 0.5 * fabs( cs->c[0] );

    if( result )
        *result = d;
    if( err )
        *err = SLC_DBL_EPSILON * e + fabs( cs->c[cs->order] );

    return PSL_OK;
}

// -------------------------------------------------------------------------
//
int psl_multiply_e( const double x, const double y, double* result, double* err )
{
    const double ax = fabs( x );
    const double ay = fabs( y );

    if( x == 0.0 || y == 0.0 ) {
        // It is necessary to eliminate this immediately
        if( result ) {
            *result = 0.0;
            if( err )
                *err = 0.0;
        }
        return PSL_SUCCESS;
    }
    else if(( ax <= 1.0 && 1.0 <= ay ) || ( ay <= 1.0 && 1.0 <= ax )) {
        // Straddling 1.0 is always safe
        if( result ) {
            *result = x * y;
            if( err )
                *err = 2.0 * SLC_DBL_EPSILON * fabs( *result );
        }
        return PSL_SUCCESS;
    }
    else {
        const double f = 1.0 - 2.0 * SLC_DBL_EPSILON;
        const double min = SLC_MIN( fabs( x ), fabs( y ));
        const double max = SLC_MAX( fabs( x ), fabs( y ));
        if( max < 0.9 * SLC_SQRT_DBL_MAX || min < ( f * SLC_DBL_MAX )/ max) {
            if( result ) {
                *result = ( x * y );
                if( err )
                    *err = 2.0 * SLC_DBL_EPSILON * fabs( *result );
                if( fabs( *result ) < SLC_DBL_MIN )
                    return PSL_ERR_UNDERFLOW;
            }
            return PSL_SUCCESS;
        }
        else {
            return PSL_ERR_OVERFLOW;
        }
    }
}


// -------------------------------------------------------------------------
//
int psl_multiply_err_e( const double x, const double dx, const double y, const double dy, double* result, double* err )
{
    int status = psl_multiply_e( x, y, result, err );
    if( err )
        *err += fabs( dx * y ) + fabs( dy * x );
    return status;
}

