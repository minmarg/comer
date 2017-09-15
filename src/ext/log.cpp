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
#include "cheb.h"
#include "log.h"

// -------------------------- Private Section ------------------------------
//
// * Chebyshev expansion for log(1 + x(t))/x(t)
// *
// * x(t) = (4t-1)/(2(4-t))
// * t(x) = (8x+1)/(2(x+2))
// * -1/2 < x < 1/2
// * -1 < t < 1
//
static double lopx_data[21] = {
    2.16647910664395270521272590407,
   -0.28565398551049742084877469679,
    0.01517767255690553732382488171,
   -0.00200215904941415466274422081,
    0.00019211375164056698287947962,
   -0.00002553258886105542567601400,
    2.9004512660400621301999384544e-06,
   -3.8873813517057343800270917900e-07,
    4.7743678729400456026672697926e-08,
   -6.4501969776090319441714445454e-09,
    8.2751976628812389601561347296e-10,
   -1.1260499376492049411710290413e-10,
    1.4844576692270934446023686322e-11,
   -2.0328515972462118942821556033e-12,
    2.7291231220549214896095654769e-13,
   -3.7581977830387938294437434651e-14,
    5.1107345870861673561462339876e-15,
   -7.0722150011433276578323272272e-16,
    9.7089758328248469219003866867e-17,
   -1.3492637457521938883731579510e-17,
    1.8657327910677296608121390705e-18
};
static Tcheb_series lopx_cs = {
    lopx_data,
    20,
    -1, 1,
    10
};

// * Chebyshev expansion for (log(1 + x(t)) - x(t))/x(t)^2
// *
// * x(t) = (4t-1)/(2(4-t))
// * t(x) = (8x+1)/(2(x+2))
// * -1/2 < x < 1/2
// * -1 < t < 1
//
static double lopxmx_data[20] = {
   -1.12100231323744103373737274541,
    0.19553462773379386241549597019,
   -0.01467470453808083971825344956,
    0.00166678250474365477643629067,
   -0.00018543356147700369785746902,
    0.00002280154021771635036301071,
   -2.8031253116633521699214134172e-06,
    3.5936568872522162983669541401e-07,
   -4.6241857041062060284381167925e-08,
    6.0822637459403991012451054971e-09,
   -8.0339824424815790302621320732e-10,
    1.0751718277499375044851551587e-10,
   -1.4445310914224613448759230882e-11,
    1.9573912180610336168921438426e-12,
   -2.6614436796793061741564104510e-13,
    3.6402634315269586532158344584e-14,
   -4.9937495922755006545809120531e-15,
    6.8802890218846809524646902703e-16,
   -9.5034129794804273611403251480e-17,
    1.3170135013050997157326965813e-17
};
static Tcheb_series lopxmx_cs = {
    lopxmx_data,
    19,
    -1, 1,
    9
};


// --------------------- Functions with Error Codes ------------------------
// log(x) for x>0
//
int psl_log_e( const double x, double* result, double* err )
{
    if( x <= 0.0 ) {
        return PSL_ERR_DOMAIN;
    }
    else {
        if( result ) {
            *result = log( x );
            if( err )
                *err = 2.0 * SLC_DBL_EPSILON * fabs( *result );
        }
        return PSL_SUCCESS;
    }
}


// -------------------------------------------------------------------------
// log(|x|) for x!=0
//
int psl_log_abs_e( const double x, double* result, double* err )
{
    if( x == 0.0 ) {
        return PSL_ERR_DOMAIN;
    }
    else {
        if( result ) {
            *result = log( fabs( x ));
            if( err )
                *err = 2.0 * SLC_DBL_EPSILON * fabs( *result );
        }
        return PSL_SUCCESS;
    }
}


// -------------------------------------------------------------------------
// log(1+x) for x>-1; accurate for small x
//
int psl_log_1plusx_e( const double x, double* result, double* err )
{
    if( x <= -1.0 ) {
        return PSL_ERR_DOMAIN;
    }
    else if( fabs( x ) < SLC_ROOT6_DBL_EPSILON ) {
        const double c1 = -0.5;
        const double c2 =  1.0 / 3.0;
        const double c3 = -1.0 / 4.0;
        const double c4 =  1.0 / 5.0;
        const double c5 = -1.0 / 6.0;
        const double c6 =  1.0 / 7.0;
        const double c7 = -1.0 / 8.0;
        const double c8 =  1.0 / 9.0;
        const double c9 = -1.0 / 10.0;
        const double t  =  c5 + x *( c6 + x *( c7 + x *( c8 + x * c9 )));
        if( result ) {
            *result = x * ( 1.0 + x *( c1 + x *( c2 + x *( c3 + x *( c4 + x * t )))));
            if( err )
                *err = SLC_DBL_EPSILON * fabs( *result );
        }
        return PSL_SUCCESS;
    }
    else if( fabs( x ) < 0.5 ) {
        double  t = 0.5 *( 8.0 * x + 1.0 )/( x + 2.0 );
        double  c, cerr;
        cheb_eval_e( &lopx_cs, t, &c, &cerr );
        if( result ) {
            *result = x * c;
            if( err )
                *err = fabs( x * cerr );
        }
        return PSL_SUCCESS;
    }
    else {
        if( result ) {
            *result = log( 1.0 + x );
            if( err ) {
                *err = SLC_DBL_EPSILON * fabs( *result );
            }
        }
        return PSL_SUCCESS;
    }
}


// -------------------------------------------------------------------------
// log(1+x)-x for x>-1; accurate for small x
//
int psl_log_1plusx_mx_e( const double x, double* result, double* err )
{
    if( x <= -1.0 ) {
        return PSL_ERR_DOMAIN;
    }
    else if( fabs( x ) < SLC_ROOT5_DBL_EPSILON ) {
        const double c1 = -0.5;
        const double c2 =  1.0 / 3.0;
        const double c3 = -1.0 / 4.0;
        const double c4 =  1.0 / 5.0;
        const double c5 = -1.0 / 6.0;
        const double c6 =  1.0 / 7.0;
        const double c7 = -1.0 / 8.0;
        const double c8 =  1.0 / 9.0;
        const double c9 = -1.0 / 10.0;
        const double t  =  c5 + x *( c6 + x *( c7 + x *( c8 + x * c9 )));
        if( result ) {
            *result = x * x *( c1 + x *( c2 + x *( c3 + x *( c4 + x * t ))));
            if( err )
                *err = SLC_DBL_EPSILON * fabs( *result );
        }
        return PSL_SUCCESS;
    }
    else if( fabs( x ) < 0.5 ) {
        double t = 0.5 *( 8.0 * x + 1.0 )/( x + 2.0 );
        double  c, cerr ;
        cheb_eval_e( &lopxmx_cs, t, &c, &cerr );
        if( result ) {
            *result = x * x * c;
            if( err )
                *err = x * x * cerr;
        }
        return PSL_SUCCESS;
    }
    else {
        const double lterm = log( 1.0 + x );
        if( result ) {
            *result = lterm - x;
            if( err )
                *err = SLC_DBL_EPSILON * ( fabs( lterm ) + fabs( x ));
        }
        return PSL_SUCCESS;
    }
}


// -------------------------------------------------------------------------
