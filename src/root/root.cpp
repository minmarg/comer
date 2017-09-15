/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "rcodes.h"
#include "root.h"

// -------------------------------------------------------------------------
// rootByNRandBisection: finds the root of a function bracketed between x1
//     and x2 using a combination of Newton-Raphson and bisection methods.
//     A Bisection step is taken whenever Newton-Raphson would take the solution
//     out of bounds, or whenever it is not reducing the size of the brackets
//     rapidly enough.
//     The root will be refined until its accuracy is reached within ±xacc or
//     maximum number of iteration, maxit, has been passed. The method must be
//     supplied with a supplied routine that returns both the function value
//     and the first derivative of the function.
// Computation logic is as in
//     Numerical Recipes in C by W.H.Press et al.
//  fdfunction -- evaluator of function value and its derivative
//  x1, x2 -- interval containing root
//  root -- variable to contain the root value on output
//  xacc -- required accuracy of result
//  maxit -- maximum number of iterations to apply
//  params -- user parameters to be passed to `fdfunction'
// -------------------------------------------------------------------------

int rootByNRandBisection(
        TRFunction fdfunction,
        double x1,
        double x2,
        double* root,
        double xacc,
        int maxit,
        void* params )
{
    double  df, dx, dxold, f, fh, fl;
    double  temp, xh, xl, rts;
    int     n;

    ( *fdfunction )( x1, &fl, &df, params );
    ( *fdfunction )( x2, &fh, &df, params );

    if(( fl > 0.0 && fh > 0.0 ) || ( fl < 0.0 && fh < 0.0 ))
        return PRT_ERR_DOMAIN;

    if( fl < 0.0 ) { //orient the search so that f (xl) < 0
        xl = x1;
        xh = x2;
    } else {
        xh = x1;
        xl = x2;
    }
    rts = 0.5 * ( x1 + x2 ); //Initialize the guess for root,
    dxold = fabs( x2 - x1 ); //the "stepsize before last,"
    dx = dxold;              //and the last step.

    ( *fdfunction )( rts, &f, &df, params );

    for( n = 0; n < maxit; n++ )
    {   //Bisect if Newton out of range, or not decreasing fast enough
        if( ((( rts - xh ) * df - f ) * (( rts - xl ) * df - f ) > 0.0 ) ||
               ( fabs( 2.0 * f ) > fabs( dxold * df ))) {
//         if(1){//purely bisection, for TESTING only
            dxold = dx;
            dx = 0.5 * ( xh - xl );
            rts = xl + dx;
            if( xl == rts ) break; //return rts; //change in root is negligible
        } else { //Newton step acceptable
            dxold = dx;
            dx = f / df;
            temp = rts;
            rts -= dx;
            if( temp == rts ) break; //return rts;
        }
        if( fabs( dx ) < xacc ) break; //return rts; //convergence criterion
        ( *fdfunction )( rts, &f, &df, params );
        if( f < 0.0 )
            xl = rts;
        else
            xh = rts;
    }
    if( root )
        *root = rts;
    if( n == maxit )
        return PRT_MAXITERATS;

    return PRT_OK;
}

