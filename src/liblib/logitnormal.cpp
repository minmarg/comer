/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "myexcept.h"
#include "logitnormal.h"

// -------------------------------------------------------------------------
// LogitNormalErrCorr: correct error of logit-normal variable to sum 
//     exactly to 1;
// NOTE: result is written back to `f'
//
void LogitNormalErrCorr( double* f, size_t fsz, double tol )
{
    double          vv, vsum, corr;
    size_t          r, non;
    char            strbuf[1024];
    const double    accuracy = 1.0e-6;

    if( f == NULL )
        return;

    if( fsz < 2 )
        throw myruntime_error("LogitNormalErrCorr: Invalid dimensionality.");

    vsum = 0.0;
    non = 0;
    for( r = 0; r < fsz; r++ ) {
        vv = f[r];
        if( vv < 0.0 )
            throw myruntime_error("LogitNormalErrCorr: Invalid logit-normal value.");
        if( vv != 0.0 )
            non++;
        vsum += vv;
    }

    corr = vsum - 1.0;
    if( corr == 0.0 )
        return;

    if( tol < corr ) {
        sprintf( strbuf, "LogitNormalErrCorr: Error tolerance exceeded: %g.", corr );
        throw myruntime_error( strbuf );
    }

    if( fabs( corr ) <= accuracy ) {
        for( r = 0; r < fsz; r++ ) {
            vv = f[r];
            //correct small error with the first value met
            if(( 0.0 < corr && corr <= vv ) ||
               ( corr < 0.0 && (( non < fsz && vv == 0.0 ) || 
                                ( fsz <= non && vv - corr <= 1.0 ))))
            {
                vv -= corr;
                f[r] = vv;//write
                break;
            }
        }
    }
    else if( accuracy < fabs( corr ) && non ) {
        corr /= ( double )non;
        for( r = 0; r < fsz; r++ ) {
            vv = f[r];
            if( vv == 0.0 )
                continue;
            vv -= corr;
            non--;
            //correct overflow/undeflow error
            if( vv < 0.0 ) {
                corr -= vv /( double )non;
                vv = 0.0;
            }
            else if( 1.0 < vv ) {
                corr += ( 1.0 - vv )/( double )non;
                vv = 1.0;
            }
            f[r] = vv;//write
        }
    }
}

// -------------------------------------------------------------------------
// LogitNormal2Normal: tranform logit-normally distributed variable to 
//      normally distributed variable;
//  NOTE: tranformed variable is written back to `f' and 
//  NOTE: dimensionality is reduced by 1
//
void LogitNormal2Normal( double* f, size_t fsz, double tol, bool correct )
{
    double          vv, vlast, vsum, corr;
    size_t          r, noz, non;
    double          fake = 1.0e-4;//fake value
    const double    accuracy = tol;//1.0e-6;

    if( f == NULL )
        return;

    if( fsz < 2 )
        throw myruntime_error("LogitNormal2Normal: Invalid dimensionality.");

    for( r = 0, vsum = 0.0; r < fsz; r++ )
        vsum += f[r];
    if( vsum < 1.0 - accuracy || 1.0 + accuracy < vsum ) {
        if( correct )
            LogitNormalErrCorr( f, fsz, tol );
        else
            throw myruntime_error("LogitNormal2Normal: Non-conserved logit-normal variable.");
    }

    vsum = 0.0;
    non = noz = 0;
    for( r = 0; r < fsz; r++ ) {
        vv = f[r];
        if( vv < 0.0 )
            throw myruntime_error("LogitNormal2Normal: Invalid logit-normal value.");
        if( vv == 0.0 ) {
            noz++;
            continue;
        }
        non++;
        vsum += vv;
    }

    if( vsum < 1.0 - accuracy || 1.0 + accuracy < vsum )
        throw myruntime_error("LogitNormal2Normal: Not a logit-normal variable.");
    if( noz < 0 || non <= 0 || noz + non != fsz )
        throw myruntime_error("LogitNormal2Normal: Invalid logit-normal vector.");

    vlast = f[fsz-1];

    if( correct ) {
        //make sure there are no null frequencies
        corr = fake * ( double )noz /( double )non;
        if( 0.0 < corr ) {
            vsum = 0.0;
            for( r = 0; r < fsz; r++ ) {
                vv = f[r];
                if( vv == 0.0 ) {
                    vsum += fake;
                    //null frequency, add fake value
                    f[r] = fake;//write
                    continue;
                }
                vv -= corr;
                non--;
                if( vv <= 0.0 ) {
                    vsum += fake;
                    //less than correction value; recalc. `corr'
                    f[r] = fake;//write
                    corr += ( fake - vv ) /( double )non;
                    continue;
                }
                vsum += vv;
                f[r] = vv;//write
            }
        }
        if( vsum < 1.0 - accuracy || 1.0 + accuracy < vsum || vlast <= 0.0 )
            throw myruntime_error("LogitNormal2Normal: Not a logit-normal variable: correct. 2.");
    }
    else
        if( 0 < noz )
            throw myruntime_error("LogitNormal2Normal: Zero values.");


    //apply transformation to normal distribution
    for( r = 0; r < fsz - 1; r++ ) {
        vv = f[r];
        if( vv <= 0.0 )
            throw myruntime_error("LogitNormal2Normal: Invalid logit-normal value: correct. 2.");
        vv = log( vv / vlast );
        f[r] = vv;//write
    }
    f[fsz-1] = 0.0;
}

// -------------------------------------------------------------------------
// Normal2LogitNormal: apply inverse logistic transformation to tranform 
//      normal random variable to logistic-normal random variable;
//  NOTE: result is written to lnv whose size lnvsz should be 1+ the 
//  size nvsz of the source nv
//
void Normal2LogitNormal( const double* nv, size_t nvsz, double* lnv, size_t lnvsz )
{
    if( nv == NULL || lnv == NULL )
        throw myruntime_error("Normal2LogitNormal: Null arguments.");
    if( nvsz < 1 || nvsz+1 != lnvsz )
        throw myruntime_error("Normal2LogitNormal: Invalid sizes.");

    int n;
    double  val, sum, tsum;
    sum = 0.0;
    for( n = 0; n < nvsz; n++ ) {
        val = exp( nv[n]);
        lnv[n] = val;
        sum += val;
    }
    tsum = 0.0;
    for( n = 0; n < nvsz; n++ ) {
        lnv[n] /= 1.0 + sum;
        tsum += lnv[n];
    }
    lnv[n] = 1.0 - tsum;
}
