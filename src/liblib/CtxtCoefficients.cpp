/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "CtxtCoefficients.h"

// number of allocated positions by default
const size_t    cAllocPoss = 100;

////////////////////////////////////////////////////////////////////////////
// CLASS CtxtCoefficients
//
// Default constructor
//

CtxtCoefficients::CtxtCoefficients( size_t size )
:   length_( 0 ),
    cweight_( 0.0 ),
    multip_( 0.0 ),
    coeffs_( NULL ),
    logcoeffs_( NULL )
{
    Init( size );
}

// Constructor
//

CtxtCoefficients::CtxtCoefficients( size_t size, double weight )
:   length_( 0 ),
    cweight_( 0.0 ),
    multip_( 0.0 ),
    coeffs_( NULL ),
    logcoeffs_( NULL )
{
    Init( size );
    SetCentralWeight( weight );
}

// Default constructor
//
CtxtCoefficients::CtxtCoefficients()
{
    throw myruntime_error( mystring( "CtxtCoefficients: Default construction is not allowed." ));
}

// -------------------------------------------------------------------------
// Destructor
// -------------------------------------------------------------------------

CtxtCoefficients::~CtxtCoefficients()
{
    Destroy();
}

// -------------------------------------------------------------------------
// reallocate: allocate necessary memory
// -------------------------------------------------------------------------

void CtxtCoefficients::Init( size_t size )
{
    Destroy();

    coeffs_ = ( double* )malloc( size * sizeof( double ));
    logcoeffs_ = ( double* )malloc( size * sizeof( double ));

    if( !coeffs_ || !logcoeffs_ )
        throw myruntime_error( mystring( "CtxtCoefficients: Not enough memory." ));

    memset( coeffs_, 0, size * sizeof( double ));
    memset( logcoeffs_, 0, size * sizeof( double ));

    SetLength( size );
}

// -------------------------------------------------------------------------
// Destroy: destroy data
//
void CtxtCoefficients::Destroy()
{
    if( coeffs_ ) {
        free( coeffs_ );
        coeffs_ = NULL;
    }
    if( logcoeffs_ ) {
        free( logcoeffs_ );
        logcoeffs_ = NULL;
    }
}

// -------------------------------------------------------------------------
// fdffunction: function and its derivative to find root of multiplier
//     f = 2w SUM x^i + w - 1
//    df = 2w SUM ix^(i-1), w is central weight
//
void CtxtCoefficients::fdfunction( double x, double* f, double* df, void* params )
{
    if( !params )
        return;

    double  w = (( double* )params )[0];
    double  w2 = TIMES2( w );
    double  ldf = 0.0;
    double  lf = 1.0;
    size_t  length = ( size_t )(( double* )params )[1];
    size_t  nd2 = length / 2;
    size_t  n;

    for( n = 0; n < nd2; n++ ) {
        ldf = ldf * x + ( nd2 - n );
        lf  = lf  * x + 1.0;
    }
    ldf = w2 * ldf;
    lf = w2 * ( lf - 1.0 ) + w - 1.0;

    if( df ) *df = ldf;
    if( f )  *f  = lf;
    return;
}

// -------------------------------------------------------------------------
// ExpandCoefficients: expand coefficients for each position
//
double CtxtCoefficients::ExpandCoefficients()
{
    double      cweight = GetCentralWeight();
    double      multipl = GetMultiplier();
    double      value = cweight;
    double      sum = cweight;
    size_t      length = GetLength();
    size_t      nd2 = length / 2;
    size_t      n;

    if(!( length % 2 ))
        return 0.0;

    SetCoefficientAt( nd2, value );
    SetLogCoefficientAt( nd2 );

    for( n = 1; n <= nd2; n++ ) {
        value *= multipl;
        SetCoefficientAt( nd2 + n, value );
        SetCoefficientAt( nd2 - n, value );
        SetLogCoefficientAt( nd2 + n );//log
        SetLogCoefficientAt( nd2 - n );//log
        sum += TIMES2( value );
    }

    return sum;
}

// -------------------------------------------------------------------------
// FindCoefficients: find missing coefficients given central weight
//
void CtxtCoefficients::FindCoefficients()
{
    double      limit = 0.001;
    double      accuracy = limit * 0.1;
    double      x1 = 0.0 + accuracy;
    double      x2 = 5.0;
    double      lolim = 0.0 + limit;
    double      uplim = 1.0 - limit;
    double      cweight = GetCentralWeight();
    double      params[2] = { cweight, GetLength() };
    double      multiplier = 0.0;
    double      consv = 0.0;
    const char* emsg = NULL;
    char        strbuf[BUF_MAX];

    if(!( GetLength() % 2 ))
        throw myruntime_error( "CtxtCoefficients: Number of coefficients should be odd." );

    if( cweight != 1. ) {
        if( cweight <= lolim || uplim <= cweight )
            throw myruntime_error( "CtxtCoefficients: Invalid central weight." );

        emsg =  root_by_NR_and_bisection(
            &CtxtCoefficients::fdfunction,
            x1, x2, accuracy, MAXIT, params, &multiplier
        );
    }
    if( emsg != NULL )
        throw myruntime_error( emsg );

    if( multiplier < 0.0 )
        throw myruntime_error( "CtxtCoefficients: Failed to find root of multiplier." );

    SetMultiplier( multiplier );

    consv = ExpandCoefficients();
    if( consv <= 1.0 - accuracy || 1.0 + accuracy <= consv ) {
        sprintf( strbuf, "CtxtCoefficients: Conservation of coefficients failed: %g.", consv );
        throw myruntime_error( strbuf );
    }
    return;
}

// =========================================================================
// Write: write coefficients to file
//
void CtxtCoefficients::Write( FILE* fp ) const
{
    size_t  n;
    if( fp == NULL )
        return;

    for( n = 0; n < GetLength(); n++ )
        fprintf( fp, " %8g", GetCoefficientAt( n ));
    fprintf( fp, "\n" );
}

// Write: write logs of coefficients to file
//
void CtxtCoefficients::WriteLogs( FILE* fp ) const
{
    size_t  n;
    if( fp == NULL )
        return;

    for( n = 0; n < GetLength(); n++ )
        fprintf( fp, " %8g", GetLogCoefficientAt( n ));
    fprintf( fp, "\n" );
}
