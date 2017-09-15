/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <stdlib.h>
#include <math.h>

#include "myexcept.h"
#include "pcmath.h"
#include "segdata.h"

// Variable definitions

const _TSQUARES             PRESQUARES;
const _TLOGARITHMS          LOGARITHMS;
const _TPARTIAL_ENTROPIES   PRT_ENTROPIES;
_TLOG_GAMMA                 LOG_GAMMA;

// -------------------------------------------------------------------------
// _TPRECOMPUTED: precomputed values
// -------------------------------------------------------------------------

_TPRECOMPUTED::_TPRECOMPUTED( size_t size )
:
    no_vals( size ),
    data( NULL )
{
    data = ( double* )malloc( sizeof( double ) * no_vals );

    if( data == NULL )
        throw myruntime_error( mystring( "_TPRECOMPUTED: Not enough memory." ));
}

_TPRECOMPUTED::~_TPRECOMPUTED()
{
    if( data )
        free( data );
}

void _TPRECOMPUTED::Init()
{
    if( data == NULL )
        throw myruntime_error( mystring( "_TPRECOMPUTED: Memory access error." ));

    for( size_t v = 0; v < no_vals; v++ )
        data[v] = Compute( v );
}

// -------------------------------------------------------------------------
// _TLOGARITHMS: precomputes logarithms
// -------------------------------------------------------------------------

_TSQUARES::_TSQUARES( size_t size )
:
    _TPRECOMPUTED( size )
{
    Init();
}

_TSQUARES::~_TSQUARES()
{
}

double _TSQUARES::Compute( size_t v ) const
{
    return SQUARE(( double )v );
}

// -------------------------------------------------------------------------
// _TLOGARITHMS: precomputes logarithms
// -------------------------------------------------------------------------

_TLOGARITHMS::_TLOGARITHMS( size_t size )
:
    _TPRECOMPUTED( size )
{
    Init();
}

_TLOGARITHMS::~_TLOGARITHMS()
{
}

double _TLOGARITHMS::Compute( size_t v ) const
{
    if( v == 0 )
        return -1.0e27;

    return log(( double )v );
}

// -------------------------------------------------------------------------
// _TPARTIAL_ENTROPIES: precomputes entropy values
// -------------------------------------------------------------------------

_TPARTIAL_ENTROPIES::_TPARTIAL_ENTROPIES( size_t size )
:
    _TPRECOMPUTED( size )
{
    Init();
}

_TPARTIAL_ENTROPIES::~_TPARTIAL_ENTROPIES()
{
}

double _TPARTIAL_ENTROPIES::Compute( size_t v ) const
{
    if( v == 0 )
        return 0.0;

    return -( double )v * log(( double )v ) / LN2;
}

// -------------------------------------------------------------------------
// _TLOG_GAMMA: precomputes log-gamma function values
// -------------------------------------------------------------------------

_TLOG_GAMMA::_TLOG_GAMMA()
:
    no_vals( 0 ),
    data( NULL )
{
    Precompute( s_SZPRLGM );
}


_TLOG_GAMMA::~_TLOG_GAMMA()
{
    if( data )
        free( data );
}


void _TLOG_GAMMA::Precompute( size_t newsize )
{
    size_t  currentsize = GetSize();

    if( newsize <= currentsize )
        return;

    if( data == NULL )
        data = ( double* )malloc( sizeof( double ) * newsize );
    else
        data = ( double* )realloc( data, sizeof( double ) * newsize );

    if( data == NULL )
        throw myruntime_error( mystring( "_TLOG_GAMMA: Not enough memory." ));


    double  gammap1 = 0.0;
    size_t  portion = 50;
    size_t  zerosps = 2;
    size_t  beginwith = ( currentsize < zerosps )? zerosps: currentsize;

    for( size_t v = currentsize; v < zerosps && v < no_vals; v++ )
        data[v] = 0.0;

    for( size_t v = beginwith, pv = v; v < newsize; v = pv ) {
        gammap1 = 1.0;

        //compute log-gamma by portions to avoid overflows
        for( pv = v; pv < v + portion && pv < newsize; pv++ ) {
            gammap1 *= ( double )pv;
            data[pv] = log( gammap1 ) + data[v-1];
        }
    }

    SetSize( newsize );
}


double _TLOG_GAMMA::GetValueOf( size_t value )
{
    if( GetSize() == 0 )
        return 0.0;

    if( value < GetSize())
        return data[value];

    if( value < GetSize() + GetSize() && value < MAX_LOGGAMMA ) {
        Precompute( GetSize() + GetSize());
        return data[value];
    }

    double  lngammap1 = 0.0;

    for( ; value >= GetSize(); value-- )
        lngammap1 += log(( double )value );

    return lngammap1 + data[ GetSize() - 1 ];
}

