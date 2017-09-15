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
#include <string.h>
#include <ctype.h>

#include "CtxtFrequencies.h"

// number of allocated positions by default
const size_t    cAllocPoss = 100;

////////////////////////////////////////////////////////////////////////////
// CLASS CtxtFrequencies
//
// Default constructor
//

CtxtFrequencies::CtxtFrequencies()
:   allocated_( 0 ),
    length_( 0 ),
    sequence_( NULL ),
    obsfreqs_( NULL ),
    minval_( 1.0e6 ),
    maxval_(-1.0e6 )
{
}

// Constructor
//

CtxtFrequencies::CtxtFrequencies( size_t size )
:   allocated_( 0 ),
    length_( 0 ),
    sequence_( NULL ),
    obsfreqs_( NULL ),
    minval_( 1.0e6 ),
    maxval_(-1.0e6 )
{
    reallocate( size );
}

// -------------------------------------------------------------------------
// Destructor
// -------------------------------------------------------------------------

CtxtFrequencies::~CtxtFrequencies()
{
    destroy();
}

// -------------------------------------------------------------------------
// clear: erase all information
// -------------------------------------------------------------------------

void CtxtFrequencies::clear()
{
    size_t  alloc = GetAllocated();
    if( alloc ) {
        if( sequence_ )
            memset( sequence_, 0, sizeof( char ) * alloc );
        if( obsfreqs_ )
            memset( obsfreqs_, 0, sizeof( double ) * NUMALPH * alloc );
        if( !sequence_ || !obsfreqs_ )
            SetAllocated( 0 );
    }
    SetLength( 0 );
}

// -------------------------------------------------------------------------
// destroy: destroy information vectors
// -------------------------------------------------------------------------

void CtxtFrequencies::destroy()
{
    if( sequence_ ) { free( sequence_ ); sequence_ = NULL; }
    if( obsfreqs_ ) { free( obsfreqs_ ); obsfreqs_ = NULL; }
    SetLength( 0 );
    SetAllocated( 0 );
}

// -------------------------------------------------------------------------
// reallocate: allocate necessary memory
// -------------------------------------------------------------------------

void CtxtFrequencies::reallocate( size_t howmuch )
{
    size_t      alloc = GetAllocated();
    char*       tmpseqn;
    double   ( *tmpfreqs )[NUMALPH];

    if( howmuch <= alloc )
        return;

    if( alloc ) {
        tmpseqn = ( char* )realloc( sequence_, sizeof( char ) * howmuch );
        tmpfreqs = ( double(*)[NUMALPH] )realloc( obsfreqs_, sizeof( double ) * NUMALPH * howmuch );
    } else {
        tmpseqn = ( char* )malloc( sizeof( char ) * howmuch );
        tmpfreqs = ( double(*)[NUMALPH] )malloc( sizeof( double ) * NUMALPH * howmuch );
    }

    if( !tmpseqn || !tmpfreqs )
        throw myruntime_error( mystring( "CtxtFrequencies: Not enough memory." ));

    sequence_ = tmpseqn;
    obsfreqs_ = tmpfreqs;

    // fill uninitialized memory with zeros
    memset( sequence_ + alloc, 0, sizeof( char ) * ( howmuch - alloc ));
    memset( obsfreqs_ + alloc, 0, sizeof( double ) * NUMALPH * ( howmuch - alloc ));

    SetAllocated( howmuch );
}

// -------------------------------------------------------------------------
// FindMinMaxValues: find min and max values
//
void CtxtFrequencies::FindMinMaxValues( double* min, double* max )
{
    double  llmin =  1.0e6;
    double  llmax = -1.0e6;
    double  value;
    size_t  rescount = GetEffectiveNoResids();
    size_t  p, r;

    for( p = 0; p < GetLength(); p++ )
        for( r = 0; r < rescount; r++ ) {
            value = GetObsFrequency( p, r );
            if( value < llmin ) llmin = value;
            if( llmax < value ) llmax = value;
        }
    SetMinValue( llmin );
    SetMaxValue( llmax );

    if( min )   *min = llmin;
    if( max )   *max = llmax;
    return;
}

// =========================================================================
// Read: read observed frequencies from file
//
bool CtxtFrequencies::Read( FILE* fp )
{
    const size_t    rescount = GetEffectiveNoResids();

    if( fp == NULL )
        return false;

    mystring        errstr;
    int             emsg;
    size_t          len, rbts, read = 0;
    const size_t    locsize = KBYTE;
    char            locbuffer[locsize] = {0};
    char            r, res;
    double          val;
    const char*     endfrq = "*";
    const int       lenendfrq = strlen( endfrq );
    const bool      lookfor = !GetAllocated();

    clear();

    if(( emsg = skip_comments( fp, locbuffer, locsize, &len )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( feof( fp ))
        return false;

    for( r = 0; r < rescount; r++ ) {
        res = 0xff;
        if( len < read )
            throw myruntime_error( mystring( "CtxtFrequencies: Residues missing." ));

        if(( emsg = read_symbol( locbuffer + read, len - read, &res, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));
        res = ( char )HashAlphSymbol( res );
        read += rbts;

        if( res != r )
            throw myruntime_error( mystring( "CtxtFrequencies: Wrong format of frequencies." ));
    }

    if( !GetAllocated())
        reallocate( TIMES2( GetAllocated() + 1 ));

    for( size_t n = 0; n < GetAllocated() && !feof( fp ); n++ )
    {
        read = 0;
        if(( emsg = skip_comments( fp, locbuffer, locsize, &len )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( lookfor ) {
            if( len && strncmp( locbuffer, endfrq, lenendfrq ) == 0 )
                break;
            if( GetAllocated() <= GetLength() + 1 )
                reallocate( TIMES2( GetAllocated() + 1 ));
        }

        for( r = 0; r < rescount; r++ )
        {
            if( len < read )
                throw myruntime_error( mystring( "CtxtFrequencies: Frequencies missing." ));

            //read frequency
            if(( emsg = read_double( locbuffer + read, len - read, &val, &rbts )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));
            read += rbts;

            if( val < 0.0 )
                throw myruntime_error( mystring( "CtxtFrequencies: Negative frequency." ));

            SetObsFrequency( n, r, val );
        }
    }
    return true;
}

// -------------------------------------------------------------------------
// Write: write observed frequencies to file
//
void CtxtFrequencies::Write( FILE* fp ) const
{
    const size_t    rescount = GetEffectiveNoResids();
    unsigned char   r;

    if( fp == NULL )
        return;

    for( r = 0; r < rescount; r++ ) fprintf( fp, "%5c ", DehashCode( r ));

    for( size_t p = 0; p < GetLength(); p++ ) {
        fprintf( fp, "\n" );
        for( r = 0; r < rescount; r++ )
            fprintf( fp, " %5d", ( int )rint( GetObsFrequency( p, r )));
    }

    fprintf( fp, "\n" );
}

// -------------------------------------------------------------------------
// WriteDoubles: write double values to file
//
void CtxtFrequencies::WriteDoubles( FILE* fp ) const
{
    const size_t    rescount = GetEffectiveNoResids();
    unsigned char   r;

    if( fp == NULL )
        return;

    for( r = 0; r < rescount; r++ ) fprintf( fp, "%12c ", DehashCode( r ));

    for( size_t p = 0; p < GetLength(); p++ ) {
        fprintf( fp, "\n" );
        for( r = 0; r < rescount; r++ )
            fprintf( fp, " %12.5g", GetObsFrequency( p, r ));
    }

    fprintf( fp, "\n" );
}

// -------------------------------------------------------------------------
// WriteDoubles1: alternative writing of double values; excludes the last
//     residue (assuming all frequencies sum to 1)
//
void CtxtFrequencies::WriteDoubles1( FILE* fp, bool title ) const
{
    const int   rescount = GetEffectiveNoResids() - 1;
    unsigned char   r;
    int     p;

    if( fp == NULL )
        return;

    if( title ) {
        for( r = 0; r < rescount && r < 1; r++ ) fprintf( fp, "## %9c ", DehashCode( r ));
        for( ; r < rescount; r++ ) fprintf( fp, "%12c ", DehashCode( r ));
        fprintf( fp, "\n" );
    }
    for( p = 0; p < GetLength(); p++ ) {
        for( r = 0; r < rescount; r++ )
            fprintf( fp, " %12.6g", GetObsFrequency( p, r ));
        fprintf( fp, "\n" );
    }
}


////////////////////////////////////////////////////////////////////////////
// CLASS CtxtProScores
//
// Constructor
//

CtxtProScores::CtxtProScores()
:   effthickness_( 0 ),
    scores_( NULL )
{
}

// -------------------------------------------------------------------------
// Destructor
// -------------------------------------------------------------------------

CtxtProScores::~CtxtProScores()
{
    destroy();
}

// -------------------------------------------------------------------------
// clear: erase all information
// -------------------------------------------------------------------------

void CtxtProScores::clear()
{
    size_t  alloc = GetAllocated();
    if( alloc ) {
        if( scores_ )
            memset( scores_, 0, sizeof( double ) * NUMALPH * alloc );
    }

    CtxtFrequencies::clear();
}

// -------------------------------------------------------------------------
// destroy: destroy information vectors
// -------------------------------------------------------------------------

void CtxtProScores::destroy()
{
    if( scores_ ) { free( scores_ ); scores_ = NULL; }

    CtxtFrequencies::destroy();
}

// -------------------------------------------------------------------------
// reallocate: allocate necessary memory
// -------------------------------------------------------------------------

void CtxtProScores::reallocate( size_t howmuch )
{
    size_t      alloc = GetAllocated();
    double   ( *tmpscors )[NUMALPH];

    if( howmuch <= alloc )
        return;

    if( alloc ) {
        tmpscors = ( double(*)[NUMALPH] )realloc( scores_, sizeof( double ) * NUMALPH * howmuch );
    } else {
        tmpscors = ( double(*)[NUMALPH] )malloc( sizeof( double ) * NUMALPH * howmuch );
    }

    if( !tmpscors )
        throw myruntime_error( mystring( "CtxtFrequencies: Not enough memory." ));

    scores_   = tmpscors;

    // fill uninitialized memory with zeros
    memset( scores_ + alloc, 0, sizeof( double ) * NUMALPH * ( howmuch - alloc ));

    CtxtFrequencies::reallocate( howmuch );
}



// =========================================================================
// Read: read profile in text format
// =========================================================================

void CtxtProScores::Read( const char* filename )
{
    FILE*   fp = NULL;

    if( !filename || !strlen( filename ))
        throw myruntime_error( mystring( "CtxtFrequencies: Wrong filename." ));

    fp = fopen( filename, "r" );

    if( fp == NULL )
        throw myruntime_error( mystring( "CtxtFrequencies: Failed to open file for reading." ));

    mystring        errstr;
    size_t          len, read = 0;
    const size_t    locsize = TIMES2( KBYTE );
    char            locbuffer[locsize] = {0};
    int             emsg;

    clear();
    reallocate( cAllocPoss );

    try {
        if(( emsg = skip_comments( fp, locbuffer, locsize, &len )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ))
            throw myruntime_error( mystring( "CtxtFrequencies: Invalid format of profile." ));

        ReadThickness( locbuffer, len, GetEffThicknessAddr());

        if(( emsg = skip_comments( fp, locbuffer, locsize, &len )) != 0 ) throw myruntime_error( TranslateReadError( emsg ));
        if(( emsg = skip_comments( fp, locbuffer, locsize, &len )) != 0 ) throw myruntime_error( TranslateReadError( emsg ));
        if(( emsg = skip_comments( fp, locbuffer, locsize, &len )) != 0 ) throw myruntime_error( TranslateReadError( emsg ));

        for( size_t n = 0; ; n++ )
        {
            if(( emsg = skip_comments( fp, locbuffer, locsize, &len )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));
            if( GetAllocated() <= n )
                reallocate( TIMES2( GetAllocated() + 1 ));

            if( !( read = ReadResidue( locbuffer, len, GetResidueAddrAt( n ))))
                break;
            read += ReadValues( locbuffer + read, len, GetScoresAddrAt( n ), true );  //scores
            read += ReadValues( locbuffer + read, len, GetObsFreqsAddrAt( n ), false);//frequencies

            IncLength();
        }

    } catch( myexception const& ex )
    {
        errstr = ex.what();
    }

    fclose( fp );

    if( !errstr.empty())
        throw myruntime_error( errstr );
}

// -------------------------------------------------------------------------
// ReadThickness: read effective thickness value
// -------------------------------------------------------------------------

void CtxtProScores::ReadThickness( const char* readfrom, size_t readlen, int* membuf )
{
    const char* p = readfrom;
    const char* pbeg = NULL;
    const char* pend = NULL;

    mystring    number;
    int         tmpval;
    char*       paux;
    const char*     pattern = "effective thickness";
    const size_t    lenpat = strlen( pattern );

#ifdef __DEBUG__
    if( !readfrom || !membuf )
        throw myruntime_error( mystring( "CtxtFrequencies: Memory access error." ));
#endif

    while( *p && *p != '\n' && *p != '\r') {
        for( ; *p == ' ' || *p == '\t'; p++ );

        if( strncmp( p, pattern, lenpat )) {
            for( ; *p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r'; p++ );
            continue;
        }

        p += lenpat;
        for( ; *p == ' ' || *p == '\t'; p++ );
        for( pbeg = p; isdigit( *p ); p++ );

        if( pbeg == p )
            throw myruntime_error( mystring( "CtxtFrequencies: Effective thickness not found in profile." ));

        number = mystring( pbeg, size_t( p - pbeg ));

        if( !number.empty()) {
            errno = 0;  //NOTE: Thread unsafe
            tmpval = strtol( number.c_str(), &paux, 10 );

            if( errno || *paux )
                throw myruntime_error( mystring( "CtxtFrequencies: Invalid thickness value in profile." ));
            if( tmpval < 0 )
                throw myruntime_error( mystring( "CtxtFrequencies: Thickness found is negative." ));

            if( membuf )
                *membuf = tmpval;

            break;
        }
    }
}

// -------------------------------------------------------------------------
// ReadResidue: read residue at the position
// -------------------------------------------------------------------------

size_t CtxtProScores::ReadResidue( const char* readfrom, size_t readlen, char* membuf )
{
    const char* p = readfrom;
    const char* pbeg = NULL;
    const char* pend = NULL;

    mystring    number;
    double      tmpval;
    char*       paux;

#ifdef __DEBUG__
    if( !readfrom || !membuf )
        throw myruntime_error( mystring( "CtxtFrequencies: Memory access error." ));
#endif

    for( ; *p == ' ' || *p == '\t'; p++ );
    for( pbeg = p; *p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r'; p++ );

    for( ; pbeg != p; pbeg++ )
        if( !isdigit( *pbeg ))
            return 0;
    //omit verification whether the number is valid

    for( ; *p == ' ' || *p == '\t'; p++ );
    for( pbeg = p; *p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r'; p++ );

    if( size_t( p - pbeg ) != 1 )
        throw myruntime_error( mystring( "CtxtFrequencies: Residue symbol not found in profile." ));

    if( membuf )
        *membuf = ( char )HashAlphSymbol( *pbeg );

    return size_t( p - readfrom );
}

// -------------------------------------------------------------------------
// ReadValues: read scores or observed frequencies at the profile position
// -------------------------------------------------------------------------

size_t CtxtProScores::ReadValues(
    const char* readfrom,
    size_t      readlen,
    double*     membuf,
    bool        allownegs )
{
    const char* p = readfrom;
    const char* pbeg = NULL;
    const char* pend = NULL;

    mystring    number;
    double      tmpval;
    char*       paux;

    size_t      current = 0;    //serial number of values
    size_t      maxvals = NUMALPH - 1 ; //gap symbol excluded

#ifdef __DEBUG__
    if( !readfrom || !membuf )
        throw myruntime_error( mystring( "CtxtFrequencies: Memory access error." ));
#endif

    while( *p && *p != '\n' && *p != '\r' && current < maxvals ) {
        for( ; *p == ' ' || *p == '\t'; p++ );
        for( pbeg = p; *p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r'; p++ );

        if( pbeg == p )
            throw myruntime_error( mystring( "CtxtFrequencies: Invalid profile format." ));

        number = mystring( pbeg, size_t( p - pbeg ));

        if( !number.empty()) {
            errno = 0;  //NOTE: Thread unsafe
            tmpval = strtod( number.c_str(), &paux );

            if( errno || *paux )
                throw myruntime_error( mystring( "CtxtFrequencies: Invalid values in profile." ));
            if( !allownegs && tmpval < 0.0 )
                throw myruntime_error( mystring( "CtxtFrequencies: Invalid values in profile." ));

            membuf[current++] = tmpval;
        }
    }

    if( current < maxvals )
        throw myruntime_error( mystring( "CtxtFrequencies: Too few values read from profile." ));

    return size_t( p - readfrom );
}

// =========================================================================
// Write: write scores and observed frequencies to file
// =========================================================================

void CtxtProScores::Write( FILE* fp ) const
{
    const size_t    rescount = NUMALPH - 1; //gap symbol excluded
    unsigned char   r;

    if( fp == NULL )
        return;

    fprintf( fp, "(effective thickness %d)\n\n", GetEffThickness());

    fprintf( fp,"%28c Position-specific scoring matrix "
                "%53c Weighted observed frequencies %30c\n",
                32, 32, 32 );

    fprintf( fp, "%9c", 32 );

    for( r = 0; r < rescount; r++ ) fprintf( fp, "%3c", DehashCode( r ));
    for( r = 0; r < rescount; r++ ) fprintf( fp, "%4c", DehashCode( r ));

    for( size_t p = 0; p < GetLength(); p++ )
    {
        fprintf( fp, "\n%5d %c   ", p + 1, DehashCode( GetResidue( p )));

        for( r = 0; r < rescount; r++ )
            fprintf( fp, "%2d ", ( int )rint( GetScore( p, r )));

        for( r = 0; r < rescount; r++ )
            fprintf( fp, "%4d", ( int )rint( GetObsFrequency( p, r )));
    }

    fprintf( fp, "\n\n" );
}

