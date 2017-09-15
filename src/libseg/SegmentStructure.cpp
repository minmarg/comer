/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <string.h>
#include <stdlib.h>

#include "mystring.h"
#include "myexcept.h"

#include "SegmentStructure.h"


// -------------------------------------------------------------------------
// constructor: initialization
// -------------------------------------------------------------------------

SegmentStructure::SegmentStructure( size_t size )
:   segments( NULL ),
    length_( 0 ),
    capacity_( 0 )
{
    Realloc( size );
}

// -------------------------------------------------------------------------
// constructor: default
// -------------------------------------------------------------------------

SegmentStructure::SegmentStructure()
:   segments( NULL ),
    length_( 0 ),
    capacity_( 0 )
{
    throw myruntime_error(
            mystring( "SegmentStructure: Default initialization prohibited." ));
}

// -------------------------------------------------------------------------
// destructor:
// -------------------------------------------------------------------------

SegmentStructure::~SegmentStructure()
{
    if( segments )
        free( segments );
}

// -------------------------------------------------------------------------
// Clear: clears all elements
// -------------------------------------------------------------------------

void SegmentStructure::Clear()
{
    memset( segments, 0, sizeof( size_t ) * TIMES2( capacity_ ));
    length_ = 0;
}

// -------------------------------------------------------------------------
// Realloc: memory allocation
// -------------------------------------------------------------------------

void SegmentStructure::Realloc( size_t newcap )
{
    size_t*     tmp_segments;

    if( newcap <= capacity_ )
        return;

    if( capacity_ == 0 ) {
        tmp_segments = ( size_t* )malloc( sizeof( size_t ) * TIMES2( newcap ));

    } else {
        tmp_segments = ( size_t* )realloc( segments, sizeof( size_t ) * TIMES2( newcap ));
    }

    if( !tmp_segments )
        throw myruntime_error( mystring( "SegmentStructure: Not enough memory." ));

    segments = tmp_segments;

    // fill uninitialized memory with zeros
    memset( segments + TIMES2( capacity_ ), 0, sizeof( size_t ) * TIMES2( newcap - capacity_ ));

    capacity_ = newcap;
}

// -------------------------------------------------------------------------
// Push: pushes segment with boundaries left and right into the structure
// -------------------------------------------------------------------------

bool SegmentStructure::Push( size_t left, size_t right )
{
    int     location = -1;

    if( right < left )
        throw myruntime_error( mystring( "SegmentStructure: Invalid segment." ));

    if( Find( left, right, &location )) {
        //intersecting segment found, adjust by expanding it if needed
        Expand( location, left, right );

        RemoveAdjacentRight( location );
        RemoveAdjacentLeft( location );

        return false;
    }

    if( location < 0 )
        return false;

    if( capacity_ < length_ + 1 ) {
        Realloc( capacity_ + capacity_ );
    }

    for( int n = length_++; n > location; n-- ) {
        SetLeftAt( n, GetLeftAt( n - 1 ));
        SetRightAt( n, GetRightAt( n - 1 ));
    }

    SetLeftAt( location, left );
    SetRightAt( location, right );

    return true;
}

// -------------------------------------------------------------------------
// Find: finds a segment that intersects with the one given; if not found
//     then set the location of the nearest but greater segment is set
// -------------------------------------------------------------------------

bool SegmentStructure::Find( size_t segleft, size_t segright, int* loc ) const
{
    int     left = 0;
    int     right = ( int )GetSize() - 1;
    int     middle = 0;
    int     comp = 1;
    size_t  midleft;        //left value of the middle segment
    size_t  midright;       //right value of the middle segment

    while( left <= right )
    {
        middle = ( left + right ) >> 1;
        midleft = GetLeftAt( middle );
        midright = GetRightAt( middle );
        comp = Compare( midleft, midright, segleft, segright );

        if( comp < 0 )      //if segment is greater than the middle element
            left  = middle + 1;
        else if( comp > 0 ) //if segment is less than the middle element
            right = middle - 1;
        else {
            if( loc )
                *loc = middle;
            return true;
        }
    }
    if( loc )
        if( comp < 0 )
            *loc = left;
        else if( comp > 0 )
            *loc = middle;

    return false;
}

// -------------------------------------------------------------------------
// RemoveAdjacentLeft: removes equal segment at left
// -------------------------------------------------------------------------

void SegmentStructure::RemoveAdjacentLeft( size_t location )
{
    if( location == 0 || GetSize() <= location )
        return;

    size_t  prvleft = GetLeftAt( location - 1 );
    size_t  prvright = GetRightAt( location - 1 );
    size_t  locleft = GetLeftAt( location );
    size_t  locright = GetRightAt( location );

    if( Compare( prvleft, prvright, locleft, locright ) == 0 ) {
        Expand( location - 1, locleft, locright );
        Remove( location );
    }
}

// -------------------------------------------------------------------------
// RemoveAdjacentRight: removes equal segment at right
// -------------------------------------------------------------------------

void SegmentStructure::RemoveAdjacentRight( size_t location )
{
    if( GetSize() <= location + 1 )
        return;

    size_t  locleft = GetLeftAt( location );
    size_t  locright = GetRightAt( location );
    size_t  nxtleft = GetLeftAt( location + 1 );
    size_t  nxtright = GetRightAt( location + 1 );

    if( Compare( locleft, locright, nxtleft, nxtright ) == 0 ) {
        Expand( location, nxtleft, nxtright );
        Remove( location + 1 );
    }
}

// -------------------------------------------------------------------------
// Expand: expand segment at the location if needed
// -------------------------------------------------------------------------

void SegmentStructure::Expand( size_t location, size_t left, size_t right )
{
    if( GetSize() <= location )
        return;

    size_t  fndleft = GetLeftAt( location );
    size_t  fndright = GetRightAt( location );

    if( left < fndleft )
        SetLeftAt( location, left );
    if( fndright < right )
        SetRightAt( location, right );
}

// -------------------------------------------------------------------------
// Remove: removes segment at the given position
// -------------------------------------------------------------------------

void SegmentStructure::Remove( size_t loc )
{
#ifdef __DEBUG__
    if( GetSize() <= loc )
        throw myruntime_error( mystring( "SegmentStructure: Memory access error." ));
#endif
    for( ; loc < GetSize() - 1; loc++ ) {
        SetLeftAt( loc, GetLeftAt( loc + 1 ));
        SetRightAt( loc, GetRightAt( loc + 1 ));
    }

    length_--;
}

// =========================================================================
// Methods for testing
// Print: prints vector of segments
//

void SegmentStructure::Print( FILE* fp )
{
    if( fp == NULL )
        return;

    fprintf( fp, "\nSegmentStructure\n Segments:\n " );
    for( size_t n = 0; n < GetSize(); n++ )
        fprintf( fp, " [%d-%d]", GetLeftAt( n ), GetRightAt( n ));

    fprintf( fp, "\n\n" );
}

