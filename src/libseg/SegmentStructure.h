/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __SegmentStructure__
#define __SegmentStructure__

#include <stdio.h>

#include "debug.h"
#include "pcmath.h"
#include "mystring.h"
#include "myexcept.h"


// _________________________________________________________________________
// Class SegmentStructure
//

class SegmentStructure
{
public:
    SegmentStructure( size_t size );
    virtual ~SegmentStructure();

    size_t          GetSize() const         { return length_; }

    size_t          GetLeftAt( size_t loc ) const;                  //obtain left boundary value
    size_t          GetRightAt( size_t loc ) const;                 //obtain right boundary value

    bool            Push( size_t left, size_t right );              //save segment boundaries in structure
    void            Clear();                                        //clear all elements

    void            Print( FILE* );         //for testing

protected:
    explicit SegmentStructure();

    void            Realloc( size_t newcap );
    size_t          GetCapacity() const     { return capacity_; }
                                                                    //find segment in structure
    bool            Find( size_t left, size_t right, int* loc = NULL ) const;

    void            SetLeftAt( size_t loc, size_t left );           //set left boundary value
    void            SetRightAt( size_t loc, size_t right );         //set right boundary value

    void            Remove( size_t );                               //remove segment at the position
    void            RemoveAdjacentLeft( size_t );                   //remove equal segment at left
    void            RemoveAdjacentRight( size_t );                  //remove equal segment at right
    void            Expand( size_t loc, size_t left, size_t right );//expand segment at the location if needed

    static int      Compare( size_t oneleft, size_t oneright, size_t anoleft, size_t anoright );

protected:
    size_t*         segments;           //array of segments
    size_t          length_;            //length of series
    size_t          capacity_;          //current capacity

};

// INLINES ...

// -------------------------------------------------------------------------
// GetLeftAt: returns left boundary value at the position

inline
size_t SegmentStructure::GetLeftAt( size_t loc ) const
{
#ifdef __DEBUG__
    if( length_ <= loc )
        throw myruntime_error( mystring( "SegmentStructure: Memory access error." ));
#endif
    return segments[ TIMES2( loc ) ];
}

// -------------------------------------------------------------------------
// GetRightAt: returns right boundary value at the position

inline
size_t SegmentStructure::GetRightAt( size_t loc ) const
{
#ifdef __DEBUG__
    if( length_ <= loc )
        throw myruntime_error( mystring( "SegmentStructure: Memory access error." ));
#endif
    return segments[ TIMES2( loc ) + 1 ];
}

// -------------------------------------------------------------------------
// SetLeftAt: sets left boundary value at the position

inline
void SegmentStructure::SetLeftAt( size_t loc, size_t left )
{
#ifdef __DEBUG__
    if( length_ <= loc )
        throw myruntime_error( mystring( "SegmentStructure: Memory access error." ));
#endif
    segments[ TIMES2( loc ) ] = left;
}

// -------------------------------------------------------------------------
// SetRightAt: sets right boundary value at the position

inline
void SegmentStructure::SetRightAt( size_t loc, size_t right )
{
#ifdef __DEBUG__
    if( length_ <= loc )
        throw myruntime_error( mystring( "SegmentStructure: Memory access error." ));
#endif
    segments[ TIMES2( loc ) + 1 ] = right;
}

// -------------------------------------------------------------------------
// Compare: segment comparison function; returns 0 if two segments
//     intersect, 1 if one is greater another, -1 otherwise
// NOTE: assumes that boundaries are correct
// -------------------------------------------------------------------------

inline
int SegmentStructure::Compare( size_t oneleft, size_t oneright, size_t anoleft, size_t anoright )
{
    if( oneright < anoleft )
        return -1;

    if( anoright < oneleft )
        return 1;

    return 0;
}


#endif//__SegmentStructure__
