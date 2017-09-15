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

#include "BinarySearchStructure.h"


// -------------------------------------------------------------------------
// constructor: initialization
// -------------------------------------------------------------------------

SimpleVector::SimpleVector( size_t size )
:   values( NULL ),
    length_( 0 ),
    capacity_( 0 )
{
    Realloc( size );
}

// -------------------------------------------------------------------------
// constructor: default
// -------------------------------------------------------------------------

SimpleVector::SimpleVector()
:   values( NULL ),
    length_( 0 ),
    capacity_( 0 )
{
    throw myruntime_error(
            mystring( "SimpleVector: Default initialization prohibited." ));
}

// -------------------------------------------------------------------------
// destructor:
// -------------------------------------------------------------------------

SimpleVector::~SimpleVector()
{
    if( values )
        free( values );
}

// -------------------------------------------------------------------------
// Realloc: memory allocation
// -------------------------------------------------------------------------

void SimpleVector::Realloc( size_t newcap )
{
    const void**    tmp_values;

    if( newcap <= capacity_ )
        return;

    if( capacity_ == 0 ) {
        tmp_values = ( const void** )malloc( sizeof( void* ) * newcap );

    } else {
        tmp_values = ( const void** )realloc( values, sizeof( void* ) * newcap );
    }

    if( !tmp_values )
        throw myruntime_error( mystring( "SimpleVector: Not enough memory." ));

    values = tmp_values;

    // fill uninitialized memory with zeros
    memset( values + capacity_, 0, sizeof( void* ) * ( newcap - capacity_ ));

    capacity_ = newcap;
}

// -------------------------------------------------------------------------
// Push: pushes back value in the vector; returns flag indicating whether or
//     not a key has been inserted
// -------------------------------------------------------------------------

bool SimpleVector::Push( const void* key, int* /*not used since an element will always be pushed*/ )
{
    if( capacity_ < length_ + 1 ) {
        Realloc(  capacity_ + capacity_ + 1 );
    }

    values[length_] = key;

    length_++;

    return true;
}

// -------------------------------------------------------------------------
// InsertValueAt: inserts value at the position by shifting elements at the
//     right to the right
// -------------------------------------------------------------------------

void SimpleVector::InsertValueAt( size_t loc, const void* key )
{
    if( capacity_ < length_ + 1 ) {
        Realloc(  capacity_ + capacity_ + 1 );
    }

    if( length_ < loc )
        throw myruntime_error( mystring( "SimpleVector: Unable to insert value." ));

    for( int n = length_; n > loc; n-- )
        values[n] = values[n-1];

    values[loc] = key;

    length_++;
}

// -------------------------------------------------------------------------
// Clear: clears all elements
// -------------------------------------------------------------------------

void SimpleVector::Clear()
{
    memset( values, 0, sizeof( void* ) * capacity_ );
    length_ = 0;
}

// /////////////////////////////////////////////////////////////////////////
// CLASS BinarySearchStructure
//
// constructor: initialization
// -------------------------------------------------------------------------

BinarySearchStructure::BinarySearchStructure( TComparator comp, size_t size, bool keep )
:   SimpleVector( size ),
    comparator( comp ),
    comparator1( NULL ),
    keep_duplicates( keep ),
    params_( NULL )
{
}


BinarySearchStructure::BinarySearchStructure( 
        TComparator1 comp, size_t size, bool keep, void* pars )
:   SimpleVector( size ),
    comparator( NULL ),
    comparator1( comp ),
    keep_duplicates( keep ),
    params_( pars )
{
}

// -------------------------------------------------------------------------
// constructor: default
// -------------------------------------------------------------------------

BinarySearchStructure::BinarySearchStructure()
:   comparator( NULL ),
    comparator1( NULL ),
    keep_duplicates( false ),
    params_( NULL )
{
    throw myruntime_error(
            mystring( "BinarySearchStructure: Default initialization prohibited." ));
}

// -------------------------------------------------------------------------
// destructor:
// -------------------------------------------------------------------------

BinarySearchStructure::~BinarySearchStructure()
{
}

// -------------------------------------------------------------------------
// Push: inserts value in the structure so that its sorted order is kept;
//     returns flag whether or not a key has been inserted; a key is not
//     inserted if a duplicate exists and the structure is not to keep them;
//     loc will be set to point to a duplicate element
// -------------------------------------------------------------------------

bool BinarySearchStructure::Push( const void* key, int* loc )
{
#ifdef __DEBUG__
    if( !comparator && !comparator1 )
        throw myruntime_error( mystring( "BinarySearchStructure: Unable to insert an element." ));
#endif

    int location;

    if( Find( key, &location ))
        //found key
        if( !KeepDuplicates()) {
            if( loc )
                *loc = location;
            return false; //nothing to do if not to keep duplicates
        }

    if( capacity_ < length_ + 1 ) {
        Realloc(  capacity_ + capacity_ + 1 );
    }

    for( int n = length_; n > location; n-- )
        values[n] = values[n-1];

    values[location] = key;

    length_++;

    if( loc )
        *loc = location;

    return true;
}

// -------------------------------------------------------------------------
// Find: find a key in the structure and set the location of the found key;
//     if no key is found then set the location of the nearest but greater
//     element
// -------------------------------------------------------------------------

bool BinarySearchStructure::Find( const void* key, int* loc ) const
{
    int     left = 0;
    int     right = ( int )GetSize() - 1;
    int     middle = 0;
    int     comp = 1;

    while( left <= right )
    {
        middle = ( left + right ) >> 1;
        if( comparator )
            comp = ( *comparator )( values[middle], key );
        else if( comparator1 )
            comp = ( *comparator1 )( values[middle], key, GetParams());
        else
            throw myruntime_error("BinarySearchStructure: Null comparator.");

        if( comp < 0 )      //if key is greater than the middle element
            left  = middle + 1;
        else if( comp > 0 ) //if key is less than the middle element
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

