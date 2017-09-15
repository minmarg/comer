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

#include "HashTable.h"



// Prime numbers...
const unsigned long HashTable::prime_list[numPrimes] = {
    53ul,           97ul,           193ul,          389ul,          769ul,      // 0-4 
    1543ul,         3079ul,         6151ul,         12289ul,        24593ul,    // 5-9 
    49157ul,        98317ul,        196613ul,       393241ul,       786433ul,   //10-14
    1572869ul,      3145739ul,      6291469ul,      12582917ul,     25165843ul, //15-19
    50331653ul,     100663319ul,    201326611ul,    402653189ul,    805306457ul,//20-24
    1610612741ul,   3221225473ul,   4294967291ul                                //25-27
};

// -------------------------------------------------------------------------
// constructor: hash table allocation
// -------------------------------------------------------------------------

HashTable::HashTable( THashFunction func, size_t size )
:   function( func ),
    values( NULL ),
    hashsize( 0 )
{
    Reserve( size );
}

// -------------------------------------------------------------------------
// constructor: default
// -------------------------------------------------------------------------

HashTable::HashTable()
:   function( NULL ),
    values( NULL ),
    hashsize( 0 )
{
    throw myruntime_error(
            mystring( "HashTable: Default initialization prohibited." ));
}

// -------------------------------------------------------------------------
// destructor:
// -------------------------------------------------------------------------

HashTable::~HashTable()
{
    if( values )
        free( values );
}

// -------------------------------------------------------------------------
// Reserve: memory allocation
// -------------------------------------------------------------------------

void HashTable::Reserve( size_t size )
{
    if( hashsize != 0 || values != NULL )
        throw myruntime_error( mystring( "HashTable: Trying to allocate space for allocated table." ));

    values = ( const void** )malloc( sizeof( void* ) * size );
    if( values == NULL )
        throw myruntime_error( mystring( "HashTable: Not enough memory." ));

    memset( values, 0, sizeof( void* ) * size );
    hashsize = size;
}

