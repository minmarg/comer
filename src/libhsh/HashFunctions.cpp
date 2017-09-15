/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#include <stdio.h>
#include <stdlib.h>

#include "HashFunctions.h"

// -------------------------------------------------------------------------
// oneatatime: method proposed by Jenkins 
// -------------------------------------------------------------------------

Uint32 oneatatime( const void *key, size_t length )
{
    Uint32                  hash = 0;
    const unsigned char*    k = ( const unsigned char* )key;

    for( size_t n = 0; n < length; n++ )
    {
        hash += k[n];
        hash += hash << 10;
        hash ^= hash >> 6;
    }

    hash += hash << 3;
    hash ^= hash >> 11;
    hash += hash << 15;

    return hash;
}

