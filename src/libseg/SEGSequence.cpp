/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <stdlib.h>
#include <string.h>

#include "SEGSequence.h"

//alphabet size for sequences
const size_t    SEGSequence::sc_sizealphabet_ = NUMAA;

// /////////////////////////////////////////////////////////////////////////
// CLASS SEGSequence
//
// constructor:
//
// NOTE: address is supposed to point at a vector of pointers void*
//

SEGSequence::SEGSequence(
    const unsigned char*    residues,
    size_t          seqlength,
    bool            hashed,
    size_t          winlen,
    double          lowent,
    double          highent,
    size_t          maxdiff )
:
    //ATTENTION!
    SEGAbstract(
        &SEGSequence::SeqValidator,
        &SEGSequence::SeqVerifier,
        &SEGSequence::SeqComparer,
        NULL,

        addresses = AllocateAddresses( seqlength ),
        seqlength,

        winlen,
        lowent,
        highent,
        maxdiff,
        GetSeqAlphabetSize() /*alphabet size*/
    ),
    length_( seqlength ),
    hashed_( hashed )
{
    if( addresses == NULL )
        throw myruntime_error( mystring( "SEGSequence: Not enough memory." ));

    Translate( residues, seqlength );
}

// destructor:
//

SEGSequence::~SEGSequence()
{
    Destroy( addresses );
}

// -------------------------------------------------------------------------
// AllocateAddresses: allocates memory required by addresses
// -------------------------------------------------------------------------

size_t* SEGSequence::AllocateAddresses( size_t newlen )
{
    size_t* localvector = ( size_t* )malloc( sizeof( size_t ) * newlen );

    if( localvector == NULL )
        return localvector;

    memset( localvector, 0, sizeof( size_t ) * newlen );

    return localvector;
}

// -------------------------------------------------------------------------
// Destroy: destroys addresses allocated previously
// -------------------------------------------------------------------------

void SEGSequence::Destroy( size_t* addrvector )
{
    if( addrvector )
        free( addrvector );
}

// -------------------------------------------------------------------------
// Translate: translates residue codes into vector of addresses suitable for
//     abstract SEG
// -------------------------------------------------------------------------

void SEGSequence::Translate( const unsigned char* residues, size_t len )
{
    if( GetLocalLength() < len )
        throw myruntime_error( mystring( "SEGSequence: Memory access error." ));

    for( size_t n = 0; n < len; n++ )
        if( GetHashed())
            SetAddressAt( n, residues[n] );
        else
            SetAddressAt( n, HashAlphSymbol( residues[n] ));
}

// -------------------------------------------------------------------------
// PrintSequence: prints formatted sequence
// -------------------------------------------------------------------------

void SEGSequence::PrintSequence( FILE* fp, size_t width )
{
    SEGAbstract::PrintSequence( fp, &SEGSequence::GetSeqResidue, width );
}

// -------------------------------------------------------------------------
// PrintSeggedSequence: prints formatted sequence with segments found by
//     running the algorithm and masked with Xs
// -------------------------------------------------------------------------

void SEGSequence::PrintSeggedSequence( FILE* fp, size_t width )
{
    SEGAbstract::PrintSeggedSequence( fp, &SEGSequence::GetSeqResidue, width );
}

