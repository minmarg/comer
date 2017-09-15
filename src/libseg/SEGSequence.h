/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __SEGSequence__
#define __SEGSequence__

#include <stdio.h>
#include "SEGAbstract.h"

#define DEFAULT_SEGSEQ_WIN_LENGTH       12
#define DEFAULT_SEGSEQ_LOW_ENTROPY      2.2
#define DEFAULT_SEGSEQ_HIGH_ENTROPY     2.5
#define DEFAULT_SEGSEQ_MAX_DIFFERENCE   100

// _________________________________________________________________________
// Class SEGSequence
//

class SEGSequence: public SEGAbstract
{
public:
    SEGSequence(
        const unsigned char*    residues,
        size_t          seqlength,
        bool            hashed = false,
        size_t          winlen = DEFAULT_SEGSEQ_WIN_LENGTH,
        double          lowent = DEFAULT_SEGSEQ_LOW_ENTROPY,
        double          highent = DEFAULT_SEGSEQ_HIGH_ENTROPY,
        size_t          maxdiff = DEFAULT_SEGSEQ_MAX_DIFFERENCE
    );

    virtual ~SEGSequence();

    void            PrintSequence( FILE*, size_t width = DEFAULT_SEGSEQ_PRINT_WIDTH );
    void            PrintSeggedSequence( FILE*, size_t width = DEFAULT_SEGSEQ_PRINT_WIDTH );

    static bool     SeqValidator( const void* );
    static bool     SeqVerifier( const void* );
    static int      SeqComparer( const void* , const void* );

protected:
    size_t          GetLocalLength() const  { return length_; }

    void            Translate( const unsigned char*, size_t len );

    size_t          GetAddressAt( size_t pos ) const;
    void            SetAddressAt( size_t pos, size_t value );

    static size_t*  AllocateAddresses( size_t newlen );
    static void     Destroy( size_t* );

    static size_t   GetSeqAlphabetSize()    { return sc_sizealphabet_; }
    bool            GetHashed() const       { return hashed_; }

    static size_t   GetSeqResidue( const void* );

private:
    size_t*     addresses;
    size_t      length_;
    bool        hashed_;

    static const size_t sc_sizealphabet_;
};

////////////////////////////////////////////////////////////////////////////
// Class SEGSequence INLINES
//
// SeqValidator: returns an indicator of letter to be valid

inline
bool SEGSequence::SeqValidator( const void* letter )
{
#ifdef __DEBUG__
    if( !letter )
        throw myruntime_error( mystring( "SEGSequence: SeqValidator: Memory access error." ));
#endif
    return *( size_t* )letter != GAP  &&  *( size_t* )letter != ASTERISK;
}

// -------------------------------------------------------------------------
// SeqVerifier: returns an indicator of letter to be considered

inline
bool SEGSequence::SeqVerifier( const void* letter )
{
#ifdef __DEBUG__
    if( !letter )
        throw myruntime_error( mystring( "SEGSequence: SeqVerifier: Memory access error." ));
#endif
    return *( size_t* )letter < GetSeqAlphabetSize();
}

// -------------------------------------------------------------------------
// SeqComparer: returns 0 if two letters are equal, 1 if the first is
//     greater than the second, -1 otherwise
//
inline
int SEGSequence::SeqComparer( const void* one, const void* another )
{
#ifdef __DEBUG__
    if( !one || !another )
        throw myruntime_error( mystring( "SEGSequence: SeqComparer: Memory access error." ));
#endif
    return *( size_t* )one - *( size_t* )another;
}

// -------------------------------------------------------------------------
// GetAddressAt: returns address (translated letter) at the position

inline
size_t SEGSequence::GetAddressAt( size_t pos ) const
{
#ifdef __DEBUG__
    if( !addresses || length_ <= pos )
        throw myruntime_error( mystring( "SEGSequence: Memory access error." ));
#endif
    return addresses[pos];
}

// -------------------------------------------------------------------------
// SetAddressAt: sets address (translated letter) at the position

inline
void SEGSequence::SetAddressAt( size_t pos, size_t value )
{
#ifdef __DEBUG__
    if( !addresses || length_ <= pos )
        throw myruntime_error( mystring( "SEGSequence: Memory access error." ));
#endif
    addresses[pos] = value;
}

// -------------------------------------------------------------------------
// GetResidue: extracts residue given its address

inline
size_t SEGSequence::GetSeqResidue( const void* letter )
{
#ifdef __DEBUG__
    if( !letter )
        throw myruntime_error( mystring( "SEGSequence: Memory access error." ));
#endif
    return *( size_t* )letter;
}


#endif//__SEGSequence__
