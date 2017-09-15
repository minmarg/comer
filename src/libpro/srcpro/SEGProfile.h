/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __SEGProfile__
#define __SEGProfile__

#include <stdio.h>
#include "libpro/srcpro/GapScheme.h"
#include "libpro/srcpro/DistributionMatrix.h"
#include "libseg/SEGAbstract.h"

#define DEFAULT_SEGPRO_WIN_LENGTH        12
#define DEFAULT_SEGPRO_LOW_ENTROPY      2.2
#define DEFAULT_SEGPRO_HIGH_ENTROPY     2.5
#define DEFAULT_SEGPRO_MAX_DIFFERENCE   100
#define DEFAULT_SEGPRO_PRINT_WIDTH       60

// #define DEFAULT_SEGPRO_EQUALITY_DISTANCE    6.480740698407860231
#define DEFAULT_SEGPRO_EQUALITY_DISTANCE    12.96

typedef ssize_t     provval_t;

// _________________________________________________________________________
// Class SEGProfile
//

class SEGProfile: public SEGAbstract
{
    enum {
        No_frequencies          = NUMALPH
    };
    enum {
        start_of_residue        = 0,
        start_of_frequencies    = start_of_residue + 1,
        start_of_weight         = start_of_frequencies + No_frequencies,
        Overall_size            = start_of_weight + 1
    };

public:
    SEGProfile(
        const FrequencyMatrix&  freq,
        const LogOddsMatrix&    logo,
        size_t          winlen = DEFAULT_SEGPRO_WIN_LENGTH,
        double          lowent = DEFAULT_SEGPRO_LOW_ENTROPY,
        double          highent = DEFAULT_SEGPRO_HIGH_ENTROPY,
        size_t          maxdiff = DEFAULT_SEGPRO_MAX_DIFFERENCE
    );

    SEGProfile(
        const FrequencyMatrix&  freq1,
        const FrequencyMatrix&  freq2,
        size_t          winlen = DEFAULT_SEGPRO_WIN_LENGTH,
        double          lowent = DEFAULT_SEGPRO_LOW_ENTROPY,
        double          highent = DEFAULT_SEGPRO_HIGH_ENTROPY,
        size_t          maxdiff = DEFAULT_SEGPRO_MAX_DIFFERENCE
    );

    virtual ~SEGProfile();

    static double   GetDistance()                   { return equalitydist_; }
    static double   GetDistanceTo2()                { return distsquared_; }

    static void     SetDistance( double value )     { distsquared_ = SQUARE( equalitydist_ = value ); }

    void            MaskSeggedPositions( FrequencyMatrix&, LogOddsMatrix&, GapScheme& ) const;

    void            PrintSequence( FILE*, size_t width = DEFAULT_SEGPRO_PRINT_WIDTH );
    void            PrintSeggedSequence( FILE*, size_t width = DEFAULT_SEGPRO_PRINT_WIDTH );

protected:
    static bool     ProValidator( const void* );
    static bool     ProVerifier( const void* );
    static int      ProComparer( const void* , const void* );
    static bool     ProEquality( const void* , const void* );

protected:
    size_t          GetLocalLength() const  { return length_; }
    provval_t**     GetAddresses() const    { return addresses; }

    provval_t       GetResidueAt( size_t pos ) const;
    void            SetResidueAt( size_t pos, provval_t value );

    provval_t       GetFrequencyAt( size_t pos, size_t r ) const;
    void            SetFrequencyAt( size_t pos, size_t r, provval_t value );

    provval_t       GetFrequencyWeightAt( size_t pos ) const;
    void            SetFrequencyWeightAt( size_t pos, provval_t value );

    static provval_t    GetResidue( const provval_t* );
    static provval_t    GetFrequency( const provval_t*, size_t r );
    static provval_t    GetFrequencyWeight( const provval_t* );

    static size_t   GetProResidue( const void* pvectaddr );

    void            AllocateAddresses( size_t newlen );
    void            Destroy();


    static size_t   GetSeqAlphabetSize()    { return sc_sizeproalphabet_; }
    static size_t   GetVectorSize()         { return Overall_size; }

    void            Translate( const FrequencyMatrix& freq, const LogOddsMatrix& logo );
    void            Translate( const FrequencyMatrix& freq1, const FrequencyMatrix& freq2 );

private:
    provval_t**     addresses;
    size_t          length_;

    static double   equalitydist_;      //threshold of Euclidean distance between two profile vectors
    static double   distsquared_;       //equalitydist_ squared

    static const size_t sc_sizeproalphabet_;
};

////////////////////////////////////////////////////////////////////////////
// Class SEGProfile INLINES
//
// ProValidator: returns an indicator of vector to be valid

inline
bool SEGProfile::ProValidator( const void* pvectaddr )
{
#ifdef __DEBUG__
    if( !pvectaddr )
        throw myruntime_error( mystring( "SEGProfile: ProValidator: Memory access error." ));
#endif
    provval_t   res = GetResidue( *( const provval_t** )pvectaddr );
    return res != GAP  &&  res != ASTERISK;
}

// -------------------------------------------------------------------------
// ProVerifier: returns an indicator of vector to be considered

inline
bool SEGProfile::ProVerifier( const void* pvectaddr )
{
#ifdef __DEBUG__
    if( !pvectaddr )
        throw myruntime_error( mystring( "SEGProfile: ProVerifier: Memory access error." ));
#endif
    provval_t   res = GetResidue( *( const provval_t** )pvectaddr );
    return res != GAP  &&  res != ASTERISK  &&  res != X;
}

// =========================================================================
// GetProResidue: extracts residue given profile vector address
//
inline
size_t SEGProfile::GetProResidue( const void* pvectaddr )
{
#ifdef __DEBUG__
    if( !pvectaddr )
        throw myruntime_error( mystring( "SEGProfile: GetProResidue: Memory access error." ));
#endif
    return ( size_t )GetResidue( *( const provval_t** )pvectaddr );
}

// GetResidue: accesses residue element in profile vector 
//
inline
provval_t SEGProfile::GetResidue( const provval_t* vectaddr )
{
#ifdef __DEBUG__
    if( !vectaddr )
        throw myruntime_error( mystring( "SEGProfile: Memory access error." ));
#endif
    return vectaddr[start_of_residue];
}

// GetFrequency: accesses frequency value of residue type in profile vector
//
inline
provval_t SEGProfile::GetFrequency( const provval_t* vectaddr, size_t r )
{
#ifdef __DEBUG__
    if( !vectaddr || GetVectorSize() <= r + start_of_frequencies )
        throw myruntime_error( mystring( "SEGProfile: Memory access error." ));
#endif
    return vectaddr[ r + start_of_frequencies ];
}

// GetFrequencyWeight: accesses frequency weight value in profile vector
//
inline
provval_t SEGProfile::GetFrequencyWeight( const provval_t* vectaddr )
{
#ifdef __DEBUG__
    if( !vectaddr )
        throw myruntime_error( mystring( "SEGProfile: Memory access error." ));
#endif
    return vectaddr[start_of_weight];
}

// =========================================================================
// GetResidueAt: returns residue at the position specified

inline
provval_t SEGProfile::GetResidueAt( size_t pos ) const
{
#ifdef __DEBUG__
    if( !addresses || length_ <= pos )
        throw myruntime_error( mystring( "SEGProfile: Memory access error." ));
#endif
    return addresses[pos][start_of_residue];
}

// -------------------------------------------------------------------------
// SetResidueAt: sets residue at the profile position specified

inline
void SEGProfile::SetResidueAt( size_t pos, provval_t value )
{
#ifdef __DEBUG__
    if( !addresses || length_ <= pos )
        throw myruntime_error( mystring( "SEGProfile: Memory access error." ));
#endif
    addresses[pos][start_of_residue] = value;
}


// -------------------------------------------------------------------------
// GetFrequencyAt: returns frequency value at the position and for residue
//     number specified

inline
provval_t SEGProfile::GetFrequencyAt( size_t pos, size_t r ) const
{
#ifdef __DEBUG__
    if( !addresses || length_ <= pos )
        throw myruntime_error( mystring( "SEGProfile: Memory access error." ));

    if( !addresses[pos] || GetVectorSize() <= r + start_of_frequencies )
        throw myruntime_error( mystring( "SEGProfile: Memory access error." ));
#endif
    return addresses[pos][ r + start_of_frequencies ];
}

// -------------------------------------------------------------------------
// SetFrequencyAt: sets frequency value at the profile position and residue
//     number specified

inline
void SEGProfile::SetFrequencyAt( size_t pos, size_t r, provval_t value )
{
#ifdef __DEBUG__
    if( !addresses || length_ <= pos )
        throw myruntime_error( mystring( "SEGProfile: Memory access error." ));

    if( !addresses[pos] || GetVectorSize() <= r + start_of_frequencies )
        throw myruntime_error( mystring( "SEGProfile: Memory access error." ));
#endif
    addresses[pos][ r + start_of_frequencies ] = value;
}


// -------------------------------------------------------------------------
// GetFrequencyWeightAt: returns weight of frequencies at the position

inline
provval_t SEGProfile::GetFrequencyWeightAt( size_t pos ) const
{
#ifdef __DEBUG__
    if( !addresses || length_ <= pos )
        throw myruntime_error( mystring( "SEGProfile: Memory access error." ));
#endif
    return addresses[pos][start_of_weight];
}

// -------------------------------------------------------------------------
// SetFrequencyWeightAt: sets weight of frequencies at the position

inline
void SEGProfile::SetFrequencyWeightAt( size_t pos, provval_t value )
{
#ifdef __DEBUG__
    if( !addresses || length_ <= pos )
        throw myruntime_error( mystring( "SEGProfile: Memory access error." ));
#endif
    addresses[pos][start_of_weight] = value;
}



#endif//__SEGProfile__
