  /***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "libseg/segdata.h"
#include "SEGProfile.h"

//thresholds of Euclidean distance between two profile vectors
double      SEGProfile::equalitydist_ = DEFAULT_SEGPRO_EQUALITY_DISTANCE;
double      SEGProfile::distsquared_  = SQUARE( DEFAULT_SEGPRO_EQUALITY_DISTANCE );

//alphabet size for profiles (number of possible columns) --
// don't specify
const size_t    SEGProfile::sc_sizeproalphabet_ = 0;

// /////////////////////////////////////////////////////////////////////////
// CLASS SEGProfile
//
// constructor:
//

SEGProfile::SEGProfile(
    const FrequencyMatrix&  freq,
    const LogOddsMatrix&    logo,
    size_t          winlen,
    double          lowent,
    double          highent,
    size_t          maxdiff )
:
    SEGAbstract(
        &SEGProfile::ProValidator,
        &SEGProfile::ProVerifier,
        &SEGProfile::ProComparer,
        &SEGProfile::ProEquality,

        winlen,
        lowent,
        highent,
        maxdiff,
        GetSeqAlphabetSize() /*alphabet size*/
    ),

    addresses( NULL ),
    length_( 0 )
{
    AllocateAddresses( freq.GetColumns());
    Translate( freq, logo );

    if( addresses == NULL )
        throw myruntime_error( mystring( "SEGProfile: Not enough memory." ));

    SetRunAddress(( void* )GetAddresses(), GetLocalLength());
}

// Alternative constructor: Vectors are constructed as differences between
//     corresponding frequency vector values
//

SEGProfile::SEGProfile(
    const FrequencyMatrix&  freq1,
    const FrequencyMatrix&  freq2,
    size_t          winlen,
    double          lowent,
    double          highent,
    size_t          maxdiff )
:
    SEGAbstract(
        &SEGProfile::ProValidator,
        &SEGProfile::ProVerifier,
        &SEGProfile::ProComparer,
        &SEGProfile::ProEquality,

        winlen,
        lowent,
        highent,
        maxdiff,
        GetSeqAlphabetSize() /*alphabet size*/
    ),

    addresses( NULL ),
    length_( 0 )
{
    if( freq1.GetColumns() != freq2.GetColumns())
        throw myruntime_error( mystring( "SEGProfile: Sizes of frequencies are invalid." ));

    AllocateAddresses( freq1.GetColumns());
    Translate( freq1, freq2 );

    if( addresses == NULL )
        throw myruntime_error( mystring( "SEGProfile: Not enough memory." ));

    SetRunAddress(( void* )GetAddresses(), GetLocalLength());
    SetHighCSearch();
}

// destructor:
//

SEGProfile::~SEGProfile()
{
    Destroy();
}

// -------------------------------------------------------------------------
// AllocateAddresses: allocates memory required by addresses
// -------------------------------------------------------------------------

void SEGProfile::AllocateAddresses( size_t newlen )
{
    void*   p = NULL;
    addresses = ( provval_t** )malloc( sizeof( provval_t* ) * newlen );

    if( addresses == NULL )
        throw myruntime_error( mystring( "SEGProfile: Not enough memory." ));

    memset( addresses, 0, sizeof( provval_t* ) * newlen );

    length_ = newlen;

    for( size_t n = 0; n < length_; n++ ) {
        p = addresses[n] = ( provval_t* )malloc( sizeof( provval_t ) * GetVectorSize());

        if( p == NULL )
            throw myruntime_error( mystring( "SEGProfile: Not enough memory." ));

        memset( p, 0, sizeof( provval_t ) * GetVectorSize());
    }
}

// -------------------------------------------------------------------------
// Destroy: destroys addresses allocated previously
// -------------------------------------------------------------------------

void SEGProfile::Destroy()
{
    if( addresses ) {
        for( size_t n = 0; n < length_; n++ )
            free( addresses[n] );

        free( addresses );
        addresses = NULL;
        length_ = 0;
    }
}

// -------------------------------------------------------------------------
// Translate: translates profile columns into vector of addresses
//     suitable for abstract SEG
// -------------------------------------------------------------------------

void SEGProfile::Translate( const FrequencyMatrix& freq, const LogOddsMatrix& logo )
{
    if( GetLocalLength() < freq.GetColumns())
        throw myruntime_error( mystring( "SEGProfile: Memory access error." ));


    for( int n = 0; n < freq.GetColumns(); n++ )
    {
        SetResidueAt( n, ( provval_t )freq.GetResidueAt( n ));

        for( size_t r = 0; r < NUMALPH; r++ )
            SetFrequencyAt( n, r, ( provval_t )rint( FREQUENCY_SUM * freq.GetValueAt( n, r )));

        SetFrequencyWeightAt( n, ( provval_t )rint( logo.GetFrequencyWeightAt( n )));
    }
}

// -------------------------------------------------------------------------
// Translate: translates differences of frequency vectorsinto vector of
//     addresses suitable for abstract SEG
// -------------------------------------------------------------------------

void SEGProfile::Translate( const FrequencyMatrix& freq1, const FrequencyMatrix& freq2 )
{
    if( GetLocalLength() < freq1.GetColumns() || freq1.GetColumns() != freq2.GetColumns())
        throw myruntime_error( mystring( "SEGProfile: Memory access error." ));


    for( int n = 0; n < freq1.GetColumns(); n++ )
    {
        SetResidueAt( n, ( provval_t )freq1.GetResidueAt( n ));

        for( size_t r = 0; r < NUMALPH; r++ )
            SetFrequencyAt( n, r, ( provval_t )rint( FREQUENCY_SUM * ( freq1.GetValueAt( n, r ) - freq2.GetValueAt( n, r ))));

        SetFrequencyWeightAt( n, ( provval_t ) 0 );
    }
}

// -------------------------------------------------------------------------
// ProEquality: returns true if two vectors are supposed to be equal, false
//     otherwise
// -------------------------------------------------------------------------

bool SEGProfile::ProEquality( const void* pone, const void* panother )
{
#ifdef __DEBUG__
    if( !pone || !panother )
        throw myruntime_error( mystring( "SEGProfile: ProComparer: Memory access error." ));
#endif

    if( pone == panother )
        return true;

    double              sqdist = 0.0;   //squared distance >0
    const provval_t*    vectone = *( const provval_t** )pone;
    const provval_t*    vectanother = *( const provval_t** )panother;
    int                 lastdiff;       //last difference of components of vectors

    //compute squared euclidean distance
    for( size_t r = 0; r < No_frequencies; r++ ) {
        lastdiff = ( int )GetFrequency( vectone, r ) - ( int )GetFrequency( vectanother, r );
        sqdist += PRESQUARES.GetValueOf( abs( lastdiff ));
    }

    //take frequency weights into account as well
//     lastdiff = ( int )GetFrequencyWeight( vectone ) - ( int )GetFrequencyWeight( vectanother );
//     sqdist += PRESQUARES.GetValueOf( abs( lastdiff ));

    if( sqdist < GetDistanceTo2())
        return true;

    return false;
}

// -------------------------------------------------------------------------
// SeqComparer: returns 0 if serial number of two vectors coinsides,
//     pos. if the number of the first is greater, neg. otherwise
// -------------------------------------------------------------------------

int SEGProfile::ProComparer( const void* pone, const void* panother )
{
#ifdef __DEBUG__
    if( !pone || !panother )
        throw myruntime_error( mystring( "SEGProfile: ProComparer: Memory access error." ));
#endif

    return ( int )(( size_t )pone - ( size_t )panother );
//     return memcmp(  vectone + start_of_frequencies,
//                     vectanother + start_of_frequencies,
//                 ( No_frequencies + 1 ) * sizeof( provval_t ));
}

// -------------------------------------------------------------------------
// MaskSeggedPositions: masks profile's positions segged by the algorithm
// -------------------------------------------------------------------------

void SEGProfile::MaskSeggedPositions( FrequencyMatrix& freq, LogOddsMatrix& logo, GapScheme& gaps ) const
{
    size_t  n = 0;
    size_t  p = 0;
    size_t  r = 0;
    size_t  left, right;
    const SegmentStructure&     loc_segments = GetSegments();

    char    maskX = X;
    double  freqweight =  0.0;
    double  information = 0.0;
    double  expMIDs[PS_NSTATES] = {1.0,0.0,0.0};
    double  freqvalues[NUMALPH];
    double  logovalues[NUMALPH];
    double  trnprobs[P_NSTATES];

    int     efective_number = NUMAA;    // effective number of residues
    double  eqpart = rint(( double )FREQUENCY_SUM / efective_number ) / FREQUENCY_SUM;

    for( r = 0; r < efective_number; r++ ) {
        freqvalues[r] = LOSCORES.PROBABility( r );//eqpart;
        logovalues[r] = LOSCORES.PrecomputedEntry( maskX, r );
    }

    for( ; r < NUMALPH; r++ ) {
        freqvalues[r] = 0.0;
        logovalues[r] = LOSCORES.PrecomputedEntry( maskX, r );
    }

    memset( trnprobs, 0, sizeof( double ) * P_NSTATES );
    if( TRANSPROBS.GetPriors())
        memcpy( trnprobs, *TRANSPROBS.GetPriors(), sizeof( double ) * P_NSTATES );

    for( n = 0; n < loc_segments.GetSize(); n++ ) {
        left = loc_segments.GetLeftAt( n );
        right = loc_segments.GetRightAt( n );

        for( p = left; p <= right; p++ ) {
            freq.PushAt( freqvalues, maskX, p );
            logo.PushAt( logovalues, maskX, freqweight, information, expMIDs, p );
//             gaps.PushAt( gaps.GetWeightsAt( p ),
//                          gaps.GetDeletesBegAt( p ),
//                          gaps.GetDeletesEndAt( p ),
//                          gaps.GetDeletesIntervalAt( p ),
//                         maskX, p );
            gaps.PushAt( &trnprobs, 0.0, 0.0, 0.0, 0, maskX, p );
        }
    }
}

// -------------------------------------------------------------------------
// PrintSequence: prints formatted profile's sequence 
// -------------------------------------------------------------------------

void SEGProfile::PrintSequence( FILE* fp, size_t width )
{
    SEGAbstract::PrintSequence( fp, &SEGProfile::GetProResidue, width );
}

// -------------------------------------------------------------------------
// PrintSeggedSequence: prints formatted profile's sequence with segments
//     found by running the algorithm and masked with Xs
// -------------------------------------------------------------------------

void SEGProfile::PrintSeggedSequence( FILE* fp, size_t width )
{
    SEGAbstract::PrintSeggedSequence( fp, &SEGProfile::GetProResidue, width );
}

