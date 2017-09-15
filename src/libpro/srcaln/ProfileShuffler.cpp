/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mystring.h"
#include "myexcept.h"
#include "libpro/srcpro/Database.h"
#include "libpro/srcpro/Serializer.h"
#include "ProfileShuffler.h"


//file extension for profiles in the text format
const char* ProfileShuffler::s_text_ext = "txt";

// -------------------------------------------------------------------------
// CLASS PositionProbabilities
//
// Constructor
//
PositionProbabilities::PositionProbabilities()
:   distFunction( NULL ),
    no_entries( 0 ),
    counter( 0 )
{
}

// Destructor
//
PositionProbabilities::~PositionProbabilities()
{
    Destroy();
}

// -------------------------------------------------------------------------
// Destroy: destroy entries
//
void PositionProbabilities::Destroy()
{
    if( distFunction )
        free( distFunction );
    no_entries = 0;
    counter = 0;
}

// -------------------------------------------------------------------------
// Allocate: allocates space for the probability and index entries
//
void PositionProbabilities::Allocate( size_t size )
{
    Destroy();
    distFunction = ( double* )malloc( sizeof( double ) * size );

    if( distFunction == NULL )
        throw myruntime_error( mystring( "PositionProbabilities: Not enough memory." ));

    no_entries = size;
}

// -------------------------------------------------------------------------
// GetIndex: gets an index of position vector given value of probability
//     distribution function by applying binary search
//
size_t PositionProbabilities::GetIndex( double prob )
{

    if( prob < 0.0 || 1.0 < prob )
        throw myruntime_error( mystring( "PositionProbabilities: Wrong probability argument." ));

    int     left = 0;
    int     right = ( int )GetCounter() - 1;
    int     middle = 0;

    double  lower = 0.0;
    double  upper = 0.0;

    while( left <= right )
    {
        middle = ( left + right ) >> 1;
        lower = ( 0 < middle ) ? GetDistFunctionAt( middle - 1 ) : 0.0;
        upper = GetDistFunctionAt( middle );

        if( upper < prob )      //if probability is greater than the upper bound of the interval
            left  = middle + 1;
        else if( lower > prob ) //if probability is less than the lower bound of the interval
            right = middle - 1;
        else
            return ( size_t )middle;
    }

    throw myruntime_error( mystring( "PositionProbabilities: Unable to process random selection." ));
    return 0;
}

// -------------------------------------------------------------------------
// Push: pushes values given the probability of the positional vector and
//     its index
//
void PositionProbabilities::Push( double prob )
{
    size_t  loc_counter = GetCounter();

    if( GetNoEntries() <= loc_counter )
        throw myruntime_error( mystring( "PositionProbabilities: Memory access error." ));

    double  pdfunction = ( loc_counter ? GetDistFunctionAt( loc_counter - 1 ) + prob: prob );

    SetDistFunctionAt( loc_counter, pdfunction );
    IncCounter();
}

// -------------------------------------------------------------------------
// Check: checks if probability distribution function is valid
//
bool PositionProbabilities::Check()
{
    if( GetCounter() &&
        GetCounter() <= GetNoEntries() &&
        GetDistFunctionAt( GetCounter() - 1 ) > 0.9999 &&
        GetDistFunctionAt( GetCounter() - 1 ) < 1.0001 )
        return true;//ok

    return false;
}

////////////////////////////////////////////////////////////////////////////
// CLASS ProfileShuffler
//
// Constructor
//
ProfileShuffler::ProfileShuffler(
        const char* database,
        const char* output,
        int noprofs,
        int length,
        int thickness,
        bool textformat,
        long seed
    )
:
    database_name( database ),
    pattern( output ),
    no_profiles( noprofs ),
    profile_length( length ),
    min_thickness( thickness ),
    print_text( textformat ),
    residues( NULL ),
    rng_seed( seed ),

    profile_db( database ),

    no_db_sequences( 0 ),
    db_size( 0 )
{
    InitProbs( residueprobs );

    if( rng_seed == -1 )
        rng_seed = -rand();
    if( 0 < rng_seed )
        rng_seed = -rng_seed;

    residues = ( char* )malloc( profile_length );

    if( !residues )
        throw myruntime_error( mystring( "ProfileShuffler: Not enough memory." ));

    memset( residues, 0, profile_length );
}

// Default constructor
//
ProfileShuffler::ProfileShuffler()
:
    database_name( NULL ),
    pattern( NULL ),
    no_profiles( 0 ),
    profile_length( 0 ),
    min_thickness( 0 ),
    print_text( false ),
    residues( NULL ),
    rng_seed( -1 ),
    profile_db(( const char* )NULL ),
    no_db_sequences( 0 ),
    db_size( 0 )
{
    throw myruntime_error(
            mystring( "ProfileShuffler: Default initialization is not allowed." ));
}

// Destructor
//
ProfileShuffler::~ProfileShuffler()
{
    if( residues )
        free( residues );
}

// -------------------------------------------------------------------------
// InitProbs: initialize interval of probabilities used to generate residues
//  probs, probability vector
//  gap_freq, gap frequency
// -------------------------------------------------------------------------

void ProfileShuffler::InitProbs( double* probs, double gap_freq )
{
    if( !probs )
        return;

    double  res_freq = 1.0 - gap_freq;  //frequency for all other residues
    double  recip_fr = 1.0 / res_freq;  //inverse of res_freq;
    int     n;

    probs[0] = 0.0;

    for( n = 1; n < NUMALPH + 1 && probs[n-1] < ( res_freq - 0.0001 ); n++ )
        probs[n] = ( probs[n-1] * recip_fr + LOSCORES.PROBABility( n-1 )) * res_freq;

    for( ; n < NUMALPH + 1; n++ )
        probs[n] = 1.0;
}

// -------------------------------------------------------------------------
// Run: start generating of profiles by shuffling and joining profile
//     positional vectors
// -------------------------------------------------------------------------

long ProfileShuffler::Run()
{
    if( !GetPattern())
        throw myruntime_error( mystring( "ProfileShuffler: Unable to start generating of profiles." ));

    const int   size = KBYTE;
    char        name[size];
    char        profname[BUF_MAX];
    char        profdesc[BUF_MAX];
    int         num = 0;
    int         len = strlen( GetPattern());
    const char* ext = GetPattern() + len - 1;

    if( size < len + 17 )
        throw myruntime_error( mystring( "ProfileShuffler: Pattern of file names given is too long." ));


    //read database header
    profile_db.Open();

    message( "Reading frequencies..." );
    //read the distinct frequency vectors
    profile_db.ReadInFrequencies();

    SetNoSequences( profile_db.GetNoSequences());
    SetDbSize( profile_db.GetDbSize());

    FillPositionProbabilities();

    message( "Generating profiles..." );

    Serializer      serializer;

    FrequencyMatrix frequencies;
    LogOddsMatrix   pssm;
    GapScheme       gaps;

    frequencies.Reserve( GetProfileLength());
    pssm.Reserve( GetProfileLength());
    gaps.Reserve( GetProfileLength());


    for( ; ext != GetPattern() && *ext != DIRSEP && *ext != '.'; ext-- );
    if( ext == GetPattern() || *ext == DIRSEP )
        ext =  GetPattern() + len;

    size_t  base = ext - GetPattern();
    strncpy( name, GetPattern(), base );


    for( int n = 0; n < GetNoProfiles(); n++, num++ )
    {
        sprintf( name + base, "_%d%s", num, ext );

        SimulateProfile( frequencies, pssm, gaps );

        sprintf( profname, "simulated_%ld", GetCurrentSeedValue());
        sprintf( profdesc, "serial %d, length %d", num, GetProfileLength());

        pssm.SetName( profname );
        pssm.SetDescription( profdesc );

//         serializer.SerializeProfile( frequencies, pssm, gaps, name );//obsolete
        serializer.WriteProfile( name, frequencies, pssm, gaps );

        sprintf( name + base, "_%d%s.%s", num, ext, s_text_ext );

        if( GetPrintText())
            OutputProfile( name, frequencies, pssm, gaps );
    }


    profile_db.Close();
    message( "Finished." );
    return rng_seed;
}

// -------------------------------------------------------------------------
// FillPositionProbabilities: fills position probabilities with values
// -------------------------------------------------------------------------

void ProfileShuffler::FillPositionProbabilities()
{
    if( profile_db.GetStore() == NULL )
        throw myruntime_error( mystring( "ProfileShuffler: Unable to compute position probability distribution function." ));
        
    const SimpleVector* frequencies = profile_db.GetStore()->GetFrequencies();

    if( frequencies == NULL )
        throw myruntime_error( mystring( "ProfileShuffler: No frequency vectors in the database." ));

    const size_t        size = frequencies->GetSize();

    if( size == 0 )
        throw myruntime_error( mystring( "ProfileShuffler: No frequency vectors in the database." ));


    positionProbs.Allocate( size );

    for( size_t n = 0; n < size; n++ ) {
        FrequencyVector vector(( const char* )frequencies->GetValueAt( n ));
        positionProbs.Push( vector.GetProbability());
    }

    //check probability distribution function for validness
    if( ! positionProbs.Check())
        throw myruntime_error( mystring( "ProfileShuffler: Probability distribution function is not valid." ));
}

// -------------------------------------------------------------------------
// SimulateProfile: simulate one profile by shuffling positional vectors
//     from the database
// -------------------------------------------------------------------------

void ProfileShuffler::SimulateProfile( FrequencyMatrix& freq, LogOddsMatrix& pssm, GapScheme& gaps )
{
    static double   posfreqvalues[NUMALPH];
    static double   pospssmvalues[NUMALPH];
    double          trnprobs[P_NSTATES];
    static double   gapweight = 0.0;
    static double   delbegwht = 0.0;
    static double   delendwht = 0.0;
    static int      deldelta = 0;
    static size_t   thick = 30; //thickness of the alignments profiles were constructed from (fake)
    double          prob;   //probability in (0,1)

    if( profile_db.GetStore() == NULL )
        throw myruntime_error( mystring( "ProfileShuffler: Unable to simulate profile." ));

    const SimpleVector* frequencies = profile_db.GetStore()->GetFrequencies();

    if( frequencies == NULL )
        throw myruntime_error( mystring( "ProfileShuffler: No frequency vectors in the database." ));


    memset( posfreqvalues, 0, sizeof( double ) * NUMALPH );
    memset( pospssmvalues, 0, sizeof( double ) * NUMALPH );
    memset( trnprobs, 0, sizeof( double ) * P_NSTATES );

    if( TRANSPROBS.GetPriors())
        memcpy( trnprobs, *TRANSPROBS.GetPriors(), sizeof( double ) * P_NSTATES );


    char*   loc_residues = GetResidues();

    if( !loc_residues )
        throw myruntime_error( mystring( "ProfileShuffler: No buffer for residues allocated." ));

    //generate residues at first
    for( int r = 0; r < GetProfileLength(); r++ )
        loc_residues[r] = GenerateResidue();

    pssm.SetNoSequences( thick );
    pssm.SetEffNoSequences( thick );

    //now generate profile consisting of existing positional columns
    //
    for( int n = 0; n < GetProfileLength(); n++ )
    {
        size_t  pos = positionProbs.GetIndex( prob = LecuyerRand( &rng_seed ) );
        const FrequencyVector   vector(( const char* )frequencies->GetValueAt( pos ));
        double  expMIDs[PS_NSTATES];

        expMIDs[PS_M] = ( double )vector.GetThickness();
        expMIDs[PS_I] = 0.0;
        expMIDs[PS_D] = 0.0;

        for( int it = 0; it < FrequencyVector::GetNoElems(); it++ ) {
            posfreqvalues[it] = ( double )vector.GetValueAt( it ) / FREQUENCY_SUM;
            pospssmvalues[it] = ( double )vector.GetScoreAt( it ) / SCALE_CONSTANT;
        }

        freq.PushAt( posfreqvalues, loc_residues[n], n );
        pssm.PushAt( pospssmvalues, loc_residues[n], vector.GetWeight(), vector.GetInfContent(), expMIDs, n );
        gaps.PushAt(&trnprobs, gapweight, delbegwht, delendwht, deldelta, loc_residues[n], n );
    }
}

// -------------------------------------------------------------------------
// GenerateResidue: generate single residue
// -------------------------------------------------------------------------

unsigned char ProfileShuffler::GenerateResidue()
{
    int             n;
    double          prob = LecuyerRand( &rng_seed );//probability in (0,1)
    const double*   prob_intervals = GetResidueProbs();

    for( n = 1; n < NUMALPH + 1; n++ )
        if( prob_intervals[n-1] <= prob && prob < prob_intervals[n] )
            break;

    return n - 1;
}
