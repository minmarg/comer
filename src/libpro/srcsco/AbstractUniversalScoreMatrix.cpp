/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rc.h"
#include "data.h"
#include "libpro/srcpro/FrequencyStore.h"
#include "AbstractUniversalScoreMatrix.h"



//fake scale factor used to indicate the termination
const double AbstractUniversalScoreMatrix::fake_scale_factor = -1.0;


// =========================================================================
// CLASS ScoreProbability
//
// constructor: initialization of values
// -------------------------------------------------------------------------

ScoreProbability::ScoreProbability( double sc, double p )
:   probb( p  )
{
    SetScores( sc );
}

// -------------------------------------------------------------------------
// constructor: default
//

ScoreProbability::ScoreProbability()
:   score( 0    ),
    double_score( 0.0 ),
    probb( 0.0  )
{
}

// destructor:
//
ScoreProbability::~ScoreProbability()
{
}


// =========================================================================
// CLASS AbstractUniversalScoreMatrix
//
// constructor: initialization
// -------------------------------------------------------------------------

AbstractUniversalScoreMatrix::AbstractUniversalScoreMatrix(
    const FrequencyStore*   store,
    TType           type,
    Configuration   config[NoSchemes],
    TBehaviour      beh,
    TScaling        scl,
    TMask           msk,
    TFVectorProbabilities   distrib,
    bool            cpu )
:
    AbstractScoreMatrix( type, config, beh, scl, msk ),
    freqstore( store ),
    prob_calculator( NULL ),
    distribution( distrib ),
    cpu_preference( cpu )
{
    if( !store || !freqstore->GetFrequencies())
            throw myruntime_error( mystring ( "AbstractUniversalScoreMatrix: No frequency vectors provided." ));

    prob_calculator = new BinarySearchStructure( ScoreProbability::ScoreComparer, GetMaxElemsOfProbCalculator());

    if( prob_calculator == NULL )
            throw myruntime_error( mystring ( "AbstractUniversalScoreMatrix: Not enough memory." ));
}

// -------------------------------------------------------------------------
// default constructor is invalid for initialization
// -------------------------------------------------------------------------

AbstractUniversalScoreMatrix::AbstractUniversalScoreMatrix()
:
    AbstractScoreMatrix(),

    freqstore( NULL ),
    prob_calculator( NULL ),
    distribution( DISCRETE ),
    cpu_preference( true )
{
}

// -------------------------------------------------------------------------
// destructor: deallocation of resources
// -------------------------------------------------------------------------

AbstractUniversalScoreMatrix::~AbstractUniversalScoreMatrix()
{
    DestroyProbCalculator( prob_calculator );
}

// -------------------------------------------------------------------------
// ClearProbCalculator: destroys members of ProbCalculator and
//     re-initializes it
// -------------------------------------------------------------------------

void AbstractUniversalScoreMatrix::ClearProbCalculator( BinarySearchStructure* pc )
{
    if( pc ) {
        for( size_t n = 0; n < pc->GetSize(); n++ )
        {
            const ScoreProbability* scp = ( const ScoreProbability* )( pc->GetValueAt( n ));
            ScoreProbability::Destroy( scp );
        }
        pc->Clear();
    }
}

// -------------------------------------------------------------------------
// DestroyProbCalculator: destroys totally ProbCalculator
// -------------------------------------------------------------------------

void AbstractUniversalScoreMatrix::DestroyProbCalculator( BinarySearchStructure* pc )
{
    if( pc ) {
        ClearProbCalculator( pc );
        delete pc;
    }
}

// -------------------------------------------------------------------------
// ComputePositionalScoreProbs: should not be called 
// -------------------------------------------------------------------------

void AbstractUniversalScoreMatrix::ComputePositionalScoreProbs( AttributableScores* )
{
    USM_THROW( "AbstractUniversalScoreMatrix: "
               "ComputePositionalScoreProbs should not be called from the class descendants." );
}

// -------------------------------------------------------------------------
// ComputeScoreProbabilities: computes probabilities to observe scores at
//     each position (i,j). The implementation depends on the cpu-mem
//     preference
// -------------------------------------------------------------------------

void AbstractUniversalScoreMatrix::ComputeScoreProbabilities( AttributableScores* PATTR_SCORES )
{
    if( PreferCPU()) {
        ComputeScoreProbabilitiesCPU( PATTR_SCORES );
    } else {
        ComputeScoreProbabilitiesMem( PATTR_SCORES );
    }
}

// -------------------------------------------------------------------------
//  Push: pushes probability structure to the private store of
//     probabilities; 
//     returns true if structure has been inserted and false if a
//     probability for the same score is observed in the private store;
//     in the latter case a new structure is not inserted but probability
//     value is adjusted instead
// -------------------------------------------------------------------------

bool AbstractUniversalScoreMatrix::Push( const ScoreProbability* scp )
{
    //get storage for scores and probabilities...
    BinarySearchStructure*  loc_probabilities = GetProbCalculator();

    if( !loc_probabilities || !scp )
        USM_THROW( "AbstractUniversalScoreMatrix: Unable to save probabilities." );

    int     location = -1;
    //push score and probability to the ordered structure
    bool    inserted = loc_probabilities->Push( scp, &location );

    if( inserted ) {
        //score with probability has been inserted
    } else {
        //the pair has not been inserted: score has already been observed,
        //update probability
        ScoreProbability*   found;

        if( location < 0 ||
            ( found = ( ScoreProbability* )loc_probabilities->GetValueAt( location ) ) == NULL )
            USM_THROW( "AbstractUniversalScoreMatrix: Fatal: Unable to access probability structure." );

        //increase probability for this particular score value
        found->IncProbability( scp->GetProbability());
    }

    return inserted;
}

// -------------------------------------------------------------------------
//  ProcessScoreProbabilities: copies probabilities stored in
//     prob_calculator, performs final computations related to
//     probabilities, and computes expected score per position pair
//
//  This method allocates required space for probabilities!
// -------------------------------------------------------------------------

void AbstractUniversalScoreMatrix::ProcessScoreProbabilities( AttributableScores* PATTR_SCORES, double norm )
{
    if( PATTR_SCORES == NULL )
        throw myruntime_error( mystring( "AbstractUniversalScoreMatrix: Unable to compute probabilities: Wrong argument." ));

    double          avgscore = 0.0;     //average score
    double          sumprobs = 0.0;     //sum of probabilities
    size_t          no_scores = 0;      //number of distinct scores
    static double   normterm = 0.0;     //normalizing term for probabilities

    if( norm )
        normterm = norm;

    //get storage for scores and probabilities...
    BinarySearchStructure*  loc_probabilities = GetProbCalculator();

#ifdef __DEBUG__
    if( !loc_probabilities || normterm <= 0.0 )
        USM_THROW( "AbstractUniversalScoreMatrix: Unable to process probabilities." );
#endif

    no_scores = loc_probabilities->GetSize();

    if( !no_scores )
        USM_THROW( "AbstractUniversalScoreMatrix: No scores to compute probabilities for." );

    PATTR_SCORES->SetMinMaxScores(
        //minimum score is at the begining; maximum score is at the end
        (( const ScoreProbability* )loc_probabilities->GetValueAt( 0 ))->GetScore(),
        (( const ScoreProbability* )loc_probabilities->GetValueAt( no_scores - 1 ))->GetScore()
    );

    //re-write probabilities to the newly allocated vector
    for( int sc = PATTR_SCORES->GetMinScore(), loc_sc = 0; sc <= PATTR_SCORES->GetMaxScore(); sc++ )
    {
        const ScoreProbability* scp = ( const ScoreProbability* )( loc_probabilities->GetValueAt( loc_sc ));
        //ensure the scores are equal
        if( scp->GetScore() != sc )
            continue;
        else
            loc_sc++;

        double  prob = scp->GetProbability();
        //destroy the ScoreProbability structure, we don't need it any more
        ScoreProbability::Destroy( scp );

        //at first probability is zero; we increase it for score sc
        PATTR_SCORES->IncProbabilityOf( prob, sc );
        PATTR_SCORES->DivideProbabilityOf( normterm, sc );    //normalizing term
        avgscore += sc * PATTR_SCORES->GetProbabilityOf( sc );
        sumprobs += PATTR_SCORES->GetProbabilityOf( sc );
    }

    loc_probabilities->Clear();

// fprintf( stderr, "sumprobs=%f\n", sumprobs );
    if( sumprobs < 0.999 || sumprobs > 1.001 )
        USM_THROW( "AbstractUniversalScoreMatrix: Score probabilities are incorrect." );

    PATTR_SCORES->SetExpectedScore( avgscore );
}

// -------------------------------------------------------------------------
//  FormatTermMessage: commpose termination message using the scale factor
//     format; returns number of bytes written to the string stream
//
//  NOTE: it is assumed that space enough to keep formatted data is
//     pre-allocated
// -------------------------------------------------------------------------

size_t AbstractUniversalScoreMatrix::FormatTermMessage( char* strstream ) const
{
    return FormatScaleFactor( strstream, fake_scale_factor );
}

// -------------------------------------------------------------------------
//  IsThatTermMessage: verifies whether the message is termination message
// -------------------------------------------------------------------------

bool AbstractUniversalScoreMatrix::IsThatTermMessage( const char* strstream ) const
{
    if( !strstream )
        USM_THROW( "AbstractUniversalScoreMatrix: Wrong argument." );

    const char* p = strstream;
    int         crc = 0;    //value of simple crc computation
    double      tmpval = 0.0;

    if( *( int* )p != ScaleHeader )
        return false;

    crc += *( int* )p; p += sizeof( int );

    tmpval = *( double* )p; crc += (( int* )p )[0]; crc += (( int* )p )[1]; p += sizeof( double );

    if( tmpval != fake_scale_factor )
        return false;

    return true;
}

// -------------------------------------------------------------------------
//  FormatScaleFactor: format message with a scale facor in it;
//     returns number of bytes written to the string stream
//
//  NOTE: it is assumed that space enough to keep formatted data is
//     pre-allocated
// -------------------------------------------------------------------------

size_t AbstractUniversalScoreMatrix::FormatScaleFactor( char* strstream ) const
{
    return FormatScaleFactor( strstream, GetMultiplier());
}

// -------------------------------------------------------------------------
//  FormatScaleFactor: format message with a scale facor in it;
//     returns number of bytes written to the string stream
//
//  NOTE: it is assumed that space enough to keep formatted data is
//     pre-allocated
// -------------------------------------------------------------------------

size_t AbstractUniversalScoreMatrix::FormatScaleFactor( char* strstream, double value ) const
{
    char*   p = strstream;
    size_t  written = 0;
    int     crc = 0;  //value of simple crc computation

    if( !strstream )
        USM_THROW( "AbstractUniversalScoreMatrix: Unable to format scale factor." );

    *( int* )p = ScaleHeader;
    crc += *( int* )p; p += sizeof( int ); written += sizeof( int );

    *( double* )p = value;
    crc += (( int* )p )[0]; crc += (( int* )p )[1]; p += sizeof( double ); written += sizeof( double );

    return written;
}

// -------------------------------------------------------------------------
//  DeformatScaleFactor: deformat message to extract scale factor in it
//
//  NOTE: scale factor will acquire a new value
// -------------------------------------------------------------------------

void AbstractUniversalScoreMatrix::DeformatScaleFactor( const char* strstream )
{
    if( !strstream )
        USM_THROW( "AbstractUniversalScoreMatrix: Unable to extract scale factor." );

    const char* p = strstream;
    int         crc = 0;    //value of simple crc computation
    double      tmpval = 0.0;

    if( *( int* )p != ScaleHeader )
        USM_THROW( "AbstractUniversalScoreMatrix: Format of data received for scale factor is invalid." );

    crc += *( int* )p; p += sizeof( int );

    tmpval = *( double* )p; crc += (( int* )p )[0]; crc += (( int* )p )[1]; p += sizeof( double );

    SetMultiplier( tmpval );
}

// -------------------------------------------------------------------------
//  FormatProbCalculator: format score probabilities and write them to the
//     given string stream;
//     returns number of bytes written to the string stream
//
//  NOTE: it is assumed that space enough to keep formatted data is
//     pre-allocated
// -------------------------------------------------------------------------

size_t AbstractUniversalScoreMatrix::FormatProbCalculator( char* strstream ) const
{
    //get storage for scores and probabilities...
    const BinarySearchStructure*    loc_probabilities = GetProbCalculator();

    if( !loc_probabilities || !strstream )
        USM_THROW( "AbstractUniversalScoreMatrix: Unable to format probabilities." );

    size_t  no_scores = loc_probabilities->GetSize();

    if( !no_scores || GetMaxElemsOfProbCalculator() < no_scores )
        USM_THROW( "AbstractUniversalScoreMatrix: Range of scores to be formatted is invalid." );

    char*   p = strstream;
    size_t  written = 0;
    int     crc = 0;  //value of simple crc computation

    *( int* )p = ProbsHeader;   crc += *( int* )p;      p += sizeof( int );      written += sizeof( int );
    *( size_t* )p = no_scores;  crc += *( size_t* )p;   p += sizeof( size_t );   written += sizeof( size_t );

    for( size_t n = 0; n < no_scores; n++ ) {
        //write values: score of double and probability
        const ScoreProbability* scp = ( const ScoreProbability* )loc_probabilities->GetValueAt( n );

        *( double* )p = scp->GetDoubleScore();
        for( int i = 0; i < sizeof( double ); i++ ) crc += p[i];
        p += sizeof( double ); written += sizeof( double );

        *( double* )p = scp->GetProbability();
        for( int i = 0; i < sizeof( double ); i++ ) crc += p[i];
        p += sizeof( double ); written += sizeof( double );
    }

    *( int* )p = crc; written += sizeof( int );

//     PrintStringinHex( strstream, written );

    return written;
}

// -------------------------------------------------------------------------
//  DeformatProbCalculator: deformat score probabilities and create
//     auxiliary binary search structure
//
//  NOTE: the method allocates new amount of memory and returns newly
//     allocated resources
// -------------------------------------------------------------------------

BinarySearchStructure* AbstractUniversalScoreMatrix::DeformatProbCalculator( const char* strstream )
{
    if( !strstream )
        USM_THROW( "AbstractUniversalScoreMatrix: Unable to extract probabilities." );

    //create local storage for scores and probabilities to be returned
    BinarySearchStructure*  loc_probabilities =
            new BinarySearchStructure( ScoreProbability::ScoreComparer, GetMaxElemsOfProbCalculator());

    if( !loc_probabilities )
        USM_THROW( "AbstractUniversalScoreMatrix: Not enough memory." );

    const char* p = strstream;
    size_t      bread = 0;
    int         crc = 0;    //value of simple crc computation
    double      sc = 0.0;   //score
    double      prob = 0.0; //probability

    if( *( int* )p != ProbsHeader ) {
        DestroyProbCalculator( loc_probabilities );
        USM_THROW( "AbstractUniversalScoreMatrix: Format of data received is invalid." );
    }

    crc += *( int* )p; p += sizeof( int ); bread += sizeof( int );

    size_t  no_scores = *( size_t* )p; crc += *( size_t* )p; p += sizeof( size_t ); bread += sizeof( size_t );

    if( GetMaxElemsOfProbCalculator() < no_scores ) {
        DestroyProbCalculator( loc_probabilities );
        USM_THROW( "AbstractUniversalScoreMatrix: Range of scores received is invalid." );
    }

    for( size_t n = 0; n < no_scores; n++ ) {
        //write values: score of double and probability

        for( int i = 0; i < sizeof( double ); i++ ) crc += p[i];
        sc   = *( double* )p; p += sizeof( double ); bread += sizeof( double );

        for( int i = 0; i < sizeof( double ); i++ ) crc += p[i];
        prob = *( double* )p; p += sizeof( double ); bread += sizeof( double );

        ScoreProbability* scp = new ScoreProbability( sc, prob );
        if( !scp ) {
            DestroyProbCalculator( loc_probabilities );
            USM_THROW( "AbstractUniversalScoreMatrix: Not enough memory." );
        }

        int     location = -1;
        bool    inserted = loc_probabilities->Push( scp, &location );

        if( !inserted ) {
            DestroyProbCalculator( loc_probabilities );
            USM_THROW( "AbstractUniversalScoreMatrix: Duplicate scores received." );
        }
    }

    bread += sizeof( int );
//     PrintStringinHex( strstream, bread );

    if( *( int* )p != crc ) {
        DestroyProbCalculator( loc_probabilities );
        USM_THROW( "AbstractUniversalScoreMatrix: CRC of the received message is invalid." );
    }

    return loc_probabilities;
}

// -------------------------------------------------------------------------
//  InfuseProbCalculator: infuse probabilities with new ones; this method
//     changes contents of private member prob_calculator;
//
//     returns flag of all scores being zero
// -------------------------------------------------------------------------

bool AbstractUniversalScoreMatrix::InfuseProbCalculator(
    const BinarySearchStructure* probs_obtained,
    double* psum,
    bool delete_duplicates )
{
    bool    all_negat = true;   //all negative scores

    if( !probs_obtained )
        USM_THROW( "AbstractUniversalScoreMatrix: Unable to infuse with probabilities received." );

    size_t  no_scores_probs_obtained = probs_obtained->GetSize();

    if( psum )
        *psum = 0.0;

    for( size_t n = 0; n < no_scores_probs_obtained; n++ ) {
        const ScoreProbability* scp = ( const ScoreProbability* )probs_obtained->GetValueAt( n );

        if( psum && scp )
            *psum += scp->GetProbability();

        if( ! Push( scp ))
            //if hasn't  been pushed, destroy if needed
            if( delete_duplicates )
                delete scp;

        if( all_negat && 0 < scp->GetScore())
            all_negat = false;
    }

    return all_negat;
}

// =========================================================================
// PRINT ROUTINES
//
// PrintStringinHex: printing of string in hexadecimal format
//
void AbstractUniversalScoreMatrix::PrintStringinHex( const char* str, size_t len ) const
{
    if( !str )
        return;

    const int   cszbuf = 5 * KBYTE;
    char        buffer[cszbuf] = {0};
    char*       p = buffer;
    int         szsym = 3;

    for( size_t n = 0; n < len && strlen( buffer ) + szsym < cszbuf; n++ ) {
        sprintf( p, "%02X ", ( unsigned char )str[n] );
        p = buffer + strlen( buffer );
    }

    *--p = 0;
    message( buffer );
}

// -------------------------------------------------------------------------
// PrintProbCalculator: print prob_calculator structure (scores and
//     probabilities)
// -------------------------------------------------------------------------

void AbstractUniversalScoreMatrix::PrintProbCalculator( FILE* fp )
{
    if( fp == NULL )
        return;

    //get storage for scores and probabilities...
    const BinarySearchStructure*    loc_probabilities = GetProbCalculator();

    if( !loc_probabilities )
        return;

    size_t  no_scores = loc_probabilities->GetSize();

    fprintf( fp,"\n%5c Score probabilities\n", 32 );

    fprintf( fp, "%9c", 32 );

    for( int n = 0; n < no_scores; n++ ) {
        const ScoreProbability* scp = ( const ScoreProbability* )loc_probabilities->GetValueAt( n );
        fprintf( fp, "\n%5d %6.4g", scp->GetScore(), scp->GetProbability());
    }
    fprintf( fp, "\n" );
}

// =========================================================================
// SERIALIZATION ROUTINES
//
// Serialize: serialization of the statistical parameters computed
//
void AbstractUniversalScoreMatrix::Serialize( Serializer& ) const
{
/*    serializer.Write(( char* )&columns, sizeof( columns ), 1 );

    for( int n = 0; n < columns; n++ )
        serializer.Write(( char* )values[n], sizeof( double ), NUMALPH );

    if( columns > 0 )
        serializer.Write( aacids, sizeof( char ), columns );*/
}

// Deserialize: deserialization of the statistical parameters
// NOTE: the method changes values of the class's own parameters
//
void AbstractUniversalScoreMatrix::Deserialize( Serializer& )
{
}

