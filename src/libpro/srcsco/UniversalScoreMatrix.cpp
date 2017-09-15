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
#include "UniversalScoreMatrix.h"


static const FrequencyMatrix    usm_dummyfreq;
static const LogOddsMatrix      usm_dummylogo;

void UniversalScoreMatrix::ResetSbjctFreq() { sbjctfreq = &usm_dummyfreq; }
void UniversalScoreMatrix::ResetSbjctLogo() { sbjctpssm = &usm_dummylogo; }

// =========================================================================
// CLASS UniversalScoreMatrix
//
// constructor: frequency matrices, log-odds matrices for the first and
//   second profiles are given with the parameters
// -------------------------------------------------------------------------

UniversalScoreMatrix::UniversalScoreMatrix(
    const FrequencyMatrix&  freq,
    const LogOddsMatrix&    pssm,
    const FrequencyStore*   store,
    double                  infrm_threshold,
    int                     thick_number,
    double                  thick_percnt,
    Configuration           config[NoSchemes],
    TBehaviour              beh,
    TScaling                a_scaling,
    TMask                   c_masking,
    TFVectorProbabilities   distrib,
    bool                    cpu )
:
    AbstractUniversalScoreMatrix(
            store,
            Universal,
            config,
            beh,
            a_scaling,
            c_masking,
            distrib,
            cpu ),

    queryfreq(  freq ),
    querypssm(  pssm ),

    sbjctfreq( &usm_dummyfreq ),
    sbjctpssm( &usm_dummylogo ),

    thickness_number( thick_number ),
    thickness_percnt( thick_percnt ),

    qspair_scores( NULL ),
    qspair_mask( NULL ),
    sbjct_length( 0 ),
    sbjct_reserved( 0 ),

    queryprob( NULL )
{
    if( !querypssm.IsCompatible( queryfreq ))
            USM_THROW( "Query profile corrupted." );

    if(  queryfreq.GetColumns() < 1 || queryfreq.GetColumns() > MAXCOLUMNS )
            USM_THROW( "Query profile has wrong number of positions." );


    size_t  no_freqs = GetStore()->GetFrequencies()->GetSize();

    if( PreferCPU()) {
        Init( 0, 0 );   //zeros to indicate that no large matrix in memory will be used
    } else {
        //usually this allocates a large amount of memory;
        //there should be enough of that if avoiding more intense computation
        Init( queryfreq.GetColumns(), no_freqs );
    }

    AllocateQueryProbabilities();
    PrecomputeQueryProbabilities();
    SetInformationThreshold( infrm_threshold );
}

// -------------------------------------------------------------------------
// default constructor is invalid for initialization
// -------------------------------------------------------------------------

UniversalScoreMatrix::UniversalScoreMatrix()
:
    AbstractUniversalScoreMatrix(),

    queryfreq( usm_dummyfreq ),
    querypssm( usm_dummylogo ),

    sbjctfreq( &usm_dummyfreq ),
    sbjctpssm( &usm_dummylogo ),

    thickness_number( 0 ),
    thickness_percnt( 0.0 ),

    qspair_scores( NULL ),
    qspair_mask( NULL ),
    sbjct_length( 0 ),
    sbjct_reserved( 0 ),

    queryprob( NULL )
{
}

// -------------------------------------------------------------------------
// destructor: deallocation of resources
// -------------------------------------------------------------------------

UniversalScoreMatrix::~UniversalScoreMatrix()
{
    DestroyPairScores();
    DestroyQueryProbabilities();
}

// -------------------------------------------------------------------------
// ReallocatePairScores: reallocates space for query-subject pair scores;
//     if sizes are greater than allocated space previously, needed space is
//     reallocated, otherwise no action is taken
// -------------------------------------------------------------------------

void UniversalScoreMatrix::ReallocatePairScores( int s_size )
{
    if( s_size < 1 || s_size > MAXCOLUMNS )
        USM_THROW( "Number of subject profile positions is wrong." );

    if( s_size <= sbjct_reserved ) {
        //we have allocated enough memory for subject of this size
        sbjct_length = s_size;
        return;
    }

    if( qspair_scores == NULL ) {
        qspair_scores = ( double** )malloc( sizeof( double* ) * queryfreq.GetColumns());

        if( qspair_scores == NULL )
            USM_THROW( "UniversalScoreMatrix: Not enough memory." );

        memset( qspair_scores, 0, sizeof( double* ) * queryfreq.GetColumns());
    }

    if( qspair_mask == NULL ) {
        qspair_mask = ( TMask** )malloc( sizeof( TMask* ) * queryfreq.GetColumns());

        if( qspair_mask == NULL )
            USM_THROW( "UniversalScoreMatrix: Not enough memory." );

        memset( qspair_mask, 0, sizeof( TMask* ) * queryfreq.GetColumns());
    }



    for( int n = 0; n < queryfreq.GetColumns(); n++ ) {
        if( qspair_scores[n] )  free( qspair_scores[n] );
        if( qspair_mask[n] )    free( qspair_mask[n] );

        qspair_scores[n] = ( double* )malloc( sizeof( double ) * s_size );
        qspair_mask[n]   = ( TMask* )malloc( sizeof( TMask ) * s_size );

        if( qspair_scores[n] == NULL || qspair_mask[n] == NULL )
            USM_THROW( "UniversalScoreMatrix: Not enough memory." );

        memset( qspair_scores[n], 0, sizeof( double ) * s_size );
        memset( qspair_mask[n], 0, sizeof( TMask ) * s_size );
    }

    sbjct_reserved = s_size;
    sbjct_length = s_size;
}

// -------------------------------------------------------------------------
// DestroyPairScores: deallocates memory allocated previously for
//     query-subject pair scores
// -------------------------------------------------------------------------

void UniversalScoreMatrix::DestroyPairScores()
{
    if( qspair_scores ) {
        for( int n = 0; n < queryfreq.GetColumns(); n++ )
            if( qspair_scores[n] )
                free( qspair_scores[n] );

        free( qspair_scores );
        qspair_scores = NULL;
    }

    if( qspair_mask ) {
        for( int n = 0; n < queryfreq.GetColumns(); n++ )
            if( qspair_mask[n] )
                free( qspair_mask[n] );

        free( qspair_mask );
        qspair_mask = NULL;
    }
}

// -------------------------------------------------------------------------
// AllocateQueryProbabilities: pre-allocate memory for query probabilities
// -------------------------------------------------------------------------

void UniversalScoreMatrix::AllocateQueryProbabilities()
{
    DestroyQueryProbabilities();

    queryprob = ( double* )malloc( sizeof( double ) * GetQuerySize());

    if( queryprob == NULL )
        USM_THROW( "UniversalScoreMatrix: Not enough memory." );

    memset( queryprob, 0, sizeof( double ) * GetQuerySize());
}

// -------------------------------------------------------------------------
// DestroyQueryProbabilities: destroys precomputed query probabilities
// -------------------------------------------------------------------------

void UniversalScoreMatrix::DestroyQueryProbabilities()
{
    if( queryprob ) {
        free( queryprob );
        queryprob = NULL;
    }
}

// -------------------------------------------------------------------------
// PrecomputeQueryProbabilities: precomputes probabilities at each
//     position of query
// -------------------------------------------------------------------------

void UniversalScoreMatrix::PrecomputeQueryProbabilities()
{
    if( ! GetQueryProbabilities()) 
        USM_THROW( "UniversalScoreMatrix: Unable to precompute query probabilities." );

    double  sum = 0.0;

    for( int n = 0; n < GetQuerySize(); n++ )
    {
        sum = 0.0;

        for( int r = 0; r < NUMALPH; r++ ) {
            if( LOSCORES.PROBABility( r ) <= 0.0 )
                continue;
            sum += GetQueryFreq()( n, r ) * LOSCORES.LogPROBABility( r );
        }

        SetQueryProbabilityAt( n, exp( sum ));
    }
}

// -------------------------------------------------------------------------
// ComputeScore: compute score given query position and subject position
//     assumed subject is valid
// -------------------------------------------------------------------------

double UniversalScoreMatrix::ComputeScore( size_t querypos, size_t sbjctpos ) const
{
    double  score = 0.0;

    for( int a = 0; a < NUMAA; a++ ) {
        score +=    GetSbjctLogo()( sbjctpos, a ) * queryfreq( querypos, a ) +
                    GetSbjctFreq()( sbjctpos, a ) * querypssm( querypos, a );
    }

    return score * GetMultiplier();
}

// -------------------------------------------------------------------------
// ComputeScore: compute score given query position and charactersitic
//     vector of frequencies
// -------------------------------------------------------------------------

double UniversalScoreMatrix::ComputeScore( size_t querypos, const FrequencyVector& vect ) const
{
    double  score = 0.0;


    for( int a = 0; a < NUMAA; a++ ) {
        //computed as frequencies of query multiplied by scores of subject plus
        //  frequencies of subject multiplied by scores of query
        //lambda is implicitly incorporated into the log-odds values of the profiles
        score +=    ( double )vect.GetScoreAt( a ) / SCALE_CONSTANT * queryfreq( querypos, a ) +
                    ( double )vect.GetValueAt( a ) / FREQUENCY_SUM  * querypssm( querypos, a );
    }

    return score * GetMultiplier();
}

// -------------------------------------------------------------------------
// VectorScoreProbability: compute score probability given query position
//     and charactersitic vector of frequencies; computation depends on
//     vector distribution type chosen
// -------------------------------------------------------------------------

double UniversalScoreMatrix::VectorScoreProbability( size_t querypos, const FrequencyVector& colvect ) const
{
    //in any case we explicitly have profile to compare database profiles with;
    //so according to probability theory, each position occurs really
    switch( GetDistributionType()) {
        case DISCRETE:
            //return column vector observed frequency
            return colvect.GetProbability();

        case PROVECTOR:
//             break;
        case MULTINOMIAL:
            //return column vector multinomial probability
            return colvect.GetProbability();

        default:
            throw myruntime_error( mystring( "UniversalScoreMatrix: Unknown vector distribution type." ));
    }

    //compute probability given profile vector distribution type;
    //this is a product of two vector probabilities
    return GetQueryProbabilityAt( querypos ) * colvect.GetProbability();
}

// -------------------------------------------------------------------------
// ComputeProfileScoringMatrix: computes scoring matrix that is to be used
//     for aligning two profiles
// -------------------------------------------------------------------------

void UniversalScoreMatrix::ComputeProfileScoringMatrix( bool final )
{
    if( PreferCPU())
        //this method should be called when one is confident about the amount of memory required
        return;


    if( !IsValid() || !GetStore() || !GetStore()->GetFrequencies())
        USM_THROW( "Failed to compute matrix of scores." );

    int                 no_elems = 0;
    int                 negatives = 0;
    TMask               currentmask = Unmasked;
    bool                all_negat = true;
    double              score = 0.0;
    const double        scoreX = -1.0;  //score of profile positions with representative Xs
    const SimpleVector* frequencies = GetStore()->GetFrequencies();

    //fill matrix with values
    for( int n = 0; n < GetPrivateQuerySize(); n++ ) {
        no_elems = 0;
        negatives = 0;

        for( int m = 0; m < GetPrivateSecndSize(); m++ )
        {
            if( GetQueryFreq().GetResidueAt( n ) == X ) {
                SetScore( scoreX, m, n );
                continue;
            }

            const FrequencyVector   vector(( const char* )frequencies->GetValueAt( m ));

            currentmask = Unmasked;

            if( false )
            {
                SetMasked( currentmask = GetMaskingApproach(), m, n );
            }

            if( currentmask != MaskToIgnore )
            {
                score = ComputeScore( n, vector );
                no_elems++;
                if( all_negat && score < 0.0 ) negatives++;   //update counter for each negative score found
            }
            else
                //give penalty to score at the positions;
                score = scoreX;

            SetScore( score, m, n );
        }

        if( all_negat && negatives != no_elems )
            all_negat = false;
    }

    SetAllNegatives( all_negat );
}

// -------------------------------------------------------------------------
// PreserveSubject: prepares scoring system for alignment of query and
//     subject; subject profile is given by the arguments to this method;
// NOTE: method is invoked AFTER necessary scaling has been performed
// -------------------------------------------------------------------------

void UniversalScoreMatrix::PreserveSubject( const FrequencyMatrix& sfreq, const LogOddsMatrix& spssm )
{
    if( ! PreferCPU())
        //we must create private score matrix in any case, thus
        //proceed computations ignoring this flag and pre-computing
        //scores to preserve effectiveness of alignment process
        ;

#ifdef __DEBUG__
    if( !GetStore())
        USM_THROW( "UniversalScoreMatrix: Unable to construct private score matrix." );
#endif

    const double        scoreX = -1.0;      //score of profile positions with representative Xs
    bool                all_negat = true;   //all negative scores
    FrequencyVector     vector( NULL );

    ReallocatePairScores( sfreq.GetColumns());

    SetSbjctFreq( sfreq );
    SetSbjctLogo( spssm );

//     if( GetCorresScores()) GetCorresScores()->Init( GetQuerySize(), GetSubjectSize());
//     if( GetScaledScores()) GetScaledScores()->Init( GetQuerySize(), GetSubjectSize());

    for( int n = 0; n < GetQuerySize(); n++ ) {
        for( int m = 0; m < GetSubjectSize(); m++ )
        {
            //do not try to find a vector where subject's corresponding residue is X
            if( GetQueryFreq().GetResidueAt( n ) == X ||
                GetSbjctFreq().GetResidueAt( m ) == X )
            {
                PreservePairScore( scoreX, n, m );
                continue;
            }

            //change the vector with values of other column
            vector.SetVector(   *GetSbjctFreq().GetVectorAt( m ),
                                 GetSbjctLogo().GetFrequencyWeightAt( m ),
                                 GetSbjctLogo().GetMIDExpNoObservationsAt( m, PS_M ),
                                *GetSbjctLogo().GetVectorAt( m ),
                                 GetSbjctLogo().GetInformationAt( m )
            );

            const char* found = ( const char* )GetStore()->Find( vector );

            if( found == NULL ) {
                vector.Destroy();
                USM_THROW( "UniversalScoreMatrix: No frequency vector found in the database." );
            }

            double                  score = 0.0;
            const FrequencyVector   quest_vector( found );

            try {
                if( GetQueryLogo().GetInformationAt( n ) < GetInformationThreshold() ||
                    GetSbjctLogo().GetInformationAt( m ) < GetInformationThreshold() ||
                    ! PairThicknessConstraintsMet( m, n ))
                {
                    //give score penalty at the positions and...
                    score = -1.0;
                    //masking approach does not play any role here since statistics has been already computed 
                    SetPairMasked( GetMaskingApproach(), n, m );
                }
                else {
                    score = ComputeScore( n, quest_vector );

                    if( all_negat && 0 < score )
                        all_negat = false;
                }

                PreservePairScore( score, n, m );

            } catch( myexception const& ex ) {
                vector.Destroy();
                USM_THROW( ex.what(), ex.eclass());
            }
        }
    }

    vector.Destroy();

    ResetSbjctFreq();
    ResetSbjctLogo();

    //if all scores are negative throw, profile searching continues with the
    //next profile by omitting this one
    if( all_negat ) {
        warning( "UniversalScoreMatrix: All scores are negative." );
        throw myruntime_error( mystring( "Unable to align: Negative scores." ));
    }
}

// -------------------------------------------------------------------------
// ComputeScoreProbabilitiesCPU: computes probabilities to observe scores at
//     each position (i,j). Probability is computed as
//
//      1     _
//    ----   \   Pfj
//   length  /_ 
//           i,j:
//         sij=sk
//
//  where Pfj are probabilities of frequency vector fj to occur; sij are
//  scores at (i,j), sk is a discrete score value; length is for query.
//  Method iterates through the large table of scores only once saving
//  probability distribution on the fly
//
//  This method allocates required space for probabilities!
// -------------------------------------------------------------------------

bool UniversalScoreMatrix::ComputeScoreProbabilitiesCPU( AttributableScores* PATTR_SCORES )
{
    size_t  length = 0;         //effective length of query
    double  probsum = 0.0;      //sum of all score probabilities
    bool    all_negat = true;   //all negative scores

    int     loc_multiplier = PATTR_SCORES->GetAutoScalingFactor();

    //get storage for scores and probabilities...
    BinarySearchStructure*  loc_probabilities = GetProbCalculator();

#ifdef __DEBUG__
    if( !GetStore() || !GetStore()->GetFrequencies() ||
        !loc_probabilities )
        USM_THROW( "UniversalScoreMatrix: Unable to compute probabilities." );
#endif

    if( loc_probabilities->GetSize())
        USM_THROW( "UniversalScoreMatrix: Error in computation of probabilities." );

    const SimpleVector*     frequencies = GetStore()->GetFrequencies();
    size_t                  no_colms = GetQuerySize();
    size_t                  no_freqs = frequencies->GetSize();

    double                  score = 0.0;
    double                  scoreprob = 0.0;
    ScoreProbability*       scp = ScoreProbability::NewScoreProbability();

    for( size_t n = 0; n < no_colms; n++ ) {
        if( GetQueryFreq().GetResidueAt( n ) == X )
            continue;

        if( GetMaskingApproach() == MaskToIgnore )
            //if positions with information content less than the threshold
            //should be excluded from statistics, continue without it
            if( GetQueryLogo().GetInformationAt( n ) < GetInformationThreshold())
                continue;

        length++;
        for( size_t m = 0; m < no_freqs; m++ )
        {
            const   FrequencyVector   vector(( const char* )frequencies->GetValueAt( m ));

            if( GetMaskingApproach() == MaskToIgnore )
                //the same...
                if(( double ) vector.GetInfContent() / INFO_SCALE_CONSTANT < GetInformationThreshold())
                    continue;

            score = ComputeScore( n, vector );
            scoreprob = VectorScoreProbability( n, vector );

            scp->SetScores( score * loc_multiplier );
            scp->SetProbability( scoreprob );
            probsum += scoreprob;

            if( all_negat && 0 < scp->GetScore())
                all_negat = false;

            if( scp->GetScore() <= SCORE_MIN || scp->GetProbability() <= 0.0 )
                continue;

            if( Push( scp )) {
                //score with probability has been inserted: construct new object
                scp = ScoreProbability::NewScoreProbability();
            }
        }
    }

    ScoreProbability::Destroy( scp );

    SetAllNegatives( all_negat );
    ProcessScoreProbabilities( PATTR_SCORES, probsum/*length*/ );

    return true;
}

// -------------------------------------------------------------------------
// ComputeScoreProbabilitiesMem: computes probabilities to observe scores at
//     each position (i,j). Probability computation is as above:
//     ComputeScoreProbabilitiesCPU
//
//  NOTE: This method allocates required space for probabilities!
// -------------------------------------------------------------------------

bool UniversalScoreMatrix::ComputeScoreProbabilitiesMem( AttributableScores* )
{
    USM_THROW( "UniversalScoreMatrix: ComputeScoreProbabilitiesMem's not implemented yet." );
    return false;
}

// =========================================================================
// PRINT ROUTINES
//
// PrintParameterTable: nothing to do since output of the parameter
//     table is in the final
//

void UniversalScoreMatrix::PrintParameterTable( TPrintFunction, void* ) const
{
}

// PrintFinal: print table of computed statistical parameter
//     values to a stream
//

void UniversalScoreMatrix::PrintFinal( TPrintFunction print_func, void* vpn ) const
{
    if( vpn == NULL )
        return;

    if( 0.0 <= GetExpectedScore()) {
        print_func( vpn, "Expected score per position is non-negative, %.4f!\n\n", GetExpectedScore());
        return;
    }

    char    _sk[BUF_MAX], _slambda[BUF_MAX];
    char    _sdk[BUF_MAX], _sdlambda[BUF_MAX];
    if( 0.0 < GetK())   sprintf( _sk, "%6.4f", GetK());
        else            sprintf( _sk, "n/a" );
    if( 0.0 < GetLambda())  sprintf( _slambda, "%6.4f", GetLambda());
        else                sprintf( _slambda, "n/a" );
    if( 0.0 < GetDerivedGappedK())  sprintf( _sdk, "%6.4f", GetDerivedGappedK());
        else                        sprintf( _sdk, "n/a" );
    if( 0.0 < GetDerivedGappedLambda()) sprintf( _sdlambda, "%6.4f", GetDerivedGappedLambda());
        else                            sprintf( _sdlambda, "n/a" );

    print_func( vpn, "%-25s  %-6s   %-6s\n", " ", "K", "Lambda" );
    print_func( vpn, "%-25s  %6.4f   %6.4f\n", "Reference ungapped,", GetRefK(),            GetRefLambda());
    print_func( vpn, "%-25s  %6.4f   %6.4f\n", "Reference gapped,",   GetExpGappedK(),      GetExpGappedLambda());

    if( GetBehaviour() != StatisticsGiven ) {
        print_func( vpn, "%-25s  %6s   %6s\n", "Computed  ungapped,", _sk, _slambda );
        print_func( vpn, "%-25s  %6s   %6s\n", "Estimated gapped,", _sdk, _sdlambda );
        print_func( vpn, "Entropy, %6.4f; Expected, %6.4f\n\n", GetEntropy(), GetExpectedScore());
    } else
        print_func( vpn, "Entropy, %6.4f\n\n", GetEntropy());
}

// -------------------------------------------------------------------------
// PrintScoringMatrix: output the computed scoring system
// -------------------------------------------------------------------------

void UniversalScoreMatrix::PrintScoringMatrix( FILE* fp )
{
    if( !GetStore() || !GetStore()->GetFrequencies())
        return;

    if( fp == NULL )
        return;

    fprintf( fp,"%12c Position-specific profile scoring matrix\n", 32 );

    fprintf( fp, "%9c", 32 );

    const SimpleVector*     frequencies = GetStore()->GetFrequencies();
    size_t                  no_colms = GetQuerySize();
    size_t                  no_freqs = frequencies->GetSize();
    int                     l = 0;

    for( size_t m = 0; m < no_freqs; m++ )
        fprintf( fp, "%4d", m );


    for( size_t n = 0; n < no_colms; n++ ) {
        fprintf( fp, "\n%5d %c   ", ++l, DehashCode( GetQueryFreq().GetResidueAt( n ) ));

        for( size_t m = 0; m < no_freqs; m++ )
        {
            const FrequencyVector   vector(( const char* )frequencies->GetValueAt( m ));
            fprintf( fp, "%3d ", ( int )rint( ComputeScore( n, vector )));
        }
    }

    fprintf( fp, "\n" );
}

