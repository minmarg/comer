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
#include "AdjustedScoreMatrix.h"



// =========================================================================
// CLASS AdjustedScoreMatrix
//
// constructor: frequency matrices, log-odds matrices for the first and
//   second profiles are given along with the frequency pool (store)
// -------------------------------------------------------------------------

AdjustedScoreMatrix::AdjustedScoreMatrix(
    const FrequencyMatrix& freq_fst, const LogOddsMatrix& logo_fst,
    const FrequencyMatrix& freq_sec, const LogOddsMatrix& logo_sec,
    const FrequencyStore*   store,
    double                  infrm_threshold,
    int                     thick_number,
    double                  thick_percnt,
    double                  mask_percnt,
    Configuration           config[NoSchemes],
    TBehaviour              s_behaviour,
    TScaling                a_scaling,
    TMask                   c_masking )
:
    ScoringMatrix(
            freq_fst, logo_fst,
            freq_sec, logo_sec,
            infrm_threshold,
            thick_number,
            thick_percnt,
            mask_percnt,
            config,
            s_behaviour,
            a_scaling,
            c_masking,
            AdjustedProfileSpecific ),

    freqstore( store ),
    vector_probabilities( NULL )
{
    if( !store || !freqstore->GetFrequencies())
            throw myruntime_error( mystring ( "AdjustedScoreMatrix: No frequency vectors provided." ));

    vector_probabilities = ( double* )malloc( sizeof( double ) * GetSubjectSize());

    if( !vector_probabilities )
            throw myruntime_error( mystring ( "AdjustedScoreMatrix: Not enough memory." ));

    memset( vector_probabilities, 0, sizeof( double ) * GetSubjectSize());
}

// -------------------------------------------------------------------------
// default constructor is invalid for initialization
// -------------------------------------------------------------------------

AdjustedScoreMatrix::AdjustedScoreMatrix()
:
    ScoringMatrix(),
    freqstore( NULL ),
    vector_probabilities( NULL )
{
}

// -------------------------------------------------------------------------
// destructor: deallocation of resources
// -------------------------------------------------------------------------

AdjustedScoreMatrix::~AdjustedScoreMatrix()
{
    if( vector_probabilities )
        free( vector_probabilities );
}

// -------------------------------------------------------------------------
// ComputeScore: compute score given query position and charactersitic
//     vector of frequencies
// -------------------------------------------------------------------------

double AdjustedScoreMatrix::ComputeScore( size_t querypos, const FrequencyVector& vect ) const
{
    double  score = 0.0;


    for( int a = 0; a < NUMAA; a++ ) {
        //computed as frequencies of query multiplied by scores of subject plus
        //  frequencies of subject multiplied by scores of query
        //lambda is implicitly incorporated into the log-odds values of the profiles
        score +=    ( double )vect.GetScoreAt( a ) / SCALE_CONSTANT * GetQueryFreq()( querypos, a ) +
                    ( double )vect.GetValueAt( a ) / FREQUENCY_SUM  * GetQueryLogo()( querypos, a );
    }

    return score;
}

// -------------------------------------------------------------------------
// ComputeProfileScoringMatrix: computes scoring matrix that is to be used
//     for aligning two profiles
// -------------------------------------------------------------------------

void AdjustedScoreMatrix::ComputeProfileScoringMatrix( bool final )
{
    if( GetQuerySize() <= 0 || GetSubjectSize() <= 0 )
        throw myruntime_error( mystring( "AdjustedScoreMatrix: No matrix." ));

    if( !IsValid() || !GetStore() || !GetStore()->GetFrequencies())
        throw myruntime_error( mystring( "AdjustedScoreMatrix: Unable to compute scores." ));

    int                 no_elems = 0;
    int                 negatives = 0;
    TMask               currentmask = Unmasked;
    bool                all_negat = true;
    double              score = 0.0;
    const double        scoreX = -1.0;
    const SimpleVector* frequencies = GetStore()->GetFrequencies();
    const char*         found = NULL;

    FrequencyVector     vector( NULL );

    //fill matrix with values
    for( int m = 0; m < GetSubjectSize(); m++ )
    {
        //omit positions of Xs
        if( GetSbjctFreq().GetResidueAt( m ) == X ) {
            for( int ns = 0; ns < GetQuerySize(); ns++ )
                SetScore( scoreX, m, ns );
            continue;
        }

        no_elems = 0;
        negatives = 0;

        //initialize vector with values of the column
        vector.SetVector(   *GetSbjctFreq().GetVectorAt( m ),
                             GetSbjctLogo().GetFrequencyWeightAt( m ),
                             GetSbjctLogo().GetMIDExpNoObservationsAt( m, PS_M ),
                            *GetSbjctLogo().GetVectorAt( m ),
                             GetSbjctLogo().GetInformationAt( m )
        );

        found = ( const char* )GetStore()->Find( vector );

        if( found == NULL ) {
            vector.Destroy();
            throw myruntime_error( mystring( "AdjustedScoreMatrix: No frequency vector found in the database." ));
        }

        const FrequencyVector   quest_vector( found );
        SaveVectorProbabilityAt( m, quest_vector.GetProbability());


        for( int n = 0; n < GetQuerySize(); n++ )
        {
            score = 0.0;
            currentmask = Unmasked;

            try {
                //give constant penalties at the positions of X
                if( GetQueryFreq().GetResidueAt( n ) == X ||
                    GetSbjctFreq().GetResidueAt( m ) == X )
                {
                    SetScore( scoreX, m, n );
                    continue;
                }

                if( GetQueryLogo().GetInformationAt( n ) < GetInformationThreshold() ||
                    GetSbjctLogo().GetInformationAt( m ) < GetInformationThreshold() ||
                    ! ThicknessConstraintsMet( m, n ))
                {
                    //mask the score possibly excluding it from the ongoing computational statistics
                    SetMasked( currentmask = GetMaskingApproach(), m, n );
                }

                if( currentmask != MaskToIgnore )
                {
                    score = ComputeScore( n, quest_vector );

//                     for( int a = 0; a < NUMAA; a++ )
//                         // lambda is implicitly incorporated into the log-odds values of the profiles
//                         score += GetQueryFreq()( n, a ) * GetSbjctLogo()( m, a ) + GetSbjctFreq()( m, a ) * GetQueryLogo()( n, a );

                    no_elems++;
                    if( all_negat && score < 0.0 ) negatives++;   //increase counter for each negative score found
                }
                else
                    //give penalty to score at the positions;
                    score = scoreX;

                SetScore( score, m, n );

            } catch( myexception const& ex ) {
                vector.Destroy();
                throw myruntime_error( ex.what(), ex.eclass());
            }
        }

        if( all_negat && negatives != no_elems )
            all_negat = false;
    }

    vector.Destroy();

    SetAllNegatives( all_negat );
}

// -------------------------------------------------------------------------
// ComputePositionalScoreProbs: computes probabilities of scores at each
//     query and subject position
// -------------------------------------------------------------------------

void AdjustedScoreMatrix::ComputePositionalScoreProbs( AttributableScores* PATTR_SCORES )
{
    if( PATTR_SCORES == NULL )
        throw myruntime_error( mystring( "ScoringMatrix: Unable to compute probabilities: Wrong argument." ));

    double  prob = 0.0;         //vector probability
    double  normterm = 0.0;     //normalizing term
    double* querynorm = NULL;   //array of normalizing terms
    TScore  loc_score, sc;
    int     m, n;

    if( 0 < GetQuerySize()) {
        querynorm = ( double* )malloc( sizeof( double ) * GetQuerySize());
        if( querynorm == NULL )
            throw myruntime_error( mystring( "ScoringMatrix: ComputePositionalScoreProbs: Not enough memory." ));
        for( n = 0; n < GetQuerySize(); n++ )
            querynorm[n] = 0.0;
    }
    try {
        for( m = 0; m < GetSubjectSize(); m++ ) {
            if( GetSbjctFreq()[m] == X )
                continue;

            normterm = 0.0;
            for( n = 0; n < GetQuerySize(); n++ )
            {
                if( GetQueryFreq()[n] == X )
                    continue;

                if( PATTR_SCORES->GetMaskedToIgnore( m, n ))
                    continue;

                loc_score = ( TScore )rint( PATTR_SCORES->GetScore( m, n ));

                if( loc_score <= SCORE_MIN )
                    continue;

                prob = GetVectorProbabilityAt( m );
                PATTR_SCORES->IncQueryInfProbOf( prob, loc_score, n );
                PATTR_SCORES->IncSbjctInfProbOf( prob, loc_score, m );
                normterm += prob;
                querynorm[n] += prob;
            }
            if( normterm ) {
                //normalize probabilities for query positions
                for( sc = PATTR_SCORES->GetSbjctMinScoreAt( m ); sc <= PATTR_SCORES->GetSbjctMaxScoreAt( m ); sc++ )
                    PATTR_SCORES->DivSbjctInfProbOf( normterm, sc, m );
            }
        }

        for( n = 0; n < GetQuerySize(); n++ )
        {   //normalize probabilities
            if( querynorm[n] ) {
                for( sc = PATTR_SCORES->GetQueryMinScoreAt( n ); sc <= PATTR_SCORES->GetQueryMaxScoreAt( n ); sc++ )
                    PATTR_SCORES->DivQueryInfProbOf( querynorm[n], sc, n );
            }
        }
    } catch( myexception const& ex ) {
        if( querynorm ) free( querynorm );
        throw myruntime_error( ex.what(), ex.eclass());
    }

    if( querynorm ) free( querynorm );
}

// -------------------------------------------------------------------------
// ComputeScoreProbabilities: computes probabilities to observe scores at
//     each position (i,j). Probability is computed as
//
//      1     _
//    ----   \   Pfj
//     SUM   /_
//           i,j:
//         sij=sk
//
//  where Pfj are probabilities of frequency vector fj to occur; sij are
//  scores at (i,j), sk is a discrete score value; SUM is normalizing term.
//  Method iterates through the table of scores once
//
//  This method allocates required space for probabilities!
// -------------------------------------------------------------------------

void AdjustedScoreMatrix::ComputeScoreProbabilities( AttributableScores* PATTR_SCORES )
{
// ScoringMatrix::ComputeScoreProbabilities( PATTR_SCORES );
// return;

    if( PATTR_SCORES == NULL )
        throw myruntime_error( mystring( "AdjustedScoreMatrix: Unable to compute probabilities: Wrong argument." ));

    double  prob = 0.0;     //probability
    double  normterm = 0.0; //normalizing term
    double  avgscore = 0.0;
    TScore  loc_score  = 0;


    for( int m = 0; m < GetSubjectSize(); m++ ) {
        if( GetSbjctFreq()[m] == X )
            continue;

        for( int n = 0; n < GetQuerySize(); n++ )
        {
            if( GetQueryFreq()[n] == X )
                continue;

            //do not consider appropriately masked scores
            if( PATTR_SCORES->GetMaskedToIgnore( m, n ))
                continue;

            loc_score = ( TScore )rint( PATTR_SCORES->GetScore( m, n ));

            if( loc_score <= SCORE_MIN )
                continue;

            prob = GetVectorProbabilityAt( m );
            PATTR_SCORES->IncProbabilityOf( prob, loc_score );
            normterm += prob;
        }
    }

    //normalize probabilities
    for( int sc = PATTR_SCORES->GetMinScore(); sc <= PATTR_SCORES->GetMaxScore(); sc++ ) {
        PATTR_SCORES->DivideProbabilityOf( normterm, sc );
        avgscore += sc * PATTR_SCORES->GetProbabilityOf( sc );
    }

    PATTR_SCORES->SetExpectedScore( avgscore );
}

