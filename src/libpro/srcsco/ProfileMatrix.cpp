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
#include "mystring.h"
#include "myexcept.h"
#include "ext/pslcodes.h"
#include "libpro/tfopt/TargetFreqOptimizerH.h"
#include "ProfileMatrix.h"



double ProfileMatrix::position[NUMALPH];

//score for masked positions
const double    g_scoreX = 0.0;

// -------------------------------------------------------------------------
// constructor:pointer initialization
// -------------------------------------------------------------------------

ProfileMatrix::ProfileMatrix(
    const double ( *pssmscores )[NUMALPH],
    const char*     resids,
    int             length,
    Configuration   config[NoSchemes],
    TScaling        a_scaling,
    TMask           c_masking )
:
    AbstractScoreMatrix( PositionSpecific, config, ComputeStatistics, a_scaling, c_masking ),
    residues( resids ),

    scores_( length * NUMALPH ),
    rprobs_( length ),
    cprobs_( NUMALPH )
{
    Init( length, NUMALPH );
    FillMatrix( pssmscores );
    SetSupportOptimFreq( false );
}

// -------------------------------------------------------------------------
// default constructor is invalid for initialization
// -------------------------------------------------------------------------

ProfileMatrix::ProfileMatrix()
:
    AbstractScoreMatrix(),
    scores_( 1 ),
    rprobs_( 1 ),
    cprobs_( 1 )
{
}

// -------------------------------------------------------------------------
// destructor: deallocation of resources
// -------------------------------------------------------------------------

ProfileMatrix::~ProfileMatrix()
{
}

// -------------------------------------------------------------------------
// FillMatrix: fills matrix with the given scores
// -------------------------------------------------------------------------

void ProfileMatrix::FillMatrix(
        const double ( *pssmscores )[NUMALPH] )
{
    if( !IsValid() || GetQuerySize() <= 0 || GetSubjectSize() != NUMALPH )
        throw myruntime_error( mystring( "ProfileMatrix: Failed to fill matrix with scores." ));

    int     negatives = 0;
    bool    all_negat = true;

    //fill matrix with values
    for( int m = 0; m < GetSubjectSize(); m++ ) {
        negatives = 0;

        for( int n = 0; n < GetQuerySize(); n++ )
        {
            double  score = pssmscores[n][m];

            SetScore( score, m, n );

            if( all_negat && score < 0 ) negatives++;   //increase counter for each negative score found
        }

        if( all_negat && negatives != GetQuerySize())
            all_negat = false;
    }

    SetAllNegatives( all_negat );
}

// -------------------------------------------------------------------------
// GetVectorAt: Returns vecor of values from the matrix given position
// -------------------------------------------------------------------------

const double ( *ProfileMatrix::GetVectorAt( int n ) const )[NUMALPH]
{
#ifdef __DEBUG__
    if( GetQuerySize() <= n || GetSubjectSize() != NUMALPH )
        throw myruntime_error( mystring( "ProfileMatrix: Memory access error." ));
#endif

    for( int r = 0; r < NUMALPH; r++ )
        position[r] = GetFinalScore( GetScore( r, n ));

    return &position;//( CONST_SCORE_ARRAY )( GetScores() + n );
}

// -------------------------------------------------------------------------
// ComputePositionalScoreProbs: should not be called
// -------------------------------------------------------------------------

void ProfileMatrix::ComputePositionalScoreProbs( AttributableScores* )
{
    throw myruntime_error( mystring( "ProfileMatrix: ComputePositionalScoreProbs should not be called from this class." ));
}

// -------------------------------------------------------------------------
// ComputeScoreProbabilities: computes probabilities to observe scores at
//     each position (i,j). Probability is computed as
//
//      1     _
//    ----   \   Pa
//   length  /_ 
//           i,j:
//         sij=sk
//
//  where Pa are background probabilities; sij are score at (i,j), sk is a
//  discrete score value; length is for query
//
//  This method allocates required space for probabilities!
// -------------------------------------------------------------------------

void ProfileMatrix::ComputeScoreProbabilities( AttributableScores* PATTR_SCORES )
{
    if( !GetResidues() || GetSubjectSize() != NUMALPH )
        throw myruntime_error( mystring( "ProfileMatrix: Unable to compute score probabilities." ));

    if( PATTR_SCORES == NULL )
        throw myruntime_error( mystring( "ProfileMatrix: Unable to compute probabilities: Wrong argument." ));

    const double    accuracy = 1.e-4;
    char    strbuf[KBYTE];
    double  avgscore = 0.0;
    double  consv = 0.0;
    double  score = 0.0;
    double  prob = 0.0;
    size_t  length = 0;
    char    res;
    int     n, m, nn, mm;

    static int  symB = HashAlphSymbol('B');
    static int  symZ = HashAlphSymbol('Z');
    static int  resN = HashAlphSymbol('N');
    static int  resD = HashAlphSymbol('D');
    static int  resQ = HashAlphSymbol('Q');
    static int  resE = HashAlphSymbol('E');

    scores_.Clear();
    rprobs_.Clear();
    cprobs_.Clear();


    for( n = 0, nn = 0; n < GetQuerySize(); n++ ) {
        res = GetResidueAt( n );
//         if( res == X )
//             continue;
        prob = LOSCORES.PROBABility( res );

        if( res == X )
            prob = 0.05;
        else if( res == symB )
            prob = ( LOSCORES.PROBABility( resN ) + LOSCORES.PROBABility( resD )) * 0.5;
        else if( res == symZ )
            prob = ( LOSCORES.PROBABility( resQ ) + LOSCORES.PROBABility( resE )) * 0.5;

        if( prob <= 0.0 )
            continue;
        rprobs_.AddValueAt( nn++, prob );
    }

    for( m = 0, mm = 0; m < GetSubjectSize(); m++ ) {
        prob = LOSCORES.PROBABility( m );
        if( prob <= 0.0 )
            continue;
        cprobs_.AddValueAt( mm++, prob );
    }



    for( n = 0; n < GetQuerySize(); n++ ) {
//         if( GetResidueAt( n ) == X )
//             continue;

//         if( LOSCORES.PROBABility( n ) <= 0.0 )
//             continue;

        length++;
        for( m = 0; m < GetSubjectSize(); m++ )
        {
            score = PATTR_SCORES->GetScore( m, n );

            if( LOSCORES.PROBABility( m ) <= 0.0 )
                continue;

            if( score <= SCORE_MIN )
                    scores_.Push( g_scoreX );
            else    scores_.Push( score );

            if( score <= SCORE_MIN )
                continue;

            PATTR_SCORES->IncProbabilityOf( LOSCORES.PROBABility( m ), score );
        }
    }

    //normalize probabilities
    for( int sc = PATTR_SCORES->GetMinScore(); sc <= PATTR_SCORES->GetMaxScore(); sc++ ) {
        PATTR_SCORES->DivideProbabilityOf( length, sc );
        avgscore += sc * PATTR_SCORES->GetProbabilityOf( sc );
    }


    consv = rprobs_.Sum();
    if( 0.0 < consv ) {
        rprobs_.MultiplyBy( 1.0 / consv );
        consv = rprobs_.Sum();
        if( consv < 1.0 - accuracy || consv > 1.0 + accuracy ) {
            sprintf( strbuf, "Background probabilities not conserved: %f", consv );
            throw myruntime_error( strbuf );
        }
    }
    consv = cprobs_.Sum();
    if( consv < 1.0 - accuracy || consv > 1.0 + accuracy ) {
        sprintf( strbuf, "Background probabilities not conserved: %f", consv );
        throw myruntime_error( strbuf );
    }


   PATTR_SCORES->SetExpectedScore( avgscore );
}


// =========================================================================
// OPTIMIZATION OF TARGET FREQUENCIES
//
// OptimizeTargetFrequencies: Optimize target frequencies
//
void ProfileMatrix::OptimizeTargetFrequencies()
{
    double  lambda;

    if( GetAllNegatives()) {
        warning( "Unable to optimize target frequencies." );
        return;
    }

    ComputeStatisticalParameters();

    if( scores_.GetSize() < 1 || rprobs_.GetSize() < 1 || cprobs_.GetSize() < 1 )
        return;

    if( 0.0 <= GetExpectedScore() || GetLambda() <= 0.0 ) {
        warning( "Unable to optimize target frequencies." );
        return;
    }


    lambda = GetLambda();
    if( GetAutoScaling())
        lambda = GetScaledScores()->GetLambda();

    scores_.MultiplyBy( lambda );


    TargetFreqOptimizerH    optimizer( scores_, rprobs_, cprobs_ );
    int                     status;
    int                     n, m, ind = 0;

//{{TESTING..
// TestPrintScoringMatrix( stderr );
// fprintf( stderr, "\n\n" );
// for( n = 0; n < rprobs_.GetSize(); n++ ) {
//     for( m = 0; m < cprobs_.GetSize(); m++ )
//         fprintf( stderr, "%12f ", scores_.GetValueAt( ind++ ));
//     fprintf( stderr, "\n");
// }
//}}

    //no need to set lambda if scores were previously multiplied by
//     optimizer.SetLambda( lambda );
    optimizer.SetConstrainedH( GetRefH());
    optimizer.SetNMResidTol( TARGET_FREQ_OPT_TOLERANCE );
    optimizer.SetMaxNoNMIterations( TARGET_FREQ_OPT_MAXITERATIONS );

    status = optimizer.Optimize();

    if( status != PSL_SUCCESS /*&& status != PSL_MAXITERATS */) {
#ifdef TFOPTESTPRINT
        warning( TranslatePSLError( status ));
#else
        warning( "Target frequencies not optimized." );
#endif
        return;
    }

    if( status != PSL_SUCCESS )
        optimizer.Normalize();

    if( optimizer.Negatives()) {
        warning( "Invalid target frequencies." );
        return;
    }

    IncorporateTargetFrequencies( optimizer.GetOptTargetFreqs());
}

// -------------------------------------------------------------------------
// IncorporateTargetFrequencies: Calculate new scores based on optimized
//     target frequencies
//
void ProfileMatrix::IncorporateTargetFrequencies( const Pslvector& optfreqs )
{
    if( rprobs_.GetSize() < 1 || cprobs_.GetSize() < 1 ||
        optfreqs.GetSize() != rprobs_.GetSize() * cprobs_.GetSize())
        return;

    const double    accuracy = 1.0e-6;
    int     n, m, nn, mm, ind;
    double  tfval = 0.0;    //target frequency value
    double  proquery = 0.0; //background probability of query position
    double  prosbjct = 0.0; //background probability of subject position
    double  loc_score = 0.0;
    double  score = 0.0;
    double  consv;
    char    res;

    consv = optfreqs.Sum();
    if( consv < 1.0 - accuracy || consv > 1.0 + accuracy ) {
        warning( "Target frequencies not conserved." );
        return;
    }

    ind = 0;
    for( n = 0, nn = 0; n < GetQuerySize(); n++ )
    {
        res = GetResidueAt( n );
//         if( res == X )
//             continue;
        if( ASTERISK <= res )
            continue;

        proquery = rprobs_.GetValueAt( nn++ );

        for( m = 0, mm = 0; m < GetSubjectSize(); m++ ) {

            if( LOSCORES.PROBABility( m ) <= 0.0 )
                continue;

            prosbjct = cprobs_.GetValueAt( mm++ );

            tfval = optfreqs.GetValueAt( ind++ );

            if( proquery <= 0.0 || prosbjct <= 0.0 )
                    throw myruntime_error( "Invalid probabilities after optimizing of target frequencies." );

            score = log( tfval / ( proquery * prosbjct ));

            loc_score = GetScore( m, n );

            if( loc_score <= SCORE_MIN )
                continue;

            SetScore( score, m, n );
        }
    }
}



// =========================================================================
// PRINT ROUTINES
//
// PrintParameterTable: print table of computed statistical parameter
//     values to a stream
//

void ProfileMatrix::PrintParameterTable( TPrintFunction print_func, void* vpn ) const
{
    if( vpn == NULL )
        return;

    if( 0.0 <= GetExpectedScore()) {
        print_func( vpn, "Expected score per position is non-negative, %.4f!\n\n", GetExpectedScore());
        return;
    }

    char    _sk[BUF_MAX], _slambda[BUF_MAX];
    if( 0.0 < GetK())   sprintf( _sk, "%6.4f", GetK());
        else            sprintf( _sk, "n/a" );
    if( 0.0 < GetLambda())  sprintf( _slambda, "%6.4f", GetLambda());
        else                sprintf( _slambda, "n/a" );

    print_func( vpn, "%-25s  %-6s   %-6s\n", " ", "K", "Lambda" );
    print_func( vpn, "%-25s  %6s   %6s\n", "Computed  ungapped,", _sk, _slambda );
    print_func( vpn, "%-25s  %6.4f   %6.4f\n", "Reference ungapped,", GetRefK(), GetRefLambda());

    print_func( vpn, "Entropy, %6.4f; Expected, %6.4f\n\n", GetEntropy(), GetExpectedScore());
}

// PrintFinal: nothing to print in the final
//

void ProfileMatrix::PrintFinal( TPrintFunction, void* ) const
{
}

// -------------------------------------------------------------------------
// PrintScoringMatrix: output the computed scoring system
// -------------------------------------------------------------------------

void ProfileMatrix::PrintScoringMatrix( FILE* fp )
{
    int l = 0;
//     const size_t    effective_nr = 20;

    if( fp == NULL )
        return;

    fprintf( fp,"%12c Position-specific scoring matrix\n", 32 );

    fprintf( fp, "%9c", 32 );

    for( int m = 0; m < GetSubjectSize(); m++ )
        fprintf( fp, "%4c", DehashCode( m ));


    for( int n = 0; n < GetQuerySize(); n++ ) {
        fprintf( fp, "\n%5d %c   ", ++l, DehashCode( GetResidueAt( n )));

        for( int m = 0; m < GetSubjectSize(); m++ )
            fprintf( fp, "%3d ", ( int )GetScore( m, n ));
    }
    fprintf( fp, "\n" );
}

// ==== Functions ==========================================================
// ComputedSubMatrixWithParams: output computed substitution matrix with
//     calculated statistical parameters for it
//
void ComputedSubMatrixWithParams()
{
    Configuration   config[NoSchemes];  //not setting filename, so that if read or write attempt is, an error occurs
    SetUngappedParams( config[ProcomUngapped] );//Fill values with parameter values of ungapped configuration
    SetUngappedParams( config[ProcomGapped] );  //not using gapped configuration, make it identical to ungapped configuration

      const char      resids[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
    ProfileMatrix   submat( COMPUTED_GONNET.data, resids, NUMAA, config, AbstractScoreMatrix::NoScaling );
    submat.ComputeStatisticalParameters();
//     submat.ScaleScoringMatrix(); //if scaling needed
    submat.PrintScoringMatrix( stderr );
    dynamic_cast<AbstractScoreMatrix&>( submat ).PrintParameterTable( stderr );
}

// -------------------------------------------------------------------------
// TestPrintScoringMatrix: output scoring matrix
//
void ProfileMatrix::TestPrintScoringMatrix( FILE* fp ) const
{
    int l = 0;

    if( fp == NULL )
        return;

    fprintf( fp,"%12c Position-specific scoring matrix\n", 32 );

    fprintf( fp, "%9c", 32 );

    for( int m = 0; m < GetSubjectSize(); m++ )
        fprintf( fp, "%8c", DehashCode( m ));


    for( int n = 0; n < GetQuerySize(); n++ ) {
        fprintf( fp, "\n%5d %c   ", ++l, DehashCode( GetResidueAt( n )));

        for( int m = 0; m < GetSubjectSize(); m++ )
            fprintf( fp, "%7f ", GetScore( m, n ));
    }
    fprintf( fp, "\n" );
}
