/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "mystring.h"
#include "myexcept.h"

#include "rc.h"
#include "data.h"
#include "ext/pslcodes.h"
#include "libpro/tfopt/TargetFreqOptimizerH.h"
#include "ScoringMatrix.h"


FrequencyMatrix smp_dummyfreq;
LogOddsMatrix   smp_dummylogo;

const double    g_scoreX = 0.0; //-1.0;

// -------------------------------------------------------------------------
// constructor: frequency matrices, log-odds matrices for the first and
//   second profiles are given with the parameters
// -------------------------------------------------------------------------

ScoringMatrix::ScoringMatrix(
    const FrequencyMatrix&  freq_fst, const LogOddsMatrix& logo_fst,
    const FrequencyMatrix&  freq_sec, const LogOddsMatrix& logo_sec,
    double          infrm_threshold,
    int             thick_number,
    double          thick_percnt,
    double          mask_percnt,
    Configuration   config[NoSchemes],
    TBehaviour      s_behaviour,
    TScaling        a_scaling,
    TMask           c_masking )
:
    AbstractScoreMatrix( ProfileSpecific, config, s_behaviour, a_scaling, c_masking ),

    freq_fst_( freq_fst ),
    logo_fst_( logo_fst ),

    freq_sec_( freq_sec ),
    logo_sec_( logo_sec ),

    thickness_number( thick_number ),
    thickness_percnt( thick_percnt ),
    maskscale_percnt( mask_percnt ),

    scores_( freq_fst_.GetColumns() * freq_sec_.GetColumns()),
    rprobs_( freq_fst_.GetColumns()),
    cprobs_( freq_sec_.GetColumns())
{
    PrivateInit();
    SetInformationThreshold( infrm_threshold );
    SetSupportOptimFreq( false );
}

// -------------------------------------------------------------------------
// constructor: overloaded
//

ScoringMatrix::ScoringMatrix(
    const FrequencyMatrix& freq_fst, const LogOddsMatrix& logo_fst,
    const FrequencyMatrix& freq_sec, const LogOddsMatrix& logo_sec,
    double          infrm_threshold,
    int             thick_number,
    double          thick_percnt,
    double          mask_percnt,
    Configuration   config[NoSchemes],
    TBehaviour      s_behaviour,
    TScaling        a_scaling,
    TMask           c_masking,
    TType           type )
:
    AbstractScoreMatrix( type, config, s_behaviour, a_scaling, c_masking ),

    freq_fst_( freq_fst ),
    logo_fst_( logo_fst ),

    freq_sec_( freq_sec ),
    logo_sec_( logo_sec ),

    thickness_number( thick_number ),
    thickness_percnt( thick_percnt ),
    maskscale_percnt( mask_percnt ),

    scores_( freq_fst_.GetColumns() * freq_sec_.GetColumns()),
    rprobs_( freq_fst_.GetColumns()),
    cprobs_( freq_sec_.GetColumns())
{
    PrivateInit();
    SetInformationThreshold( infrm_threshold );
}

// -------------------------------------------------------------------------
// default constructor is invalid for initialization
// -------------------------------------------------------------------------

ScoringMatrix::ScoringMatrix()
:
    AbstractScoreMatrix(),

    freq_fst_( smp_dummyfreq ),
    logo_fst_( smp_dummylogo ),

    freq_sec_( smp_dummyfreq ),
    logo_sec_( smp_dummylogo ),

    thickness_number( 0 ),
    thickness_percnt( 0.0 ),
    maskscale_percnt( 0.0 ),

    scores_( 1 ),
    rprobs_( 1 ),
    cprobs_( 1 )
{
}

// -------------------------------------------------------------------------
// destructor: deallocation of resources
// -------------------------------------------------------------------------

ScoringMatrix::~ScoringMatrix()
{
}

// -------------------------------------------------------------------------
// PrivateInit: initialization method private for this class
// -------------------------------------------------------------------------

void ScoringMatrix::PrivateInit()
{
    if( !logo_fst_.IsCompatible( freq_fst_ ) || !logo_sec_.IsCompatible( freq_sec_ ))
            throw myruntime_error( mystring( "ScoringMatrix: Profile matrices are incompatible." ));

    if(  freq_fst_.GetColumns() < 1 || freq_fst_.GetColumns() > MAXCOLUMNS ||
         freq_sec_.GetColumns() < 1 || freq_sec_.GetColumns() > MAXCOLUMNS )
            throw myruntime_error( mystring( "ScoringMatrix: Wrong profile matrices." ));

    if( GetMaskscalePercents() < 0.0 )
            throw myruntime_error( mystring( "ScoringMatrix: Negative scaling factor of masked positions." ));

    Init( freq_fst_.GetColumns(), freq_sec_.GetColumns());

#ifdef USEPROFBACKPROBS
    LOSCORES.StoreProbabilities_1( logo_fst_.GetBackProbs());
    LOSCORES.StoreProbabilities_2( logo_sec_.GetBackProbs());
#endif
}

// -------------------------------------------------------------------------
// ComputeScore: compute score given position of query and subject
//
double ScoringMatrix::ComputeScore( int m, int n, bool final ) const
{
    double  score = 0.0;
    double  lambda = LOSCORES.StatisParam( Ungapped, Lambda ); //reference lambda

    size_t  fst_thckn = 0;
    size_t  sec_thckn = 0;

    double  log_fst_thckn = 0.0;
    double  log_sec_thckn = 0.0;
    double  sum_thckn = 0.0;

    double  fst_weight = 1.0;
    double  sec_weight = 1.0;

    if( final ) {
        fst_thckn = GetQueryLogo().GetMIDExpNoObservationsAt( n, PS_M );
        sec_thckn = GetSbjctLogo().GetMIDExpNoObservationsAt( m, PS_M );


//         if( GetQueryLogo().GetEffNoSequences() &&
//             GetSbjctLogo().GetEffNoSequences())
//         {
//             log_fst_thckn = ( double )fst_thckn / ( double )GetQueryLogo().GetEffNoSequences();
//             log_sec_thckn = ( double )sec_thckn / ( double )GetSbjctLogo().GetEffNoSequences();
//         }
//         else {
//             if( fst_thckn ) log_fst_thckn = 1.0 - sqrt( fst_thckn ) / ( double ) fst_thckn; //log( fst_thckn );
//             if( sec_thckn ) log_sec_thckn = 1.0 - sqrt( sec_thckn ) / ( double ) sec_thckn; //log( sec_thckn );

            if( fst_thckn ) log_fst_thckn = 1.0 - ( 1.0 + log( fst_thckn )) / ( double ) fst_thckn; //log( fst_thckn );
            if( sec_thckn ) log_sec_thckn = 1.0 - ( 1.0 + log( sec_thckn )) / ( double ) sec_thckn; //log( sec_thckn );
//         }

        if( log_fst_thckn <= 0.0 && log_sec_thckn <= 0.0 ) {
            log_fst_thckn = 1.0;
            log_sec_thckn = 1.0;
        }

        sum_thckn = log_fst_thckn + log_sec_thckn;

        if( sec_thckn == 0.0 )
            sum_thckn = 1.0;

        fst_weight = ( log_fst_thckn + log_fst_thckn ) / sum_thckn;
        sec_weight = ( log_sec_thckn + log_sec_thckn ) / sum_thckn;
    }

    for( int a = 0; a < NUMAA; a++ )
        // lambda is implicitly incorporated into the log-odds values of the profiles
        score += fst_weight * GetQueryFreq()( n, a ) * GetSbjctLogo()( m, a ) +
                 sec_weight * GetSbjctFreq()( m, a ) * GetQueryLogo()( n, a );

    return ( score - 0.0 ) * GetMultiplier();
}

// -------------------------------------------------------------------------
// ComputeProfileScoringMatrix: compute profile-profile scoring matrix
//
void ScoringMatrix::ComputeProfileScoringMatrix( bool final )
{
    if( GetQuerySize() <= 0 || GetSubjectSize() <= 0 )
        throw myruntime_error( mystring( "ScoringMatrix: Unable to compute matrix of scores." ));

    if( !IsValid())
        throw myruntime_error( mystring( "ScoringMatrix: Failed to compute matrix of scores." ));

    int     m, n;
    int     no_elems = 0;
    int     negatives = 0;
    TMask   currentmask = Unmasked;
    double  score = 0.0;


    //fill matrix with values
    for( m = 0; m < GetSubjectSize(); m++ ) {
        for( n = 0; n < GetQuerySize(); n++ )
        {
            //give constant penalties at the positions of X
//             if( GetQueryFreq().GetResidueAt( n ) == X ||
//                 GetSbjctFreq().GetResidueAt( m ) == X )
//             {
//                 SetScore( g_scoreX, m, n );
//                 continue;
//             }

            score = 0.0;
            currentmask = Unmasked;

            if( logo_fst_.GetInformationAt( n ) < GetInformationThreshold() ||
                logo_sec_.GetInformationAt( m ) < GetInformationThreshold() ||
                ! ThicknessConstraintsMet( m, n ))
            {
                //mask the score possibly excluding it from the ongoing computational statistics
                SetMasked( currentmask = GetMaskingApproach(), m, n );
            }

            if( 1 || currentmask != MaskToIgnore )
            {
                score = ComputeScore( m, n, final );

                if( currentmask != MaskToIgnore ) {
                    no_elems++;
                    if( score < 0.0 ) negatives++;    //increase counter for each negative score found
                }
                if( currentmask == MaskToConsider || currentmask == MaskToIgnore )
                    score *= GetMaskscalePercents();
            }
            else
                //give penalty to score at the positions;
                score = g_scoreX;

            SetScore( score, m, n );
        }
    }


    if( negatives == no_elems ) {
        SetAllNegatives( true );
        return;
    }

    PreliminaryVerification();
}

// #include "alnsrc/FastAlignment.h"

// -------------------------------------------------------------------------
// PreliminaryVerification: apply fast alignment algorithm before scaling
//
void ScoringMatrix::PreliminaryVerification()
{
    double  expect;

    if( 0.0 < GetLambda())
        //the method has been called earlier
        return;

//     fprintf( stderr, "\n" );
//     for( int n = 0; n < GetQueryFreq().GetColumns(); n++ )
//         fprintf( stderr, "%c", DehashCode( GetQueryFreq().GetResidueAt( n )));
//     fprintf( stderr, "\n\n" );
// 
//     for( int m = 0; m < GetSbjctFreq().GetColumns(); m++ )
//         fprintf( stderr, "%c", DehashCode( GetSbjctFreq().GetResidueAt( m )));
//     fprintf( stderr, "\n" );

//     FastAlignment   fasta( this, GetQueryFreq().GetResidues(), GetSbjctFreq().GetResidues());
//     fasta.Run();

//     ComputeExpectation( fasta.GetScore(), &expect );
//     SetPrelimExpect( expect );
//     SetPrelimScore( fasta.GetScore());
// fprintf( stderr, "S %9.3f \tE %12.3g \t%s\n", fasta.GetScore(), expect, GetSbjctLogo().GetName());

//     fasta.PrintDPMatrix( stderr );
//     fasta.Print( stderr );
}

// -------------------------------------------------------------------------
// ComputeExpectation: computes e-value given previously computed
//     statistical parameters
// -------------------------------------------------------------------------

double ScoringMatrix::ComputeExpectation(
        double score,
        double* ref_expect,
        double* pure_expect,
        double* pair_expect,
        double* bitscore ) const
{
    double  expect = -1.0;
    return AbstractScoreMatrix::ComputeExpectation( score, ref_expect, pure_expect, pair_expect, bitscore );
}

// -------------------------------------------------------------------------
// ComputePositionalScoreProbs: computes probabilities of scores at each
//     query and subject position
// -------------------------------------------------------------------------

void ScoringMatrix::ComputePositionalScoreProbs( AttributableScores* PATTR_SCORES )
{
    if( PATTR_SCORES == NULL )
        throw myruntime_error( mystring( "ScoringMatrix: Unable to compute probabilities: Wrong argument." ));

    double  sum = 0.0;          //initial sum
    double  consval = 0.0;      //probability conservation value
    double  exp2sum = 0.0;      //e to sum
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
    try{
        for( m = 0; m < GetSubjectSize(); m++ ) {
//             if( freq_sec_[m] == X )
//                 continue;

            normterm = 0.0;
            for( n = 0; n < GetQuerySize(); n++ )
            {
//                 if( freq_fst_[n] == X )
//                     continue;
                sum = 0.0;

                if( PATTR_SCORES->GetMaskedToIgnore( m, n ))
                    continue;

                loc_score = ( TScore )rint( PATTR_SCORES->GetScore( m, n ));

                if( loc_score <= SCORE_MIN )
                    continue;

                if( 0 ) {
                    for( int r = 0; r < NUMALPH; r++ ) {
                        if( LOSCORES.PROBABility( r ) <= 0.0 )
                            continue;
//                         sum += ( freq_fst_( n, r ) + freq_sec_( m, r )) * LOSCORES.LogPROBABility( r );
                        sum +=  freq_fst_( n, r ) * LOSCORES.LogPROBABILITY_2( r ) +
                                freq_sec_( m, r ) * LOSCORES.LogPROBABILITY_1( r );
                    }
                    exp2sum = exp( sum );
                } else
                    exp2sum = 1.0;
                PATTR_SCORES->IncQueryInfProbOf( exp2sum, loc_score, n );
                PATTR_SCORES->IncSbjctInfProbOf( exp2sum, loc_score, m );
                normterm += exp2sum;
                querynorm[n] += exp2sum;
            }
            if( normterm ) {
                consval = 0.0;
                PATTR_SCORES->SetSbjctInfProbNorm( normterm, m );
                //normalize probabilities for query positions
                for( sc = PATTR_SCORES->GetSbjctMinScoreAt( m ); sc <= PATTR_SCORES->GetSbjctMaxScoreAt( m ); sc++ )
                    consval += PATTR_SCORES->DivSbjctInfProbOf( normterm, sc, m );
                if( consval < 0.999 || consval > 1.001 )
                    warning( "ScoringMatrix: ComputePositionalScoreProbs: Probabilities not conserved." );
            }
        }

        for( n = 0; n < GetQuerySize(); n++ )
        {   //normalize probabilities
            if( querynorm[n] ) {
                consval = 0.0;
                PATTR_SCORES->SetQueryInfProbNorm( querynorm[n], n );
                for( sc = PATTR_SCORES->GetQueryMinScoreAt( n ); sc <= PATTR_SCORES->GetQueryMaxScoreAt( n ); sc++ )
                    consval += PATTR_SCORES->DivQueryInfProbOf( querynorm[n], sc, n );
                if( consval < 0.999 || consval > 1.001 )
                    warning( "ScoringMatrix: ComputePositionalScoreProbs: Probabilities not conserved." );
            }
        }
    } catch( myexception const& ex ) {
        if( querynorm ) free( querynorm );
        throw myruntime_error( ex.what(), ex.eclass());
    }

    if( querynorm ) free( querynorm );
}

#if 1
// -------------------------------------------------------------------------
// ComputeScoreProbabilities: computes probabilities to observe scores at
//     each position (i,j). Probability is computed as
//
//     1    _   __    (Fia + Fja)
//   ----  \   |  | Pa
//   Norm  /_   a
//         i,j:
//       sij=sk
//
//  where Pa are background probabilities, Fia and Fja are observed
//  frequencies from profiles i and j, respectively; sij are score
//  at (i,j), sk is a discrete score value; Norm is normalizing term
//  which is the sum over all i,j
//
//  This method allocates required space for probabilities!
//  (Older version)
// -------------------------------------------------------------------------

void ScoringMatrix::ComputeScoreProbabilities( AttributableScores* PATTR_SCORES )
{
    if( PATTR_SCORES == NULL )
        throw myruntime_error( mystring( "ScoringMatrix: Unable to compute probabilities: Wrong argument." ));

    const double    accuracy = 1.0e-4;

    int     n, m, nn, mm, r;
    char    strbuf[KBYTE];
    double  sum = 0.0;      //initial sum
    double  proquery = 0.0; //non-normalized probability of query position
    double  prosbjct = 0.0; //non-normalized probability of subject position
    double  exp2sum = 0.0;  //e to sum
    double  normprob = 0.0; //normalizing term
    double  normterm = 0.0; //normalizing term
    double  avgscore = 0.0;
    double  consv;
    TScore  loc_score  = 0;

    scores_.Clear();
    rprobs_.Clear();
    cprobs_.Clear();


    for( n = 0, nn = 0; n < GetQuerySize(); n++ ) {
//         if( freq_fst_[n] == X )
//             continue;
        if( 0 ) {
            sum = 0.0;
            for( r = 0; r < NUMALPH; r++ ) {
                if( LOSCORES.PROBABility( r ) <= 0.0 )
                    continue;
                sum += freq_fst_( n, r ) * LOSCORES.LogPROBABILITY_1( r );
            }
            rprobs_.AddValueAt( nn++, exp( sum ));
        } else
            rprobs_.SetValueAt( nn++, 1.0 );
    }

    for( m = 0, mm = 0; m < GetSubjectSize(); m++ ) {
//         if( freq_sec_[m] == X )
//             continue;
        if( 0 ) {
            sum = 0.0;
            for( r = 0; r < NUMALPH; r++ ) {
                if( LOSCORES.PROBABility( r ) <= 0.0 )
                    continue;
                sum += freq_sec_( m, r ) * LOSCORES.LogPROBABILITY_2( r );
            }
            cprobs_.AddValueAt( mm++, exp( sum ));
        } else
            cprobs_.SetValueAt( mm++, 1.0 );
    }

    avgscore = normterm = 0.0;
    for( n = 0, nn = 0; n < GetQuerySize(); n++ )
    {
//         if( freq_fst_[n] == X )
//             continue;
        proquery = rprobs_.GetValueAt( nn++ );

        for( m = 0, mm = 0; m < GetSubjectSize(); m++ ) {
//             if( freq_sec_[m] == X )
//                 continue;
            prosbjct = cprobs_.GetValueAt( mm++ );

            loc_score = ( TScore )rint( PATTR_SCORES->GetScore( m, n ));

            exp2sum = proquery * prosbjct;

            if( loc_score <= SCORE_MIN )
                    scores_.Push( g_scoreX );
            else    scores_.Push( loc_score );

            if( loc_score <= SCORE_MIN )
                continue;

            if( PATTR_SCORES->GetMaskedToIgnore( m, n ))
                continue;

            PATTR_SCORES->IncProbabilityOf( exp2sum, loc_score );
            normterm += exp2sum;
        }
    }
    //normalize probabilities
    if( normterm )
        for( int sc = PATTR_SCORES->GetMinScore(); sc <= PATTR_SCORES->GetMaxScore(); sc++ ) {
            PATTR_SCORES->DivideProbabilityOf( normterm, sc );
            avgscore += sc * PATTR_SCORES->GetProbabilityOf( sc );
        }

    normprob = rprobs_.Sum();
    if( 0.0 < normprob ) {
        rprobs_.MultiplyBy( 1.0 / normprob );
        consv = rprobs_.Sum();
        if( consv < 1.0 - accuracy || consv > 1.0 + accuracy ) {
            sprintf( strbuf, "Background probabilities not conserved: %f", consv );
            throw myruntime_error( strbuf );
        }
    }
    normprob = cprobs_.Sum();
    if( 0.0 < normprob ) {
        cprobs_.MultiplyBy( 1.0 / normprob );
        consv = cprobs_.Sum();
        if( consv < 1.0 - accuracy || consv > 1.0 + accuracy ) {
            sprintf( strbuf, "Background probabilities not conserved: %f", consv );
            throw myruntime_error( strbuf );
        }
    }

    PATTR_SCORES->SetExpectedScore( avgscore );
}
#else//versions of ComputeScoreProbabilities
// -------------------------------------------------------------------------
// ComputeScoreProbabilities: computes probabilities to observe scores at
//     each position (i,j). Probability is computed according to multinomial
//     distribution of frequencies
//
//                   _
//                ( \  Fja ) !
//     1    _       /_ a       __    (Fja)
//   ----  \      ----------- |  | Pa
//   Norm  /_       __         a
//         i,j:    |  | Fja!
//       sij=sk     a
//
//  where Pa are background probabilities, Fja are observed frequencies from
//  profile j (subject profile); sij are score at (i,j), sk is a discrete
//  score value; Norm is equal to the length of query
//
//  This method allocates required space for probabilities!
// -------------------------------------------------------------------------

void ScoringMatrix::ComputeScoreProbabilities( AttributableScores* PATTR_SCORES )
{
    if( PATTR_SCORES == NULL )
        throw myruntime_error( mystring( "ScoringMatrix: Unable to compute probabilities: Wrong argument." ));

    double          sum = 0.0;      //initial sum
    double          exp2sum = 0.0;  //exp to sum
    double          avgscore = 0.0;
    unsigned int    freq = 0;   //particular frequency value in interval [0-100]
    TScore          loc_score  = 0;


    for( int n = 0; n < GetQuerySize(); n++ ) {
        if( freq_fst_[n] == X )
            continue;

        for( int m = 0; m < GetSubjectSize(); m++ )
        {
            if( freq_sec_[m] == X )
                continue;

            sum = 0.0;

            if( PATTR_SCORES->GetMaskedToIgnore( m, n ))
                continue;

            loc_score = ( TScore )rint( PATTR_SCORES->GetScore( m, n ));

            if( loc_score <= SCORE_MIN )
                continue;

            for( int r = 0; r < NUMALPH; r++ ) {
                if( LOSCORES.PROBABility( r ) <= 0.0 )
                    continue;//consider only residues which have background probabilities

                freq = ( unsigned )rint( 100.0 * freq_sec_( m, r ));
                sum += freq * LOSCORES.LogPROBABility( r ) - LOG_FREQUENCIES.SumOf( freq );
            }
            sum += LOG_FREQUENCIES.Total();//log gamma(sum freq)

            exp2sum = exp( sum );
            PATTR_SCORES->IncProbabilityOf( exp2sum, loc_score );
        }
    }

    //normalize probabilities
    for( int sc = PATTR_SCORES->GetMinScore(); sc <= PATTR_SCORES->GetMaxScore(); sc++ ) {
        DivideProbabilityOf( GetQuerySize(), sc );    //normalizing term
        avgscore += sc * PATTR_SCORES->GetProbabilityOf( sc );
    }

    PATTR_SCORES->SetExpectedScore( avgscore );
}
#endif//versions of ComputeScoreProbabilities

// -------------------------------------------------------------------------
// OptimizeTargetFrequencies: Optimize target frequencies
//
void ScoringMatrix::OptimizeTargetFrequencies()
{
    double  lambda;

    if( GetAllNegatives()) {
        //warning( "Unable to optimize target frequencies." );
        return;
    }

    ComputeStatisticalParameters();

    if( scores_.GetSize() < 1 || rprobs_.GetSize() < 1 || cprobs_.GetSize() < 1 )
        return;

    if( 0.0 <= GetExpectedScore() || GetLambda() <= 0.0 ) {
        //warning( "Unable to optimize target frequencies." );
        return;
    }


    lambda = GetLambda();
    if( GetAutoScaling())
        lambda = GetScaledScores()->GetLambda();

    scores_.MultiplyBy( lambda );


    TargetFreqOptimizerH    optimizer( scores_, rprobs_, cprobs_ );
    int                     status;

    //no need to set lambda if scores were previously multiplied by
//     optimizer.SetLambda( GetLambda());
    optimizer.SetConstrainedH( GetRefH());
    optimizer.SetNMResidTol( TARGET_FREQ_OPT_TOLERANCE );
    optimizer.SetMaxNoNMIterations( TARGET_FREQ_OPT_MAXITERATIONS );

    status = optimizer.Optimize();

    if( status != PSL_SUCCESS /*&& status != PSL_MAXITERATS*/ ) {
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
        warning( "Invalid target frequencies after optimization." );
        return;
    }

    IncorporateTargetFrequencies( optimizer.GetOptTargetFreqs());
}

// -------------------------------------------------------------------------
// IncorporateTargetFrequencies: Calculate new scores based on optimized
//     target frequencies
//
void ScoringMatrix::IncorporateTargetFrequencies( const Pslvector& optfreqs )
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

    consv = optfreqs.Sum();
    if( consv < 1.0 - accuracy || consv > 1.0 + accuracy ) {
        warning( "Target frequencies after optimization not conserved." );
        return;
    }

    ind = 0;
    for( n = 0, nn = 0; n < GetQuerySize(); n++ )
    {
//         if( freq_fst_[n] == X )
//             continue;

        proquery = rprobs_.GetValueAt( nn++ );

        for( m = 0, mm = 0; m < GetSubjectSize(); m++ ) {
//             if( freq_sec_[m] == X )
//                 continue;

            prosbjct = cprobs_.GetValueAt( mm++ );

            tfval = optfreqs.GetValueAt( ind++ );

            if( proquery <= 0.0 || prosbjct <= 0.0 )
                throw myruntime_error( "Invalid probabilities after optimizing of target frequencies." );

            score = log( tfval / ( proquery * prosbjct ));

            loc_score = GetScore( m, n );

            if( loc_score <= SCORE_MIN )
                continue;

//             if( GetMaskedToIgnore( m, n ))
//                 continue;

            SetScore( score, m, n );
        }
    }
}



// =========================================================================
// PRINT ROUTINES
//
// PrintParameterTable: print table of computed statistical parameter
//     values to stream
//
void ScoringMatrix::PrintParameterTable( TPrintFunction print_func, void* vpn ) const
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
//    print_func( vpn, "%-25s  %6.4f   %6.4f\n", "Reference ungapped,", GetRefK(),            GetRefLambda());
//    print_func( vpn, "%-25s  %6.4f   %6.4f\n", "Reference gapped,",   GetExpGappedK(),      GetExpGappedLambda());
    print_func( vpn, "%-25s  %6s   %6s\n", "Computed  ungapped,", _sk, _slambda );
    print_func( vpn, "%-25s  %6s   %6s\n", "Estimated gapped,", _sdk, _sdlambda );

    print_func( vpn, "Entropy, %6.4f; Expected, %6.4f; Min/Max, %d/%-d\n\n",
                GetEntropy(), GetExpectedScore(), GetMinScore(), GetMaxScore());
}

// PrintFinal: nothing to print in the final
//
void ScoringMatrix::PrintFinal( TPrintFunction print_func, void* vpn ) const
{
    if( vpn == NULL )
        return;

    PrintReferenceParameterTable( print_func, vpn );
}

// -------------------------------------------------------------------------
// PrintScoringMatrix: output computed scoring system
//
void ScoringMatrix::PrintScoringMatrix( FILE* fp )
{
    int l = 0;
    int m, n;
    const double* const queryentrops = GetQueryInfContent();
    const double* const sbjctentrops = GetSbjctInfContent();

    if( fp == NULL )
        return;

    fprintf( fp,"%12c Position-specific profile scoring matrix\n", 32 );

    fprintf( fp, "%9c", 32 );
    for( m = 0; m < GetSubjectSize(); m++ )
        fprintf( fp, " %4d", m + 1 );

    fprintf( fp, "\n%9c", 32 );
    for( m = 0; m < GetSubjectSize(); m++ )
        fprintf( fp, " %4c", DehashCode( logo_sec_[m] ));

    for( n = 0; n < GetQuerySize(); n++ ) {

        fprintf( fp, "\n%5d %c   ", ++l, DehashCode( logo_fst_[n] ));

        for( m = 0; m < GetSubjectSize(); m++ )
            fprintf( fp, "%4d ", ( int )rint( GetImageScore/*GetScoreBase*/( m, n ) ));

        if( queryentrops )
            fprintf( fp, " %4.2f", queryentrops[n] );
    }

    if( sbjctentrops ) {
        fprintf( fp, "\n\n          " );
        for( m = 0; m < GetSubjectSize(); m++ )
            fprintf( fp, "%4.2f ", sbjctentrops[m] );
    }

    fprintf( fp, "\n" );
}

