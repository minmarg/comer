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
#include "AbstractScoreMatrix.h"
#include "AttributableScoresII.h"
#include "AttributableScoresFPI.h"


size_t AbstractScoreMatrix::deltaLength = 0;    //initialize length adjustment to zero at first
Uint64 AbstractScoreMatrix::searchSpace = 0;    //equal search space to zero at initialization
Uint64 AbstractScoreMatrix::raw_search_space = 0;//raw search space

// -------------------------------------------------------------------------
// constructor: initializing members to the default values
// -------------------------------------------------------------------------

AbstractScoreMatrix::AbstractScoreMatrix(
        TType           type,
        Configuration   config[NoSchemes],
        TBehaviour      s_behaviour,
        TScaling        a_scaling,
        TMask           c_masking )
:
    matrixType( type ),
    behaviour( s_behaviour ),
    commonmasking( c_masking ),
    name_( NULL ),

    image( NULL ),
    modimage_( NULL ),
    usemodimage_( false ),
    mask( NULL ),
    corres_scores( NULL ),
    scaled_scores( NULL ),

    scadjrelent_( false ),
    scadjcorrsc_( false ),

    deletioncoeff( 1.0 ),

    information_content( 0.0 ),
    ainfoupper2nd( 0.0 ),
    ainfonumerator2nd( 1.0 ),
    ainfoscale2nd( 1.0 ),
    ainfonumeratoralt( 1.0 ),
    ainfoscalealt( 1.0 ),

    score_multiplier( 1.0 ),

    prescore( 0.0 ),
    prexpect( 10.0 ),

    min_score( 0 ),
    max_score( 0 ),

    auto_scaling( a_scaling ),
    auto_scaling_constant( PP_SCALE_CONSTANT ),

    expscore( 0.0 ),
    allnegatives( false ),
    derivedGappedLambda( -1.0 ),
    derivedGappedK( -1.0 ),

    lambda( -1.0 ),
    entropy( 0.0 ),
    parameterK( -1.0 ),

    querySize( 0 ),
    subjectSize( 0 ),

    supportoptfreq_( false ),

    configuration( config )
{
    InitializeSSParameters();
}

// -------------------------------------------------------------------------
// Init: data initialization
// -------------------------------------------------------------------------

void AbstractScoreMatrix::Init( int sz_query, int sz_sbjct )
{
    SetQuerySize( sz_query );
    SetSubjectSize( sz_sbjct );

    bool    keep_in_memory = ( 0 < GetQuerySize() ) && ( 0 < GetSubjectSize());

    if( keep_in_memory ) {
        image = ( double** )malloc( sizeof( double* ) * subjectSize );
        mask = ( TMask** )malloc( sizeof( TMask* ) * subjectSize );

        if( GetUseModImage()) {
            modimage_ = ( double** )malloc( sizeof( double* ) * subjectSize );
            if( !modimage_ )
                throw myruntime_error("AbstractScoreMatrix: Not enough memory.");
        }

        if( !image || !mask )
            throw myruntime_error("AbstractScoreMatrix: Not enough memory.");

        for( int m = 0; m < subjectSize; m++ ) {
            image[m] = ( double* )malloc( sizeof( double ) * querySize );
            mask[m] = ( TMask* )malloc( sizeof( TMask ) * querySize );

            if( GetUseModImage()) {
                modimage_[m] = ( double* )malloc( sizeof( double ) * querySize );
                if( !modimage_[m])
                    throw myruntime_error("AbstractScoreMatrix: Not enough memory.");
                memset( modimage_[m], 0, sizeof( double ) * querySize );
            }

            if( !image[m] || !mask[m])
                throw myruntime_error("AbstractScoreMatrix: Not enough memory.");

            memset( image[m], 0, sizeof( double ) * querySize );
            memset( mask[m], 0, sizeof( TMask ) * querySize );
        }
    }

    if( GetFPScaling())
        corres_scores = new AttributableScoresFPI(
                this,
                &AbstractScoreMatrix::ComputeScoreProbabilities,
                &AbstractScoreMatrix::ReduceScores,
                &AbstractScoreMatrix::MultiplierProgress,
                GetRefLambda(),
                GetRefH(),
                GetRefK(),
                GetImage(),
                GetMask(),
                keep_in_memory
        );
    else
        corres_scores = new AttributableScoresII(
                this,
                &AbstractScoreMatrix::ComputeScoreProbabilities,
                &AbstractScoreMatrix::ReduceScores,
                &AbstractScoreMatrix::MultiplierProgress,
                GetRefLambda(),
                GetRefH(),
                GetRefK(),
                GetImage(),
                GetMask(),
                keep_in_memory
        );

    if( !corres_scores )
        throw myruntime_error( mystring( "AbstractScoreMatrix: Not enough memory." ));

    corres_scores->Init( GetQuerySize(), GetSubjectSize());

    if( GetAutoScaling()) {
        scaled_scores = new AttributableScoresII(
            this,
            &AbstractScoreMatrix::ComputeScoreProbabilities,
            &AbstractScoreMatrix::ReduceScores,
            &AbstractScoreMatrix::MultiplierProgress,
            GetRefLambda() / GetAutoScalingFactor(),    //divide lambda by a factor and...
            GetRefH(),
            GetRefK(),
            GetImage(),
            GetMask(),
            keep_in_memory,
            GetAutoScalingFactor()  //use the factor to scale scores in order to gain precision
        );

        if( !scaled_scores )
            throw myruntime_error( mystring( "AbstractScoreMatrix: Not enough memory." ));

        scaled_scores->Init( GetQuerySize(), GetSubjectSize());
    }
}

// -------------------------------------------------------------------------
// default constructor is invalid for initialization
// -------------------------------------------------------------------------

AbstractScoreMatrix::AbstractScoreMatrix()
:
    matrixType( Dummy ),
    behaviour( StatisticsGiven ),
    commonmasking( Unmasked ),
    name_( NULL ),

    image( NULL ),
    modimage_( NULL ),
    usemodimage_( false ),
    mask( NULL ),
    corres_scores( NULL ),
    scaled_scores( NULL ),

    scadjrelent_( false ),
    scadjcorrsc_( false ),

    information_content( 0.0 ),
    ainfoupper2nd( 0.0 ),
    ainfonumerator2nd( 1.0 ),
    ainfoscale2nd( 1.0 ),
    ainfonumeratoralt( 1.0 ),
    ainfoscalealt( 1.0 ),

    score_multiplier( 1.0 ),

    min_score( 0 ),
    max_score( 0 ),

    auto_scaling( NoScaling ),
    auto_scaling_constant( PP_SCALE_CONSTANT ),

    expscore( 0.0 ),
    allnegatives( false ),
    derivedGappedLambda( -1.0 ),
    derivedGappedK( -1.0 ),

    lambda( -1.0 ),
    entropy( 0.0 ),
    parameterK( -1.0 ),

    querySize( 0 ),
    subjectSize( 0 ),

    supportoptfreq_( false ),

    configuration( NULL )
{
    throw myruntime_error(
        mystring( "AbstractScoreMatrix: Default initialization is not allowed." ));
}

// -------------------------------------------------------------------------
// destructor: deallocate memory used by this class
// -------------------------------------------------------------------------

AbstractScoreMatrix::~AbstractScoreMatrix()
{
    if( corres_scores )
        delete corres_scores;

    if( scaled_scores )
        delete scaled_scores;

    if( image ){
        for( int m = 0; m < subjectSize; m++ )
            if( image[m] )
                free( image[m] );
        free( image );
        image = NULL;
    }

    if( modimage_ ){
        for( int m = 0; m < subjectSize; m++ )
            if( modimage_[m] )
                free( modimage_[m] );
        free( modimage_ );
        modimage_ = NULL;
    }

    if( mask ){
        for( int m = 0; m < subjectSize; m++ )
            if( mask[m] )
                free( mask[m] );
        free( mask );
        mask = NULL;
    }
}

// -------------------------------------------------------------------------
// InitializeStatisticalSignificanceParameters: initialize statistical
//     significance parameters given scoring scheme
// -------------------------------------------------------------------------

void AbstractScoreMatrix::InitializeSSParameters()
{
    const Configuration&    ungapped_config = GetConfiguration( ProcomUngapped );
    const Configuration&    gapped_config = GetConfiguration( ProcomGapped );

    SetRefLambda(   ungapped_config.GetLambda());
    SetRefH(        ungapped_config.GetH());
    SetRefK(        ungapped_config.GetK());
    SetExpGappedLambda( gapped_config.GetLambda());
    SetExpGappedK(      gapped_config.GetK());

    if( GetBehaviour() == StatisticsGiven ) {
        //all necessary parameters are assumed to have been
        //precomputed earlier and saved to configuration
        SetLambda(  ungapped_config.GetLambda());
        SetEntropy( ungapped_config.GetH());
        SetK(       ungapped_config.GetK());
        //artificially set expected score to indicate valid statistics
        SetExpectedScore( -1.0 );

        SetDerivedGappedLambda( gapped_config.GetLambda());
        SetDerivedGappedK(      gapped_config.GetK());

        SetMultiplier( ungapped_config.GetScaleFactor());
    }
}

// -------------------------------------------------------------------------
// ObtainReferenceParameters: obtain statistical reference parameters given
//     gap cost scheme
// -------------------------------------------------------------------------

void AbstractScoreMatrix::ObtainReferenceParameters( int ss,
        double* refLambda, double* refK, double* refH, double* refAlpha, double* refBeta,
        double* expLambda, double* expK, double* expH, double* expAlpha, double* expBeta )
{
    if( refLambda ) *refLambda  = LOSCORES.StatisParam( Ungapped, Lambda );
    if( refK )      *refK       = LOSCORES.StatisParam( Ungapped, K );
    if( refH )      *refH       = LOSCORES.StatisParam( Ungapped, H );
    if( refAlpha )  *refAlpha   = LOSCORES.StatisParam( Ungapped, alpha );
    if( refBeta )   *refBeta    = LOSCORES.StatisParam( Ungapped, beta );

    if( expLambda ) *expLambda  = LOSCORES.StatisParam( ss, Lambda );
    if( expK )      *expK       = LOSCORES.StatisParam( ss, K );
    if( expH )      *expH       = LOSCORES.StatisParam( ss, H );
    if( expAlpha )  *expAlpha   = LOSCORES.StatisParam( ss, alpha );
    if( expBeta )   *expBeta    = LOSCORES.StatisParam( ss, beta );
}


// -------------------------------------------------------------------------
// ScanForHSPs: performs verification of multiple high-scoring pairs (hsp)
//     in the same diagonal of the score system
// minhspscore,- minimum score of hsp
// hsplen,- length of hsp
// nohsps,- minimum number of hsps required
// maxdist,- max distance between hsps
// possbjct, posquery,- subject and query positions of hsps found in the
//     same diagonal
// returns true if such hsps have been found
// -------------------------------------------------------------------------

bool AbstractScoreMatrix::ScanForHSPs(
    double minhspscore, int hsplen, int nohsps, int maxdist,
    int* possbjct, int* posquery )
{
    bool    bret = false;

    if( !IsValid())
        throw myruntime_error( mystring( "Scan of HSPs failed: no score matrix." ));
    if( !GetScaledScores() || !GetCorresScores())
        throw myruntime_error( mystring( "ScanForHSPs: Memory access error." ));

    if( GetAutoScaling()) {
        GetScaledScores()->MultiplyScoresByPrivateMultiplier();
        bret = GetScaledScores()->SearchForHSPs(
            minhspscore * GetScaledScores()->GetAutoScalingFactor(),
            hsplen, nohsps, maxdist,
            possbjct, posquery );
    } else {
//         GetCorresScores()->MultiplyScoresByPrivateMultiplier();
        bret = GetCorresScores()->SearchForHSPs( minhspscore, hsplen, nohsps, maxdist, possbjct, posquery );
    }

    return bret;
}

// -------------------------------------------------------------------------
// ComputeStatisticalParameters: implies computation of all the required
//     statistical parameters
// -------------------------------------------------------------------------

void AbstractScoreMatrix::ComputeStatisticalParameters()
{
    if( !IsValid())
        throw myruntime_error( mystring( "AbstractScoreMatrix: Unable to compute statistical parameters." ));

    GetCorresScores()->ComputeStatisticalParameters( true, !GetAutoScaling());

    SetLambda(          GetCorresScores()->GetLambda());
    SetEntropy(         GetCorresScores()->GetH());
    SetK(               GetCorresScores()->GetK());
    SetExpectedScore(   GetCorresScores()->GetExpectedScore());

    if( GetAutoScaling()) {
        GetScaledScores()->ComputeProbabilitiesAndLambda();

        if( 0.0 < GetScaledScores()->GetLambda())
            SetLambda( GetScaledScores()->GetLambda() * GetScaledScores()->GetAutoScalingFactor());
        SetExpectedScore( GetScaledScores()->GetExpectedScore() / GetScaledScores()->GetAutoScalingFactor());
    }

    SetMinScore( GetCorresScores()->GetMinScore());
    SetMaxScore( GetCorresScores()->GetMaxScore());

    DeriveGappedParameters();
}

// -------------------------------------------------------------------------
// ScaleScoringMatrix: scales scoring matrix so that its scores to be
//     comparable to a reference scoring system.
// -------------------------------------------------------------------------

void AbstractScoreMatrix::ScaleScoringMatrix()
{
    if( GetBehaviour() == StatisticsGiven ) {
        ComputeStatisticalParameters();//NOTE!!
        //Nothing doing since all the required statistical
        //parameters are assumed to be known from the configuration
        return;
    }

    if( !IsValid())
        throw myruntime_error( mystring( "Failed to scale matrix: Score matrix is not initialized." ));

//     TransformScoresToAchieveH();

    if( GetAutoScaling()) {
        //lambda is scaled according to the scaled scores
        //other parameters may be computed either on the scaled scores or on the corres. scores
        //
        GetScaledScores()->ScaleScoringMatrix();
//        GetScaledScores()->ComputeStatisticalParameters( false /*computelambda*/);  ///***
        GetCorresScores()->SetPrivateMultiplier( GetScaledScores()->GetPrivateMultiplier());
        GetCorresScores()->MultiplyScoresByPrivateMultiplier();
        GetCorresScores()->ComputeStatisticalParameters( true, false/*wrn*/);

        if( GetScaledScores()->GetLambda() < 0.0 )
            SetLambda(  GetCorresScores()->GetLambda());
        else
            SetLambda(  GetScaledScores()->GetLambda() * GetScaledScores()->GetAutoScalingFactor()); ///This...
///        SetLambda(          GetCorresScores()->GetLambda());    ///OR MAYBE THIS ONE

        SetEntropy(     GetCorresScores()->GetH());  //alternative
        SetK(           GetCorresScores()->GetK());  //alternative
///        //more accurate, but may take much more of computational time...
///        SetEntropy(         GetScaledScores()->GetH()); //H slightly differs from that for corres. scores (sys. error)
///        SetK(               GetScaledScores()->GetK()); //consequenctly, so does K

        SetExpectedScore(   GetScaledScores()->GetExpectedScore() / GetScaledScores()->GetAutoScalingFactor()); ///This...
///        SetExpectedScore(   GetCorresScores()->GetExpectedScore()); ///OR MAYBE THIS ONE

        SetMultiplier( GetScaledScores()->GetPrivateMultiplier());

    } else {
        GetCorresScores()->ScaleScoringMatrix();
        GetCorresScores()->ComputeStatisticalParameters( false/*computelambda*/ );

        SetLambda(          GetCorresScores()->GetLambda());
        SetEntropy(         GetCorresScores()->GetH());
        SetK(               GetCorresScores()->GetK());
        SetExpectedScore(   GetCorresScores()->GetExpectedScore());

        SetMultiplier( GetCorresScores()->GetPrivateMultiplier());
    }

    //We keep min and max scores from the representative system only
    SetMinScore( GetCorresScores()->GetMinScore());
    SetMaxScore( GetCorresScores()->GetMaxScore());

// PrintScoringMatrix( stderr );//***TEST***
    DeriveGappedParameters();
}

// -------------------------------------------------------------------------
// MultiplierProgress: method used to save (and show) progress of the values of score_multiplier
// -------------------------------------------------------------------------

void AbstractScoreMatrix::MultiplierProgress( AttributableScores* PATTR_SCORES )
{
    SetLambda(      PATTR_SCORES->GetLambda());
    SetMultiplier(  PATTR_SCORES->GetPrivateMultiplier());
}

// -------------------------------------------------------------------------
// SetInfoThresholdByEval: sets information content threshold by e-value
//     given
// -------------------------------------------------------------------------

void AbstractScoreMatrix::SetInfoThresholdByEval( double pair_eval )
{
    if( pair_eval < 0.0 )
        return;

    if( GetInfoCorrectionNumerator2nd() < 0.0 ||
        GetInfoCorrectionUpperBound2nd() < 0.0 )
        throw myruntime_error( mystring( "AbstractScoreMatrix: Invalid parameters of information correction." ));

    if( GetInfoCorrectionNumerator2nd() == 0.0 ||
        GetInfoCorrectionUpperBound2nd() == 0.0 )
        return;


    double  upper = GetInfoCorrectionUpperBound2nd();
    double  numerator = GetInfoCorrectionNumerator2nd();
    double  scale = GetInfoCorrectionScale2nd();
    double  numerator_alt = GetInfoCorrectionNumeratorAlt();
    double  scale_alt = GetInfoCorrectionScaleAlt();
    double  thresh = upper;
    double  thresh_alt = upper;
    double  value;

    if( 0.0 < pair_eval ) {
        value = log( pair_eval );
        //manage phase transition cases
        if( pair_eval == 1.0 || value == scale )
            thresh = upper;
        else if( exp( -scale_alt ) <= pair_eval )
            thresh = upper;
        else if( exp( scale ) <= pair_eval )
            thresh = upper;
        else {
            thresh += numerator /( value - scale );
            thresh_alt = -numerator_alt /( value + scale_alt );
            if( thresh < thresh_alt )
                thresh = thresh_alt;
            if( thresh < 0.0 )
                thresh = upper;//thresh = 0.0;
            if( upper < thresh )
                thresh = upper;
        }
    }

    SetInformationThreshold( thresh );
}

// -------------------------------------------------------------------------
// TransformScoresToAchieveH: first of all computes adjustment constant for
//     relative entropy so that transformation of scores by using it ensures
//     relative entropy desired to achieve; afterwards transforms scores
// -------------------------------------------------------------------------

void AbstractScoreMatrix::TransformScoresToAchieveH()
{
    if( GetCorresScores() == NULL )
        return;

    GetCorresScores()->ComputeHConstant();
    double  constant = GetCorresScores()->GetConstantForHAdjustment();

    if( constant < 0.0 ) {
        warning( "AbstractScoreMatrix: Unable to transform scores to get entropy desired." );
        return;
    }

    SetLambda(          GetCorresScores()->GetLambda());
    SetExpectedScore(   GetCorresScores()->GetExpectedScore());

    TransformScoresInMatrix( constant );
}

// -------------------------------------------------------------------------
// TransformScoresInMatrix: transform scores in the matrix given previously
//     computed constant c used in the adjusted relative entropy expression.
// Transformed scores sn(k) are related to the original ones s(k) by the
// relation:
//
//   lambda s(k) = 2 log(c) + lambda sn(k)
// -------------------------------------------------------------------------

void AbstractScoreMatrix::TransformScoresInMatrix( double c )
{
#ifdef __DEBUG__
    if( !image || !mask )
        throw myruntime_error( mystring( "AbstractScoreMatrix: Unable to transform scores." ));
#endif

    if( c <= 0.0 )
        return;

    double  scaling = GetLambda();

    if( scaling <= 0.0 )
        throw myruntime_error( mystring( "AbstractScoreMatrix: Failed to transform scores." ), SCALING );

    double  term = ( log( c ) + log( c )) / scaling;

    for( int m = 0; m < subjectSize; m++ ) {
        for( int n = 0; n < querySize; n++ )
        {
            if( SCORE_MIN < image[m][n] /*&& ! GetMaskedToIgnore( m, n )*/) {
                image[m][n] -= term;
            }
        }
    }
}

// -------------------------------------------------------------------------
// DeriveGappedParameters: estimate statistical gapped parameters given
//     experimentally defined values and computed ungapped ones
// -------------------------------------------------------------------------

void AbstractScoreMatrix::DeriveGappedParameters()
{
    //(computed ungapped for the scoring system) * (experimentally defined gapped) / (analitically computed ungapped)
    //
    if( 0.0 < GetLambda())
//         SetDerivedGappedLambda( GetExpGappedLambda());
        SetDerivedGappedLambda( GetLambda() * GetExpGappedLambda() / GetRefLambda());
    if( 0.0 < GetK())
        SetDerivedGappedK( GetK() * GetExpGappedK() / GetRefK());
}

// -------------------------------------------------------------------------
// ComputeExpectation: computes e-value given previously computed
//     statistical parameters
// -------------------------------------------------------------------------

double AbstractScoreMatrix::ComputeExpectation(
        double score,
        double* ref_expect,
        double* pure_expect,
        double* pair_expect,
        double* bitscore ) const
{
    double  expect = -1.0;
    double  pairsspace = 0.0;

    if( GetLambda() <= 0.0 || GetK() < 0.0 ) {
        if( ref_expect )
            //compute expectation for ungapped alignment (trying to predict approx. e-value
            //  assuming the score of ungapped alignment is given: expected score per position is positive)
            *ref_expect = GetExpGappedK() * GetSearchSpace() * exp( -GetExpGappedLambda() * score );

        if( pair_expect ) {
            if( GetDeltaLength() < GetQuerySize() && GetDeltaLength() < GetSubjectSize())
                pairsspace = ( GetQuerySize() - GetDeltaLength()) * ( GetSubjectSize() - GetDeltaLength());
            else
                pairsspace = GetQuerySize() * GetSubjectSize();
            *pair_expect = GetRefK() * pairsspace * exp( -GetRefLambda() * score );
        }

        if( pure_expect )
            *pure_expect = GetRefK() * GetSearchSpace() * exp( -GetRefLambda() * score );

        if( bitscore && 0.0 < GetExpGappedK())
            *bitscore = ( GetExpGappedLambda() * score - log( GetExpGappedK())) / LN2;

        return expect;
    }

    if( GetBehaviour() == StatisticsGiven )
        expect = GetExpGappedK() * GetSearchSpace() * exp( -GetExpGappedLambda() * score );
    else
        expect = GetDerivedGappedK() * GetSearchSpace() * exp( -GetDerivedGappedLambda() * score );

    if( bitscore && 0.0 < GetExpGappedK())
        *bitscore = ( GetDerivedGappedLambda() * score - log( GetDerivedGappedK())) / LN2;

    if( pair_expect ) {
        if( GetDeltaLength() < GetQuerySize() && GetDeltaLength() < GetSubjectSize())
            pairsspace = ( GetQuerySize() - GetDeltaLength()) * ( GetSubjectSize() - GetDeltaLength());
        else
            pairsspace = GetQuerySize() * GetSubjectSize();
        *pair_expect = GetK() * pairsspace * exp( -GetLambda() * score );
    }
    if( pure_expect )
        *pure_expect = GetK() * GetSearchSpace() * exp( -GetLambda() * score );

    return expect;
}

// -------------------------------------------------------------------------
// GetMeanLengthGivenExpectation: computes expected mean alignment length
//     for ungapped alignment given e-value
// -------------------------------------------------------------------------

double AbstractScoreMatrix::GetMeanLengthGivenExpectation( double expect, double* score ) const
{
    double  l_sc;   //lambda times score which is computed according to the given expect
    double  mlen;   //expected mean length

    if( expect < 0.0 || GetDerivedGappedK() < 0.0 || GetEntropy() < 0.0 )
        return 0.0;

    if( expect <= 0.0 || GetDerivedGappedK() <= 0.0 || GetEntropy() <= 0.0 ||
        GetSearchSpace() == 0 || GetRawSearchSpace() == 0 )
        return 1.0e5;

//     l_sc = log( GetDerivedGappedK() * GetSearchSpace() / expect );
    //compute mean length for ungapped alignments to make it longer
    //NOTE: most appropriately to be used here is effective search space for ungapped alignment searching...
    //since it is nearly equal to the raw search space, it is used here
    l_sc = log( GetK() * GetRawSearchSpace() / expect );
    mlen = l_sc / GetEntropy();

    if( score && 0.0 < GetLambda())
        *score = l_sc / GetLambda();

    return mlen;
}

// -------------------------------------------------------------------------
// ComputeLengthAdjustment: helper method for computation of length
//     adjustment to edge-effect correction
// -------------------------------------------------------------------------

bool AbstractScoreMatrix::ComputeLengthAdjustment(
    const Configuration& config, size_t query_len, Uint64 db_len, size_t no_sequences )
{
    return ComputeLengthAdjustment(
        config.GetLambda(), config.GetK(), config.GetAlpha(), config.GetBeta(),
        query_len, db_len, no_sequences );
}

// -------------------------------------------------------------------------
// ComputeLengthAdjustment: compute length adjustment for edge-effect
//     correction
// Altschul et al. in their paper in Nucleic Acids Res. 29 (2001) showed
// that expected alignment length can be expressed in linear regression of
// scores s: l(s) = alpha s + beta. Substituting this expression into the
// e-value expression K (m-l) (n-Nl) exp( -lambda s ), and solving it for l
// given e-value equals to 1, one obtains
//
//  l = alpha / lambda ln (K (m-l)(n-Nl)) + beta ,
//
// where m is query length, n is DB length and N is number of
// sequences in the DB. The equation is to be solved with constraint
//  K (m-l)(n-Nl) >= max{m,n}, and solution has to be integer value.
//
// E-value of 1 is not chosen by chance. The greater this number is,
// the smaller value of l is obtained and consequently, the greater
// e-value in result is obtained. Since length adjustment is computed for
// unrelated sequences, setting here e-value to one, that corresponds to
// p-value of 0.6, is logical.
// The computation algorithm is as given by E. Gertz (2005)
// -------------------------------------------------------------------------

bool AbstractScoreMatrix::ComputeLengthAdjustment(
        double lambda, double K,    double alpha,       double beta,
        size_t m/*query_len*/,      Uint64 n/*db_len*/, size_t N/*no_sequences*/ )
{
    SetRawSearchSpace( m * n );
    SetSearchSpace( m * n );
    SetDeltaLength( 0 );

    if( lambda <= 0.0 || K <= 0.0 || N == 0 ) {
        return false;
    }

    const int       maxit = LENGTH_ADJUSTMENT_MAXIT;
    const double    logK = log( K );
    const double    alratio = alpha / lambda;
    const double    max_mn = ( m > n )? m: n;

    double      length = 0.0;       //expected alignment length
    double      min_length = 0.0;   //lower bound of interval of alignment length
    double      max_length = 0.0;   //upper bound of interval of alignment length
    double      nxt_length = 0.0;   //next iterated value of alignment length
    double      sspace = 0.0;       //effective search space
    double      space2 = 0.0;       //effective search space
    bool        converged = false;  //whether an iteration converged to the solution of length

    //compute upper bound of interval, length can take a value from satisfying the constraint
    //we have quadratic equation of l: mn - max{m,n}/K - (n + mN)l + Nl*l
    double      a = N;                      //coefficient a of the quadratic equation
    double      b = n + m * N;              //coefficient -b
    double      c = m * n - max_mn / K;

    if( c < 0.0 ) {
        return false;
    }
    //take the smaller root of the equation
    //max_length = ( b - sqrt( SQUARE( b ) - 4.0 * a * c )) / ( 2.0 * a );
    //since plausible are the cases when b*b>>4ac, use alternative form for the solution
    max_length = 2.0 * c / ( b + sqrt( SQUARE( b ) - 4.0 * a * c ));


    for( int j = 0; j < maxit; j++ ) {
        length = nxt_length;
        sspace = ( m - length ) * ( n - N * length );
        nxt_length = alratio * ( logK + log( sspace )) + beta;

        if( length <= nxt_length ) {
            min_length = length;
            if( nxt_length - length <= 1.0 ) {
                converged = true;
                break;
            }
            if( min_length == max_length )
                break;
        } else
            max_length = length;

        if( nxt_length < min_length || max_length < nxt_length )
            //outside the range
            nxt_length = ( !j )? max_length : ( min_length + max_length ) / 2.0;
    }

    if( converged ) {
        //make sure that floor(min_length) + 1 != floor(min_length)
        length = ceil( min_length );
        if( length <= max_length ) {
            space2 = ( m - length ) * ( n - N * length );
            nxt_length = alratio * ( logK + log( space2 )) + beta;

            if( length <= nxt_length ) {
                min_length = length;
                sspace = space2;
            }
        }
    }
    //if not converged, save the closest value to the solution
    SetSearchSpace(( Uint64 )sspace );
    SetDeltaLength(( size_t )min_length );

    return converged;
}

// =========================================================================
// PRINT ROUTINES
//
// PrintReferenceParameterTable: print refrence parameter table to string
//     stream; space for stream must be PRE-ALLOCATED before!
//

void AbstractScoreMatrix::PrintReferenceParameterTable( char* sp ) const
{
    if( !sp )
        return;
    *sp = 0;  //to ensure printing to the end of the stream
    PrintReferenceParameterTable( &string_print, sp );
}

// PrintReferenceParameterTable: print reference parameter table to file
//

void AbstractScoreMatrix::PrintReferenceParameterTable( FILE* fp ) const
{
    PrintReferenceParameterTable( &file_print, fp );
}

// PrintReferenceParameterTable: print reference values of statistical
//     pareameters
//
void AbstractScoreMatrix::PrintReferenceParameterTable( TPrintFunction print_func, void* vpn ) const
{
    if( vpn == NULL )
        return;

    print_func( vpn, "%-20s  %-6s   %-6s\n", "Reference values of", "K", "Lambda" );
    print_func( vpn, "%-20s  %6.4f   %6.4f\n", "Ungapped", GetRefK(), GetRefLambda());
    print_func( vpn, "%-20s  %6.4f   %6.4f\n", "Gapped",   GetExpGappedK(), GetExpGappedLambda());
}

// PrintParameterTable: print table of computed statistical parameter values
// -------------------------------------------------------------------------

// PrintParameterTable: print attributes of the scoring system to string
//     stream; space for stream must be PRE-ALLOCATED before!
//

void AbstractScoreMatrix::PrintParameterTable( char* sp ) const
{
    if( !sp )
        return;
    *sp = 0;  //to ensure printing to the end of the stream
    PrintParameterTable( &string_print, sp );
}

// PrintParameterTable: print attributes of the scoring system to file
//

void AbstractScoreMatrix::PrintParameterTable( FILE* fp ) const
{
    PrintParameterTable( &file_print, fp );
}

// -------------------------------------------------------------------------
// PrintProbabilities: print score probabilities if available
//
// PrintFinal: final print to string stream; space for stream must be
//     PRE-ALLOCATED before!
//

void AbstractScoreMatrix::PrintFinal( char* sp ) const
{
    if( !sp )
        return;
    *sp = 0;  //to ensure printing to the end of the stream
    PrintFinal( &string_print, sp );
}

// PrintFinal: final print to file
//

void AbstractScoreMatrix::PrintFinal( FILE* fp ) const
{
    PrintFinal( &file_print, fp );
}

