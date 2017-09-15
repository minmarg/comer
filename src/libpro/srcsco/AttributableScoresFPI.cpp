/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "ext/psl.h"
#include "liblib/rc.h"
#include "liblib/data.h"
#include "libpro/srcsco/AbstractScoreMatrix.h"
#include "AttributableScoresFPI.h"


// -------------------------------------------------------------------------
// constructor: invoke parent constructor and initialize the class members
// -------------------------------------------------------------------------

AttributableScoresFPI::AttributableScoresFPI(
        AbstractScoreMatrix*    prnt,
        TProbabilityFunction    probfun,
        TScoreReductFunction    redtfun,
        TrogressFunction        progressfun,
        double                  reflambda,
        double                  refH,
        double                  refK,
        const double**          im,
        const TMask**           msk,
        bool                    keep,
        int                     as_factor )
:
    AttributableScores(
        prnt,
        probfun,
        redtfun,
        progressfun,
        reflambda,
        refH,
        refK,
        im,
        msk,
        keep,
        as_factor
    ),
    scores( NULL )
{
}

// -------------------------------------------------------------------------
// default construction
// -------------------------------------------------------------------------

AttributableScoresFPI::AttributableScoresFPI()
:
    AttributableScores(),
    scores( NULL )
{
}

// -------------------------------------------------------------------------
// Init: initialization of integer score table
// -------------------------------------------------------------------------

void AttributableScoresFPI::Init( int sz_query, int sz_sbjct )
{
    if( GetKeepInMemory()) {

        SetQuerySize( sz_query );
        SetSubjectSize( sz_sbjct );

        AttributableScores::Init(); //query and subject size have to be initialized before

        scores = ( double** )malloc( sizeof( double* ) * GetSubjectSize());

        if( !scores ) 
            throw myruntime_error( mystring( "AttributableScoresFPI: Not enough memory." ));


        for( int m = 0; m < GetSubjectSize(); m++ )
        {
            scores[m] = ( double* )malloc( sizeof( double ) * GetQuerySize());

            if( !scores[m] )
                throw myruntime_error( mystring( "AttributableScoresFPI: Not enough memory." ));

            memset( scores[m], 0, sizeof( double ) * GetQuerySize());
        }
    }
}

// -------------------------------------------------------------------------
// destructor: deallocate memory used by this class
// -------------------------------------------------------------------------

AttributableScoresFPI::~AttributableScoresFPI()
{
    if( scores ){
        for( int m = 0; m < GetSubjectSize(); m++ )
            free( scores[m] );
        free( scores );
        scores = NULL;
    }
}

// -------------------------------------------------------------------------
// ScaleToAttainLambda: perform iterative scaling until required lambda is
//     attained; binary search is much more precise since otherwise
//     operating with integer scores carries errors and exact solution with
//     adjusting scores by multiplying them by ratio of lambdas is rarely
//     possible. However, binary search is much slower
// -------------------------------------------------------------------------

void AttributableScoresFPI::ScaleToAttainLambda()
{
    double  multiplier = 1.0;
    ComputeProbabilitiesAndLambda();
    AdjustScoresInMatrix( multiplier = GetLambda() / GetRefLambda(), ImageTable );
    //here it is possible to explicitly set lambda keeping in mind that
    //lambda almost always will be equal to given ratio, but we also
    //need to compute probabilities of integer scores (for K)
    ComputeProbabilitiesAndLambda();
    SetPrivateMultiplier( multiplier );
}

// -------------------------------------------------------------------------
// ComputeScoreProbabilities: compute probabilities to observe scores at
//     each position of the scoring system
//
void AttributableScoresFPI::ComputeScoreProbabilities()
{
    const double    accuracy = 1.e-6;
    char            strbuf[KBYTE];
    double          prob, consv;

    if( GetParent() == NULL || GetParentProbFunction() == NULL )
        throw myruntime_error( mystring( "AttributableScores: Unable to compute score probabilities." ));

    if( GetKeepInMemory()) {
        if( GetQuerySize() <= 0 || GetSubjectSize() <= 0 )
            throw myruntime_error( mystring( "AttributableScores: Unable to compute score probabilities." ));

        if( !GetScores() || !GetMask())
            throw myruntime_error( mystring( "AttributableScores: Unable to compute score probabilities." ));

        SetMinMaxScores( 0, GetSubjectSize() * GetQuerySize());
    }

    ( GetParent()->*GetParentProbFunction())( this );

    consv = 0.0;
    for( int sc = GetMinScore(); sc <= GetMaxScore(); sc ++ ) {
        prob = GetProbabilityOf( sc );
        if( !prob )
            continue;
        consv += prob;
    }
    if( consv < 1.0 - accuracy || consv > 1.0 + accuracy ) {
        sprintf( strbuf, "Score probabilities not conserved: %f", consv );
        throw myruntime_error( strbuf );
    }
}

// -------------------------------------------------------------------------
// ComputeEntropyGivenLambda: compute relative entropy given previously
//     computed lambda.
// Relative entropy or information per aligned pair of positions:
//
//           _                   lambda s(k)
//   lambda \  s(k) * p(s(k)) * e
//          /_
//          k
//
// where p(s(k)) is a probability of the score s(k)
// -------------------------------------------------------------------------

void AttributableScoresFPI::ComputeEntropyGivenLambda()
{
    double  lmbda = GetLambda();

    if( lmbda < 0.0 ) {
        warning( "AttributableScoresFPI: Unable to compute relative entropy." );
        return;
    }

    double  y = exp( lmbda );
    double  r = 0.0;
    double  H = 0.0;
    size_t  p = 0;
    double  loc_score, prob;
    bool    term = false;
    int n, m, ind;

    for( n = 0; n < GetQuerySize() && !term; n++ ) {
        for( m = 0; m < GetSubjectSize(); m++ )
        {
            if( GetMaskedToIgnore( m, n ))
                continue;

            loc_score = GetScore( m, n );
            if( loc_score <= SCORE_MIN )
                continue;

            ind = n * GetSubjectSize() + m;
            prob = GetProbabilityOf(( TScore )ind );

//             r = loc_score * prob * exp( lmbda * loc_score );

            if( prob < 0.0 )
                throw myruntime_error("AttributableScoresFPI: Unable to calculate lambda: Domain error.");
            if( prob == 0.0 )
                r = SLC_LOG_DBL_MIN - 1.0;
            else {
                r = log( prob );
                r += lmbda * loc_score;
            }

            if( SLC_LOG_DBL_MAX < r ) {
                H = SLC_DBL_MAX;
                term = true;
                break;
            }
            r = loc_score * exp( r );
            if( SLC_DBL_MAX - H < r ) {
                H = SLC_DBL_MAX;
                term = true;
                break;
            }
            H += r;

            p++;
        }
    }

    if( p == 0 )
        throw myruntime_error( mystring( "AttributableScoresFPI: Failed to compute relative entropy: Illegal operation." ));

    if( H && lmbda < SLC_DBL_MAX / H )
        H *= lmbda;
    SetH( H );
}

// -------------------------------------------------------------------------
// ComputePositionalEntropy: compute entropy at specified position of score
//     system; to compute entropy, use all needed data given by parameters
// -------------------------------------------------------------------------

double AttributableScoresFPI::ComputePositionalEntropy(
    double lambda,
    int minscore,
    int maxscore,
    int pos,
    TGetProbFunction getProbfunction,
    TGetProbNormFunction getProbNormfunction ) const
{
    if( getProbfunction == NULL )
        throw myruntime_error( mystring( "AttributableScoresFPI: ComputePositionalEntropy: Memory access error." ));

    if( lambda < 0.0 )
        throw myruntime_error( mystring( "AttributableScoresFPI: ComputePositionalEntropy: Unable to compute entropy." ));

    double  N = ( this->*getProbNormfunction )( pos );
    double  y = exp( -lambda );
    double  m = pow( y, maxscore ); // e to -lambda * maxscore
    double  H = 0.0;

    if( N <= 0.0 )
        N = 0.0;
    else
        N = log( N ) / lambda;

    for( int sc = minscore; sc <= maxscore; sc++ )
        H = H * y + sc * ( this->*getProbfunction )( sc, pos );

    if( 0.0 < m )
        H /= m;
    else
        if( 0.0 < H )
            H = exp( lambda * maxscore + log( H ));

    return( lambda * H );
}

// -------------------------------------------------------------------------
// ComputeKarlinsK: calculate statistical parameter K
//
void AttributableScoresFPI::ComputeKarlinsK()
{
    double  loc_lambda = GetLambda();
    double  loc_expected_score = GetExpectedScore();

    if( 0.0 <= loc_expected_score ) {
        warning( "Unable to compute statistical parameter K." );
        return;
    }
    if( loc_lambda < 0.0 || GetH() <= 0.0 ) {
        warning( "Unable to compute statistical parameter K." );
        return;
    }

    double  K = 1.0;

    SetK( K );
}

// -------------------------------------------------------------------------
// FindLambdaRootAtPosition: helper method to find root of lambda given
//     parameters of positional scores
// -------------------------------------------------------------------------

double AttributableScoresFPI::FindLambdaRootAtPosition(
        TConservationFunction fdfunction,
        int gcd,
        double x1,
        double x2,
        double xacc,
        int maxit,
        int minscore,
        int maxscore,
        int pos,
        TGetProbFunction getProbfunction )
{
    char    data[BUF_MAX];
    char*   privpar = data;
    double  rts = -1.0;

    *(( int* )privpar ) = minscore; privpar += sizeof( int );
    *(( int* )privpar ) = maxscore; privpar += sizeof( int );
    *(( int* )privpar ) = pos;      privpar += sizeof( int );
    *(( TGetProbFunction* )privpar ) = getProbfunction;

    x1 = exp( -x1 );
    x2 = exp( -x2 );

    rts = rootByNRandBisection( fdfunction, gcd, x1, x2, xacc, maxit, data, false/*warn*/ );

    if( rts < 0.0 )
        return rts;

    return -log( rts ) / gcd;
}

// -------------------------------------------------------------------------
// findLambdaRoot: helper method to find root of lambda
// -------------------------------------------------------------------------

double AttributableScoresFPI::findLambdaRoot( TConservationFunction fdfunction, int gcd,
                            double x1, double x2, double xacc, int maxit, void* params )
{
    double  rts = -1.0;

    rts = rootByNRandBisection( fdfunction, gcd, x1, x2, xacc, maxit, params, false );

    if( rts < 0.0 )
        return rts;

    return rts;
}

// -------------------------------------------------------------------------
// findConstantRoot: helper method to find root of constant to adjust
//     relative entropy
// -------------------------------------------------------------------------

double AttributableScoresFPI::findConstantRoot( TConservationFunction fdfunction, int gcd,
                            double x1, double x2, double xacc, int maxit, void* params )
{
    double  rts = -1.0;

    rts = rootByNRandBisection( fdfunction, gcd, x1, x2, xacc, maxit, params );

    if( rts < 0.0 )
        return rts;

    return rts;
}

// -------------------------------------------------------------------------
// conservation: function describing probability conservation equation
//     which is used to find the root of lambda
// When double is used as type of scores, then the conservation function is
// computed as
//
//     1     __  lambda s(m,n)
//   ------ \   e             - 1 = 0
//   l1 l2  /__
//          m,n
//
// -------------------------------------------------------------------------

void AttributableScoresFPI::conservation( double x, double* f, double* df, int /*step*/, void* )
{
    double  y = exp( x );
    double  r = 0.0;
    double  h = 0.0;
    double  g = 0.0;
    size_t  p = 0;
    double  loc_score, prob;
    bool    term = false;
    int m, n, ind;

    for( n = 0; n < GetQuerySize() && !term; n++ ) {
        for( m = 0; m < GetSubjectSize(); m++ )
        {
            if( GetMaskedToIgnore( m, n ))
                continue;

            loc_score = GetScore( m, n );
            if( loc_score <= SCORE_MIN )
                continue;

            ind = n * GetSubjectSize() + m;
            prob = GetProbabilityOf(( TScore )ind );

//             r = prob * pow( y, loc_score );
//             h += r;
//             g += r * loc_score;

            if( prob < 0.0 )
                throw myruntime_error("AttributableScoresFPI: Unable to calculate lambda: Domain error.");
            if( prob == 0.0 )
                r = SLC_LOG_DBL_MIN - 1.0;
            else {
                r = log( prob );
                r += x * loc_score;
            }

            if( SLC_LOG_DBL_MAX < r + loc_score ) {
                h = g = SLC_DBL_MAX;
                term = true;
                break;
            }
            r = exp( r );
            if( SLC_DBL_MAX - g < r * loc_score ) {
                h = g = SLC_DBL_MAX;
                term = true;
                break;
            }
            h += r;
            g += r * loc_score;
            p++;
        }
    }

    if( p == 0 )
        throw myruntime_error( mystring( "AttributableScoresFPI: Failed to compute lambda: Illegal operation." ));

    *f = h - 1.0;     //function value
    *df= g;           //derivative value

}

// -------------------------------------------------------------------------
// lambdaconspos: alternative lambda conservation function which uses
//     data given by params to positional properties of scores
//
void AttributableScoresFPI::lambdaconspos( double x, double* f, double* df, int step, void* params )
{
    if( params == NULL )
        throw myruntime_error( mystring( "AttributableScoresII: lambdaconspos: Wrong parameters." ));

    char*   privpar = ( char* )params;
    int     minscore = *(( int* )privpar ); privpar += sizeof( int );
    int     maxscore = *(( int* )privpar ); privpar += sizeof( int );
    int     pos      = *(( int* )privpar ); privpar += sizeof( int );
    TGetProbFunction    getProbfunction = *(( TGetProbFunction* )privpar );

    double  y = x;      //initially x is supposed to be exp(-lambda)
    double  ldf = 0.0;
    double  lf = 0.0;

    for( int sc = minscore; sc <= maxscore; sc += step ) {
        ldf = ldf * y + lf;
        lf  = lf  * y + ( this->*getProbfunction )( sc, pos );
        //zero will always be reached
        if( sc == 0 )
            lf -= 1.0;
    }
    *df = ldf;
    *f  = lf;
}

// -------------------------------------------------------------------------
// relentropy_conservation: this method describes equation for relative
//     entropy, however modified so that scores to be transformed to
//     achieve a specific value of relative entropy.
// It is known that relative entropy originally is computed as follows:
//
//           _                   lambda s(k)
//   lambda \  s(k) * p(s(k)) * e
//          /_
//          k
//
// where p(s(k)) is a probability of the score s(k).
// If, on the other hand, we allow transformation of scores, relative
// entropy becomes:
//
//    1   _                                      lambda s(k)
//   --- \  (lambda s(k) - 2 log c) * p(s(k)) * e
//     2 /_
//   c   k
//
// where c is a constant to be found in order to obtain the desired relative
// entropy value. Derivative of this expression with respect to c is
//
//   -2   _                                          lambda s(k)
//   --- \  (lambda s(k) - 2 log c + 1) * p(s(k)) * e
//     3 /_
//   c   k
//
// Having c, new scores sn(k) are related to s(k) by relation
//
//   s(k) = 2 / lambda log(c) + sn(k)
// -------------------------------------------------------------------------

void AttributableScoresFPI::relentropy_conservation( double /*x*/, double* /*f*/, double* /*df*/, int /*step*/, void* )
{
    throw myruntime_error( mystring(
        "AttributableScoresFPI: "
        "Floating-point computations for adjustment of relative entropy is not implemented yet." ));
}

