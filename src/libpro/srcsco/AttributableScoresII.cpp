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


// -------------------------------------------------------------------------
// constructor: invoke parent constructor and initialize the class members
// -------------------------------------------------------------------------

AttributableScoresII::AttributableScoresII(
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

AttributableScoresII::AttributableScoresII()
:
    AttributableScores(),
    scores( NULL )
{
}

// -------------------------------------------------------------------------
// Init: initialization of integer score table
// -------------------------------------------------------------------------

void AttributableScoresII::Init( int sz_query, int sz_sbjct )
{
    if( GetKeepInMemory()) {

        SetQuerySize( sz_query );
        SetSubjectSize( sz_sbjct );

        AttributableScores::Init(); //query and subject size have to be initialized before

        scores = ( TScore** )malloc( sizeof( TScore* ) * GetSubjectSize());

        if( !scores ) 
            throw myruntime_error( mystring( "AttributableScoresII: Not enough memory." ));


        for( int m = 0; m < GetSubjectSize(); m++ )
        {
            scores[m] = ( TScore* )malloc( sizeof( TScore ) * GetQuerySize());

            if( !scores[m] )
                throw myruntime_error( mystring( "AttributableScoresII: Not enough memory." ));

            memset( scores[m], 0, sizeof( TScore ) * GetQuerySize());
        }
    }
}

// -------------------------------------------------------------------------
// destructor: deallocate memory used by this class
// -------------------------------------------------------------------------

AttributableScoresII::~AttributableScoresII()
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

void AttributableScoresII::ScaleToAttainLambda()
{
#ifdef SCALE_BINARY_SEARCH
    bool    USE_BOUND_PRECOMPUTATION = false;//true;

    double  prev_lambda;
    double  diff_target = 1.0;//difference of lambdas
    double  min_diff_target = 10.0;//minimum value of the difference
    double  multiplier = 1.0;

    ComputeProbabilitiesAndLambda();
    if( GetAutoScalingFactor() < 1 )
        throw myruntime_error("AttributableScoresII: ScaleToAttainLambda: Invalid scaling factor.");

    double          right = 1.0;    //right interval bound
    double          left  = 1.0;    //left interval bound
    double          step = 0.05;    //initial step to narrow interval by
    double          max_multiplier = 1.0;
    const double    accuracy = 1.e-4 /( double )GetAutoScalingFactor();
    int i = 0;

    bool    greater = GetRefLambda() < GetLambda();

    if( !USE_BOUND_PRECOMPUTATION ) {

//        if( greater )   right = 10.0;
//        else            left = 0.0;
        double lmb_l = GetRefLambda()/1.1;
        double lmb_u = GetRefLambda()/1.;//0.98;
        if( greater ) {
            if( lmb_u<GetLambda())
                return;
            right = 1.02;
        } else {
            //if( GetLambda()<lmb_l)
            //    return;
            left = 0.9;//.2;
        }
//        if( greater )   return;
//        else if(1.0<fct) right = fct;
//             else        left = fct;

    } else {//USE_BOUND_PRECOMPUTATION
// fprintf(stderr,"init: fct=%f, left=%f, right=%f l=%g ref=%g\n", multiplier,left,right,GetLambda(),GetRefLambda());
        for( i = 0; 0.0 < GetLambda() && i < MAX_SCALE_ITERATIONS; i++, step += step ) {
            if( greater ) {
                    if( GetLambda() <= GetRefLambda())   break;
            } else  if( GetLambda() >= GetRefLambda())   break;

            if( greater ) {
                left   = right;
                right += step;
                multiplier = right;
            } else {
                right = left;
                left -= step;
                multiplier = left;
                if( left < 0.0 ) {
                    left = 0.0;
                    break;
                }
            }

            prev_lambda = GetLambda();
            AdjustScoresInMatrix( multiplier, ImageTable );
            ComputeProbabilitiesAndLambda( false/*wrn*/);
// fprintf(stderr,"bound: fct=%f, left=%f, right=%f l=%g ref=%g\n", multiplier,left,right,GetLambda(),GetRefLambda());

            if( GetLambda() < 0.0 )
                break;

            diff_target = fabs( GetRefLambda() - GetLambda());

            if( diff_target < min_diff_target ) {
                min_diff_target = diff_target;
                max_multiplier = multiplier;
            }
            if( diff_target < accuracy )
                break;
            if( fabs( GetLambda() - prev_lambda ) < accuracy )
                break;
        }//for
        if( 0.0 < GetLambda()) {
            if( greater ) {
                    if( GetLambda() > GetRefLambda())   right += step + step;
            } else  if( GetLambda() < GetRefLambda())   left = 0.0;
        }
    }//if USE_BOUND_PRECOMPUTATION

    if( 0.0 < GetLambda())
        for( ; i < MAX_SCALE_ITERATIONS && accuracy < diff_target; i++ ) {
            multiplier = ( left + right ) * 0.5;

            prev_lambda = GetLambda();
            AdjustScoresInMatrix( multiplier, ImageTable );
            ComputeProbabilitiesAndLambda( false/*wrn*/);
// fprintf(stderr,"fct=%f, left=%f, right=%f l=%g ref=%g\n", multiplier,left,right,GetLambda(),GetRefLambda());

            if( GetLambda() < 0.0 )
                break;

            diff_target = fabs( GetRefLambda() - GetLambda());

            if( diff_target < min_diff_target ) {
                min_diff_target = diff_target;
                max_multiplier = multiplier;
            }
            if( fabs( GetLambda() - prev_lambda ) < accuracy )
                break;

            if( GetRefLambda() < GetLambda())
                left  = multiplier;
            else
                right = multiplier;
        }//for

    if( GetLambda() < 0.0 || min_diff_target < diff_target ) {
        AdjustScoresInMatrix( multiplier = max_multiplier, ImageTable );
        ComputeProbabilitiesAndLambda( false/*wrn*/);
// fprintf(stderr,"adjust: fct=%f, l=%g ref=%g\n", multiplier,GetLambda(),GetRefLambda());
    }

#else //if not defined SCALE_BINARY_SEARCH

//     ComputeStatisticalParameters(); //TEST!!! {{
//     AdjustScoresInMatrix();
//     ComputeStatisticalParameters(); //}}

    //sometimes scaling with floating point numbers and rounding then to nearest
    //integer can produce positive expected score per position; this may happen
    //when expected score is nearly zero. In this case, approximate integer scaling
    //is used
    multiplier = IntegerScaleMatrix();

#endif

    SetPrivateMultiplier( multiplier );
}

// -------------------------------------------------------------------------
// IntegerScaleMatrix: perform iterative scaling until required lambda is
//     attained; integer scaling is used 
// -------------------------------------------------------------------------

double AttributableScoresII::IntegerScaleMatrix()
{
    double  prev_lambda;
    double  diff_lambda;
    double  multiplier = 1.0;
    int     iter = 0;
    const double    least_difference = GetRefLambda() * 0.001;

    if( ! GetKeepInMemory())
        return multiplier;

    ComputeProbabilitiesAndLambda();

    while( 0.0 < GetLambda() && iter++ < MAX_SCALE_ITERATIONS )
    {
        prev_lambda = GetLambda();

        //we must necessary use the 'scores' table;
        //otherwise, no proper value will be obtained
        AdjustScoresInMatrix( multiplier = GetLambda() / GetRefLambda(), ScoreTable );

        ComputeProbabilitiesAndLambda( false/*wrn*/);
        diff_lambda = GetLambda() - prev_lambda;

        if( -least_difference < diff_lambda && diff_lambda < least_difference )
            break;
    }

    return multiplier;
}

// -------------------------------------------------------------------------
// ComputeEntropyGivenLambda: compute relative entropy given previously
//     computed lambda; We compute entropy before adjusting the scores.
//     Mathematically the result would be obtained the same if computing
//     entropy after the score adjustment, since the probabilities of scores
//     retain the same. The only difference is that the reference lambda
//     instead of computed lambda should be used.
// Relative entropy or information per aligned pair of positions is computed
// as follows:
//
//           _                   lambda s(k)
//   lambda \  s(k) * p(s(k)) * e
//          /_
//          k
//
// where p(s(k)) is a probability of the score s(k)
// -------------------------------------------------------------------------

void AttributableScoresII::ComputeEntropyGivenLambda()
{
    double  scaling = GetLambda();

    if( scaling < 0.0 ) {
        warning( "AttributableScoresII: Unable to compute relative entropy." );
        return;
    }

    double  y = exp( -scaling );
    double  m = pow( y, GetMaxScore()); // e to -lambda * max_score
    double  H = 0.0;

    for( int sc = GetMinScore(); sc <= GetMaxScore(); sc ++ )
        H = H * y + sc * GetProbabilityOf( sc );

    if( 0.0 < m )
        H /= m;
    else
        if( 0.0 < H )
            H = exp( scaling * GetMaxScore() + log( H ));

    SetH( scaling * H );
}

// -------------------------------------------------------------------------
// ComputePositionalEntropy: compute entropy at specified position of score
//     system; to compute entropy, use all needed data given by parameters
// -------------------------------------------------------------------------

double AttributableScoresII::ComputePositionalEntropy(
    double lambda,
    int minscore,
    int maxscore,
    int pos,
    TGetProbFunction getProbfunction,
    TGetProbNormFunction getProbNormfunction ) const
{
    if( getProbfunction == NULL )
        throw myruntime_error( mystring( "AttributableScoresII: ComputePositionalEntropy: Memory access error." ));

    if( lambda < 0.0 )
        throw myruntime_error( mystring( "AttributableScoresII: ComputePositionalEntropy: Unable to compute entropy." ));

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
// FindLambdaRootAtPosition: helper method to find root of lambda given
//     parameters of positional scores
// -------------------------------------------------------------------------

double AttributableScoresII::FindLambdaRootAtPosition(
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

double AttributableScoresII::findLambdaRoot( TConservationFunction fdfunction, int gcd,
                            double x1, double x2, double xacc, int maxit, void* params )
{
    double  rts = -1.0;

    x1 = exp( -x1 );
    x2 = exp( -x2 );

    rts = rootByNRandBisection( fdfunction, gcd, x1, x2, xacc, maxit, params, false/*warn*/);

    if( rts < 0.0 )
        return rts;

    return -log( rts ) / gcd;
}

// -------------------------------------------------------------------------
// findConstantRoot: helper method to find root of constant to adjust
//     relative entropy
// -------------------------------------------------------------------------

double AttributableScoresII::findConstantRoot( TConservationFunction fdfunction, int gcd,
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
// However such computation is very slow and is not accurate because of
// false probability estimates.
//     Therefore using scores as integers is encouraged. When this is the
// case, function whose roots to be found is
//
//    _   lambda s(k)
//   \  e            p(s(k)) - 1
//   /_
//    k
//
// where p(s(k)) are probabilities of scores s(k). To find a solution, a
// procedure as in psi-blast is used: the equation is transformed into
// the polynomial form,
//     exp(max s(k) lambda) * poly(exp(-lambda))
// where poly is a polynomial of exp(-lambda) and has exactly two zeros:
// 1 and another one in interval (0,1)
// -------------------------------------------------------------------------

void AttributableScoresII::conservation( double x, double* f, double* df, int step, void* )
{
    double  y = x;      //initially x is supposed to be exp(-lambda)
    double  ldf = 0.0;
    double  lf = 0.0;

    for( int sc = GetMinScore(); sc <= GetMaxScore(); sc += step ) {
        ldf = ldf * y + lf;
        lf  = lf  * y + GetProbabilityOf( sc );
        //zero will always be reached
        if( sc == 0 )
            lf -= 1.0;
    }
    *df = ldf;
    *f  = lf;
}

// -------------------------------------------------------------------------
// lambdaconspos: alternative lambda conservation function which differs
//     from the one above by how data are used
//
void AttributableScoresII::lambdaconspos( double x, double* f, double* df, int step, void* params )
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

void AttributableScoresII::relentropy_conservation( double x, double* f, double* df, int step, void* )
{
    if( x <= 0.0 )
        throw myruntime_error(
            mystring( "AttributableScoresII: Computation of adjusted relative entropy failed." ));

    double  c = x;          //constant c to be found
    double  lf = 0.0;       //functin of relative entropy
    double  ldf = 0.0;      //derivative function of relative entropy with respect to c
    double  scp = 0.0;      //score multiplied by its probability

    double  scaling = GetLambda();
    //y should be e to -lambda; however in the beginning of computation of
    // statistical parameters, scores implicitly include parameter lambda by
    // implicit product of the two
    double  y = exp( -scaling );
    double  m = pow( y, GetMaxScore()); // e to -max_score

    double  logc = log( c );    //logarithm of c
    double  dlogc = logc + logc;//two logarithms of c
    double  dlogcdl = dlogc / scaling;//2 log c / lambda


    for( int sc = GetMinScore(); sc <= GetMaxScore(); sc += step ) {
        scp = ( sc - dlogcdl ) * GetProbabilityOf( sc );
        lf  = lf  * y + scp;
        ldf = ldf * y + scp + GetProbabilityOf( sc );
        //zero will always be reached
        if( sc == 0 )
            //substract value of relative entropy to be achieved
            lf -= GetRefH();
    }

    if( 0.0 < m ) { lf  /= m;
                    ldf /= m;
    } else {
            if( 0.0 < lf  ) lf  = exp( GetMaxScore() + log( lf  ));
            if( 0.0 < ldf ) ldf = exp( GetMaxScore() + log( ldf ));
        }

    *f  = lf * scaling / SQUARE( c );
    *df = -2.0 * ldf * scaling / CUBE( c );
}

