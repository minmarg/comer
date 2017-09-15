/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __AttributableScoresFPI__
#define __AttributableScoresFPI__

#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "data.h"
#include "mystring.h"
#include "myexcept.h"

#include "AttributableScores.h"

class AbstractScoreMatrix;

////////////////////////////////////////////////////////////////////////////
// CLASS AttributableScoresFPI
// Floating-point implementation of AttributableScores
//
class AttributableScoresFPI: public AttributableScores
{
public:
    AttributableScoresFPI(
            AbstractScoreMatrix* prnt,
            TProbabilityFunction probfun,
            TScoreReductFunction redtfun,
            TrogressFunction progressfun,
            double reflambda,
            double refH,
            double refK,
            const double** im,
            const TMask** msk,
            bool keep,
            int as_factor = 1 );

    virtual ~AttributableScoresFPI();

    virtual double      GetScore( int m, int n ) const;             //returns score at specified profile positions
    virtual void        SetScore( double, int m, int n );           //set score at specified positions

    virtual void        Init( int sz_query, int sz_sbjct );

protected:

    explicit AttributableScoresFPI();

    virtual double      FindLambdaRootAtPosition(                   //find lambda at specified position given probs.
        TConservationFunction,
        int gcd,
        double x1,
        double x2,
        double xacc,
        int maxit,
        int minscore,
        int maxscore,
        int pos,
        TGetProbFunction );

    virtual double      findLambdaRoot( TConservationFunction, int, double, double, double, int, void* = NULL );
    virtual double      findConstantRoot( TConservationFunction, int, double, double, double, int, void* = NULL );

    //functions describing conservation equations
    virtual void        lambdaconspos( double x, double* f, double* df, int step, void* params );
    virtual void        conservation( double, double*, double*, int step = 1, void* = NULL );
    virtual void        relentropy_conservation( double x, double* f, double* df, int step = 1, void* = NULL );

    virtual void        ScaleToAttainLambda();                      //perform iterative scaling until required lambda is attained

    virtual void        ComputeScoreProbabilities();                //compute probabilities to observe scores at each position
    virtual void        ComputeEntropyGivenLambda();                //compute relative entropy given previously computed lambda
    virtual double      ComputePositionalEntropy(                   //compute entropy at specified position of score system
                            double lambda, int minscore, int maxscore, int pos,
                            TGetProbFunction, TGetProbNormFunction ) const;
    virtual void        ComputeKarlinsK();                          //compute parameter K

    virtual void        MultiplyScoreBy( double, int m, int n );    //multiply score at the specified position by a factor

protected:
    virtual void** const    GetScores() const               { return ( void** )scores; }

private:
    double**                scores;             //matrix of floating-point scores
};

// INLINES ...

// -------------------------------------------------------------------------
// GetScore: used to access score at the profile positions
// -------------------------------------------------------------------------

inline
double AttributableScoresFPI::GetScore( int m, int n ) const
{
#ifdef __DEBUG__
    if( !scores )
        throw myruntime_error(
            mystring( "AttributableScoresFPI: Memory access error." ));

    if( GetSubjectSize() <= m || m < 0 )
        throw myruntime_error(
            mystring( "AttributableScoresFPI: Memory access error." ));

    if( GetQuerySize() <= n || n < 0 )
        throw myruntime_error(
            mystring( "AttributableScoresFPI: Memory access error." ));
#endif

    return scores[m][n];
}

// -------------------------------------------------------------------------
// SetScore: change the score at the specified location of the matrix
// -------------------------------------------------------------------------

inline
void AttributableScoresFPI::SetScore( double value, int m, int n )
{
#ifdef __DEBUG__
    if( !scores || GetSubjectSize() <= m || m < 0 || GetQuerySize() <= n || n < 0 )
        throw myruntime_error( mystring( "AttributableScoresFPI: Memory access error." ));
#endif

    if( SCORE_MIN < value )
        value *= GetAutoScalingFactor();
    scores[m][n] = value;
}

// -------------------------------------------------------------------------
// MultiplyScoreBy: multiply score at m and n by the value given
// -------------------------------------------------------------------------

inline
void AttributableScoresFPI::MultiplyScoreBy( double value, int m, int n )
{
#ifdef __DEBUG__
    if( !scores || GetSubjectSize() <= m || m < 0 || GetQuerySize() <= n || n < 0 )
        throw myruntime_error( mystring( "AttributableScoresFPI: Memory access error." ));
#endif

    scores[m][n] *= value;
}


#endif//__AttributableScoresFPI__
