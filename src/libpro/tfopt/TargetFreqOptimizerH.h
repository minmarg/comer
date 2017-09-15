/***************************************************************************
 *   Copyright (C) 2010 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __TargetFreqOptimizerH__
#define __TargetFreqOptimizerH__

#include "mystring.h"
#include "myexcept.h"

#include "ext/pslvector.h"

// #ifndef TFOPTESTPRINT
// #define TFOPTESTPRINT
// #endif

// _________________________________________________________________________
// Class TargetFreqOptimizerH
//
// Optimization of target frequencies of score system to keep relative
// entropy (statistical parameter) of the frequencies constant.
// The method and the algorithm is based on
//  Yu and Altschul (2005) Bioinformatics 21(7), and
//  Altschul et al. (2005) FEBS J. 272(20)
//
class TargetFreqOptimizerH
{
public:
    TargetFreqOptimizerH( Pslvector& scores, int norows, int nocols );
    TargetFreqOptimizerH( Pslvector& scores, const Pslvector& rprobs, const Pslvector& cprobs );
    ~TargetFreqOptimizerH();

    double              GetLambda() const                   { return lambda_; }
    void                SetLambda( double value )           { lambda_ = value; }

    double              GetConstrainedH() const             { return relentropy_; }
    void                SetConstrainedH( double value )     { relentropy_ = value; }

    double              GetNMResidTol() const               { return nmtolerance_; }
    void                SetNMResidTol( double value )       { nmtolerance_ = value; }

    int                 GetMaxNoNMIterations() const        { return maxnmiters_; }
    void                SetMaxNoNMIterations( int value )   { maxnmiters_ = value; }

    int                 Optimize();
    void                Normalize();
    bool                Negatives();

    const Pslvector&    GetOptTargetFreqs() const           { return optfreqs_; }
    int                 GetNoIterations() const             { return iters_; }

protected:
    void                InitData();
    int                 ComputeInitialTargetFreqs();
    int                 ComputeProbsUsingTargetFreqs();

    int                 Iterate();
    int                 CalculateGradient();
    int                 CalculateResiduals();

    int                 MultiplyByA_T( Pslvector& u, const Pslvector& v, bool sub = false );
    int                 MultiplyByA( Pslvector& v, const Pslvector& u, bool sub = false );
    int                 FindStep( const Pslvector& x, const Pslvector& dx, double* alpha );

    int                 ReduceJacobian();
    int                 CholeskyDecompose();
    int                 SolveNewtonSystem();
    int                 SolveDecomposed();


    void                ResetNoIterations()                 { iters_ = 0; }
    void                IncNoIterations()                   { iters_++; }

    double              GetTarFunValue() const              { return tarfunval_; }
    void                SetTarFunValue( double value )      { tarfunval_ = value; }
    void                IncTarFunValue( double value )      { tarfunval_ += value; }

    double              GetRelentValue() const              { return relentval_; }
    void                SetRelentValue( double value )      { relentval_ = value; }
    void                IncRelentValue( double value )      { relentval_ += value; }

    bool                GetProbsGiven() const               { return givenps_; }

private:
    Pslvector&          scores_;//scores
    Pslvector           rowprobs_;//row background probabilities
    Pslvector           colprobs_;//column background probabilities
    Pslvector           initfreqs_;//initial target frequencies derived from scores
    Pslvector           optfreqs_;//optimal frequencies when done
    double              lambda_;//statistical parameter lambda, scale of scores
    double              relentropy_;//negative if not used
    double              nmtolerance_;//residual tolerance (precision) of the Newton's method
    int                 maxnmiters_;//maximum no. iterations allowed
    int                 iters_;//no. iterations passed so far

    Pslvector           duals_;//dual variables (multipliers) of Lagrangian
    Pslvector           obsduals_;//previously computed dual variables
    Pslvector           obsfreqs_;//previously computed target frequencies
    Pslvector           lineresids_;//residuals of linear constraints
    Pslvector           dualresids_;//dual residuals of gradient of Lagrangian
    Pslvector           primsteps_;//steps of primary variables
    Pslvector           dualsteps_;//steps of dual variables
    Pslvector           tarfungrad_;//gradient of target function to optimize
    Pslvector           relentgrad_;//gradient of relative entropy constraint
    double              tarfunval_;//value of target function
    double              relentval_;//value of relative entropy
    Pslvector           dglinverse_;//inverse of diagonal matrix (just diagonal elements)
    Pslvector           blkproduct_;//lower triangle of the symmetric block-product matrix: A D^(-1) A<T>
    Pslvector           tempvector_;//vector of temporary usage

    bool                givenps_;//probabilities given
};


#endif//__TargetFreqOptimizerH__
