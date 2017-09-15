/***************************************************************************
 *   Copyright (C) 2010 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "pcmath.h"
#include "pslcodes.h"
#include "TargetFreqOptimizerH.h"



// =========================================================================
// CLASS TargetFreqOptimizerH
//
// constructor: Initialization
//
TargetFreqOptimizerH::TargetFreqOptimizerH( Pslvector& scores, int norows, int nocols )
:   scores_( scores ),
    rowprobs_( norows ),
    colprobs_( nocols ),

    initfreqs_( scores.GetSize()),
    optfreqs_( scores.GetSize()),
    lambda_( 1.0 ),
    relentropy_( -1.0 ),
    nmtolerance_( 1.0 ),
    maxnmiters_( 1 ),
    iters_( 0 ),

    //number of dual variables is as large as a number of linear constraints
    duals_( colprobs_.GetSize() + rowprobs_.GetSize()),
    obsduals_( colprobs_.GetSize() + rowprobs_.GetSize()),
    obsfreqs_( scores.GetSize()),
    lineresids_( colprobs_.GetSize() + rowprobs_.GetSize()),
    dualresids_( scores.GetSize()),
    primsteps_( scores.GetSize()),
    dualsteps_( colprobs_.GetSize() + rowprobs_.GetSize()),
    tarfungrad_( scores.GetSize()),
    relentgrad_( scores.GetSize()),
    tarfunval_( 0.0 ),
    relentval_( 0.0 ),
    dglinverse_( scores.GetSize()),
    blkproduct_( duals_.GetSize() * ( duals_.GetSize() + 1 ) / 2 ),
    tempvector_( scores.GetSize()),
    givenps_( 0 )
{
    if( rowprobs_.GetSize() < 1 || rowprobs_.GetSize() < 1 || colprobs_.GetSize() < 1 ||
        scores.GetSize() != rowprobs_.GetSize() * colprobs_.GetSize())
        throw myruntime_error( "TargetFreqOptimizerH: Inconsistent scores and background probabilities." );
    throw myruntime_error( "TargetFreqOptimizerH: Derivation of probabilities from scores alone not implemented yet." );
}

// constructor: Initialization given background probabilities
//
TargetFreqOptimizerH::TargetFreqOptimizerH( Pslvector& scores, const Pslvector& rprobs, const Pslvector& cprobs )
:   scores_( scores ),
    rowprobs_( rprobs ),
    colprobs_( cprobs ),
    initfreqs_( scores.GetSize()),
    optfreqs_( scores.GetSize()),
    lambda_( 1.0 ),
    relentropy_( -1.0 ),
    nmtolerance_( 1.0 ),
    maxnmiters_( 1 ),
    iters_( 0 ),

    //number of dual variables is as large as a number of linear constraints
    duals_( cprobs.GetSize() + rprobs.GetSize()),
    obsduals_( cprobs.GetSize() + rprobs.GetSize()),
    obsfreqs_( scores.GetSize()),
    lineresids_( cprobs.GetSize() + rprobs.GetSize()),
    dualresids_( scores.GetSize()),
    primsteps_( scores.GetSize()),
    dualsteps_( cprobs.GetSize() + rprobs.GetSize()),
    tarfungrad_( scores.GetSize()),
    relentgrad_( scores.GetSize()),
    tarfunval_( 0.0 ),
    relentval_( 0.0 ),
    dglinverse_( scores.GetSize()),
    blkproduct_( duals_.GetSize() * ( duals_.GetSize() + 1 ) / 2 ),
    tempvector_( scores.GetSize()),

    givenps_( 1 )
{
    if( rprobs.GetSize() < 1 || rprobs.GetSize() < 1 || cprobs.GetSize() < 1 ||
        scores.GetSize() != rprobs.GetSize() * cprobs.GetSize())
        throw myruntime_error( "TargetFreqOptimizerH: Inconsistent scores and background probabilities." );
//     InitData();
}

// -------------------------------------------------------------------------
// destructor: deallocation of resources
//
TargetFreqOptimizerH::~TargetFreqOptimizerH()
{
}

// -------------------------------------------------------------------------
// InitData: Initialization of private data
//
void TargetFreqOptimizerH::InitData()
{
    duals_.Zero();
}

// -------------------------------------------------------------------------
// Normalize: Normalize optimized target frequencies
//
void TargetFreqOptimizerH::Normalize()
{
    double  sum;

    sum = optfreqs_.Sum();
    if( !sum )
        return;

    optfreqs_.MultiplyBy( 1.0 / sum );
}

// -------------------------------------------------------------------------
// Negatives: Check whether exists negative frequencies
//
bool TargetFreqOptimizerH::Negatives()
{
    int n;

    for( n = 0; n < optfreqs_.GetSize(); n++ )
        if( optfreqs_.GetValueAt( n ) < 0.0 )
            return true;
    return false;
}

// -------------------------------------------------------------------------
// ComputeInitialTargetFreqs: Compute initial target frequencies given
//     log-odds scores and background probabilities
//
int TargetFreqOptimizerH::ComputeInitialTargetFreqs()
{
    int     n, m, ind;
    double  consv;
    double  value;

    if( initfreqs_.GetSize() != scores_.GetSize())
        return PSL_ERR_DIM;


    for( n = 0; n < rowprobs_.GetSize(); n++ )
        if( rowprobs_.GetValueAt( n ) <= 0.0 ) {
#ifdef TFOPTESTPRINT
            fprintf( stderr, "Invalid background probabilities\n" );
#endif
            return PSL_ERR_ILLEGAL;
        }
    for( m = 0; m < colprobs_.GetSize(); m++ )
        if( colprobs_.GetValueAt( m ) <= 0.0 ) {
#ifdef TFOPTESTPRINT
            fprintf( stderr, "Invalid background probabilities\n" );
#endif
            return PSL_ERR_ILLEGAL;
        }



    for( n = 0; n < rowprobs_.GetSize(); n++ ) {
        for( m = 0; m < colprobs_.GetSize(); m++ ) {
            ind = n * colprobs_.GetSize() + m;

            if( GetLambda() == 1.0 )
                value = rowprobs_.GetValueAt( n ) * colprobs_.GetValueAt( m ) *
                    exp( scores_.GetValueAt( ind ));
            else
                value = rowprobs_.GetValueAt( n ) * colprobs_.GetValueAt( m ) *
                    exp( GetLambda() * scores_.GetValueAt( ind ));
            initfreqs_.SetValueAt( ind, value );
        }
    }

    consv = initfreqs_.Sum();
    if( consv < 0.9 || 1.1 < consv ) {
#ifdef TFOPTESTPRINT
        fprintf( stderr, "Initial frequencies unstable: %g\n", consv );
#endif
    }
    if( 0.0 < consv && consv != 1.0 ) {
        initfreqs_.MultiplyBy( 1.0 / consv );
        //now adjust appropriately scores
        for( n = 0; n < rowprobs_.GetSize(); n++ ) {
            for( m = 0; m < colprobs_.GetSize(); m++ ) {
                ind = n * colprobs_.GetSize() + m;

                value = log( initfreqs_.GetValueAt( ind ) / ( rowprobs_.GetValueAt( n ) * colprobs_.GetValueAt( m )));
                if( GetLambda() && GetLambda() != 1.0 )
                    value /= GetLambda();
                scores_.SetValueAt( ind, value );
            }
        }
    }

    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// ComputeProbsUsingTargetFreqs: Compute background probabilities using
//     initial target frequencies
//
int TargetFreqOptimizerH::ComputeProbsUsingTargetFreqs()
{
    const double    accuracy = 1.e-8;
    int     n, m, ind;
    double  consv;
    double  value;

    if( rowprobs_.GetSize() < 1 || colprobs_.GetSize() < 1 ||
        initfreqs_.GetSize() != rowprobs_.GetSize() * colprobs_.GetSize())
        return PSL_ERR_DIM;

    rowprobs_.Zero();
    colprobs_.Zero();

    for( n = 0; n < rowprobs_.GetSize(); n++ ) {
        for( m = 0; m < colprobs_.GetSize(); m++ ) {
            ind = n * colprobs_.GetSize() + m;

            value = initfreqs_.GetValueAt( ind );
            rowprobs_.AddValueAt( n, value );
            colprobs_.AddValueAt( m, value );
        }
    }

    consv = rowprobs_.Sum();
    if( consv < 1.0 - accuracy || consv > 1.0 + accuracy ) {
#ifdef TFOPTESTPRINT
        fprintf( stderr, "TargetFreqOptimizerH::ComputeProbsUsingTargetFreqs: "
                        "Background probabilities not conserved: %f", consv );
#endif
        return PSL_ERR_INVALID;
    }
    consv = colprobs_.Sum();
    if( consv < 1.0 - accuracy || consv > 1.0 + accuracy ) {
#ifdef TFOPTESTPRINT
        fprintf( stderr, "TargetFreqOptimizerH::ComputeProbsUsingTargetFreqs: "
                        "Background probabilities not conserved: %f", consv );
#endif
        return PSL_ERR_INVALID;
    }

    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// Optimize: Begin multidimensional root finding by the Newton's method
//  Target function is relative entropy of target frequencies:
//      SUM x[i][j] log( x[i][j]/q[i][j] )
//  with linear constraints:
//      SUM x[i][:] = P[i] and SUM x[:][j] = P'[j]
//  and optional constraint of relative entropy:
//      SUM x[i][j] log( x[i][j] /( P[i]P'{j] )) = h
// 
//
int TargetFreqOptimizerH::Optimize()
{
    const double    accuracy = 1.e-10;
    int     maxiters = GetMaxNoNMIterations();
    int     status;
    double  consv;

    if( optfreqs_.GetSize() != initfreqs_.GetSize())
        return PSL_ERR_DIM;

    InitData();

#ifdef TFOPTESTPRINT
        fprintf( stderr, "*** TargetFreqOptimizerH::Optimize ***\n" );
        fprintf( stderr, "Lambda=%g, H=%g\n", GetLambda(), GetConstrainedH());
#endif

    if(( status = ComputeInitialTargetFreqs()) != PSL_SUCCESS )
        return status;

//     if(( status = ComputeProbsUsingTargetFreqs()) != PSL_SUCCESS )
//         return status;

    optfreqs_.Copy( initfreqs_ );
    consv = optfreqs_.Sum();

    if( consv < 1.0 - accuracy || consv > 1.0 + accuracy ) {
#ifdef TFOPTESTPRINT
        fprintf( stderr, "Frequencies not conserved: %f\n", consv );
#endif
        return PSL_ERR_INVALID;
    }

    status = PSL_MAXITERATS;
    ResetNoIterations();

    for( ; GetNoIterations() < maxiters && status == PSL_MAXITERATS; IncNoIterations()) {
        status = Iterate();
#ifdef TFOPTESTPRINT
        fprintf( stderr, "* it=%d status=%d\n", GetNoIterations(), status );
#endif
    }

    return status;
}

// -------------------------------------------------------------------------
// Iterate: Perform one iteration of Newton's method
//
int TargetFreqOptimizerH::Iterate()
{
    int     status;
    double  alpha;
    double  norm, dn, ln;
    const double    camul = 0.95;
    static double   prenorm = 0.0;
    double  amul = camul;

    if( primsteps_.GetSize() != dualresids_.GetSize() ||
        dualsteps_.GetSize() != lineresids_.GetSize() ||
        obsfreqs_.GetSize() != optfreqs_.GetSize() ||
        obsduals_.GetSize() != duals_.GetSize())
        return PSL_ERR_DIM;

    if(( status = CalculateGradient()) != PSL_SUCCESS )
        return status;

    if(( status = CalculateResiduals()) != PSL_SUCCESS )
        return status;

    dn = dualresids_.Norm2();
    ln = lineresids_.Norm2();
    norm = sqrt( SQUARE( dn ) + SQUARE( ln ));

#ifdef TFOPTESTPRINT
        fprintf( stderr, "dual/res.norm=%g, linear/res.norm=%g, norm=%g\n", dn, ln, norm );
        fprintf( stderr, "rel.ent=%g, KL-distance=%g\n", GetRelentValue(), GetTarFunValue());
#endif

    if( norm <= GetNMResidTol())
        //converged
        return PSL_SUCCESS;

    if( 0 && prenorm && prenorm <= norm ) {
        amul *= 0.5;
        alpha *= amul;

        if( prenorm <= norm && alpha < 1.e-6 )
            return PSL_ERR_NOPROG;

#ifdef TFOPTESTPRINT
        fprintf( stderr, "ADJUSTED step/mult=%g\n", alpha );
#endif

        if( obsfreqs_.Superposition( alpha, primsteps_ ))
            return PSL_ERR_DIM;
        if( obsduals_.Superposition( alpha, dualsteps_ ))
            return PSL_ERR_DIM;

        optfreqs_.Copy( obsfreqs_ );
        duals_.Copy( obsduals_ );

        return PSL_MAXITERATS;
    }


    if(( status = ReduceJacobian()) != PSL_SUCCESS )
        return status;

    if(( status = CholeskyDecompose()) != PSL_SUCCESS )
        return status;

    if(( status = SolveNewtonSystem()) != PSL_SUCCESS )
        return status;

    //now, solutions in steps of dual and primal variables
    //(target frequencies) are in lineresids_ and dualresids_, respectively.

    primsteps_.Copy( dualresids_ );
    dualsteps_.Copy( lineresids_ );

    amul = camul;
    alpha = 1.0 / camul;
    if(( status = FindStep( optfreqs_, primsteps_, &alpha )) != PSL_SUCCESS )
        return status;
    alpha *= amul;//to avoid situation when optfreqs_[n]==0

    prenorm = norm;

#ifdef TFOPTESTPRINT
    fprintf( stderr, "step/mult=%g\n", alpha );
#endif

    obsfreqs_.Copy( optfreqs_ );
    obsduals_.Copy( duals_ );

    //update vectors of primal and dual variables:
    //x += a dx
    if( optfreqs_.Superposition( alpha, primsteps_ ))
        return PSL_ERR_DIM;
    if( duals_.Superposition( alpha, dualsteps_ ))
        return PSL_ERR_DIM;

    return PSL_MAXITERATS;
}

// -------------------------------------------------------------------------
// CalculateGradient: Calculate gradient of target function and
//     constraint of relative entropy
//  in case of relative entropy, its gradient is calculated as follows:
//      SUM x[i][j]  log( x[i][j] /( P[i]P'[j] )) =
//      SUM x[i][j]( log( x[i][j] / q[i][j] ) + log( q[i][j] /( P[i]P'[j] )))
//
int TargetFreqOptimizerH::CalculateGradient()
{
    double  x, q, val;
    int     n;

    if( tarfungrad_.GetSize() != scores_.GetSize() ||
        relentgrad_.GetSize() != scores_.GetSize())
        return PSL_ERR_DIM;

    SetTarFunValue( 0.0 );
    SetRelentValue( 0.0 );

    for( n = 0; n < scores_.GetSize(); n ++ )
    {
        x = optfreqs_.GetValueAt( n );
        q = initfreqs_.GetValueAt( n );

        if( x <= 0.0 || q <= 0.0 ) {
#ifdef TFOPTESTPRINT
            fprintf( stderr, "PSL_ERR_ILLEGAL: x=%g, q=%g\n", x, q );
#endif
            return PSL_ERR_ILLEGAL;
        }

        val = log( x / q );

        IncTarFunValue( x * val );
        tarfungrad_.SetValueAt( n, val + 1.0 );

        if( GetConstrainedH() < 0.0 )
            //relative entropy constraint is not in use
            continue;

        if( GetLambda() == 1.0 )
            val += scores_.GetValueAt( n );
        else
            val += scores_.GetValueAt( n ) * GetLambda();

        IncRelentValue( x * val );
        relentgrad_.SetValueAt( n, val + 1.0 );
    }

    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// CalculateResiduals: Calculate residuals of linear constraints and dual
//  residuals of gradient of Lagrangian
//
int TargetFreqOptimizerH::CalculateResiduals()
{
    double  kappa;//multiplier (variable) for constraint of relative entropy
    double  val;
    int     n, m;
    int     status;

    if( duals_.GetSize() < 1 || lineresids_.GetSize() < 1 || dualresids_.GetSize() < 1 )
        return PSL_ERR_DIM;

    if( duals_.GetSize()      != rowprobs_.GetSize() + colprobs_.GetSize() ||
        lineresids_.GetSize() != rowprobs_.GetSize() + colprobs_.GetSize() ||
        dualresids_.GetSize() != scores_.GetSize() ||
        dualresids_.GetSize() != tarfungrad_.GetSize() ||
        dualresids_.GetSize() != relentgrad_.GetSize())
        return PSL_ERR_DIM;

    //1.Calculate DUAL RESIDUALS
    dualresids_.Zero();

    //if constraint of relative entropy is in use
    if( 0.0 < GetConstrainedH())
    {
        kappa = duals_.GetValueAt( duals_.GetSize() - 1 );
        for( n =0; n < dualresids_.GetSize(); n++ ) {
            val = kappa * relentgrad_.GetValueAt( n );
            dualresids_.SetValueAt( n, val );
        }
    }

    for( n =0; n < dualresids_.GetSize(); n++ ) {
        val = -tarfungrad_.GetValueAt( n );
        dualresids_.AddValueAt( n, val );
    }

    if(( status = MultiplyByA_T( dualresids_, duals_ )) != PSL_SUCCESS )
        return status;


    //2.Calculate RESIDUALS of LINEAR constraints
    lineresids_.Zero();

    //first, column linear constraints follow
    for( m = 0; m < colprobs_.GetSize(); m++ )
        lineresids_.SetValueAt( m, colprobs_.GetValueAt( m ));

    //start with 1, the first row linear constraint is omitted as redundant
    for( n = 1; n < rowprobs_.GetSize(); n++ )
        lineresids_.SetValueAt( n - 1 + colprobs_.GetSize(), rowprobs_.GetValueAt( n ));

    if(( status = MultiplyByA( lineresids_, optfreqs_, true )) != PSL_SUCCESS )
        return status;

    //if constraint of relative entropy is in use
    if( 0.0 < GetConstrainedH())
        lineresids_.SetValueAt( lineresids_.GetSize() - 1, GetConstrainedH() - GetRelentValue());

    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// MultiplyByA_T: Multiply vectors by A<T> (A transpose): u +/-= A<T> v
//  A is a block matrix from the Jacobian matrix (J) of the Newton system:
//      J dx = -F(x),
//      where J has the form of:( D -A<T> )
//                              ( A    0  )
//      D is a diagonal matrix with partial derivatives with respect to
//      target frequencies, A is a matrix of derivatives with respect to
//      primal (target frequencies from linear constraints and relative
//      entropy constraint) (dual variables/Lagrange multipliers in case of
//      A transpose).
//      All elements from A are 0s and 1s but the last line (column),
//      which correspond to constraint of relative entropy. The last line
//      (column) is not included in multiplication, this is accomplished
//      elsewhere.
//  Dimensions of D, A are (NM x NM), ((N+M) x NM), respectively;
//  N and M are dimensions of the score system.
//  In multiplication, A explicitly is not used, vectors are multiplied
//  using the pattern observed in A instead.
//
int TargetFreqOptimizerH::MultiplyByA_T( Pslvector& u, const Pslvector& v, bool sub )
{
    int n, m, k, ind;

    if( u.GetSize() != scores_.GetSize() ||
        v.GetSize() != rowprobs_.GetSize() + colprobs_.GetSize())
        return PSL_ERR_DIM;

    for( n = 0; n < rowprobs_.GetSize(); n++ )
        for( m = 0; m < colprobs_.GetSize(); m++ ) {
            ind = n * colprobs_.GetSize() + m;

            //first, column dual variables come
            u.AddValueAt( ind, sub? -v.GetValueAt( m ): v.GetValueAt( m ));
            if( n ) {
                //one (row) linear constraint along with corresponding
                // dual variable (multiplier) is omitted as redundant
                k = n - 1 + colprobs_.GetSize();
                u.AddValueAt( ind, sub? -v.GetValueAt( k ): v.GetValueAt( k ));
            }
        }

    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// MultiplyByA: Multiply vectors by matrix A: v +/-= A u
//  A is a block matrix from the Jacobian matrix (J) of the Newton system:
//      J dx = -F(x),
//      where J has the form of:( D -A<T> )
//                              ( A    0  )
//      All elements from A are 0s and 1s but the last line (column),
//      which correspond to constraint of relative entropy. The last line
//      (column) is not included in multiplication, this is accomplished
//      elsewhere.
//  Dimensions of D, A are (NM x NM), ((N+M) x NM), respectively;
//  N and M are dimensions of the score system.
//  In multiplication, A explicitly is not used, vectors are multiplied
//  using the pattern observed in A instead.
//
int TargetFreqOptimizerH::MultiplyByA( Pslvector& v, const Pslvector& u, bool sub )
{
    int n, m, k, ind;

    if( v.GetSize() != rowprobs_.GetSize() + colprobs_.GetSize() ||
        u.GetSize() != scores_.GetSize())
        return PSL_ERR_DIM;

    for( n = 0; n < rowprobs_.GetSize(); n++ )
        for( m = 0; m < colprobs_.GetSize(); m++ ) {
            ind = n * colprobs_.GetSize() + m;

            //first, column linear constraints come
            v.AddValueAt( m, sub? -u.GetValueAt( ind ): u.GetValueAt( ind ));
            if( n ) {
                //one (row) linear constraint along with corresponding
                // dual variable (multiplier) is omitted as redundant
                k = n - 1 + colprobs_.GetSize();
                v.AddValueAt( k, sub? -u.GetValueAt( ind ): u.GetValueAt( ind ));
            }
        }

    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// ReduceJacobian: Reduce Jacobian of the Newton system:
//      J dx = -F(x), ( J[i][j] = dF[i]/dx[j] )
//      where J has the form of:( D -A<T> )
//                              ( A    0  )
//      D is a diagonal matrix with partial derivatives with respect to
//      target frequencies, A is a matrix of derivatives with respect to
//      primal (target frequencies from linear constraints and relative
//      entropy constraint) (dual variables/Lagrange multipliers in case of
//      A transpose).
//      All elements from A are 0s and 1s but the last line (column),
//      which correspond to constraint of relative entropy.
//  Dimensions of D, A are (NM x NM), ((N+M) x NM), respectively;
//  N and M are dimensions of the score system.
//      J is block-reduced to:( D          -A<T> )
//                            ( 0  A D^(-1) A<T> )
//  A D^(-1) A<T> is symmetric and positive definite, thus it is
//  decomposed by the Cholesky decomposition so that it equals to L L<T>.
//
int TargetFreqOptimizerH::ReduceJacobian()
{
    double  kappa;//multiplier (variable) for constraint of relative entropy
    double  val, grad, dgval;
    int     n, m, k, ind;
    int     bli, blr, lrow, last;

    if( duals_.GetSize() < 1 )
        return PSL_ERR_DIM;

    if( dglinverse_.GetSize() != optfreqs_.GetSize() ||
        dglinverse_.GetSize() != relentgrad_.GetSize() ||
        duals_.GetSize()      != rowprobs_.GetSize() + colprobs_.GetSize() ||
        blkproduct_.GetSize() != duals_.GetSize() * ( duals_.GetSize() + 1 ) / 2 )
        return PSL_ERR_DIM;

    //1.Compute inverse of diagonal matrix D^(-1)
    //
    //if constraint of relative entropy is in use
    if( 0.0 < GetConstrainedH())
    {
        kappa = duals_.GetValueAt( duals_.GetSize() - 1 );
        if( kappa == 1.0 ) {
#ifdef TFOPTESTPRINT
            fprintf( stderr, "PSL_ERR_ILLEGAL: ReduceJacobian: kappa=%g\n", kappa );
#endif
            return PSL_ERR_ILLEGAL;
        }

        for( n = 0; n < dglinverse_.GetSize(); n++ ) {
            val = optfreqs_.GetValueAt( n ) / ( 1.0 - kappa );
            dglinverse_.SetValueAt( n, val );
        }
    }
    else
        dglinverse_.Copy( optfreqs_ );

    //2.Compute symmetric matrix A D^(-1) A<T>;
    // fill in the lower triangle of the matrix instead;
    // e.g. for problem of alphabet size 3x3, this matrix has the pattern of:
    //
    //  d[0]+d[3]+d[6]
    //       0      d[1]+d[4]+d[7]
    //       0           0      d[2]+d[5]+d[8]
    //       d[3]        d[4]        d[5]   d[3]+d[4]+d[5]
    //       d[6]        d[7]        d[8]        0      d[6]+d[7]+d[8]
    // h0d0+h3d3+h6d6   ...         ...         ...         ...     SUM h[i]d[i]h[i]
    //
    blkproduct_.Zero();

    lrow = rowprobs_.GetSize() + colprobs_.GetSize() - 1;
    lrow = ( lrow * ( lrow + 1 )) / 2;//blkproduct_ last row index
    last = blkproduct_.GetSize() - 1;//blkproduct_ last element index

    for( n = 0; n < rowprobs_.GetSize(); n++ )
        for( m = 0; m < colprobs_.GetSize(); m++ )
        {
            ind = n * colprobs_.GetSize() + m;
            val = dglinverse_.GetValueAt( ind );
            grad = relentgrad_.GetValueAt( ind );

            bli = (( m + 1 ) * ( m + 2 )) / 2 - 1;//0-based index of diagonal element [m][m]
            blkproduct_.AddValueAt( bli, val );

            //if constraint of relative entropy is in use
            if( 0.0 < GetConstrainedH()) {
                dgval = val * grad;
                blkproduct_.AddValueAt( last, grad * dgval );//last element
                blkproduct_.AddValueAt( lrow + m, dgval );//element at [lrow][m]
            }

            if( n ) {
                //one (row) linear constraint along with corresponding
                // dual variable (multiplier) is omitted as redundant
                k = n - 1 + colprobs_.GetSize();//row index
                blr = ( k * ( k + 1 )) / 2;//blkproduct_ index of the first element in row k
                bli = blr + m;//blkproduct_ index of regular element at [k][m]
                blkproduct_.AddValueAt( bli, val );

                //bli = (( k + 1 ) * ( k + 2 )) / 2 - 1;//diagonal element [k][k]
                bli = blr + k;//diagonal element [k][k]
                blkproduct_.AddValueAt( bli, val );

                if( 0.0 < GetConstrainedH()) {
                    blkproduct_.AddValueAt( lrow + k, dgval );//element at [lrow][k]
                }
            }
        }

    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// CholeskyDecompose: Block product matrix obtained in ReduceJacobian,
//  A D^(-1) A<T> is symmetric and positive definite and amenable to
//  the Cholesky decomposition: A D^(-1) A<T> = L L<T>.
//  On exit, blkproduct_ will contain calculated L
//
int TargetFreqOptimizerH::CholeskyDecompose()
{
    double  tmpval, tmp1, tmp2;
    int     n, m, k, norows;
    int     bli, blr, blc;

    if( duals_.GetSize() < 1 )
        return PSL_ERR_DIM;

    if( duals_.GetSize()      != rowprobs_.GetSize() + colprobs_.GetSize() ||
        blkproduct_.GetSize() != duals_.GetSize() * ( duals_.GetSize() + 1 ) / 2 )
        return PSL_ERR_DIM;

    norows = duals_.GetSize() - 1;
    //if constraint of relative entropy is in use
    if( 0.0 < GetConstrainedH())
        norows = duals_.GetSize();

    bli = blr = 0;
    for( n = 0; n < norows; n++ )
    {
        for( m = 0; m < n; m++ )
        {
            blc = ( m * ( m + 1 )) / 2;//M[m][0]
            bli = blr + m;
            tmpval = blkproduct_.GetValueAt( bli );//M[n][m]
            for( k = 0; k < m; k++ ) {
                tmp1 = blkproduct_.GetValueAt( blr + k );//M[n][k]
                tmp2 = blkproduct_.GetValueAt( blc + k );//M[m][k]
                tmpval -= tmp1 * tmp2;
            }
            tmp1 = blkproduct_.GetValueAt( blc + m );//M[m][m]
            if( tmp1 == 0.0 ) {
#ifdef TFOPTESTPRINT
                fprintf( stderr, "PSL_ERR_ILLEGAL: CholeskyDecompose: M[%d][%d]=%g\n", m, m, tmp1 );
#endif
                return PSL_ERR_ILLEGAL;
            }
            blkproduct_.SetValueAt( bli, tmpval / tmp1 );//tmpval/M[m][m]
        }

        tmpval = blkproduct_.GetValueAt( blr + n );//M[n][n]

        for( k = 0; k < n; k++ ) {
            bli = blr + k;
            tmp1 = blkproduct_.GetValueAt( bli );//M[n][k]
            tmpval -= tmp1 * tmp1;
        }
        if( tmpval < 0.0 ) {
#ifdef TFOPTESTPRINT
            fprintf( stderr, "PSL_ERR_ILLEGAL: CholeskyDecompose: M[%d][%d]-SUM M[%d][k]^2=%g\n",
                n, n, n, tmpval );
#endif
            return PSL_ERR_ILLEGAL;
        }
        blkproduct_.SetValueAt( blr + n, sqrt( tmpval ));//M[n][n]

        blr += n + 1;//M[n+1][0]
    }

    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// SolveNewtonSystem: Solve Newton system: J dx = F(x), in two stages:
//  first, solve for step (dx) of dual variables (Lagrange multipliers)
//  using relation A D^(-1) A<T> = L L<T>, second, solve for step of target
//  frequencies.
//  Solutions are back written to lineresids_ and dualresids_.
//
int TargetFreqOptimizerH::SolveNewtonSystem()
{
    double  deltakappa;//step in multiplier (variable) for constraint of relative entropy
    double  val;
    int     n, m;
    int     status;

    if( duals_.GetSize() < 1 || lineresids_.GetSize() < 1 || dualresids_.GetSize() < 1 )
        return PSL_ERR_DIM;

    if( duals_.GetSize()      != rowprobs_.GetSize() + colprobs_.GetSize() ||
        lineresids_.GetSize() != rowprobs_.GetSize() + colprobs_.GetSize() ||
        dualresids_.GetSize() != scores_.GetSize() ||
        dualresids_.GetSize() != relentgrad_.GetSize() ||
        dglinverse_.GetSize() != dualresids_.GetSize() ||
        dglinverse_.GetSize() != tempvector_.GetSize())
        return PSL_ERR_DIM;

    tempvector_.Zero();

    //1.Solve for STEP in DUAL VARIABLES
    //
    //apply block reduction to the right-hand side, i.e.
    //multiplication of dual residuals by -A D^(-1) and addition to
    //appropriate block of F(x) (residuals of linear constraints)
    for( n = 0; n < optfreqs_.GetSize(); n++ ) {
        //calculate -D^(-1) F(x) (dual residuals)
        val  = dglinverse_.GetValueAt( n );
        val *= dualresids_.GetValueAt( n );
        tempvector_.SetValueAt( n, -val );
    }

    //calculate F(x) (residuals of linear constraints) - A [ D^(-1) F(x) ]
    if(( status = MultiplyByA( lineresids_, tempvector_ )) != PSL_SUCCESS )
        return status;

    //if constraint of relative entropy is in use
    if( 0.0 < GetConstrainedH()) {
        for( n = 0; n < relentgrad_.GetSize(); n++ ) {
            //apply reduction operation to the last residual of linear constraint
            val  = relentgrad_.GetValueAt( n );
            val *= tempvector_.GetValueAt( n );
            lineresids_.AddValueAt( lineresids_.GetSize() - 1, val );
        }
    }

    //Right-hand side is filled, solve for step in dual variables
    //using decomposed A D^(-1) A<T> = L L<T> on the left-hand side
    if(( status = SolveDecomposed()) != PSL_SUCCESS )
        return status;

    //2.Solve for STEP in TARGET FREQUENCIES (Newton's method)
    //
    //use newly calculated step in dual variables dz (lineresids_):
    //dx = D^(-1)( -F(x) + A<T> dz )
    //
    //if constraint of relative entropy is in use
    if( 0.0 < GetConstrainedH())
    {
        //include step in multiplier of relative entropy
        deltakappa = lineresids_.GetValueAt( lineresids_.GetSize() - 1 );
        for( n =0; n < dualresids_.GetSize(); n++ ) {
            val = deltakappa * relentgrad_.GetValueAt( n );
            dualresids_.AddValueAt( n, val );
        }
    }

    //calculate: -F(x) + A<T> dz
    if(( status = MultiplyByA_T( dualresids_, lineresids_ )) != PSL_SUCCESS )
        return status;

    //calculate final step: D^(-1)( -F(x) + A<T> dz )
    for( n =0; n < dualresids_.GetSize(); n++ ) {
        val = dglinverse_.GetValueAt( n );
        dualresids_.MulValueAt( n, val );
    }

    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// SolveDecomposed: Solve matrix system with decomposed matrix on the
//  left-hand side (A D^(-1) A<T> = L L<T>) and vector values on the
//  right-hand side (F(x), reduced residuals of linear constraints).
//  L L<T> x = b is solved in two steps: 1) L y = b, 2) L<T> x = y.
//  Set solution of step vector in residuals of linear constraints.
//
int TargetFreqOptimizerH::SolveDecomposed()
{
    double  val, tmpval;
    int     n, m, norows;
    int     bli, blr;

    if( duals_.GetSize() < 1 || lineresids_.GetSize() < 1 )
        return PSL_ERR_DIM;

    if( duals_.GetSize()      != rowprobs_.GetSize() + colprobs_.GetSize() ||
        lineresids_.GetSize() != rowprobs_.GetSize() + colprobs_.GetSize() ||
        blkproduct_.GetSize() != duals_.GetSize() * ( duals_.GetSize() + 1 ) / 2 )
        return PSL_ERR_DIM;

    norows = duals_.GetSize() - 1;
    //if constraint of relative entropy is in use
    if( 0.0 < GetConstrainedH())
        norows = duals_.GetSize();

    //1.Solve L y = b for y, L is lower triangle
    //
    bli = blr = 0;
    for( n = 0; n < norows; n++ )
    {
        tmpval = lineresids_.GetValueAt( n );
        for( m = 0; m < n; m++ )
        {
            bli = blr + m;
            val = blkproduct_.GetValueAt( bli );//L[n][m]
            tmpval -= val * lineresids_.GetValueAt( m );//lineresids_[m] already solved
        }
        val = blkproduct_.GetValueAt( blr + n );//L[n][n]
        if( val == 0.0 ) {
#ifdef TFOPTESTPRINT
            fprintf( stderr, "PSL_ERR_ILLEGAL: SolveDecomposed: L[%d][%d]=%g\n", n, n, val );
#endif
            return PSL_ERR_ILLEGAL;
        }
        //set solution back to residuals of linear constraints
        lineresids_.SetValueAt( n, tmpval / val );

        blr += n + 1;//L[n+1][0]
    }

    //2.Solve L<T> x = y for x given y in lineresids_
    //
    blr = blkproduct_.GetSize() - duals_.GetSize() - 1;
    if( 0.0 < GetConstrainedH())
        blr = blkproduct_.GetSize() - 1;

    for( m = norows - 1; 0 <= m; m-- )
    {
        tmpval = lineresids_.GetValueAt( m );
        val = blkproduct_.GetValueAt( blr );//L[m][m]
        if( val == 0.0 ) {
#ifdef TFOPTESTPRINT
            fprintf( stderr, "PSL_ERR_ILLEGAL: SolveDecomposed: L[%d][%d]=%g\n", m, m, val );
#endif
            return PSL_ERR_ILLEGAL;
        }
        lineresids_.SetValueAt( m, tmpval / val );        

        tmpval = lineresids_.GetValueAt( m );
        for( n = 0; n < m; n++ )
        {
            bli = blr - ( m - n );
            val = blkproduct_.GetValueAt( bli );//L[m][n];
            lineresids_.AddValueAt( n, -val * tmpval );
        }

        blr -= ( m + 1 );//L[m-1][m-1]
    }

    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// FindStep: find valid, maximum available Newton's step <1 such that
//  xi + alpha * dxi >= 0 for all i.
//
int TargetFreqOptimizerH::FindStep( const Pslvector& x, const Pslvector& dx, double* alpha )
{
    int     n;
    double  aatn;//alpha at n

    if( !alpha ) {
#ifdef TFOPTESTPRINT
        fprintf( stderr, "PSL_ERR_ILLEGAL: FindStep: alpha=%g\n", alpha );
#endif
        return PSL_ERR_ILLEGAL;
    }

    if( !alpha || x.GetSize() != dx.GetSize())
        return PSL_ERR_DIM;

    for( n = 0; n < x.GetSize(); n++ ) {
        if( !dx.GetValueAt( n )) {
#ifdef TFOPTESTPRINT
            fprintf( stderr, "PSL_ERR_ILLEGAL: FindStep: dx[%d]=%g\n", n, dx.GetValueAt( n ));
#endif
            return PSL_ERR_ILLEGAL;
        }
        aatn = -x.GetValueAt( n ) / dx.GetValueAt( n );

        if( 0 <= aatn && aatn < *alpha )
            *alpha = aatn;
    }

    return PSL_SUCCESS;
}

