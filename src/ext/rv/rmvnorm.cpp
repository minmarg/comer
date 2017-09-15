/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rvnorm.h"
#include "rmvnorm.h"

// =========================================================================
// RMVNorm
//
RMVNorm::RMVNorm( Rng& rng )
:   RMV( rng ),
    mu_( NULL ),
    sigma_( NULL )
{
}

RMVNorm::~RMVNorm()
{
}

// -------------------------------------------------------------------------
// Gen: generate normally distributed random multivariate by Cholesky
//  factorization 
//  (ref. S.Asmussen. Stochastic Simulation: Algorithms and Analysis, p.49);
//  if not sigma and mu are given, generate standard normal multivariate
//
int RMVNorm::Gen( Pslvector* rmv )
{
    int i, j, err;
    double  rv1, rv2;

    if( rmv == NULL )
        return PSL_ERR_ADDRESS;

    int nn = rmv->GetSize();//dimensions

    if( GetSigma() == NULL ) {
        //assume cov. matrix is identity matrix
        if(( err = GenStd( rmv )) != PSL_OK )
            return err;
        //just return if no mean vector is given
        if( GetMu() == NULL )
            return PSL_OK;
        if(( err = rmv->Superposition( 1.0, *GetMu())) != PSL_OK )
            return err;
        return PSL_OK;
    }

    if( nn != GetSigma()->GetNoCols())
        return PSL_ERR_DIM;

    SPDmatrix   cfac( *GetSigma());

    if(( err = cfac.CholeskyDecompose()) != PSL_OK )
        return err;
    //now cfac in lower triangular contains Cholesky factor

    //vector to contain i.i.d. std. normals
    Pslvector   stdn( nn );

    if(( err = GenRVNorms2( nn, &stdn )) != PSL_OK )
        return err;

    for( i = 0; i < nn; i++ ) {
        rv2 = 0.0;
        for( j = 0; j <= i; j++ ) {
            rv1 = stdn.GetValueAt( j );
            rv1 *= cfac.GetValueAt( i, j );
            rv2 += rv1;
        }
        rmv->SetValueAt( i, rv2 );
    }

    //if null, assume origin for mean vector
    if( GetMu() == NULL )
        return PSL_OK;
    if(( err = rmv->Superposition( 1.0, *GetMu())) != PSL_OK )
        return err;
    return PSL_OK;
}

// -------------------------------------------------------------------------
// GenStd: generate standard normal multivariate
//
int RMVNorm::GenStd( Pslvector* rmv )
{
    if( rmv == NULL )
        return PSL_ERR_ADDRESS;

    int nn = rmv->GetSize();//dimensions
    int err;

    //just i.i.d. std. normals
    if(( err = GenRVNorms2( nn, rmv )) != PSL_OK )
        return err;

    return PSL_OK;
}

// -------------------------------------------------------------------------
// GenRVNorms: generate number nn of i.i.d. standard normal random variables
//  NOTE: space (nn) for resulting vector should be reallocated
//
int RMVNorm::GenRVNorms( int nn, Pslvector* rvs )
{
    if( rvs == NULL || rvs->GetSize() < nn )
        return PSL_ERR_ADDRESS;

    //use Polar method to generate two r.vs. at a time
    RVNorm  rvn( GetRng(), RVNorm::TRVN_Polar );
    double  rv1, rv2;
    int i, err;

    for( i = 0; i < nn; ) {
        if(( err = rvn.Gen( &rv1, &rv2 )) != PSL_OK )
            return err;
        rvs->SetValueAt( i++, rv1 );
        if( i < nn )
            rvs->SetValueAt( i++, rv2 );
    }

    return PSL_OK;
}

// -------------------------------------------------------------------------
// GenRVNorms2: generate number nn of i.i.d. standard normal deviates;
//  NOTE: space (nn) for resulting vector should be reallocated
//
int RMVNorm::GenRVNorms2( int nn, Pslvector* rvs )
{
    if( rvs == NULL || rvs->GetSize() < nn )
        return PSL_ERR_ADDRESS;

    RVNorm  rvn( GetRng(), RVNorm::TRVN_Ziggurat_M00 );
    double  rv1;
    int i, err;

    for( i = 0; i < nn; ) {
        if(( err = rvn.Gen( &rv1 )) != PSL_OK )
            return err;
        rvs->SetValueAt( i++, rv1 );
    }

    return PSL_OK;
}
