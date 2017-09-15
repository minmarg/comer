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
#include "ext/psl.h"
#include "rvnorm.h"
#include "rvchi2.h"
#include "rmvnorm.h"
#include "rmvt.h"

// =========================================================================
// RMVt
//
RMVt::RMVt( Rng& rng )
:   RMV( rng ),
    df_( 1.0 ),
    mu_( NULL ),
    scalemat_( NULL )
{
}

RMVt::~RMVt()
{
}

// -------------------------------------------------------------------------
// Gen: generate t-distributed random multivariates
//  (ref. S.Asmussen. Stochastic Simulation: Algorithms and Analysis, p.51;
//  Lin. (1972) J. Multivariate Anal. 2, 339-44;
//  Kotz & Nadarajah. Multivariate t Distributions and Their Applications, p.2);
//  if no scale matrix and location vector are given, generate standard t multivariates
//
//  pdf: Gamma((df+p)/2)/((pi df)^(p/2) Gamma(df/2) |S|^{1/2})
//          (1+(x-mu)'S^{-1}(x-mu)/df)^(-(df+p)/2),
//      S -- scale matrix, p -- the dimensions of vectors
//
int RMVt::Gen( Pslvector* rmv )
{
    if( rmv == NULL )
        return PSL_ERR_ADDRESS;

    double  df = GetDF();
    double  c2;
    int nn = rmv->GetSize();//dimensions
    int err;

    if( df <= 0.0 )
        return PSL_ERR_INVALID;

    RMVNorm rmvn( GetRng());
    RVChi2  rc2( GetRng());

    rmvn.SetSigma( GetScaleM());
    rc2.SetDF( df );

    if(( err = rmvn.Gen( rmv )) != PSL_OK )
        return err;
    if(( err = rc2.Gen( &c2 )) != PSL_OK )
        return err;

    if( c2 < 0.0 )
        return PSL_ERR_ILLEGAL;
    if( c2 == 0.0 ) {
        rmv->SetAllToValue( SLC_DBL_MAX );
        return PSL_OK;
    }
    c2 = sqrt( df / c2 );
//     c2 = 1.0 /sqrt( c2 / df );

    if(( err = rmv->MultiplyBy( c2 )) != PSL_OK )
        return err;

    //if null, assume origin for mean vector
    if( GetMu() == NULL )
        return PSL_OK;
    if(( err = rmv->Superposition( 1.0, *GetMu())) != PSL_OK )
        return err;
    return PSL_OK;
}
