/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include "ext/psl.h"
#include "ext/pslcodes.h"
#include "ext/rng.h"
#include "rmvdir.h"


// =========================================================================
// RMVDir
//
RMVDir::RMVDir( Rng& rng, RVGamma::TRVG gmet )
:   RMV( rng ),
    rvg_( rng, gmet ),
    alphas_( NULL )
{
    rvg_.SetScale( 1.0 );
}

RMVDir::~RMVDir()
{
}

// -------------------------------------------------------------------------
// Gen: generate r.v. distributed according to Dirichlet distribution;
//  (ref. S.Asmussen. Stochastic Simulation: Algorithms and Analysis, p.51; 
//  L.Devroye. (1986) Non Uniform Random Variate Generation, p.593);
//
//  pdf: Gamma(SUM alpha_i) / PROD Gamma(alpha_i) PROD x_i^(alpha_i-1),
//          (alpha_i>0, 0<=x_i<=1 and SUM x_i=1, dimensions>=2)
//
//  NOTE: alpha parameters must be preset
//
int RMVDir::Gen( Pslvector* rmv )
{
    if( !rmv )
        return PSL_ERR_ADDRESS;

    const Pslvector* alphs = GetAlphs();

    if( !alphs )
        return PSL_ERR_INVALID;

    double  a, val, xsum = 0.0;
    int p = alphs->GetSize();//dimensions
    int n, err;

    if( p < 2 || p != rmv->GetSize())
        return PSL_ERR_INVALID;

    for( n = 0; n < p; n++ ) {
        a = alphs->GetValueAt( n );
        rvg_.SetShape( a );
        if(( err = rvg_.Gen( &val )) != PSL_OK )
            return err;
        rmv->SetValueAt( n, val );
    }

    xsum = rmv->Sum();

    if( xsum < SLC_SQRT_DBL_MIN ) {
        //alhpa parameters are small, use scaling to increase precision
        return GenScaled( rmv );
    }

    if(( err = rmv->MultiplyBy( 1.0 / xsum )) != PSL_OK )
        return err;

    return PSL_OK;
}

// -------------------------------------------------------------------------
// GenScaled: generate Dirichlet r.v. when concentration parameters are 
//  small (<1), using prescaling of the parameters;
//  for small alphas use recursive gamma variable relation given in
//  Marsaglia & Tsang. (2000) ACM Transactions on Mathematical Software 26(3), 363-72;
//  see also RVGamma::GenMT00.
//
int RMVDir::GenScaled( Pslvector* rmv )
{
    if( !rmv )
        return PSL_ERR_ADDRESS;

    const Pslvector* alphs = GetAlphs();

    if( !alphs )
        return PSL_ERR_INVALID;

    double  u, a, val, xsum = 0.0, asum = 0.0;
    double  lmax = -SLC_DBL_MAX;
    int p = alphs->GetSize();//dimensions
    int n, err;

    if( p < 2 || p != rmv->GetSize())
        return PSL_ERR_INVALID;

    //find largest eponent
    for( n = 0; n < p; n++ ) {
        u = GetRng().GetDouble0();
        a = alphs->GetValueAt( n );
        val = log(u) / a;
        rmv->SetValueAt( n, val );
        if( lmax < val )
            lmax = val;
    }
    //scale each variable
    for( n = 0; n < p; n++ ) {
        val = rmv->GetValueAt( n );
        rmv->SetValueAt( n, exp( val - lmax ));
    }

    //apply recursive relation to generate gamma r.v.;
    //(valid for all shape parameters, as charact. functions of a 
    //  product and gamma r.v. are equal)
    for( n = 0; n < p; n++ ) {
        a = alphs->GetValueAt( n );
        rvg_.SetShape( a + 1.0 );
        if(( err = rvg_.Gen( &val )) != PSL_OK )
            return err;
        //multiply scaled u^(1/a) by gamma_(a+1)
        rmv->MulValueAt( n, val );
    }

    xsum = rmv->Sum();

    if( xsum <= 0.0 ) {
        asum = alphs->Sum();
        if( 0.0 < asum ) {
            for( n = 0; n < p; n++ ) {
                a = alphs->GetValueAt( n );
                rmv->SetValueAt( n, a / asum );
            }
            return PSL_OK;
        }
        return PSL_ERR_ILLEGAL;
    }

    if(( err = rmv->MultiplyBy( 1.0 / xsum )) != PSL_OK )
        return err;

    return PSL_OK;
}
