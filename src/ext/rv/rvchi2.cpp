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
#include "rvchi2.h"


// =========================================================================
// RVChi2
//
RVChi2::RVChi2( Rng& rng, RVGamma::TRVG met )
:   RVar( rng ),
    rvg_( rng, met ),
    df_( 1.0 )
{
    rvg_.SetScale( 2.0 );
}

RVChi2::~RVChi2()
{
}

// -------------------------------------------------------------------------
// chi-squared distributed random variable;
//
//  pdf: 1 / (Gamma(df/2) 2^(df/2)) x^(df/2-1) e^(-x/2),  (x,df>=0)
//
int RVChi2::Gen( double* result )
{
    if( !result )
        return PSL_ERR_ADDRESS;

    double df = GetDF();

    if( df <= 0.0 )
        return PSL_ERR_DOMAIN;

    rvg_.SetShape( 0.5 * df );
    return rvg_.Gen( result );
}
