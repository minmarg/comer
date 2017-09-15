/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __rmvnorm__
#define __rmvnorm__

#include "ext/pslcodes.h"
#include "ext/rng.h"
#include "ext/pslvector.h"
#include "ext/spdmatrix.h"
#include "rmv.h"

// -------------------------------------------------------------------------
// Normal random variable
//
class RMVNorm: public RMV
{
public:
    RMVNorm( Rng& rng );
    virtual ~RMVNorm();

    //if not sigma and mu are given, generate standard normal multivariates
    virtual int     Gen( Pslvector* );
    //generate standard normal multivariates
    int     GenStd( Pslvector* );

    const Pslvector*    GetMu() const { return mu_; }
    void SetMu( const Pslvector* mu ) { mu_ = mu; }

    const SPDmatrix*    GetSigma() const { return sigma_; }
    void SetSigma( const SPDmatrix* sigma ) { sigma_ = sigma; }

    //generate nn i.i.d. standard normal r.vs
    int GenRVNorms( int nn, Pslvector* rvs );
    int GenRVNorms2( int nn, Pslvector* rvs );

private:
    const Pslvector*    mu_;
    const SPDmatrix*    sigma_;
};

#endif//__rmvnorm__
