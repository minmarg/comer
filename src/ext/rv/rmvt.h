/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __rmvt__
#define __rmvt__

#include "ext/pslcodes.h"
#include "ext/rng.h"
#include "ext/pslvector.h"
#include "ext/spdmatrix.h"
#include "rmv.h"

// -------------------------------------------------------------------------
// Normal random variable
//
class RMVt: public RMV
{
public:
    RMVt( Rng& rng );
    virtual ~RMVt();

    //if not scale matrix and location vector are given, 
    //  generate standard t multivariates
    virtual int     Gen( Pslvector* );

    double  GetDF() const { return df_; }
    void    SetDF( double value ) { df_ = value; }

    const Pslvector*    GetMu() const { return mu_; }
    void SetMu( const Pslvector* mu ) { mu_ = mu; }

    const SPDmatrix*    GetScaleM() const { return scalemat_; }
    void SetScaleM( const SPDmatrix* sigma ) { scalemat_ = sigma; }

private:
    double              df_;//degrees of freedom
    const Pslvector*    mu_;//location vector
    const SPDmatrix*    scalemat_;//scale matrix
};

#endif//__rmvt__
