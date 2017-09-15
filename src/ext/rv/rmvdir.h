/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __rmvdir__
#define __rmvdir__

#include "ext/pslcodes.h"
#include "ext/rng.h"
#include "ext/pslvector.h"
#include "rmv.h"
#include "rvgamma.h"

// -------------------------------------------------------------------------
// Dirichlet random variable
//
class RMVDir: public RMV
{
public:
    RMVDir( Rng& rng, RVGamma::TRVG gmet = RVGamma::TRVG_AD82 );
    virtual ~RMVDir();

    virtual int Gen( Pslvector* rmv );

    const Pslvector*    GetAlphs() const { return alphas_; }
    void SetAlphs( const Pslvector* as ) { alphas_ = as; }

    RVGamma::TRVG   GetGMet() const { return rvg_.GetMet(); }
    void    SetGMet( RVGamma::TRVG value ) { rvg_.SetMet( value ); }

protected:
    int GenScaled( Pslvector* rmv );//scaled version to increase precision for small alphas

private:
    RVGamma rvg_;//gamma random object with its state variables
    const Pslvector* alphas_;//concentration parameters
};

#endif//__rmvdir__
