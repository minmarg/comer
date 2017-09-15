/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __rvchi2__
#define __rvchi2__

#include "ext/pslcodes.h"
#include "ext/rng.h"
#include "rv.h"
#include "rvgamma.h"

// -------------------------------------------------------------------------
// Chi-squared random variable
//
class RVChi2: public RVar
{
public:
    RVChi2( Rng& rng, RVGamma::TRVG met = RVGamma::TRVG_AD82 );
    virtual ~RVChi2();

    virtual int Gen( double* res );

    double  GetDF() const { return df_; }
    void    SetDF( double value ) { df_ = value; }

    RVGamma::TRVG   GetGMet() const { return rvg_.GetMet(); }
    void    SetGMet( RVGamma::TRVG value ) { rvg_.SetMet( value ); }

private:
    RVGamma rvg_;//gamma random object with its state variables
    double  df_;//degrees of freedom
};

#endif//__rvchi2__
