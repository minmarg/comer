/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __rvbeta__
#define __rvbeta__

#include "ext/pslcodes.h"
#include "ext/rng.h"
#include "rv.h"
#include "rvgamma.h"

// -------------------------------------------------------------------------
// Beta random variable
//
class RVBeta: public RVar
{
public:
    enum TRVB {
        TRVB_Basic,
        TRVB_Ch78//more effective than the basic algorithm when shape pars.<1
    };
    RVBeta( Rng& rng, TRVB met = TRVB_Basic, RVGamma::TRVG gmet = RVGamma::TRVG_AD82 );
    virtual ~RVBeta();

    virtual int Gen( double* res );
    int     GenBasic( double* );
    int     GenCh78( double* );

    TRVB    GetMet() const { return met_; }
    void    SetMet( TRVB value ) { met_ = value; }

    double  GetA() const { return a_; }
    void    SetA( double value ) { a_ = value; }

    double  GetB() const { return b_; }
    void    SetB( double value ) { b_ = value; }

    RVGamma::TRVG   GetGMet() const { return rvg_.GetMet(); }
    void    SetGMet( RVGamma::TRVG value ) { rvg_.SetMet( value ); }

private:
    TRVB    met_;
    RVGamma rvg_;//gamma random object with its state variables
    double  a_;//shape alpha
    double  b_;//shape beta
    //state variables of the Ch78 algorithm
    double  Cap_, Cbp_;
    double  alpha_, beta_, gamma_;
    double  delta_, kap1_, kap2_;
};

// -------------------------------------------------------------------------
//
inline
int RVBeta::Gen( double* rv )
{
    if( rv == ( double* )0 )
        return PSL_ERR_ADDRESS;
    int err;

    switch( GetMet()) {
      case TRVB_Basic: err = GenBasic( rv ); break;
      case TRVB_Ch78: err = GenCh78( rv ); break;
      default:  return PSL_ERR_INVALID;
    };

    if( err != PSL_OK )
        return err;
    return err;
}

#endif//__rvbeta__
