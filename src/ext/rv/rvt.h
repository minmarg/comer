/***************************************************************************
 *   Copyright (C) 2010 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __rvt__
#define __rvt__

#include "ext/pslcodes.h"
#include "ext/rng.h"
#include "rv.h"

// -------------------------------------------------------------------------
// Student's t random variable
//
class RVt: public RVar
{
public:
    enum TRVt {
        TRVt_Basic,//more effective
        TRVt_M80
    };
    RVt( Rng& rng, TRVt met = TRVt_Basic );
    virtual ~RVt();

    virtual int     Gen( double* res );
    int             GenBasic( double* );
    int             GenM80( double* );

    TRVt    GetMet() const { return met_; }
    void    SetMet( TRVt value ) { met_ = value; }

    double  GetDF() const { return df_; }
    void    SetDF( double value ) { df_ = value; }

    double  GetMean() const { return mean_; }
    void    SetMean( double value ) { mean_ = value; }

    double  GetScale() const { return scl_; }
    void    SetScale( double value ) { scl_ = value; }

private:
    TRVt    met_;
    double  df_;//degrees of freedom
    double  mean_;//mean parameter
    double  scl_;//scale parameter
};

// -------------------------------------------------------------------------
//
inline
int RVt::Gen( double* rv )
{
    if( rv == ( double* )0 )
        return PSL_ERR_ADDRESS;
    int err;
    double  df = GetDF();
    double  scl = GetScale();
    double  mean = GetMean();

    switch( GetMet()) {
      case TRVt_Basic: err = GenBasic( rv ); break;
      case TRVt_M80: 
          if( df <= 2.0 )
              err = GenBasic( rv );
          else
              err = GenM80( rv );
          break;
      default:  return PSL_ERR_INVALID;
    };

    if( err != PSL_OK )
        return err;
    if( scl != 1.0 ) {
        *rv *= scl;
    }
    if( mean ) {
        *rv += mean;
    }
    return err;
}

#endif//__rvt__
