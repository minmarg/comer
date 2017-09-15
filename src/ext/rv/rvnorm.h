/***************************************************************************
 *   Copyright (C) 2010 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __rvnorm__
#define __rvnorm__

#include "ext/pslcodes.h"
#include "ext/rng.h"
#include "rv.h"

// -------------------------------------------------------------------------
// Normal random variable
//
class RVNorm: public RVar
{
public:
    enum TRVN {
        TRVN_Polar,
        TRVN_Ratio,
        TRVN_Ziggurat_M00
    };
    RVNorm( Rng& rng, TRVN met = TRVN_Ziggurat_M00 );
    virtual ~RVNorm();

    virtual int     Gen( double* res ) { return Gen( res, (double*)0 ); }
    int             Gen( double*, double* );
    int             GenPolar( double*, double* = (double*)0 );
    int             GenRatio( double* );
    int             GenZigMsg00( double* );

    TRVN    GetMet() const { return met_; }
    void    SetMet( TRVN value ) { met_ = value; }

    double  GetStd() const { return std_; }
    void    SetStd( double value ) { std_ = value; }

    double  GetMean() const { return mean_; }
    void    SetMean( double value ) { mean_ = value; }

private:
    TRVN    met_;
    double  std_;
    double  mean_;
};

// -------------------------------------------------------------------------
//
inline
int RVNorm::Gen( double* rv1, double* rv2 )
{
    if( rv1 == NULL )
        return PSL_ERR_ADDRESS;
    int err;
    double  std = GetStd();
    double  mean = GetMean();
    bool    brv2 = GetMet() == TRVN_Polar && rv2 != NULL;

    switch( GetMet()) {
      case TRVN_Polar: err = GenPolar( rv1, rv2 ); break;
      case TRVN_Ratio: err = GenRatio( rv1 ); break;
      case TRVN_Ziggurat_M00: err = GenZigMsg00( rv1 ); break;
      default:  return PSL_ERR_INVALID;
    };

    if( err != PSL_OK )
        return err;
    if( std != 1.0 ) {
        *rv1 *= std;
        if( brv2 )
            *rv2 *= std;
    }
    if( mean ) {
        *rv1 += mean;
        if( brv2 )
            *rv2 += mean;
    }
    return err;
}

#endif//__rvnorm__
