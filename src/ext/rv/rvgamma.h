/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __rvgamma__
#define __rvgamma__

#include "ext/pslcodes.h"
#include "ext/rng.h"
#include "rv.h"

// -------------------------------------------------------------------------
// Gamma random variable
// NOTE: TRVG_AD82 uses state variables: DEFINE object GLOBALLY
//
class RVGamma: public RVar
{
public:
    enum TRVG {
        TRVG_AD82,//starting to be more effective than MT00 for shape par.>10
        TRVG_MT00
    };
    RVGamma( Rng& rng, TRVG met = TRVG_AD82 );
    virtual ~RVGamma();

    virtual int     Gen( double* res );
    int             GenAD82( double*, double shape );
    int             GenMT00( double*, double shape );

    TRVG    GetMet() const { return met_; }
    void    SetMet( TRVG value ) { met_ = value; }

    double  GetShape() const { return shape_; }
    void    SetShape( double value ) { shape_ = value; }

    double  GetScale() const { return scale_; }
    void    SetScale( double value ) { scale_ = value; }

private:
    TRVG    met_;
    double  shape_;
    double  scale_;
    //state variables of the AD82 algorithm
    double  Aap_;
    double  Aapp_;
    double  As_, As_2_, Ad_;//step 1
    double  Aq0_, Ab_, Asgm_, Ac_;//step 4
    //state variables of the MT00 algorithm
    double  Map_;
    double  Md_, Mc_;
};

// -------------------------------------------------------------------------
//
inline
int RVGamma::Gen( double* rv )
{
    if( rv == ( double* )0 )
        return PSL_ERR_ADDRESS;
    int err;
    double  shape = GetShape();
    double  scale = GetScale();

    if( scale <= 0.0 )
        return PSL_ERR_INVALID;

    switch( GetMet()) {
      case TRVG_AD82: err = GenAD82( rv, shape ); break;
      case TRVG_MT00: err = GenMT00( rv, shape ); break;
      default:  return PSL_ERR_INVALID;
    };

    if( err != PSL_OK )
        return err;
    if( scale != 1.0 ) {
        *rv *= scale;
    }
    return err;
}

#endif//__rvgamma__
