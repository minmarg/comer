/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __rvexp__
#define __rvexp__

#include "ext/pslcodes.h"
#include "ext/rng.h"
#include "rv.h"

// -------------------------------------------------------------------------
// Exponential random variable
//
class RVExp: public RVar
{
public:
    RVExp( Rng& rng );
    virtual ~RVExp();

    virtual int Gen( double* res ) { return rvgexp( res ); }

    double  GetScale() const { return scale_; }
    void    SetScale( double value ) { scale_ = value; }

protected:
    int     rvgexp( double* result );
    int     rvgexprand( double* result );//std. exponential distribution

private:
    double  scale_;//1/scale * exp(-x/scale)
};



#endif//__rvexp__
