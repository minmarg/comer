/***************************************************************************
 *   Copyright (C) 2010 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __rv__
#define __rv__

#include "ext/pslcodes.h"
#include "ext/rng.h"

// -------------------------------------------------------------------------
// Interface for random variable
//
class RVar
{
public:
    RVar( Rng& rng ): rng_( rng ) {};
    virtual ~RVar() {};

    virtual int     Gen( double* ) = 0;
    Rng&            GetRng() { return rng_; }

protected:
    Rng&    rng_;
};

#endif//__rv__
