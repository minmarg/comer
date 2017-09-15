/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __rmv__
#define __rmv__

#include "ext/pslcodes.h"
#include "ext/rng.h"
#include "ext/pslvector.h"

// -------------------------------------------------------------------------
// Interface for random multivariates
//
class RMV
{
public:
    RMV( Rng& rng ): rng_( rng ) {};
    virtual ~RMV() {};

    virtual int     Gen( Pslvector* ) = 0;
    Rng&            GetRng() { return rng_; }

protected:
    Rng&    rng_;
};

#endif//__rmv__
