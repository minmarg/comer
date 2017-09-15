/***************************************************************************
 *   Copyright (C) 2009 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __PrlLSOUniversalScoreMatrix__
#define __PrlLSOUniversalScoreMatrix__

#include "mystring.h"
#include "myexcept.h"
#include "debug.h"
#include "ParallelUniversalScoreMatrix.h"


// _________________________________________________________________________
// Class PrlLSOUniversalScoreMatrix
//
// Scaling of profile scoring matrix to the given lambda value.
// Parallelized via the MPI interface.
// Using LSO scores
//
class PrlLSOUniversalScoreMatrix: public ParallelUniversalScoreMatrix
{
public:
    PrlLSOUniversalScoreMatrix(
            const FrequencyStore*   store,
            Configuration           config[NoSchemes],
            bool                    no_scaling,
            bool                    using_spec_lambda,
            double                  infrm_threshold,
            TFVectorProbabilities   distrib,
            bool                    master_flag,
            int                     rank,
            int                     ring_size,
            TBcastFunction          func_bcast,
            TSendFunction           func_send,
            TRecvFunction           func_recv,
            TBlockFunction          func_block,
            TScaling                a_scaling,
            TMask                   c_masking
    );

    explicit PrlLSOUniversalScoreMatrix();
    virtual ~PrlLSOUniversalScoreMatrix();

    virtual const char* GetMethodName() const { return "Parallel universal LSO vector-based"; }

protected:
    virtual double  ComputeScore( const FrequencyVector&, const FrequencyVector&, size_t, size_t, bool = false );
    virtual double  VectorScoreProbability( const FrequencyVector&, const FrequencyVector&, 
                                            size_t, size_t, bool = false ) const;
private:
};


// =========================================================================
// PrlLSOUniversalScoreMatrix INLINES
//


#endif//__PrlLSOUniversalScoreMatrix__
