/***************************************************************************
 *   Copyright (C) 2009 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "lmpi/msgcodes.h"

#include "rc.h"
#include "data.h"
#include "libpro/srcpro/FrequencyStore.h"
#include "PrlLSOUniversalScoreMatrix.h"



// =========================================================================
// CLASS PrlLSOUniversalScoreMatrix
//
// constructor
//
PrlLSOUniversalScoreMatrix::PrlLSOUniversalScoreMatrix(
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
    )
:
    ParallelUniversalScoreMatrix(
        store,
        config,
        no_scaling,
        using_spec_lambda,
        infrm_threshold,
        distrib,
        master_flag,
        rank,
        ring_size,
        func_bcast,
        func_send,
        func_recv,
        func_block,
        a_scaling,
        c_masking,
        ParallelLSOUniversal )
{
}

// -------------------------------------------------------------------------
// default constructor is invalid
//
PrlLSOUniversalScoreMatrix::PrlLSOUniversalScoreMatrix()
:   ParallelUniversalScoreMatrix()
{
}

// -------------------------------------------------------------------------
// destructor
//
PrlLSOUniversalScoreMatrix::~PrlLSOUniversalScoreMatrix()
{
}

// -------------------------------------------------------------------------
// ComputeScore: compute score given two vectors of frequencies
//
double PrlLSOUniversalScoreMatrix::ComputeScore(
    const FrequencyVector& vect_1st,
    const FrequencyVector& vect_2nd, 
    size_t, size_t, bool )
{
    static const double lambda = LOSCORES.StatisParam( Ungapped, Lambda ); //reference lambda
    double  score = 0.0;
    double  offset = 0.0003;
    size_t  fst_enos = vect_1st.GetThickness();
    size_t  sec_enos = vect_2nd.GetThickness();
    int a;

    for( a = 0; a < NUMAA; a++ )
        score += exp(( double )vect_2nd.GetScoreAt( a )/( double )SCALE_CONSTANT * lambda ) * 
                 exp(( double )vect_1st.GetScoreAt( a )/( double )SCALE_CONSTANT * lambda ) * 
                 LOSCORES.PROBABility( a );

    score = log( score ) - offset;
    score *= GetMultiplier();
    return score;
}

// -------------------------------------------------------------------------
// VectorScoreProbability: compute score probability given two vectors of
//  frequencies
//
double PrlLSOUniversalScoreMatrix::VectorScoreProbability(
    const FrequencyVector& rowvect,
    const FrequencyVector& colvect, 
    size_t, size_t, bool ) const
{
    double  pro1 = rowvect.GetProbability();
    double  pro2 = colvect.GetProbability();

    return 1.;
    return pro1 * pro2;
}
