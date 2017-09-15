/***************************************************************************
 *   Copyright (C) 2012 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __PrlHDPUniversalScoreMatrix__
#define __PrlHDPUniversalScoreMatrix__

#include "mystring.h"
#include "myexcept.h"
#include "debug.h"
#include "ext/pslvector.h"
#include "libHDP/HDPbase.h"
#include "ParallelUniversalScoreMatrix.h"


// _________________________________________________________________________
// Class PrlHDPUniversalScoreMatrix
//
// Scaling of profile scoring matrix to the given lambda value.
// Parallelized via the MPI interface.
// Using HDP scores
//
class PrlHDPUniversalScoreMatrix: public ParallelUniversalScoreMatrix
{
public:
    PrlHDPUniversalScoreMatrix(
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

    explicit PrlHDPUniversalScoreMatrix();
    virtual ~PrlHDPUniversalScoreMatrix();

    virtual const char* GetMethodName() const { return "Parallel universal HDP vector-based"; }
    virtual void        ScaleScoringMatrix();

protected:
    virtual double  ComputeScore( const FrequencyVector&, const FrequencyVector&, size_t, size_t, bool = false );
    virtual double  VectorScoreProbability( const FrequencyVector&, const FrequencyVector&, 
                                            size_t, size_t, bool = false ) const;

    const Basin*        GetBasin() const { return basin_; }
    Basin*              GetBasin() { return basin_; }

    const HDPbase*      GetHDPbase() const { return hdpbase_; }
    HDPbase*            GetHDPbase() { return hdpbase_; }

    const Pslvector&    GetQueValues() const { return quevals_; }
    Pslvector&          GetQueValues() { return quevals_; }

    const Pslvector&    GetSubValues() const { return subvals_; }
    Pslvector&          GetSubValues() { return subvals_; }

    double              GetProbNormTerm() const { return prnormterm_; }
    void                SetProbNormTerm( double value ) { prnormterm_ = value; }

    int                 GetNoSupPoss() const { return nosupports_; }
    void                SetNoSupPoss( int nops ) { nosupports_ = nops; }

    static bool         GetUseHDPBProbs() { return hdpbprobs_; }
    static void         SetUseHDPBProbs( bool value ) { hdpbprobs_ = value; }

    void    InitHDP();
    void    InitBasin( int size );
    void    InitHDPbase();
    void    DestroyBasin() { if( basin_ ) delete basin_; basin_ = NULL; }
    void    DestroyHDPbase() { if( hdpbase_ ) delete hdpbase_; hdpbase_ = NULL; }

private:
    Basin*          basin_;//basin of DB vectors
    HDPbase*        hdpbase_;//base to calculate HDP scores
    Pslvector       quevals_;//precalculated values of P(Q|B)
    Pslvector       subvals_;//precalculated values of P(Q)
    double          prnormterm_;//prob. normalizing term
    int             nosupports_;//number of support vectors
    static bool     hdpbprobs_;//background probabilities are multivariate t-dist.
};


// =========================================================================
// PrlHDPUniversalScoreMatrix INLINES
//
// InitBasin: initialize basin of vectors
//
inline
void PrlHDPUniversalScoreMatrix::InitBasin( int size )
{
    if( size < 1 )
        USM_THROW("PrlHDPUniversalScoreMatrix: InitBasin: Invalid size.");
    DestroyBasin();
    basin_ = new Basin( size );
    if( basin_ == NULL )
        USM_THROW("PrlHDPUniversalScoreMatrix: InitBasin: Not enough memory.");
    basin_->SetDestroy( true );
}

// -------------------------------------------------------------------------
// InitHDPbase: initialize hdpbase structure
//
inline
void PrlHDPUniversalScoreMatrix::InitHDPbase()
{
    DestroyHDPbase();
    hdpbase_ = new HDPbase();
    if( hdpbase_ == NULL )
        throw myruntime_error("PrlHDPUniversalScoreMatrix: Not enough memory.");
}



#endif//__PrlHDPUniversalScoreMatrix__
