/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __ParallelUniversalScoreMatrix__
#define __ParallelUniversalScoreMatrix__

#include "mystring.h"
#include "myexcept.h"
#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "libpro/srcsco/AbstractUniversalScoreMatrix.h"

class FrequencyVector;


// _________________________________________________________________________
// Class ParallelUniversalScoreMatrix
//
// Scaling of profile scoring matrix to the given lambda value.
// Parallelized via the MPI interface
//
class ParallelUniversalScoreMatrix: public AbstractUniversalScoreMatrix
{
public:
    //typedefs...
    typedef int ( *TBcastFunction )( char* msg, size_t, bool throw_on_error );
    typedef int ( *TSendFunction )( char* msg, size_t, bool throw_on_error );
    typedef int ( *TRecvFunction )( char* msg, size_t, bool throw_on_error );
    typedef int ( *TBlockFunction )( bool throw_on_error );

public:
    ParallelUniversalScoreMatrix(
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

    explicit ParallelUniversalScoreMatrix();
    virtual ~ParallelUniversalScoreMatrix();

    virtual const char* GetMethodName() const { return "Parallel universal vector-based"; }

    void                TellSlavesToTerminate( bool throw_on_error );

    virtual void        ComputeProfileScoringMatrix( bool = false );
    virtual void        ScaleScoringMatrix();

    virtual void        PrintParameterTable( TPrintFunction, void* vpn ) const;
    virtual void        PrintScoringMatrix( FILE* );            //print the computed scoring system
    virtual void        PrintFinal( TPrintFunction, void* vpn ) const;//output parameters in final print

protected:
    ParallelUniversalScoreMatrix(
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
            TMask                   c_masking,
            TType                   type
    );
    virtual double      GetScoreBase( int m, int n ) const;     //returns score for two profile positions
                                                                //compute score given query position and vector of frequencies
    virtual double  ComputeScore( const FrequencyVector&, const FrequencyVector&, size_t, size_t, bool = false );
    virtual double  VectorScoreProbability( const FrequencyVector&, const FrequencyVector&, 
                                            size_t, size_t, bool = false ) const;

    virtual bool        ComputeScoreProbabilitiesCPU( AttributableScores* );    //an implementation of ComputeScoreProbabilities
    virtual bool        ComputeScoreProbabilitiesMem( AttributableScores* );    //an implementation of ComputeScoreProbabilities
    virtual void        MultiplierProgress( AttributableScores* );              //method invoked just before saving of the multiplier

    bool                IsMasterProcess() const                 { return master;        }
    int                 GetMyMPIRank() const                    { return mpi_rank;      }
    int                 GetMPIRingSize() const                  { return mpi_size;      }

    char*               GetSRMBuffer()                          { return srm_buffer; }
    const char*         GetSRMBuffer() const                    { return srm_buffer; }

    bool                NotToPerformScaling() const             { return not_to_scale;  }

    virtual void        USM_THROW( const char*, int = NOCLASS ) const;  //private throw method

private:
                                                                //computation logic for the master process
    void                MasterComputations( AttributableScores* );
    void                SlaveComputations();                    //computation logic for the slave processes

    bool                UsingSpecificLambda() const             { return use_spec_lambda; }

private:
    bool                    master;         //whether a process is master
    int                     mpi_rank;       //rank of this process in the MPI ring
    int                     mpi_size;       //number of processes running in the MPI ring

    TBcastFunction          bcast_function; //broadcast-function address
    TSendFunction           send_function;  //address of function to send messages
    TRecvFunction           recv_function;  //address of function to receive messages
    TBlockFunction          block_function; //address of function to block until all processes are ready

    char*                   srm_buffer;     //send-receive buffer
    bool                    not_to_scale;   //whether or not to perform scaling of matrix
    bool                    use_spec_lambda;//flag of usage of a specific user-specified lambda
};


// =========================================================================
// INLINES of CLASS ParallelUniversalScoreMatrix
//
// USM_THROW: throws an exception
//
inline
void ParallelUniversalScoreMatrix::USM_THROW( const char* errstr, int edesc ) const
{
    throw myruntime_error( mystring( errstr ), edesc );
}

// -------------------------------------------------------------------------
// GetScoreBase: should not be used by objects of this class
//
inline
double ParallelUniversalScoreMatrix::GetScoreBase( int, int ) const
{
    USM_THROW( "ParallelUniversalScoreMatrix: The class should not be used for reference of scores." );

    return 0.0;
}

// -------------------------------------------------------------------------
// ComputeProfileScoringMatrix should not be called from the class
//
inline
void ParallelUniversalScoreMatrix::ComputeProfileScoringMatrix( bool )
{
    USM_THROW( "ParallelUniversalScoreMatrix: ComputeProfileScoringMatrix should not be called." );
}



#endif//__ParallelUniversalScoreMatrix__
