/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __AdjustedScoreMatrix__
#define __AdjustedScoreMatrix__

#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "mystring.h"
#include "myexcept.h"
#include "libpro/srcpro/DistributionMatrix.h"
#include "ScoringMatrix.h"


class FrequencyStore;
class FrequencyVector;


// _________________________________________________________________________
// Class AdjustedScoreMatrix
//
// This class implements score system between two given profiles assuming
// that the statistical parameters were previously computed for the
// database by applying universal scoring method.
//
class AdjustedScoreMatrix: public ScoringMatrix
{
public:
    AdjustedScoreMatrix(
            const FrequencyMatrix& freq_fst, const LogOddsMatrix& logo_fst, 
            const FrequencyMatrix& freq_sec, const LogOddsMatrix& logo_sec,
            const FrequencyStore*   store,
            double          infrm_threshold,
            int             thick_number,
            double          thick_percnt,
            double          mask_percnt,
            Configuration   config[NoSchemes],
            TBehaviour      s_behaviour,
            TScaling        a_scaling,
            TMask           c_masking
    );

    explicit AdjustedScoreMatrix();
    virtual ~AdjustedScoreMatrix();

    virtual const char* GetMethodName() const { return "Adjusted vector-based"; }

    virtual void        ComputeProfileScoringMatrix( bool final = false );
                                                                //compute probabilities to observe scores at each position
    virtual void        ComputeScoreProbabilities( AttributableScores* );
                                                            //compute positional probabilities of scores
    virtual void        ComputePositionalScoreProbs( AttributableScores* );

protected:
                                                                //compute score given query position and vector of frequencies
    double              ComputeScore( size_t query_pos, const FrequencyVector& ) const;

    const FrequencyStore*   GetStore() const                    { return freqstore; }
                                                                //save probability of vector in subject profile
    void                SaveVectorProbabilityAt( int pos, double prob );
    double              GetVectorProbabilityAt( int pos ) const;//return probability of vector to occur

private:
    const FrequencyStore*   freqstore;              //big array of frequency vectors
    double*                 vector_probabilities;   //probabilities of vectors comprising the subject and extracted from the pool
};


// =========================================================================
// INLINES of CLASS AdjustedScoreMatrix 
//
// SaveVectorProbabilityAt: saves probability value of vector residing in
//     subject profile
//
inline
void AdjustedScoreMatrix::SaveVectorProbabilityAt( int pos, double prob )
{
#ifdef __DEBUG__
    if( !vector_probabilities )
        throw myruntime_error( mystring( "AdjustedScoreMatrix: Memory access error." ));

    if( GetSubjectSize() <= pos || pos < 0 )
        throw myruntime_error( mystring( "AbstractScoreMatrix: Memory access error." ));
#endif

    vector_probabilities[pos] = prob;
}

// -------------------------------------------------------------------------
// GetVectorProbabilityAt: returns probability of vector to occur at the
//     alignment process
//
inline
double AdjustedScoreMatrix::GetVectorProbabilityAt( int pos ) const
{
#ifdef __DEBUG__
    if( !vector_probabilities )
        throw myruntime_error( mystring( "AdjustedScoreMatrix: Memory access error." ));

    if( GetSubjectSize() <= pos || pos < 0 )
        throw myruntime_error( mystring( "AbstractScoreMatrix: Memory access error." ));
#endif

    return vector_probabilities[pos];
}


#endif//__AdjustedScoreMatrix__
