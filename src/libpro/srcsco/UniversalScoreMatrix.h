/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __UniversalScoreMatrix__
#define __UniversalScoreMatrix__

#include <math.h>

#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "mystring.h"
#include "myexcept.h"
#include "libpro/srcpro/DistributionMatrix.h"
#include "AbstractUniversalScoreMatrix.h"

// extern double rint( double x );

class FrequencyVector;


// _________________________________________________________________________
// Class UniversalScoreMatrix
//
// This class uses the frequency pool to comprise universal scoring
// system for the profiles from the database.
//
class UniversalScoreMatrix: public AbstractUniversalScoreMatrix
{
public:
    UniversalScoreMatrix(
            const FrequencyMatrix&  freq, const LogOddsMatrix& pssm,
            const FrequencyStore*   store,
            double                  infrm_threshold,
            int                     thick_number,
            double                  thick_percnt,
            Configuration           config[NoSchemes],
            TBehaviour              beh,
            TScaling                a_scaling,
            TMask                   c_masking,
            TFVectorProbabilities   distrib,
            bool                    cpu = true
    );
    explicit UniversalScoreMatrix();
    virtual ~UniversalScoreMatrix();

    virtual const char* GetMethodName() const { return "Universal vector-based"; }
                                                                //prepares scoring system for alignment of query and subject
    void                PreserveSubject( const FrequencyMatrix&, const LogOddsMatrix& );

    virtual int         GetQuerySize() const                { return GetQueryFreq().GetColumns(); }
    virtual int         GetSubjectSize() const              { return sbjct_length; }

    virtual void        ComputeProfileScoringMatrix( bool final = false );

    int                 GetThicknessNumber() const          { return thickness_number; }
    double              GetThicknessPercents() const        { return thickness_percnt; }

    virtual void        PrintParameterTable( TPrintFunction, void* vpn ) const;
    virtual void        PrintScoringMatrix( FILE* );//print scoring system
    virtual void        PrintFinal( TPrintFunction, void* vpn ) const;//output parameters in final print

protected:
    virtual double      GetScoreBase( int m, int n ) const;     //returns score for two profile positions
    TMask               GetPairMasked( int m, int n ) const;

                                                                //compute score given query position and vector of frequencies
    double              ComputeScore( size_t querypos, size_t sbjctpos ) const;
    double              ComputeScore( size_t querypos, const FrequencyVector& ) const;
    double              VectorScoreProbability( size_t querypos, const FrequencyVector& ) const;

    virtual bool        ComputeScoreProbabilitiesCPU( AttributableScores* );    //an implementation of ComputeScoreProbabilities
    virtual bool        ComputeScoreProbabilitiesMem( AttributableScores* );    //an implementation of ComputeScoreProbabilities

    void                ReallocatePairScores( int sbjct_sz );   //allocate space for pair scores
    void                DestroyPairScores();                    //destroy pair scores

    void                AllocateQueryProbabilities();           //pre-allocate space for query probabilities
    void                DestroyQueryProbabilities();            //destroy precomputed query probabilities
    void                PrecomputeQueryProbabilities();         //precompute probabilities at each position of query 
    double              GetQueryProbabilityAt( int n ) const;   //get probability at query position
    void                SetQueryProbabilityAt( int n, double ); //set probability value at query position
    const double*       GetQueryProbabilities() const       { return queryprob; }

                                                                //save score in the query-subject pair score table
    void                PreservePairScore( double score, int query_n, int sbjct_m );
    void                SetPairMasked( TMask value, int query_n, int sbjct_m );

    bool                PairThicknessConstraintsMet( int sbjct_m, int query_n ) const;

    const FrequencyMatrix&  GetQueryFreq() const    { return queryfreq; }
    const LogOddsMatrix&    GetQueryLogo() const    { return querypssm; }

    const FrequencyMatrix&  GetSbjctFreq() const;
    const LogOddsMatrix&    GetSbjctLogo() const;

    void                ResetSbjctFreq();
    void                ResetSbjctLogo();

    void                SetSbjctFreq( const FrequencyMatrix& freq ) { sbjctfreq = &freq; }
    void                SetSbjctLogo( const LogOddsMatrix& logo )   { sbjctpssm = &logo; }


    virtual void        USM_THROW( const char*, int = NOCLASS ) const;         //private throw method

private:
    const FrequencyMatrix&  queryfreq;      //reference to the query weighted frequency matrix
    const LogOddsMatrix&    querypssm;      //reference to the query log-odds matrix

    const FrequencyMatrix*  sbjctfreq;      //reference to the subject weighted frequency matrix
    const LogOddsMatrix*    sbjctpssm;      //reference to the subject log-odds matrix

    const int               thickness_number;   //thickness in number of sequences a position must have, otherwise it will be ignored
    const double            thickness_percnt;   //thickness in percentage

    double**                qspair_scores;  //query-subject pair scores
    TMask**                 qspair_mask;    //masks of scores
    int                     sbjct_length;   //length of subject
    int                     sbjct_reserved; //currently reserved length of subject
    double*                 queryprob;      //vector of probabilities of all query positions
};


// =========================================================================
// INLINES of CLASS UniversalScoreMatrix 
//
// GetScoreBase: used to access score at the profile positions
//
inline
double UniversalScoreMatrix::GetScoreBase( int sbjct_m, int query_n ) const
{
    if( ! PreferCPU())
        //since we're using private query-subject pair score table, use it anyway
        ;

#ifdef __DEBUG__
    if( qspair_scores == NULL ||
        query_n < 0 || GetQuerySize() <= query_n ||
        sbjct_m < 0 || GetSubjectSize() <= sbjct_m )
        USM_THROW( "UniversalScoreMatrix: Memory access error." );
#endif

    return qspair_scores[query_n][sbjct_m];
}

// -------------------------------------------------------------------------
// PreservePairScore: saves score in the query-subject pair score table
//
inline
void UniversalScoreMatrix::PreservePairScore( double score, int query_n, int sbjct_m )
{
#ifdef __DEBUG__
    if( qspair_scores == NULL ||
        query_n < 0 || GetQuerySize() <= query_n ||
        sbjct_m < 0 || GetSubjectSize() <= sbjct_m )
        USM_THROW( "UniversalScoreMatrix: Memory access error." );
#endif

    if( GetAutoScaling())
        score *= GetAutoScalingFactor();

    qspair_scores[query_n][sbjct_m] = score;
}

// -------------------------------------------------------------------------
// GetSbjctFreq: returns subject weighted frequency matrix
//
inline
const FrequencyMatrix& UniversalScoreMatrix::GetSbjctFreq() const
{
#ifdef __DEBUG__
    if( sbjctfreq == NULL )
        USM_THROW( "UniversalScoreMatrix: Memory access error." );
#endif
    return *sbjctfreq;
}

// GetSbjctLogo: returns subject log-odds matrix
//
inline
const LogOddsMatrix& UniversalScoreMatrix::GetSbjctLogo() const
{
#ifdef __DEBUG__
    if( sbjctpssm == NULL )
        USM_THROW( "UniversalScoreMatrix: Memory access error." );
#endif
    return *sbjctpssm;
}

// -------------------------------------------------------------------------
// GetPairMasked: returns score masking value at the position
//
inline
TMask UniversalScoreMatrix::GetPairMasked( int sbjct_m, int query_n ) const
{
#ifdef __DEBUG__
    if( qspair_mask == NULL ||
        query_n < 0 || GetQuerySize() <= query_n ||
        sbjct_m < 0 || GetSubjectSize() <= sbjct_m )
        USM_THROW( "UniversalScoreMatrix: Memory access error." );
#endif

    return qspair_mask[query_n][sbjct_m];
}

// -------------------------------------------------------------------------
// SetPairMasked: sets score masking value at the position
//
inline
void UniversalScoreMatrix::SetPairMasked( TMask value, int query_n, int sbjct_m )
{
#ifdef __DEBUG__
    if( qspair_mask == NULL ||
        query_n < 0 || GetQuerySize() <= query_n ||
        sbjct_m < 0 || GetSubjectSize() <= sbjct_m )
        USM_THROW( "UniversalScoreMatrix: Memory access error." );
#endif

    qspair_mask[query_n][sbjct_m] = value;
}

// -------------------------------------------------------------------------
// GetQueryProbabilityAt: get precomputed probability at the given query
//     position
// -------------------------------------------------------------------------

inline
double UniversalScoreMatrix::GetQueryProbabilityAt( int n ) const
{
#ifdef __DEBUG__
    if( queryprob == NULL ||
        n < 0 || GetQuerySize() <= n )
        USM_THROW( "UniversalScoreMatrix: Memory access error." );
#endif

    return queryprob[n];
}

// -------------------------------------------------------------------------
// SetQueryProbabilityAt: set probability value at query position
// -------------------------------------------------------------------------

inline
void UniversalScoreMatrix::SetQueryProbabilityAt( int n, double value )
{
#ifdef __DEBUG__
    if( queryprob == NULL ||
        n < 0 || GetQuerySize() <= n )
        USM_THROW( "UniversalScoreMatrix: Memory access error." );
#endif

    queryprob[n] = value;
}

// -------------------------------------------------------------------------
// ThicknessConstraintsBroken: returns false if thickness constraints at the
//     positions indicated are broken
// -------------------------------------------------------------------------

inline
bool UniversalScoreMatrix::PairThicknessConstraintsMet( int m, int n ) const
{
#ifdef __DEBUG__
    if( n < 0 || GetQuerySize() <= n ||
        m < 0 || GetSubjectSize() <= m )
        USM_THROW( "UniversalScoreMatrix: Memory access error." );

    if( !GetSbjctLogo().GetEffNoSequences() || !GetQueryLogo().GetEffNoSequences())
        throw myruntime_error( mystring( "UniversalScoreMatrix: Profile thickness found to be wrong." ));
#endif

    if(           GetSbjctLogo().GetMIDExpNoObservationsAt( m, PS_M ) < ( size_t )GetThicknessNumber() ||
                  GetQueryLogo().GetMIDExpNoObservationsAt( n, PS_M ) < ( size_t )GetThicknessNumber() ||
        ( double )GetSbjctLogo().GetMIDExpNoObservationsAt( m, PS_M ) / GetSbjctLogo().GetEffNoSequences() < GetThicknessPercents() ||
        ( double )GetQueryLogo().GetMIDExpNoObservationsAt( n, PS_M ) / GetQueryLogo().GetEffNoSequences() < GetThicknessPercents())
            //if number of residues within a column is less than a number specified or
            //residues ratio in a column is less than a percentage specified
            return false;

    return true;
}

// -------------------------------------------------------------------------
// USM_THROW: throws an exception
//
inline
void UniversalScoreMatrix::USM_THROW( const char* errstr, int edesc ) const
{
    throw myruntime_error( mystring( errstr ), edesc );
}


#endif//__UniversalScoreMatrix__
