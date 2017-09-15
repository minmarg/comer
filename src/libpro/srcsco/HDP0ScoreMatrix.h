/***************************************************************************
 *   Copyright (C) 2012 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __HDP0ScoreMatrix__
#define __HDP0ScoreMatrix__

#include "debug.h"
#include "types.h"
#include "compdef.h"

#include "mystring.h"
#include "myexcept.h"
#include "ext/pslvector.h"
#include "libpro/srcpro/DistributionMatrix.h"
#include "libHDP/HDPbase.h"
#include "AbstractScoreMatrix.h"


class HDP0ScoreMatrix: public AbstractScoreMatrix
{
public:
    HDP0ScoreMatrix(
            const HDPbase*  hdpparent,
            const FrequencyMatrix&  freq_fst, const LogOddsMatrix& logo_fst, 
            const FrequencyMatrix&  freq_sec, const LogOddsMatrix& logo_sec,
            double          infrm_threshold,
            int             thick_number,
            double          thick_percnt,
            double          mask_percnt,
            Configuration   config[NoSchemes],
            TBehaviour      s_behaviour,
            TScaling        a_scaling,
            TMask           c_masking
    );
    explicit HDP0ScoreMatrix();
    virtual ~HDP0ScoreMatrix();
    virtual double      GetScore( int m, int n ) const;             //get score at the positions
    virtual const char* GetMethodName() const   { return "Profile-specific HDP0-based"; }
    virtual void        ComputeProfileScoringMatrix( bool final = false );
    virtual double      ComputeExpectation( double, double* = NULL, double* = NULL, double* = NULL, double* = NULL ) const;

    int                 GetThicknessNumber() const      { return thickness_number; }
    double              GetThicknessPercents() const    { return thickness_percnt; }
    double              GetMaskscalePercents() const    { return maskscale_percnt; }

    virtual void        OptimizeTargetFrequencies();
    void                OptimizeTargetFrequenciesObs();

    virtual void        ComputeScoreProbabilities( AttributableScores* );
    virtual void        ComputePositionalScoreProbs( AttributableScores* );
    virtual bool        ReduceScores( AttributableScores* );

    virtual void        PrintParameterTable( TPrintFunction, void* vpn ) const;
    virtual void        PrintFinal( TPrintFunction, void* vpn ) const;

    virtual void        PrintScoringMatrix( FILE* );

protected:
    HDP0ScoreMatrix(
            const HDPbase*  hdpparent,
            const FrequencyMatrix& freq_fst, const LogOddsMatrix& logo_fst, 
            const FrequencyMatrix& freq_sec, const LogOddsMatrix& logo_sec,
            double          infrm_threshold,
            int             thick_number,
            double          thick_percnt,
            double          mask_percnt,
            Configuration   config[NoSchemes],
            TBehaviour      s_behaviour,
            TScaling        a_scaling,
            TMask           c_masking,
            TType           type
    );

    void                    PrivateInit();
    double                  ComputeScore( int m, int n, bool final );
    double                  ComputeScoreMtx( int m, int n, bool final );
    double                  ComputeScoreVec( int m, int n, bool final );
    double                  ComputeScoreVec( int m, int n, bool final,
                                          const LogOddsMatrix& sec, const LogOddsMatrix& fst, 
                                          Pslvector& secvec, Pslvector& fstvec,
                                          Basin& secbas, Basin& fstbas, HDPbase& basHDP );
    void                    ReduceScoresBy( double delta );

    void                    ComputeScoreProbabilitiesII( AttributableScores* );
    void                    ComputeScoreProbabilitiesFPI( AttributableScores* );
 
    const FrequencyMatrix&  GetQueryFreq() const    { return freq_fst_; }
    const LogOddsMatrix&    GetQueryLogo() const    { return logo_fst_; }

    const FrequencyMatrix&  GetSbjctFreq() const    { return freq_sec_; }
    const LogOddsMatrix&    GetSbjctLogo() const    { return logo_sec_; }

    bool                    ThicknessConstraintsMet( int m, int n ) const;

    void    InitHDP();
    void    InitHDPMtx();
    void    InitHDPVec();
    void    InitHDPVec( const LogOddsMatrix& fst, const LogOddsMatrix& sec,
                     Pslvector& fstvec, Pslvector& secvec,
                     Basin& fstbas, Basin& secbas, HDPbase& basHDP );
    void    InitHDPVecobs( const LogOddsMatrix& fst, const LogOddsMatrix& sec,
                     Pslvector& fstvec, Pslvector& secvec,
                     Basin& fstbas, Basin& secbas, HDPbase& basHDP );

    void    InitHDPbase();
    void    InitQueryBasin();
    void    InitSbjctBasin();

    void    DestroyHDPbase() { if( hdpbase_ ) delete hdpbase_; hdpbase_ = NULL; }
    void    DestroyQueryBasin() { if( quebasin_ ) delete quebasin_; quebasin_ = NULL; }
    void    DestroySbjctBasin() { if( subbasin_ ) delete subbasin_; subbasin_ = NULL; }

    const HDPbase*  GetHDPparent() const { return hdpparent_; }
    const HDPbase*  GetHDPbase() const { return hdpbase_; }
    HDPbase*        GetHDPbase() { return hdpbase_; }
    const Basin*    GetQueryBasin() const { return quebasin_; }
    Basin*          GetQueryBasin() { return quebasin_; }
    const Basin*    GetSbjctBasin() const { return subbasin_; }
    Basin*          GetSbjctBasin() { return subbasin_; }

    const Pslvector&    GetQueValues() const { return quevals_; }
    Pslvector&          GetQueValues() { return quevals_; }

    const Pslvector&    GetSubValues() const { return subvals_; }
    Pslvector&          GetSubValues() { return subvals_; }

    //helper routines for optimization of target frequencies
    void            IncorporateTargetFrequencies( const Pslvector& );

    static bool         GetUseMtxVarT() { return mtxvart_; }
    static bool         GetUseHDPBProbs() { return hdpbprobs_; }
    static void         SetUseHDPBProbs( bool value ) { hdpbprobs_ = value; }

private:
    const FrequencyMatrix&  freq_fst_;//reference to the first weighted frequency matrix
    const LogOddsMatrix&    logo_fst_;//reference to the first log-odds matrix

    const FrequencyMatrix&  freq_sec_;//reference to the second weighted frequency matrix
    const LogOddsMatrix&    logo_sec_;//reference to the second log-odds matrix

    const int       thickness_number;//thickness in number of sequences a position must have, otherwise it will be ignored
    const double    thickness_percnt;//thickness in percentage
    const double    maskscale_percnt;//scaling of masked positions in percentage

    Pslvector               scores_;//scores for optimization of target frequencies
    Pslvector               rprobs_;//row background probabilities of the score system
    Pslvector               cprobs_;//column background probabilities of the score system

    const HDPbase*  hdpparent_;//parent HDP structure
    HDPbase*        hdpbase_;//HDP substructure
    Basin*          quebasin_;//query vectors
    Basin*          subbasin_;//subject vectors
    Pslvector       quevals_;//precalculated values for query
    Pslvector       subvals_;//precalculated values for subject

    static bool     mtxvart_;//use matrix variate t-distribution
    static bool     hdpbprobs_;//background probabilities are multivariate t-dist.
    static const int    s_defbasize;//default basin size
};


// INLINES ...

// -------------------------------------------------------------------------
// GetScore: get score at the positions
// -------------------------------------------------------------------------

inline
double HDP0ScoreMatrix::GetScore( int m, int n ) const
{
    const double    scoreX = 0.0; //-1.0;

//     if( GetMaskedToConsider( m, n ))
//         return scoreX;

//     if( GetMaskedToConsider( m, n ) || GetMaskedToIgnore( m, n ))
//         return GetScoreBase( m, n ) * GetMaskscalePercents();

    return GetScoreBase( m, n );
}

// -------------------------------------------------------------------------
// ThicknessConstraintsBroken: returns false if thickness constraints at the
//     positions indicated are broken
// -------------------------------------------------------------------------

inline
bool HDP0ScoreMatrix::ThicknessConstraintsMet( int m, int n ) const
{
#ifdef __DEBUG__
    if( GetSubjectSize() <= m || m < 0 )
        throw myruntime_error( mystring( "HDP0ScoreMatrix: Memory access error." ));

    if( GetQuerySize() <= n || n < 0 )
        throw myruntime_error( mystring( "HDP0ScoreMatrix: Memory access error." ));

    if( !logo_sec_.GetEffNoSequences() || !logo_fst_.GetEffNoSequences())
        throw myruntime_error( mystring( "HDP0ScoreMatrix: Profile thickness found to be wrong." ));
#endif

    if(           logo_sec_.GetMIDExpNoObservationsAt( m, PS_M ) < ( size_t )GetThicknessNumber() ||
                  logo_fst_.GetMIDExpNoObservationsAt( n, PS_M ) < ( size_t )GetThicknessNumber() ||
        ( double )logo_sec_.GetMIDExpNoObservationsAt( m, PS_M ) / logo_sec_.GetEffNoSequences() < GetThicknessPercents() ||
        ( double )logo_fst_.GetMIDExpNoObservationsAt( n, PS_M ) / logo_fst_.GetEffNoSequences() < GetThicknessPercents())
            //if number of residues within a column is less than a number specified or
            //residues ratio in a column is less than a percentage specified
            return false;

    return true;
}


// -------------------------------------------------------------------------
// InitHDPbase: initialize hdpbase structure
//
inline
void HDP0ScoreMatrix::InitHDPbase()
{
    DestroyHDPbase();
    hdpbase_ = new HDPbase();
    if( hdpbase_ == NULL )
        throw myruntime_error("HDP0ScoreMatrix: Not enough memory.");
}

// -------------------------------------------------------------------------
// InitQueryBasin: initialize query basin of vectors
//
inline
void HDP0ScoreMatrix::InitQueryBasin()
{
    DestroyQueryBasin();
    int size = GetQuerySize();
    if( size <= 0 )
        size = s_defbasize;
    //extra one position for background prob. vector
    quebasin_ = new Basin( size + 1 );
    if( quebasin_ == NULL )
        throw myruntime_error("HDP0ScoreMatrix: Not enough memory.");
    quebasin_->SetDestroy( true );
}

// -------------------------------------------------------------------------
// InitSbjctBasin: initialize subject basin of vectors
//
inline
void HDP0ScoreMatrix::InitSbjctBasin()
{
    DestroySbjctBasin();
    int size = GetSubjectSize();
    if( size <= 0 )
        size = s_defbasize;
    //extra one position for background prob. vector
    subbasin_ = new Basin( size + 1 );
    if( subbasin_ == NULL )
        throw myruntime_error("HDP0ScoreMatrix: Not enough memory.");
    subbasin_->SetDestroy( true );
}


#endif//__HDP0ScoreMatrix__
