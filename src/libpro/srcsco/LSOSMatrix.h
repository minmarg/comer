/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __LSOSMatrix__
#define __LSOSMatrix__

#include "debug.h"
#include "types.h"
#include "compdef.h"

#include "mystring.h"
#include "myexcept.h"
#include "ext/pslvector.h"
#include "libpro/srcpro/DistributionMatrix.h"
#include "AbstractScoreMatrix.h"


class LSOSMatrix: public AbstractScoreMatrix
{
public:
    LSOSMatrix(
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
    LSOSMatrix(
            const FrequencyMatrix&  freq_fst, const LogOddsMatrix& logo_fst, 
            const FrequencyMatrix&  freq_sec, const LogOddsMatrix& logo_sec,
            const HDPbase*  hdpbase,
            const HDPbase*  hdpctbase,
            double          infrm_threshold,
            int             thick_number,
            double          thick_percnt,
            double          mask_percnt,
            Configuration   config[NoSchemes],
            TBehaviour      s_behaviour,
            TScaling        a_scaling,
            TMask           c_masking
    );
    explicit LSOSMatrix();

    virtual ~LSOSMatrix();

    virtual double      GetScore( int m, int n ) const;
    virtual const char* GetMethodName() const   { return "Profile-specific LSO"; }

    virtual void        ComputeProfileScoringMatrix( bool final = false );
    virtual double      ComputeExpectation( double, double* = NULL, double* = NULL, double* = NULL, double* = NULL ) const;

    int                 GetThicknessNumber() const      { return thickness_number; }
    double              GetThicknessPercents() const    { return thickness_percnt; }
    double              GetMaskscalePercents() const    { return maskscale_percnt; }

    const HDPbase*      GetHDPBase() const { return hdpbase_; }
    const HDPbase*      GetHDPctBase() const { return hdpctbase_; }

    virtual void        OptimizeTargetFrequencies();
    virtual void        ComputeScoreProbabilities( AttributableScores* );
    virtual void        ComputePositionalScoreProbs( AttributableScores* );

    virtual void        PrintParameterTable( TPrintFunction, void* vpn ) const;
    virtual void        PrintFinal( TPrintFunction, void* vpn ) const;
    virtual void        PrintScoringMatrix( FILE* );

    double                  ComputeScore( int m, int n, bool final, double* = NULL ) const;
    double                  ComputeSSSScore( int m, int n, double* modscore ) const;
    double                  RegulateHDPSSSS( int m, int n, double* modscore ) const;
    double                  ComputeCVS2ScoreHelper( int m, int n, bool sson ) const;
    double                  ComputeCVS2Score( int m, int n, bool sson, double* modscore ) const;
    double                  ComputeHDPctScore( int m, int n, bool sson, double* modscore ) const;
    double                  ComputeHDPctDrvdScore( int m, int n, bool sson, double* modscore ) const;
    double                  ComputeHDPScore( int m, int n, bool sson, double* modscore ) const;
    double                  ComputeHDPDrvdScore( int m, int n, bool sson, double* modscore ) const;
    double                  ComputeHDPStatScore( int m, int n, bool sson, double* modscore ) const;

protected:
    LSOSMatrix(
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

    void                    PrivateInit();                          //private initialization method

    const FrequencyMatrix&  GetQueryFreq() const    { return freq_fst_; }
    const LogOddsMatrix&    GetQueryLogo() const    { return logo_fst_; }

    const FrequencyMatrix&  GetSbjctFreq() const    { return freq_sec_; }
    const LogOddsMatrix&    GetSbjctLogo() const    { return logo_sec_; }

    double                  GetMixBPAt( int r ) const;
    void                    SetMixBPAt( int r, double value );

    bool                    ThicknessConstraintsMet( int m, int n ) const;

    //helper routines for optimization of target frequencies
    void                    IncorporateTargetFrequencies( const Pslvector& );

    void                    DestroyCVS2scores();

private:
    const FrequencyMatrix&  freq_fst_;  //reference to the first weighted frequency matrix
    const LogOddsMatrix&    logo_fst_;  //reference to the first log-odds matrix

    const FrequencyMatrix&  freq_sec_;  //reference to the second weighted frequency matrix
    const LogOddsMatrix&    logo_sec_;  //reference to the second log-odds matrix

    const int               thickness_number;   //thickness in number of sequences a position must have, otherwise it will be ignored
    const double            thickness_percnt;   //thickness in percentage
    const double            maskscale_percnt;   //scaling of masked positions in percentage

    Pslvector               scores_;    //scores for optimization of target frequencies
    Pslvector               rprobs_;    //row background probabilities of the score system
    Pslvector               cprobs_;    //column background probabilities of the score system

    const double            (*pqry_)[NUMAA];//query probabilities
    double                  (*psub_)[NUMAA];//subject probabilities
    double                  mixbp_[NUMAA];//mix of profile background probabilities

    const HDPbase*          hdpbase_;
    const HDPbase*          hdpctbase_;
    double*                 qryhdpbuf_;

    double**                cvs2scores_;//translated context vector scores
};


// INLINES ...

// -------------------------------------------------------------------------
// GetScore: get score at the position
//
inline
double LSOSMatrix::GetScore( int m, int n ) const
{
    return GetScoreBase( m, n );
}

// -------------------------------------------------------------------------
// ThicknessConstraintsBroken: ret false if thickness constraints aren't met
//
inline
bool LSOSMatrix::ThicknessConstraintsMet( int m, int n ) const
{
#ifdef __DEBUG__
    if( GetSubjectSize() <= m || m < 0 )
        throw myruntime_error( mystring( "LSOSMatrix: Memory access error." ));

    if( GetQuerySize() <= n || n < 0 )
        throw myruntime_error( mystring( "LSOSMatrix: Memory access error." ));

    if( !logo_sec_.GetEffNoSequences() || !logo_fst_.GetEffNoSequences())
        throw myruntime_error( mystring( "LSOSMatrix: Profile thickness found to be wrong." ));
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
// GetMixBPAt: get mixed profile background probability at the position
//
inline
double LSOSMatrix::GetMixBPAt( int r ) const
{
#ifdef __DEBUG__
    if( r < 0 || NUMAA <= r )
        throw myruntime_error("LSOSMatrix: GetMixBPAt: Memory access error.");
#endif
    return mixbp_[r];
}

// -------------------------------------------------------------------------
// SetMixBPAt: set mixed profile background probability at the position
//
inline
void LSOSMatrix::SetMixBPAt( int r, double value )
{
#ifdef __DEBUG__
    if( r < 0 || NUMAA <= r )
        throw myruntime_error("LSOSMatrix: GetMixBPAt: Memory access error.");
#endif
    mixbp_[r] = value;
}

// -------------------------------------------------------------------------
// DestroyCVS2scores: Destroy translated context vector scores
//
inline
void LSOSMatrix::DestroyCVS2scores()
{
    if( cvs2scores_ ) {
        for( int m = 0; m < GetSubjectSize(); m++ )
            if( cvs2scores_[m] )
                free( cvs2scores_[m] );
        free( cvs2scores_ );
        cvs2scores_ = NULL;
    }
}

#endif//__LSOSMatrix__
