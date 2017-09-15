/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __ScoringMatrix__
#define __ScoringMatrix__

#include "debug.h"
#include "types.h"
#include "compdef.h"

#include "mystring.h"
#include "myexcept.h"
#include "ext/pslvector.h"
#include "libpro/srcpro/DistributionMatrix.h"
#include "AbstractScoreMatrix.h"


class ScoringMatrix: public AbstractScoreMatrix
{
public:
    ScoringMatrix(
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
    explicit ScoringMatrix();

    virtual ~ScoringMatrix();

    virtual double      GetScore( int m, int n ) const;             //get score at the positions

    virtual const char* GetMethodName() const   { return "Profile-specific vector-based"; }

    virtual void        ComputeProfileScoringMatrix( bool final = false );
                                                                    //compute e-value
    virtual double      ComputeExpectation( double, double* = NULL, double* = NULL, double* = NULL, double* = NULL ) const;

    int                 GetThicknessNumber() const      { return thickness_number; }
    double              GetThicknessPercents() const    { return thickness_percnt; }
    double              GetMaskscalePercents() const    { return maskscale_percnt; }

    virtual void        OptimizeTargetFrequencies();

                                                            //compute probabilities to observe scores at each position
    virtual void        ComputeScoreProbabilities( AttributableScores* );
                                                            //compute positional probabilities of scores
    virtual void        ComputePositionalScoreProbs( AttributableScores* );

    virtual void        PrintParameterTable( TPrintFunction, void* vpn ) const;
    virtual void        PrintFinal( TPrintFunction, void* vpn ) const;

    virtual void        PrintScoringMatrix( FILE* );


protected:
    ScoringMatrix(
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
    double                  ComputeScore( int m, int n, bool final ) const;

    const FrequencyMatrix&  GetQueryFreq() const    { return freq_fst_; }
    const LogOddsMatrix&    GetQueryLogo() const    { return logo_fst_; }

    const FrequencyMatrix&  GetSbjctFreq() const    { return freq_sec_; }
    const LogOddsMatrix&    GetSbjctLogo() const    { return logo_sec_; }

    bool                    ThicknessConstraintsMet( int m, int n ) const;

    void                    PreliminaryVerification();

    //helper routines for optimization of target frequencies
    void                    IncorporateTargetFrequencies( const Pslvector& );

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
};


// INLINES ...

// -------------------------------------------------------------------------
// GetScore: get score at the position
//
inline
double ScoringMatrix::GetScore( int m, int n ) const
{
    const double    scoreX = 0.0; //-1.0;

//     if( GetMaskedToConsider( m, n ))
//         return scoreX;

//     if( GetMaskedToConsider( m, n ) || GetMaskedToIgnore( m, n ))
//         return GetScoreBase( m, n ) * GetMaskscalePercents();

    return GetScoreBase( m, n );
}

// -------------------------------------------------------------------------
// ThicknessConstraintsBroken: ret false if thickness constraints aren't met
//
inline
bool ScoringMatrix::ThicknessConstraintsMet( int m, int n ) const
{
#ifdef __DEBUG__
    if( GetSubjectSize() <= m || m < 0 )
        throw myruntime_error( mystring( "ScoringMatrix: Memory access error." ));

    if( GetQuerySize() <= n || n < 0 )
        throw myruntime_error( mystring( "ScoringMatrix: Memory access error." ));

    if( !logo_sec_.GetEffNoSequences() || !logo_fst_.GetEffNoSequences())
        throw myruntime_error( mystring( "ScoringMatrix: Profile thickness found to be wrong." ));
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


#endif//__ScoringMatrix__
