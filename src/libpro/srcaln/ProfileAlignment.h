/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __ProfileAlignment__
#define __ProfileAlignment__

#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "rc.h"
#include "mystring.h"
#include "myexcept.h"
#include "stat.h"
#include "libpro/srcsco/ScoringMatrix.h"
#include "libpro/srcsco/UniversalScoreMatrix.h"
#include "libpro/srcpro/GapScheme.h"
#include "libpro/srcpro/DistributionMatrix.h"
#include "libpro/srcpro/DBProfileProbs.h"

//global variable externs
extern int      gnMAXMARGIN;//maximum margin length at both ends of profile
extern double   gdMARGINPRC;//percentage of profile length, set as margin length

// typedefs -
// pairwise alignment score
typedef double  TPAScore;

// enumeration of indices
enum ProInd {
    nQuery, //query index
    nSbjct, //subject index
    noInds  //number of indices
};

// _________________________________________________________________________
// Class ProfileAlignment
//
// This class implements an algorithm for alignment of two profiles
//

class ProfileAlignment
{
public:
    enum TAlgorithm {
            Optimizational,
            Hybrid,
            noAlgorithms
    };
    enum State { //states for dynamic programming matrices
        stateMM, //match-to-match state
        stateMI, //match-to-insert
        stateIM, //insert-to-match
        stateDG, //delete-to-no-state(gap)
        stateGD, //no-state-to-delete
        noStates //number of states
    };
    enum Direction {//direction codes for back-tracing
        dNo =  0,//no valid direction
        dMM =  1,//00001
        dMI =  2,//00010
        dIM =  4,//00100
        dDG =  8,//01000
        dGD = 16,//10000

        dMM_MI = dMM | dMI, dMM_IM = dMM | dIM, dMM_DG = dMM | dDG, dMM_GD = dMM | dGD,
        dMI_MM = dMM_MI,    dIM_MM = dMM_IM,    dDG_MM = dMM_DG,    dGD_MM = dMM_GD,

        dMI_IM = dMI | dIM, dMI_DG = dMI | dDG, dMI_GD = dMI | dGD,
        dIM_MI = dMI_IM,    dDG_MI = dMI_DG,    dGD_MI = dMI_GD,

        dIM_DG = dIM | dDG, dIM_GD = dIM | dGD,
        dDG_IM = dIM_DG,    dGD_IM = dIM_GD,

        dDG_GD = dDG | dGD, dGD_DG = dDG_GD,

        dAll = dMM | dMI | dIM | dDG | dGD,
        noDirections
    };

    ProfileAlignment(
            const FrequencyMatrix& freq_fst, const LogOddsMatrix& logo_fst, const GapScheme& gaps_fst,
            const FrequencyMatrix& freq_sec, const LogOddsMatrix& logo_sec, const GapScheme& gaps_sec,
            const AbstractScoreMatrix*  usc_system,
            bool ungapped = false
    );

    virtual ~ProfileAlignment();

    virtual void                Run( size_t );
    void                        MARealign();
    void                        AdjustScore( double value );

    int                         GetQuerySize() const            { return querySize_; }
    int                         GetSubjectSize() const          { return subjectSize_; }

    double                      GetScore() const                { return finscore_; }
    int                         GetAlnSteps() const             { return alnsteps_; }
    int                         GetNoMatched() const            { return mtcsteps_; }
    int                         GetNoIdentities() const         { return idnsteps_; }
    int                         GetNoPositives() const          { return pststeps_; }
    double                      GetBitScore() const             { return bitscore; }
    double                      GetReferenceExpectation() const { return ref_expectation; }
    double                      GetExpectPerAlignment() const   { return expectperaln_; }
    double                      GetRawExpectation() const       { return raw_expect_; }
    double                      GetExpectation() const          { return expectation; }
    double                      GetPvalue() const               { return p_value; }

    int                         GetAlnLength() const;
    void                        RelEntropyAdjustment();
    void                        CorrelatedScoreAdjustment();

    void        CopyFrequencyColumn( const FrequencyMatrix& freq_from, FrequencyMatrix& freq_to, int ind ) const;
    void        CopyProfileColumn( const LogOddsMatrix& pssm_from, LogOddsMatrix& pssm_to, int ind ) const;
    void        ProcessUngappedSegment(
                    const FrequencyMatrix& f1, const LogOddsMatrix& l1,
                    const FrequencyMatrix& f2, const LogOddsMatrix& l2 );

    void                        Output( const char* filename );             //output alignment information

    void                        Print( char* sp, bool showpars = true );    //print alignment with statistical significance
    void                        Print( FILE* fp, bool showpars = true );    //print alignment with statistical significance
                                                                            //universal print method
    void                        Print( TPrintFunction, void* vpn, bool = true );

    void                        OutputScoringMatrix( const char* = NULL );  //output scoring system

    size_t                      GetMaxReqSizeForAlignment() const; //size to hold alignment information

    static void                 SetInformationThreshold( double value )     { s_information_thrsh = value; }
    static void                 SetSEGdistanceThreshold( double value )     { s_segdistance_thrsh = value; }

    static double               GetInformationThreshold()                   { return s_information_thrsh; }
    static double               GetSEGdistanceThreshold()                   { return s_segdistance_thrsh; }

    void                        SetProProbs( const DBProfileProbs* value )  { proprobs_ = value; }
    const DBProfileProbs*       GetProProbs() const { return proprobs_; }

protected:
    explicit ProfileAlignment();

    void                        Initialize();
    void                        Destroy();

    virtual void                Clear();
    virtual void                ClearF();
    virtual void                ClearPointer();
    virtual void                ClearPath();

    virtual void                AlignProfiles();
    virtual void                PostProcess();
    virtual void                SetFinalScore( double value );
    virtual void                MakeAlignmentPath();
    void                        ComputeStatistics();
    void                        ComputeStatisticsTEST( size_t );
    static double               ProbProduct( double term1, double term2, bool apply );
    static double               LogOfProduct( double term1, double term2, bool apply );

    double                      AutocorrScoreCmb( const AbstractScoreMatrix*, int sbjctpos, int querypos );
    double                      AutocorrScore( const AbstractScoreMatrix*, int sbjctpos, int querypos );
    double                      AutocorrScore( const AbstractScoreMatrix*, int sbjctpos, int querypos, 
                                               int, int, int, int, int winsz );

    State                       GetState( int dir ) const;

    void                        SetAlnSteps( int steps )        { alnsteps_ = steps; }
    void                        SetNoMatched( int steps )       { mtcsteps_ = steps; }
    void                        SetNoIdentities( int steps )    { idnsteps_ = steps; }
    void                        SetNoPositives( int steps )     { pststeps_ = steps; }

    double                      GetAlnScore() const             { return alnscore_; }
    void                        SetAlnScore( double value )     { alnscore_ = value; }
    void                        SubtractScore( double value );

    void                        SetBitScore( double value )             { bitscore = value; }
    void                        SetReferenceExpectation( double value ) { ref_expectation = value; }
    void                        SetExpectPerAlignment( double value )   { expectperaln_ = value; }
    void                        SetRawExpectation( double value )       { raw_expect_ = value; }
    void                        SetExpectation( double value )          { expectation = value; }
    void                        SetPvalue( double value )               { p_value = value; }

    const AbstractScoreMatrix*  GetScoreMatrix() const;

    const AbstractScoreMatrix*  GetScoreSystem() const { return scoreSystem; }
    bool                        IsUngapped() const { return ungapped_;  }

    static  double              ComputePvalue( double expect );

protected:
    TPAScore    ( **F )[noStates];              //dynamic programming matrices
    int         ( **pointer )[noStates];        //backtracing pointers
    int         ( *path )[noInds];      //alignment path containing indices of two profiles in each position

    double      finscore_;                      //final alignment score
    double      alnscore_;                      //alignment score
    int         alnsteps_;                      //length of alignment
    int         mtcsteps_;                      //matched positions
    int         idnsteps_;                      //number of identities
    int         pststeps_;                      //number of positives (positive identities included)

    int                         querySize_;     //length of query sequence (profile)
    int                         subjectSize_;   //length of subject sequence (profile)

    const FrequencyMatrix&      freq_fst_;      //reference to the first weighted frequency matrix
    const LogOddsMatrix&        logo_fst_;      //reference to the first log-odds matrix
    const GapScheme&            gaps_fst_;      //reference to the first gap cost vector

    const FrequencyMatrix&      freq_sec_;      //reference to the second weighted frequency matrix
    const LogOddsMatrix&        logo_sec_;      //reference to the second log-odds matrix
    const GapScheme&            gaps_sec_;      //reference to the second gap cost vector

    const DBProfileProbs*       proprobs_;      //profile pair probabilities
    const AbstractScoreMatrix*  scoreSystem;    //score system used to align profiles

private:
    StatModel   model;                      //statistical model
    double      bitscore;                   // bit score of alignment
    double      ref_expectation;            // reference expectation for the case expected score per position is positive
    double      expectation;                // E-value
    double      expectperaln_;              // expect per alignment
    double      raw_expect_;                // raw E-value without any experimental adjustment
    double      p_value;                    // P-value

    bool        ungapped_;                  //align profiles without using gaps

    static double   s_information_thrsh;    //information content threshold
    static double   s_segdistance_thrsh;    //SEG distance threshold

};


////////////////////////////////////////////////////////////////////////////
// INLINES
//
inline ProfileAlignment::State ProfileAlignment::GetState( int dir ) const
{
    State   state = noStates;
    if( dir & dMM ) state = stateMM;
    else if( dir & dMI ) state = stateMI;
    else if( dir & dIM ) state = stateIM;
    else if( dir & dDG ) state = stateDG;
    else if( dir & dGD ) state = stateGD;
    return state;
}

// -------------------------------------------------------------------------
// GetScoreMatrix: obtain one of the two possible scoring matrices
// -------------------------------------------------------------------------

inline
const AbstractScoreMatrix* ProfileAlignment::GetScoreMatrix() const
{
    if( !scoreSystem )
        throw myruntime_error( "ProfileAlignment: No score matrix." );
    return scoreSystem;
}

// -------------------------------------------------------------------------
// SetFinalScore: scale and save final alignment score
// 
inline
void ProfileAlignment::SetFinalScore( double value )
{
    finscore_ = GetScoreMatrix()->GetFinalScore( value );
}

// -------------------------------------------------------------------------
// SubtractScore: subtract score from the final sum
// -------------------------------------------------------------------------

inline
void ProfileAlignment::SubtractScore( double value )
{
    if( value <= 0.0 )
        return;
    if( alnscore_ - value < 0.0 )
        return;
    alnscore_ -= value;
}

// -------------------------------------------------------------------------
// ComputePvalue: compute p-value given expect value
//
inline
double ProfileAlignment::ComputePvalue( double expect )
{
    if( expect < 0.0 )
        throw myruntime_error( "ProfileAlignment: Expect value is negative." );
    return 1.0 - exp( -expect );
}

#endif//__ProfileAlignment__
