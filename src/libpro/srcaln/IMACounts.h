/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __IMACounts__
#define __IMACounts__

#include "compdef.h"
#include "debug.h"
#include "rc.h"
#include "data.h"

#include "mystring.h"
#include "myexcept.h"


// _________________________________________________________________________
// CLASS IMACounts
//
class IMACounts
{
public:
    IMACounts();
    ~IMACounts();

    void    Read( const char* filename );
    void    Write( FILE* fp );

    void    ReadScores( const char* filename );
    void    WriteScores( FILE* fp );

    void    SumUp( const IMACounts& );

    void    ComputeScores();
    void    ComputeParameters();

    double  GetResidCount( size_t ) const;
    void    SetResidCount( size_t, double value );
    void    IncResidCountBy( size_t, double value );

    double  GetPairCount( size_t, size_t ) const;
    void    SetPairCount( size_t, size_t, double value );
    void    IncPairCountBy( size_t, size_t, double value );

    static double   GetScoreScale()                 { return scorescale_; }
    static void     SetScoreScale( double value )   { scorescale_ = value; }

    static bool     GetAutoScoreScale()             { return autoscale_; }
    static void     SetAutoScoreScale( bool value ) { autoscale_ = value; }

    double  GetMargProbability( size_t ) const;
    double  GetBackProbability( size_t ) const;

    double  GetSubstScore( size_t, size_t ) const;
    double  GetOddsRatio( size_t, size_t ) const;

    double  GetStatLambda() const               { return statlambda_; }
    double  GetStatK() const                    { return statK_; }
    double  GetStatH() const                    { return statH_; }

protected:
    size_t  GetEffectiveNoResids() const    { return NUMAA; }

    void    ComputeAbstractScores();
    void    DetermineScale() const;
    void    DetermineScaleThroughScores() const;

    void    ReadFooter( FILE*, char* buffer, size_t bufsize, size_t* readlen );
    void    SkipComments( FILE*, char* buffer, size_t bufsize, size_t* readlen );
    void    ReadCounts( const char* readfrom, size_t readlen, double* membuffer, size_t no_cnts, bool allownegs = false );

    void    ReadStatisParams( const char* readfrom, size_t readlen );
    void    ReadAssocValue( const char** p, mystring name, double* value );

    bool    GetScoresRead() const               { return scoresread_; }
    void    SetScoresRead( bool value )         { scoresread_ = value; }

    double*     GetResidCountsAddr()            { return obsresidues_; }
    double*     GetPairCountsAddrAt( size_t res );

    void    ResetResidCounts();
    void    ResetPairCounts();
    void    ResetTargetFrequencies();
    void    ResetMargProbabilities();
    void    ResetBackProbabilities();
    void    ResetSubstScores();
    void    ResetOddsRatios();

    static double   GetMinScore()               { return minscore_; }

    double  GetTotalResidCount() const          { return totresidcount_; }
    void    SetTotalResidCount( double value )  { totresidcount_ = value; }

    double  GetTotalPairCount() const           { return totpaircount_; }
    void    SetTotalPairCount( double value )   { totpaircount_ = value; }

    double  GetTargetFrequency( size_t, size_t ) const;
    void    SetTargetFrequency( size_t, size_t, double value );
    double* GetTargetFreqsAddrAt( size_t );

    const double*   GetMargProbabilities() const { return marginprobs_; }
    void    SetMargProbability( size_t, double value );
    void    IncMargProbabilityBy( size_t, double value );
    double* GetMargProbsAddr()                  { return marginprobs_; }

    const double*   GetBackProbabilities() const { return backgrprobs_; }
    void    SetBackProbability( size_t, double value );
    double* GetBackProbsAddr()                  { return backgrprobs_; }

    void    SetSubstScore( size_t, size_t, double value );
    double* GetSubstScoresAddrAt( size_t );

    void    SetOddsRatio( size_t, size_t, double value );
    double* GetOddsRatiosAddrAt( size_t );

    double  GetEntropy() const                  { return scorentropy_; }
    void    SetEntropy( double value )          { scorentropy_ = value; }

    double  GetExpected() const                 { return expectedscore_; }
    void    SetExpected( double value )         { expectedscore_ = value; }

    double  GetMinSubstScore() const            { return minsubscore_; }
    void    SetMinSubstScore( double value )    { minsubscore_ = value; }

    double  GetMaxSubstScore() const            { return maxsubscore_; }
    void    SetMaxSubstScore( double value )    { maxsubscore_ = value; }


    mystring    GetLambdaName()                 { return mystring( "Lambda" ); }
    mystring    GetKName()                      { return mystring( "K" ); }
    mystring    GetHName()                      { return mystring( "H" ); }

    void    SetStatLambda( double value )       { statlambda_ = value; }
    double* GetLambdaAddr()                     { return &statlambda_; }

    void    SetStatK( double value )            { statK_ = value; }
    double* GetKAddr()                          { return &statK_; }

    void    SetStatH( double value )            { statH_ = value; }
    double* GetHAddr()                          { return &statH_; }

private:
    double  obsresidues_[NUMALPH];          //observed amino acid counts
    double  obscounts_[NUMALPH][NUMALPH];   //observed amino acid pair counts
    double  totresidcount_;                 //total residue count
    double  totpaircount_;                  //total substitution count
    bool    scoresread_;                    //flag of whether scores are read from file
    double  targetfreqs_[NUMALPH][NUMALPH]; //target frequencies derived from counts
    double  marginprobs_[NUMALPH];          //marginal probabilities
    double  backgrprobs_[NUMALPH];          //bacground probabilities
    double  subscores_[NUMALPH][NUMALPH];   //substitution scores derived
    double  oddratios_[NUMALPH][NUMALPH];   //odds ratios matrix (for reading)
    double  scorentropy_;                   //entropy of scores
    double  expectedscore_;                 //expected score
    double  minsubscore_;                   //minimum substitution score
    double  maxsubscore_;                   //maximum substitution score
    double  statlambda_;                    //statistical paramater lambda
    double  statK_;                         //statistical parameter K
    double  statH_;                         //statistical parameter H

    static const double minscore_;          //minimum score
    static double       scorescale_;        //scale factor of scores
    static bool         autoscale_;         //automatically computed scale factor
};

////////////////////////////////////////////////////////////////////////////
// Class IMACounts inlines
//
// -------------------------------------------------------------------------
// GetResidCount: returns observed amino acid count

inline
double IMACounts::GetResidCount( size_t res ) const
{
#ifdef __DEBUG__
    if( NUMALPH <= res )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif
    return obsresidues_[res];
}

// SetResidCount: sets observed amino acid count
//
inline
void IMACounts::SetResidCount( size_t res, double value )
{
#ifdef __DEBUG__
    if( NUMALPH <= res )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif
    obsresidues_[res] = value;
}

// IncResidCountBy: increases observed amino acid count by a value given
//
inline
void IMACounts::IncResidCountBy( size_t res, double value )
{
#ifdef __DEBUG__
    if( NUMALPH <= res )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif
    obsresidues_[res] += value;
}

// -------------------------------------------------------------------------
// GetPairCount: returns observed amino acid pair count

inline
double IMACounts::GetPairCount( size_t res1, size_t res2 ) const
{
#ifdef __DEBUG__
    if( NUMALPH <= res1 || NUMALPH <= res2 )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif
    return obscounts_[res1][res2];
}

// SetPairCount: sets observed amino acid pair count
//
inline
void IMACounts::SetPairCount( size_t res1, size_t res2, double value )
{
#ifdef __DEBUG__
    if( NUMALPH <= res1 || NUMALPH <= res2 )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif
    obscounts_[res1][res2] = value;
}

// IncPairCountBy: increases observed amino acid pair count by a value given
//
inline
void IMACounts::IncPairCountBy( size_t res1, size_t res2, double value )
{
#ifdef __DEBUG__
    if( NUMALPH <= res1 || NUMALPH <= res2 )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif
    obscounts_[res1][res2] += value;
}

// -------------------------------------------------------------------------
// GetPairCountsAddrAt: returns address of pair count vector at the given
//     position
inline
double* IMACounts::GetPairCountsAddrAt( size_t res )
{
#ifdef __DEBUG__
    if( NUMALPH <= res )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif
    return obscounts_[res];
}

// -------------------------------------------------------------------------
// GetTargetFrequency: returns target frequency

inline
double IMACounts::GetTargetFrequency( size_t res1, size_t res2 ) const
{
#ifdef __DEBUG__
    if( NUMALPH <= res1 || NUMALPH <= res2 )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif
    return targetfreqs_[res1][res2];
}

// SetTargetFrequency: sets target frequency
//
inline
void IMACounts::SetTargetFrequency( size_t res1, size_t res2, double value )
{
#ifdef __DEBUG__
    if( NUMALPH <= res1 || NUMALPH <= res2 )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif
    targetfreqs_[res1][res2] = value;
}

// GetTargetFreqsAddrAt: returns address of target frequencies vector at the
//     given position
inline
double* IMACounts::GetTargetFreqsAddrAt( size_t res )
{
#ifdef __DEBUG__
    if( NUMALPH <= res )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif
    return targetfreqs_[res];
}

// -------------------------------------------------------------------------
// GetMargProbability: returns marginal residue probability

inline
double IMACounts::GetMargProbability( size_t res ) const
{
#ifdef __DEBUG__
    if( NUMALPH <= res )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif
    return marginprobs_[res];
}

// SetMargProbability: sets marginal probability
//
inline
void IMACounts::SetMargProbability( size_t res, double value )
{
#ifdef __DEBUG__
    if( NUMALPH <= res )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif
    marginprobs_[res] = value;
}

// IncMargProbabilityBy: increases marginal probability by value given
//
inline
void IMACounts::IncMargProbabilityBy( size_t res, double value )
{
#ifdef __DEBUG__
    if( NUMALPH <= res )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif
    marginprobs_[res] += value;
}

// -------------------------------------------------------------------------
// GetBackProbability: returns background residue probability

inline
double IMACounts::GetBackProbability( size_t res ) const
{
#ifdef __DEBUG__
    if( NUMALPH <= res )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif
    return backgrprobs_[res];
}

// SetBackProbability: sets background probability
//
inline
void IMACounts::SetBackProbability( size_t res, double value )
{
#ifdef __DEBUG__
    if( NUMALPH <= res )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif
    backgrprobs_[res] = value;
}

// -------------------------------------------------------------------------
// GetSubstScore: returns substitution score

inline
double IMACounts::GetSubstScore( size_t res1, size_t res2 ) const
{
#ifdef __DEBUG__
    if( NUMALPH <= res1 || NUMALPH <= res2 )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif
    return subscores_[res1][res2];
}

// SetSubstScore: sets substitution score
//
inline
void IMACounts::SetSubstScore( size_t res1, size_t res2, double value )
{
#ifdef __DEBUG__
    if( NUMALPH <= res1 || NUMALPH <= res2 )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif
    subscores_[res1][res2] = value;
}

// GetSubstScoresAddrAt: returns address of scores vector at the given
//     position
inline
double* IMACounts::GetSubstScoresAddrAt( size_t res )
{
#ifdef __DEBUG__
    if( NUMALPH <= res )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif
    return subscores_[res];
}

// -------------------------------------------------------------------------
// GetOddsRatio: returns odds ratio value

inline
double IMACounts::GetOddsRatio( size_t res1, size_t res2 ) const
{
#ifdef __DEBUG__
    if( NUMALPH <= res1 || NUMALPH <= res2 )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif
    return oddratios_[res1][res2];
}

// SetOddsRatio: sets odds ratio value
//
inline
void IMACounts::SetOddsRatio( size_t res1, size_t res2, double value )
{
#ifdef __DEBUG__
    if( NUMALPH <= res1 || NUMALPH <= res2 )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif
    oddratios_[res1][res2] = value;
}

// GetOddsRatiosAddrAt: returns address of odds ratios vector at the given
//     position
inline
double* IMACounts::GetOddsRatiosAddrAt( size_t res )
{
#ifdef __DEBUG__
    if( NUMALPH <= res )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif
    return oddratios_[res];
}


#endif//__IMACounts__

