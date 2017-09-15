/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __ProfileMatrix__
#define __ProfileMatrix__

#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "mystring.h"
#include "myexcept.h"
#include "pslvector.h"

#include "AbstractScoreMatrix.h"


class ProfileMatrix: public AbstractScoreMatrix
{
public:
    ProfileMatrix(
            const double ( *pssmscores )[NUMALPH],
            const char*     resids,
            int             length,
            Configuration   config[NoSchemes],
            TScaling        a_scaling = AutoScalling,
            TMask           c_masking = Unmasked
    );
    virtual ~ProfileMatrix();

    virtual const char* GetMethodName() const { return "Position-specific"; }

    const double     ( *GetVectorAt( int n ) const )[NUMALPH];

    virtual void        ComputeProfileScoringMatrix( bool = false );
                                                                //compute probabilities to observe scores at each position
    virtual void        ComputeScoreProbabilities( AttributableScores* );
                                                        //compute positional probabilities of scores
    virtual void        ComputePositionalScoreProbs( AttributableScores* );

    virtual void        OptimizeTargetFrequencies();

                                                                //print statistical information
    virtual void        PrintParameterTable( TPrintFunction, void* vpn ) const;
    void                TestPrintScoringMatrix( FILE* fp ) const;
    virtual void        PrintScoringMatrix( FILE* );//print the scoring system
    virtual void        PrintFinal( TPrintFunction, void* vpn ) const;

protected:
    explicit ProfileMatrix();

    const char*         GetResidues() const { return residues; }
    char                GetResidueAt( int n ) const;
                                                                //fill matrix with scores
    void                FillMatrix( const double ( *pssmscores )[NUMALPH] );

    //helper routines for optimization of target frequencies
    void                IncorporateTargetFrequencies( const Pslvector& );

private:
    const char*         residues;                               //vector of residues
    static double       position[NUMALPH];

    Pslvector           scores_;    //scores for optimization of target frequencies
    Pslvector           rprobs_;    //row background probabilities of the score system
    Pslvector           cprobs_;    //column background probabilities of the score system
};

void ComputedSubMatrixWithParams();

// INLINES ...

// -------------------------------------------------------------------------
// Return residue at the position

inline
char ProfileMatrix::GetResidueAt( int n ) const
{
#ifdef __DEBUG__
    if( n < GetQuerySize())
#endif
        return residues[n];

    throw myruntime_error( mystring( "ProfileMatrix: Memory access error." ));
}

// -------------------------------------------------------------------------
// ComputeProfileScoringMatrix should not be called from the class

inline
void ProfileMatrix::ComputeProfileScoringMatrix( bool )
{
    throw myruntime_error( mystring( "ProfileMatrix: ComputeProfileScoringMatrix should not be called." ));
}


#endif//__ProfileMatrix__
