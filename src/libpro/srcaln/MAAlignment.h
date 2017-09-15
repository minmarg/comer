/***************************************************************************
 *   Copyright (C) 2012 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __MAAlignment__
#define __MAAlignment__

#include "debug.h"
#include "compdef.h"
#include "rc.h"
#include "mystring.h"
#include "myexcept.h"
#include "libpro/srcaln/ProfileAlignment.h"


#ifndef SCALEDFWDBWD
#define SCALEDFWDBWD
#endif

// #ifndef MAATESTPRINT
// #define MAATESTPRINT
// #endif

// _________________________________________________________________________
// Class MAAlignment
//  implementation of maximum accuracy alignment
//
class MAAlignment: public ProfileAlignment
{
public:
    MAAlignment(
            const FrequencyMatrix& freq_fst, const LogOddsMatrix& logo_fst, const GapScheme& gaps_fst,
            const FrequencyMatrix& freq_sec, const LogOddsMatrix& logo_sec, const GapScheme& gaps_sec,
            int beg_fst, int end_fst, int beg_sec, int end_sec,
            const AbstractScoreMatrix*  usc_system,
            bool ungapped = false
    );

    virtual ~MAAlignment();
    virtual void    Run();

    bool    GetRestricted() const { return restricted_; }
    void    SetRestricted( bool value ) { restricted_ = value; }

    double  GetAReg() const { return areg_; }
    void    SetAReg( double value ) { areg_ = value; }

    double  GetARegExt() const { return aregext_; }
    void    SetARegExt( double value ) { aregext_ = value; }

    static bool GetGlobal() { return global_; }
    static void SetGlobal( bool value ) { global_ = value; }

protected:
    explicit MAAlignment();
    virtual void    ClearF();
    virtual void    ClearPointer();
    virtual void    AlignProfiles();
    virtual void    PostProcess();
    virtual void    MakeAlignmentPath();
    void            MakeAlignmentPath2();
    void            MakeAlignmentPathObs();

    void    ForwardDPScaled();
    void    BackwardDPScaled();
    void    ForwardDP();
    void    BackwardDP();
    void    AccuracyDP();
    void    AccuracyDP2();
    void    AccuracyDPObs();

    int     GetQueryBeg() const { return beg_fst_; }
    int     GetQueryEnd() const { return end_fst_; }
    int     GetSbjctBeg() const { return beg_sec_; }
    int     GetSbjctEnd() const { return end_sec_; }

    TPAScore    GetFE() const { return FE_; }
    void        SetFE( TPAScore value ) { FE_ = value; }
    void        CalculateFE();

    void        CalcGlobalFE();
    void        CalcLocalFE();

    void        CalcGlobalFEScaled();
    void        CalcLocalFEScaled();

    void    MAInitialize();
    void    MADestroy();
    void    ClearB();
    void    ClearBptr();
    void    ClearA();
    void    ClearAptr();
    void    ClearS();
    void    ClearScale();

    void    ClearFMarginHorz();
    void    ClearFMarginVert();

    void    ClearFMarginVertScaledAt( int n, int sbeg );
    void    ClearFMarginHorzScaledAt( int qbeg, int m );


private:
    int     beg_fst_;//beginning position of alignment for query
    int     end_fst_;//end position of alignment for query
    int     beg_sec_;//beginning position of alignment for subject
    int     end_sec_;//end position of alignment for query
    bool    restricted_;//restricted-to-given-boundaries calculations

    TPAScore    FE_;    //finalized F at end state (N,M)
    TPAScore    ( **B_ )[noStates];         //backward DP matrices
    int         ( **Bptr_ )[noStates];      //backtracing pointers
    TPAScore    ( **A_ )[noStates];         //target accuracy DP matrices
    int         ( **Aptr_ )[noStates];      //backtracing pointers
    TPAScore    ( **S_ )[noStates];         //matrices of original scores pathed in relation to A_
    double*     scale_;     //scale factor for forward/backward matrices
    double      powscale_;  //power scale factor

    static double   areg_;  //regulator (subt. term) for accuracy matrix
    static double   aregext_;//regulator for extention of gaps
    static bool     global_;//global configuration
};


////////////////////////////////////////////////////////////////////////////
// INLINES
//



#endif//__MAAlignment__
