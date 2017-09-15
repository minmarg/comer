/***************************************************************************
 *   Copyright (C) 2012 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "data.h"
#include "mystring.h"
#include "myexcept.h"
#include "ext/psl.h"
#include "libpro/srcpro/MOptions.h"
#include "MAAlignment.h"


#ifdef SCALEDFWDBWD
const TPAScore  cgMINPVAL = 0;
const TPAScore  cgMAXPVAL = 1;
#else
const TPAScore  cgMINPVAL = LOG_PROB_MIN;
const TPAScore  cgMAXPVAL = 0;
#endif

//use image scores
const bool      gbUSEIMAGE = true;
double          gOffset = 0.07;//global score offset
double          gModScMul = 0.7;//global multiplier for mod score
double          gScBase = 1.0;//principal score base

//global configuration
bool MAAlignment::global_ = false;//true;

//accuracy matrix regulator
double MAAlignment::areg_ = 0.3;
//regulator for extention of gaps
double MAAlignment::aregext_ = MAAlignment::areg_ * 0.5;

////////////////////////////////////////////////////////////////////////////
// CLASS MAAlignment
//
// constructor
//
MAAlignment::MAAlignment(
        const FrequencyMatrix& freq_fst, const LogOddsMatrix& logo_fst, const GapScheme& gaps_fst,
        const FrequencyMatrix& freq_sec, const LogOddsMatrix& logo_sec, const GapScheme& gaps_sec,
        int beg_fst, int end_fst, int beg_sec, int end_sec,
        const AbstractScoreMatrix*  usc_system,
        bool ungapped
    )
:   ProfileAlignment( freq_fst, logo_fst, gaps_fst, 
                      freq_sec, logo_sec, gaps_sec, usc_system, ungapped ),
    beg_fst_( beg_fst ),
    end_fst_( end_fst ),
    beg_sec_( beg_sec ),
    end_sec_( end_sec ),
    restricted_( false ),
    FE_( 0 ),
    B_( NULL ),
    Bptr_( NULL ),
    A_( NULL ),
    Aptr_( NULL ),
    S_( NULL ),
    scale_( NULL ),
    powscale_( 1.0 )
{
    if( beg_fst_ < 0 || beg_sec_ < 0 ||
        end_fst_ < beg_fst_ || end_sec_ < beg_sec_ ||
        freq_fst_.GetColumns() <= end_fst_ || freq_sec_.GetColumns() <= end_sec_ )
        throw myruntime_error("MAAlignment: Invalid start-end positions.");

    SetAReg( OPTIONS.GetMINPP());
    SetARegExt( 0.5 * GetAReg());

    if( !OPTIONS.GetSSSWGT()) {
        gOffset = 0.0;
        gModScMul = 0.0;
    }

    MAInitialize();
}

// -------------------------------------------------------------------------
// default constructor is not allowed
//
MAAlignment::MAAlignment()
:   ProfileAlignment(),
    beg_fst_( 0 ),
    end_fst_( 0 ),
    beg_sec_( 0 ),
    end_sec_( 0 ),
    restricted_( false ),
    FE_( 0 ),
    B_( NULL ),
    Bptr_( NULL ),
    A_( NULL ),
    Aptr_( NULL ),
    S_( NULL ),
    scale_( NULL ),
    powscale_( 1.0 )
{
}

// -------------------------------------------------------------------------
// destructor: deallocate memory
//
MAAlignment::~MAAlignment()
{
    MADestroy();
}

// -------------------------------------------------------------------------
// MADestroy: deallocate memory used by class structures
//
void MAAlignment::MADestroy()
{
    int n;

    if( B_ ) {
        for( n = 0; n < querySize_ + 2; n++ )
            if( B_[n])
                free( B_[n]);
        free( B_ );
        B_ = NULL;
    }
    if( Bptr_ ) {
        for( n = 0; n < querySize_ + 2; n++ )
            if( Bptr_[n])
                free( Bptr_[n]);
        free( Bptr_ );
        Bptr_ = NULL;
    }
    if( A_ ) {
        for( n = 0; n < querySize_ + 1; n++ )
            if( A_[n])
                free( A_[n]);
        free( A_ );
        A_ = NULL;
    }
    if( Aptr_ ) {
        for( n = 0; n < querySize_ + 1; n++ )
            if( Aptr_[n])
                free( Aptr_[n]);
        free( Aptr_ );
        Aptr_ = NULL;
    }
    if( S_ ) {
        for( n = 0; n < querySize_ + 1; n++ )
            if( S_[n])
                free( S_[n]);
        free( S_ );
        S_ = NULL;
    }
    if( scale_ ) {
        free( scale_ );
        scale_ = NULL;
    }
}

// -------------------------------------------------------------------------
// MAInitialize: create and initialize related structures and variables
//
void MAAlignment::MAInitialize()
{
    int n;

    if( GetQuerySize() < 1 || MAXCOLUMNS < GetQuerySize() ||
        GetSubjectSize() < 1 || MAXCOLUMNS < GetSubjectSize())
        throw myruntime_error( "MAAlignment: Invalid profile size." );

    MADestroy();

    //reserve with extra size of 2: 1+1 at beginning and end
    B_ = ( TPAScore(**)[noStates] )malloc( sizeof( TPAScore* ) * ( GetQuerySize() + 2 ));
    Bptr_ = ( int(**)[noStates] )malloc( sizeof( int* ) * ( GetQuerySize() + 2 ));
    //allocate memory for accuracy matrix
    A_ = ( TPAScore(**)[noStates] )malloc( sizeof( TPAScore* ) * ( GetQuerySize() + 1 ));
    Aptr_ = ( int(**)[noStates] )malloc( sizeof( int* ) * ( GetQuerySize() + 1 ));
    //allocation of memory for matrices of scores
    S_ = ( TPAScore(**)[noStates] )malloc( sizeof( TPAScore* ) * ( GetQuerySize() + 1 ));
    scale_ = ( double* )malloc( sizeof( double ) * ( GetQuerySize() + 2 ));

    if( !B_ || !Bptr_ ||
        !A_ || !Aptr_ || !S_ || !scale_ )
        throw myruntime_error( "MAAlignment: Not enough memory." );

    for( n = 0; n < GetQuerySize() + 2; n++ )
    {
        B_[n] = ( TPAScore(*)[noStates] )malloc( sizeof( TPAScore ) * ( GetSubjectSize() + 2 ) * noStates );
        Bptr_[n] = ( int(*)[noStates] )malloc( sizeof( int ) * ( GetSubjectSize() + 2 ) * noStates );
        if( !B_[n] || !Bptr_[n] )
            throw myruntime_error( "MAAlignment: Not enough memory." );
    }
    for( n = 0; n < GetQuerySize() + 1; n++ )
    {
        A_[n] = ( TPAScore(*)[noStates] )malloc( sizeof( TPAScore ) * ( GetSubjectSize() + 1 ) * noStates );
        Aptr_[n] = ( int(*)[noStates] )malloc( sizeof( int ) * ( GetSubjectSize() + 1 ) * noStates );
        S_[n] = ( TPAScore(*)[noStates] )malloc( sizeof( TPAScore ) * ( GetSubjectSize() + 1 ) * noStates );
        if( !A_[n] || !Aptr_[n] || !S_[n])
            throw myruntime_error( "MAAlignment: Not enough memory." );
    }
}

// -------------------------------------------------------------------------
// ClearF: clear dynamic programming structure
//
void MAAlignment::ClearF()
{
    int qbeg = GetRestricted()? GetQueryBeg(): 0;
    int sbeg = GetRestricted()? GetSbjctBeg(): 0;
    int n, m;
    if( F == NULL )
        return;

    for( n = 0; n < GetQuerySize() + 1; n++ ) {
        if( F[n] == NULL )
            throw myruntime_error("MAAlignment: Null vector of DPS.");
        for( m = 0; m < GetSubjectSize() + 1; m++ ) {
            if( GetGlobal()) {
                if( n == qbeg && m == sbeg ) {
                    //beginning position
                    F[n][m][stateMM] = cgMAXPVAL;
                    F[n][m][stateMI] = cgMAXPVAL;//cgMINPVAL;
                    F[n][m][stateIM] = cgMAXPVAL;//cgMINPVAL;
                    F[n][m][stateDG] = cgMAXPVAL;//cgMINPVAL;
                    F[n][m][stateGD] = cgMAXPVAL;//cgMINPVAL;
                }
                else {
                    F[n][m][stateMM] = cgMINPVAL;
                    F[n][m][stateMI] = cgMINPVAL;
                    F[n][m][stateIM] = cgMINPVAL;
                    F[n][m][stateDG] = cgMINPVAL;
                    F[n][m][stateGD] = cgMINPVAL;
                }
            }
            else {
                if( n == qbeg || m == sbeg )
                    //beginning position
                    F[n][m][stateMM] = cgMAXPVAL;
                else
                    F[n][m][stateMM] = cgMINPVAL;
                F[n][m][stateMI] = cgMINPVAL;
                F[n][m][stateIM] = cgMINPVAL;
                F[n][m][stateDG] = cgMINPVAL;
                F[n][m][stateGD] = cgMINPVAL;
            }
        }
    }

    ClearB();
    ClearA();
    ClearS();
    ClearScale();
}

// -------------------------------------------------------------------------
// ClearB: clear backward DP matrices
//
void MAAlignment::ClearB()
{
    int qend = GetRestricted()? GetQueryEnd() + 2: GetQuerySize() + 1;
    int send = GetRestricted()? GetSbjctEnd() + 2: GetSubjectSize() + 1;
    int n, m;
    if( B_ == NULL )
        return;

    if( !GetGlobal()) {
        qend--;
        send--;
    }

    for( n = GetQuerySize() + 1; 0 <= n; n-- ) {
        if( B_[n] == NULL )
            throw myruntime_error("MAAlignment: Null vector of backward DPS.");
        for( m = GetSubjectSize() + 1; 0 <= m; m-- ) {
            if( GetGlobal()) {
                if( n == qend && m == send ) {
                    //end position
                    B_[n][m][stateMM] = cgMAXPVAL;
                    B_[n][m][stateMI] = cgMAXPVAL;
                    B_[n][m][stateIM] = cgMAXPVAL;
                    B_[n][m][stateDG] = cgMAXPVAL;
                    B_[n][m][stateGD] = cgMAXPVAL;
                }
                else {
                    B_[n][m][stateMM] = cgMINPVAL;
                    B_[n][m][stateMI] = cgMINPVAL;
                    B_[n][m][stateIM] = cgMINPVAL;
                    B_[n][m][stateDG] = cgMINPVAL;
                    B_[n][m][stateGD] = cgMINPVAL;
                }
            }
            else {
                if( n == qend || m == send )
                    B_[n][m][stateMM] = cgMAXPVAL;
                else
                    B_[n][m][stateMM] = cgMINPVAL;
                B_[n][m][stateMI] = cgMINPVAL;
                B_[n][m][stateIM] = cgMINPVAL;
                B_[n][m][stateDG] = cgMINPVAL;
                B_[n][m][stateGD] = cgMINPVAL;
            }
        }
    }
}

// -------------------------------------------------------------------------
// ClearA: clear accuracy DP matrices
//
void MAAlignment::ClearA()
{
    int n, m;
    if( A_ == NULL )
        return;

    for( n = 0; n < GetQuerySize() + 1; n++ ) {
        if( A_[n] == NULL )
            throw myruntime_error("MAAlignment: Null vector of DPS.");
        for( m = 0; m < GetSubjectSize() + 1; m++ ) {
            A_[n][m][stateMM] = 0;
            A_[n][m][stateMI] = 0;
            A_[n][m][stateIM] = 0;
            A_[n][m][stateDG] = 0;
            A_[n][m][stateGD] = 0;
        }
    }
}

// -------------------------------------------------------------------------
// ClearS: clear matrix of running scores
//
void MAAlignment::ClearS()
{
    int n, m;
    if( S_ == NULL )
        return;

    for( n = 0; n < GetQuerySize() + 1; n++ ) {
        if( S_[n] == NULL )
            throw myruntime_error("MAAlignment: Null vector of matrix of running scores.");
        for( m = 0; m < GetSubjectSize() + 1; m++ ) {
            S_[n][m][stateMM] = 0;
            S_[n][m][stateMI] = 0;
            S_[n][m][stateIM] = 0;
            S_[n][m][stateDG] = 0;
            S_[n][m][stateGD] = 0;
        }
    }
}

// -------------------------------------------------------------------------
// ClearScale: clear vector of scale factors
//
void MAAlignment::ClearScale()
{
    int n;
    powscale_ = 1.0;
    if( scale_ == NULL )
        return;
    for( n = 0; n < GetQuerySize() + 2; n++ )
        scale_[n] = 1.0;
}

// -------------------------------------------------------------------------
// ClearPointer: clear back-tracing pointer structure
//
void MAAlignment::ClearPointer()
{
    ProfileAlignment::ClearPointer();
    ClearBptr();
    ClearAptr();
}

// -------------------------------------------------------------------------
// ClearBptr: clear backward DPS' back-tracing pointer structure
//
void MAAlignment::ClearBptr()
{
    int n;
    if( Bptr_ )
        for( n = 0; n < GetQuerySize() + 2; n++ )
            if( Bptr_[n])
                memset( Bptr_[n], 0, sizeof( int )*( GetSubjectSize() + 2 ) * noStates );
}

// -------------------------------------------------------------------------
// ClearAptr: clear accuracy DPS' back-tracing pointer structure
//
void MAAlignment::ClearAptr()
{
    int n;
    if( Aptr_ )
        for( n = 0; n < GetQuerySize() + 1; n++ )
            if( Aptr_[n])
                memset( Aptr_[n], 0, sizeof( int )*( GetSubjectSize() + 1 ) * noStates );
}

// -------------------------------------------------------------------------
// Run: prepare and start alignment of profiles
//
void MAAlignment::Run()
{
    Clear();
    AlignProfiles();                //1.
#ifdef MAATESTPRINT
    MakeAlignmentPathObs();
#else
    MakeAlignmentPath();//Prev.
//     MakeAlignmentPath2();
#endif
    PostProcess();
    SetFinalScore( GetAlnScore());
}

// -------------------------------------------------------------------------
// PostProcess: post process alignment
//
void MAAlignment::PostProcess()
{
}

// -------------------------------------------------------------------------
// AlignProfiles: dynamic programming to align a pair of profiles;
//  semi-probabilistic approach to alignment
//
void MAAlignment::AlignProfiles()
{
#ifdef SCALEDFWDBWD
    ForwardDPScaled();
    BackwardDPScaled();
#else
    ForwardDP();
    BackwardDP();
#endif
    CalculateFE();
#ifdef MAATESTPRINT
    AccuracyDPObs();
#else
    AccuracyDP();//Prev.
//     AccuracyDP2();
#endif
}

// -------------------------------------------------------------------------
// CalculateFE: calculate full probability at end state of DPS F applied 
//  forward calculations;
//
void MAAlignment::CalculateFE()
{
    if( GetGlobal())
#ifdef SCALEDFWDBWD
        CalcGlobalFEScaled();
#else
        CalcGlobalFE();
#endif
    else
#ifdef SCALEDFWDBWD
        CalcLocalFEScaled();
#else
        CalcLocalFE();
#endif
}

// -------------------------------------------------------------------------
// ScBack: scale back log probabilities
//
inline double ScBack( double lprob, double scale )
{
    if( scale == 1.0 )
        return lprob;
    if( scale )
        return lprob / scale;
    throw myruntime_error( "MAAlignment: ScBack: Invalid scale factor." );
}

// -------------------------------------------------------------------------
// GetLogTransProbsAt: get log transition probability at the position
//
inline double GetLogTransProbsAt( const GapScheme& gsch, int trans, int pos )
{
    if( gbUSEIMAGE )
        return gsch.GetNSLogTransProbsAt( trans, pos );
    return gsch.GetLogTransProbsAt( trans, pos );
}



// -------------------------------------------------------------------------
// CalcGlobalFE: calculate full probability at end state of DPS F applied 
//  forward calculations;
//  NOTE: in global cfg. this full probability should be equal to the 
//  backward DPS B_ value at the beginning state: 
//  FE = B_[0][0][stateMM] (if B_ was calculated down to beginning state 
//  point (0,0) and F[0][0][stateMM]=0, F[0][0][.]=-inf)
//
void MAAlignment::CalcGlobalFE()
{
    TPAScore    bestMM;
    TPAScore    oMM;
    TPAScore    oMI, oIM;
    TPAScore    oDG, oGD;

    double      tpqMM, tpsMM;   //transition probabilities for query and subject
    double      tpqIM, tpsIM;
    double      tpqDM, tpsDM;

    double      maxlog; //maximum value of logarithms

    int n = GetRestricted()? GetQueryEnd() + 1: GetQuerySize();
    int m = GetRestricted()? GetSbjctEnd() + 1: GetSubjectSize();

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();
    if( !scmatrix )
        throw myruntime_error( "MAAlignment: Null score matrix." );

    double scs = scmatrix->GetAutoScalingFactor();
    if( gbUSEIMAGE )
        scs = 1.0;

    //-----------------
    //state MM
    tpqMM = GetLogTransProbsAt( gaps_fst_, P_MM, n-1 );//index value of -1 is ok: beg. state
    tpsMM = GetLogTransProbsAt( gaps_sec_, P_MM, m-1 );
    oMM = F[n][m][stateMM];
    oMM += LogOfProduct( tpqMM, tpsMM, !IsUngapped());

    //state MI
    tpsIM = GetLogTransProbsAt( gaps_sec_, P_IM, m-1 );
    oMI = F[n][m][stateMI];
    oMI += LogOfProduct( tpqMM, tpsIM, !IsUngapped());

    //state IM
    tpqIM = GetLogTransProbsAt( gaps_fst_, P_IM, n-1 );
    oIM = F[n][m][stateIM];
    oIM += LogOfProduct( tpqIM, tpsMM, !IsUngapped());

    //state DG
    tpqDM = GetLogTransProbsAt( gaps_fst_, P_DM, n-1 );
    oDG = F[n][m][stateDG];
    oDG += LogOfProduct( tpqDM, tpsMM, !IsUngapped());

    //state GD
    tpsDM = GetLogTransProbsAt( gaps_sec_, P_DM, m-1 );
    oGD = F[n][m][stateGD];
    oGD += LogOfProduct( tpqMM, tpsDM, !IsUngapped());

    maxlog = SLC_MAX5( oMM, oMI, oIM, oDG, oGD );
    oMM -= maxlog;
    oMI -= maxlog;
    oIM -= maxlog;
    oDG -= maxlog;
    oGD -= maxlog;
    //save log Z_MM(n,m)
    bestMM = exp(ScBack(oMM,scs)) + exp(ScBack(oMI,scs)) + exp(ScBack(oIM,scs)) + 
             exp(ScBack(oDG,scs)) + exp(ScBack(oGD,scs));
    if( 0.0 < bestMM )
        bestMM = maxlog + log( bestMM ) * scs;
    else
        bestMM = LOG_PROB_MIN;

    SetFE( bestMM );
}

// -------------------------------------------------------------------------
// CalcLocalFE: calculate full probability at end state of DPS F applied 
//  forward calculations;
//
void MAAlignment::CalcLocalFE()
{
    double      maxlog = 0;
    TPAScore    FMM = 0;
    TPAScore    BMM = 0;

    int         qbeg = GetRestricted()? GetQueryBeg() + 1: 1;
    int         qend = GetRestricted()? GetQueryEnd() + 1: GetQuerySize();
    int         sbeg = GetRestricted()? GetSbjctBeg() + 1: 1;
    int         send = GetRestricted()? GetSbjctEnd() + 1: GetSubjectSize();
    int n, m;

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();
    if( !scmatrix )
        throw myruntime_error( "MAAlignment: Null score matrix." );

    double scs = scmatrix->GetAutoScalingFactor();
    if( gbUSEIMAGE )
        scs = 1.0;

    for( n = qbeg; n <= qend; n++ )
        for( m = sbeg; m <= send; m++ )
            if( maxlog < F[n][m][stateMM])
                maxlog = F[n][m][stateMM];

    FMM += exp( ScBack(-maxlog,scs));
    for( n = qbeg; n <= qend; n++ )
        for( m = sbeg; m <= send; m++ )
            FMM += exp( ScBack( F[n][m][stateMM] - maxlog,scs));

    if( 0 < FMM )
        FMM = maxlog + log( FMM ) * scs;
    else
        FMM = LOG_PROB_MIN;

    SetFE( FMM );
}

// -------------------------------------------------------------------------
// CalcGlobalFEscaled: calculate full probability at end state of scaled F
//  NOTE: in global cfg. this full probability should be equal to the 
//  backward DPS B_ value at the beginning state: 
//  FE = B_[0][0][stateMM] (if B_ was calculated down to beginning state 
//  point (0,0) and F[0][0][stateMM]=0, F[0][0][.]=-inf)
//
void MAAlignment::CalcGlobalFEScaled()
{
    TPAScore    bestMM;
    TPAScore    oMM;
    TPAScore    oMI, oIM;
    TPAScore    oDG, oGD;

    double      tpqMM, tpsMM;   //transition probabilities
    double      tpqIM, tpsIM;
    double      tpqDM, tpsDM;

    int n = GetRestricted()? GetQueryEnd() + 1: GetQuerySize();
    int m = GetRestricted()? GetSbjctEnd() + 1: GetSubjectSize();

    //-----------------
    //state MM
    tpqMM = gaps_fst_.GetTransProbsAt( P_MM, n-1 );//index value of -1 is ok: beg. state
    tpsMM = gaps_sec_.GetTransProbsAt( P_MM, m-1 );
    oMM = F[n][m][stateMM];
    oMM *= ProbProduct( tpqMM, tpsMM, !IsUngapped());

    //state MI
    tpsIM = gaps_sec_.GetTransProbsAt( P_IM, m-1 );
    oMI = F[n][m][stateMI];
    oMI *= ProbProduct( tpqMM, tpsIM, !IsUngapped());

    //state IM
    tpqIM = gaps_fst_.GetTransProbsAt( P_IM, n-1 );
    oIM = F[n][m][stateIM];
    oIM *= ProbProduct( tpqIM, tpsMM, !IsUngapped());

    //state DG
    tpqDM = gaps_fst_.GetTransProbsAt( P_DM, n-1 );
    oDG = F[n][m][stateDG];
    oDG *= ProbProduct( tpqDM, tpsMM, !IsUngapped());

    //state GD
    tpsDM = gaps_sec_.GetTransProbsAt( P_DM, m-1 );
    oGD = F[n][m][stateGD];
    oGD *= ProbProduct( tpqMM, tpsDM, !IsUngapped());

    bestMM = scale_[n+1] * ( oMM + oMI + oIM + oDG + oGD );
    SetFE( bestMM );
}

// -------------------------------------------------------------------------
// CalcLocalFEScaled: calculate full probability at end state of scaled F
//
void MAAlignment::CalcLocalFEScaled()
{
    TPAScore    FMM = 0;
    TPAScore    BMM = 0;

    int         qbeg = GetRestricted()? GetQueryBeg() + 1: 1;
    int         qend = GetRestricted()? GetQueryEnd() + 1: GetQuerySize();
    int         sbeg = GetRestricted()? GetSbjctBeg() + 1: 1;
    int         send = GetRestricted()? GetSbjctEnd() + 1: GetSubjectSize();
    int n, m;

    FMM += 1.0;
    for( n = qbeg; n <= qend; n++ ) {
        for( m = sbeg; m <= send; m++ )
            FMM += F[n][m][stateMM];
        FMM *= scale_[n+1];
    }
    SetFE( FMM );
}





// -------------------------------------------------------------------------
// ClearFMarginVert: clear and reinitialize vertical margin of DPS F
//
void MAAlignment::ClearFMarginVert()
{
    TPAScore    currentMM;
    TPAScore    currentMI, bestMI;
    TPAScore    currentDG, bestDG;

    TPAScore    upMI_MM_MI, upMI_MM_II;
    TPAScore    upDG_MD_1, upDG_DD_1;

    double      tpqMM_1;//transition probabilities for query and subject
    double      tpqMD_1;
    double      tpqDD_1;
    double      tpsMI;
    double      tpsII;

    double      maxlog;//maximum value of logarithms
    int         ptr;//value of pointer
    int         n;
    int         qbeg = GetRestricted()? GetQueryBeg() + 1: 1;
    int         qend = GetRestricted()? GetQueryEnd() + 1: GetQuerySize();
    int         sbeg = GetRestricted()? GetSbjctBeg(): 0;

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();
    if( !scmatrix )
        throw myruntime_error( "MAAlignment: Null score matrix." );

    double scs = scmatrix->GetAutoScalingFactor();
    if( gbUSEIMAGE )
        scs = 1.0;

    //assign (.,0) probabilities
    for( n = qbeg; n <= qend; n++ )
    {
        //-----------------
        //state MI: enter
        tpqMM_1 = GetLogTransProbsAt( gaps_fst_, P_MM, n-2 );//-1 is beg. state

        currentMM = F[n-1][sbeg][stateMM];
        upMI_MM_MI = currentMM; //query

        tpsMI = GetLogTransProbsAt( gaps_sec_, P_MI, sbeg-1 );
        upMI_MM_MI += LogOfProduct( tpqMM_1, tpsMI, !IsUngapped());

        //state MI: extend
        currentMI = F[n-1][sbeg][stateMI];
        upMI_MM_II = currentMI; //query

        tpsII = GetLogTransProbsAt( gaps_sec_, P_II, sbeg-1 );
        upMI_MM_II += LogOfProduct( tpqMM_1, tpsII, !IsUngapped());

        //-----------------
        //state DG: enter
        upDG_MD_1 = currentMM;  //query

        tpqMD_1 = GetLogTransProbsAt( gaps_fst_, P_MD, n-2 );
        upDG_MD_1 += LogOfProduct( tpqMD_1, 0.0, !IsUngapped());

        //state DG: extend
        currentDG = F[n-1][sbeg][stateDG];
        upDG_DD_1 = currentDG;  //query

        tpqDD_1 = GetLogTransProbsAt( gaps_fst_, P_DD, n-2 );
        upDG_DD_1 += LogOfProduct( tpqDD_1, 0.0, !IsUngapped());

        //process state stateMI
        //.....................
        if( upMI_MM_II < upMI_MM_MI ) {
            maxlog = upMI_MM_MI;
            ptr = dMM; //state for backtracing
        } else if( upMI_MM_MI < upMI_MM_II ) {
                    maxlog = upMI_MM_II;
                    ptr = dMI; //up direction
                } else { //upMI_MM_MI == upMI_MM_II
                    maxlog = upMI_MM_MI;
                    ptr = dMM_MI; //diagonal or up
                }
        upMI_MM_MI -= maxlog;
        upMI_MM_II -= maxlog;
        //save log Z_MI(n,m)
        bestMI = exp(ScBack(upMI_MM_MI,scs)) + exp(ScBack(upMI_MM_II,scs));
        if( 0.0 < bestMI )
            bestMI = maxlog + log( bestMI ) * scs;
        else
            bestMI = LOG_PROB_MIN;
        pointer[n][sbeg][stateMI] = ptr;
        F[n][sbeg][stateMI] = bestMI;

        //process state stateDG
        //.....................
        if( upDG_DD_1 < upDG_MD_1 ) {
            maxlog = upDG_MD_1;
            ptr = dMM; //backtracing state
        } else if( upDG_MD_1 < upDG_DD_1 ) {
                    maxlog = upDG_DD_1;
                    ptr = dDG; //up direction
                } else { //upDG_MD_1 == upDG_DD_1
                    maxlog = upDG_MD_1;
                    ptr = dMM_DG; //diagonal or up
                }
        upDG_MD_1 -= maxlog;
        upDG_DD_1 -= maxlog;
        //save log Z_DG(n,m)
        bestDG = exp(ScBack(upDG_MD_1,scs)) + exp(ScBack(upDG_DD_1,scs));
        if( 0.0 < bestDG )
            bestDG = maxlog + log( bestDG ) * scs;
        else
            bestDG = LOG_PROB_MIN;
        pointer[n][sbeg][stateDG] = ptr;
        F[n][sbeg][stateDG] = bestDG;
    }
}

// -------------------------------------------------------------------------
// ClearFMarginHorz: clear and reinitialize horizontal margin of DPS F
//
void MAAlignment::ClearFMarginHorz()
{
    TPAScore    bestMM;
    TPAScore    bestIM;
    TPAScore    bestGD;

    TPAScore    leftIM_MI_MM, leftIM_II_MM;
    TPAScore    leftGD_1_MD, leftGD_1_DD;

    double      tpsMM_1;//transition probabilities for query and subject
    double      tpsMD_1;
    double      tpsDD_1;
    double      tpqMI;
    double      tpqII;

    double      maxlog;//maximum value of logarithms
    int         ptr;//value of pointer
    int         m;
    int         sbeg = GetRestricted()? GetSbjctBeg() + 1: 1;
    int         send = GetRestricted()? GetSbjctEnd() + 1: GetSubjectSize();
    int         qbeg = GetRestricted()? GetQueryBeg(): 0;

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();
    if( !scmatrix )
        throw myruntime_error( "MAAlignment: Null score matrix." );

    double scs = scmatrix->GetAutoScalingFactor();
    if( gbUSEIMAGE )
        scs = 1.0;

    //assign (0,.) probabilities
    for( m = sbeg; m <= send; m++ )
    {
        //-----------------
        //state IM: enter
        tpsMM_1 = GetLogTransProbsAt( gaps_sec_, P_MM, m-2 );

        bestMM = F[qbeg][m-1][stateMM];
        leftIM_MI_MM = bestMM;   //subject

        tpqMI = GetLogTransProbsAt( gaps_fst_, P_MI, qbeg-1 );
        leftIM_MI_MM += LogOfProduct( tpqMI, tpsMM_1, !IsUngapped());

        //state IM: extend
        bestIM = F[qbeg][m-1][stateIM];
        leftIM_II_MM = bestIM;   //subject

        tpqII = GetLogTransProbsAt( gaps_fst_, P_II, qbeg-1 );
        leftIM_II_MM += LogOfProduct( tpqII, tpsMM_1, !IsUngapped());

        //-----------------
        //state GD: enter
        leftGD_1_MD = bestMM;    //subject

        tpsMD_1 = GetLogTransProbsAt( gaps_sec_, P_MD, m-2 );
        leftGD_1_MD += LogOfProduct( 0.0, tpsMD_1, !IsUngapped());

        //state GD: extend
        bestGD = F[qbeg][m-1][stateGD];
        leftGD_1_DD = bestGD;   //subject

        tpsDD_1 = GetLogTransProbsAt( gaps_sec_, P_DD, m-2 );
        leftGD_1_DD += LogOfProduct( 0.0, tpsDD_1, !IsUngapped());

        //process state stateIM
        //.....................
        if( leftIM_II_MM < leftIM_MI_MM ) {
            maxlog = leftIM_MI_MM;
            ptr = dMM; //state for backtracing
        } else if( leftIM_MI_MM < leftIM_II_MM ) {
                    maxlog = leftIM_II_MM;
                    ptr = dIM; //left direction
                } else { //leftIM_MI_MM == leftIM_II_MM
                    maxlog = leftIM_MI_MM;
                    ptr = dMM_IM; //diagonal or left
                }
        leftIM_MI_MM -= maxlog;
        leftIM_II_MM -= maxlog;
        //save log Z_IM(n,m)
        bestIM = exp(ScBack(leftIM_MI_MM,scs)) + exp(ScBack(leftIM_II_MM,scs));
        if( 0.0 < bestIM )
            bestIM = maxlog + log( bestIM ) * scs;
        else
            bestIM = LOG_PROB_MIN;
        pointer[qbeg][m][stateIM] = ptr;
        F[qbeg][m][stateIM] = bestIM;

        //process state stateGD
        //.....................
        if( leftGD_1_DD < leftGD_1_MD ) {
            maxlog = leftGD_1_MD;
            ptr = dMM; //backtracing state
        } else if( leftGD_1_MD < leftGD_1_DD ) {
                    maxlog = leftGD_1_DD;
                    ptr = dGD; //left direction
                } else { //leftGD_1_MD == leftGD_1_DD
                    maxlog = leftGD_1_MD;
                    ptr = dMM_GD; //diagonal or left
                }
        leftGD_1_MD -= maxlog;
        leftGD_1_DD -= maxlog;
        //save log Z_GD(n,m)
        bestGD = exp(ScBack(leftGD_1_MD,scs)) + exp(ScBack(leftGD_1_DD,scs));
        if( 0.0 < bestGD )
            bestGD = maxlog + log( bestGD ) * scs;
        else
            bestGD = LOG_PROB_MIN;
        pointer[qbeg][m][stateGD] = ptr;
        F[qbeg][m][stateGD] = bestGD;
    }
}

// -------------------------------------------------------------------------
// ForwardDP: forward algorithm of dynamic programming to get forward
//  probabilities
//
void MAAlignment::ForwardDP()
{
    const TPAScore  cllprob = GetGlobal()? LOG_PROB_MIN: 0;//log of local probability

    TPAScore    bestMM, currentMM, oMM; //best, current, and operating scores for state stateMM
    TPAScore    bestMI, currentMI, oMI; //best, current, and operating scores for state stateMI
    TPAScore    bestIM, currentIM, oIM; //best, current, and operating scores for state stateIM
    TPAScore    bestDG, currentDG, oDG; //best, current, and operating scores for state stateDG
    TPAScore    bestGD, currentGD, oGD; //best, current, and operating scores for state stateGD

    TPAScore    upMI_MM_MI, upMI_MM_II;     //state MI: score of entering MI state (trans:MM-MI) and extending it (MM-II)
    TPAScore    leftIM_MI_MM, leftIM_II_MM; //state IM: score of entering IM state (trans:MI-MM) and extending it (II-MM)

    TPAScore    upDG_MD_1, upDG_DD_1;       //state DG: score of entering DG state (trans:MD-1) and extending it (DD-1)
    TPAScore    leftGD_1_MD, leftGD_1_DD;   //state GD: score of entering GD state (trans:1-MD) and extending it (1-DD)

    double      tpqMM_1, tpsMM_1;   //transition probabilities for query and subject
    double      tpqIM_1, tpsIM_1;
    double      tpqMD_1, tpsMD_1;
    double      tpqDM_1, tpsDM_1;
    double      tpqDD_1, tpsDD_1;
    double      tpqMI, tpsMI;
    double      tpqII, tpsII;

    double      maxlog; //maximum value of logarithms

    TPAScore    score;                  //profile score at certain positions
    int         ptr;                    //value of pointer
    int         n, m;
    int         qbeg = GetRestricted()? GetQueryBeg() + 1: 1;
    int         qend = GetRestricted()? GetQueryEnd() + 1: GetQuerySize();
    int         sbeg = GetRestricted()? GetSbjctBeg() + 1: 1;
    int         send = GetRestricted()? GetSbjctEnd() + 1: GetSubjectSize();

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();
    if( !scmatrix )
        throw myruntime_error( "MAAlignment: Unable to align profiles." );

    double          epa = GetExpectPerAlignment();
    double          dyno = 0.0;//dynamic offset
    double          dmul = 0.0;//dynamic multiplier
    double          dgrd = 1.0;
    double          delt = dgrd - gScBase;
    //s=y+x; s^n=ay+bx=a(s-x)+bx=as+cx, a<1,b>1
    if( SLC_SQRT_DBL_MIN < epa ) {
        if( epa < 1.0 ) {
              dyno = gOffset * epa;
              if( 0 < delt )
                  dgrd -= delt * ( 1.0 - epa );
              dmul = gModScMul * ( 1.0 - epa );
        }
        else  dyno = gOffset;
    }

    const double    pmult = scmatrix->GetMultiplier();
    double scs = scmatrix->GetAutoScalingFactor();
    if( gbUSEIMAGE )
        scs = 1.0;

// // const_cast<GapScheme&>(gaps_fst_).OutputGapScheme(); fprintf(stderr,"\n\n\n");
// // const_cast<GapScheme&>(gaps_sec_).OutputGapScheme();

    if( GetGlobal()) {
        //prepare margins of F
        ClearFMarginVert();
        ClearFMarginHorz();
    }

    //start full probabilistic (global alignment) calculations
    //iterate over query positions
    for( n = qbeg; n <= qend; n++ )
    {
        bestMM      = F[n][sbeg-1][stateMM];
        bestIM      = F[n][sbeg-1][stateIM];
        bestGD      = F[n][sbeg-1][stateGD];

        currentMM   = F[n-1][sbeg-1][stateMM];
        currentMI   = F[n-1][sbeg-1][stateMI];
        currentDG   = F[n-1][sbeg-1][stateDG];

// fprintf(stderr,"\n" );
        //iterate over subject positions
        for( m = sbeg; m <= send; m++ )
        {
            if( gbUSEIMAGE ) {
                if( scmatrix->GetUseModImage() && dmul )
                    score =(scmatrix->GetImageScore(m-1,n-1)*dgrd - dyno +
                            scmatrix->GetModImageScore(m-1,n-1)*dmul ) * pmult;
                else
//                     score = scmatrix->GetImageScore( m-1, n-1 ) * pmult;
                    score =(scmatrix->GetImageScore( m-1, n-1 ) - dyno) * pmult;
            } else {
                if( scmatrix->GetUseModImage())
                    score = scmatrix->GetModScore( m-1, n-1 );
                else
                    score = scmatrix->GetScore( m-1, n-1 ); //profile score at positions m and n
//                 score  = AutocorrScore( scmatrix, m-1, n-1 ); //TEST:autocorrelation score
// fprintf(stderr,"\n %d %d: score: %.4f\n", n, m, exp(score/scs));
            }

            //-----------------
            //state MI: enter
            tpqMM_1 = GetLogTransProbsAt( gaps_fst_, P_MM, n-2 );//index value of -1 is ok: beg. state
            tpsMM_1 = GetLogTransProbsAt( gaps_sec_, P_MM, m-2 );
            oMM = currentMM;
            oMM += LogOfProduct( tpqMM_1, tpsMM_1, !IsUngapped());
// fprintf(stderr," MM_MM_MM: %.4f tpqMM_1=%.4f tpsMM_1=%.4f\n", 
// exp(LogOfProduct(tpqMM_1,tpsMM_1,!IsUngapped())/scs),exp(tpqMM_1/scs),exp(tpsMM_1/scs));

            //F has EXTRA row and column
            currentMM = F[n-1][m][stateMM];
            upMI_MM_MI = currentMM; //query

            tpsMI = GetLogTransProbsAt( gaps_sec_, P_MI, m-1 );
            upMI_MM_MI += LogOfProduct( tpqMM_1, tpsMI, !IsUngapped());
// fprintf(stderr," MI_MM_MI: %.4f tpsMI=%.4f\n",
// exp(LogOfProduct(tpqMM_1,tpsMI,!IsUngapped())/scs),exp(tpsMI/scs));

            //state MI: extend
            tpsIM_1 = GetLogTransProbsAt( gaps_sec_, P_IM, m-2 );
            oMI = currentMI;
            oMI += LogOfProduct( tpqMM_1, tpsIM_1, !IsUngapped());
// fprintf(stderr," MM_MM_IM: %.4f tpsIM_1=%.4f\n",
// exp(LogOfProduct(tpqMM_1,tpsIM_1,!IsUngapped())/scs),exp(tpsIM_1/scs));

            currentMI = F[n-1][m][stateMI];
            upMI_MM_II = currentMI; //query

            tpsII = GetLogTransProbsAt( gaps_sec_, P_II, m-1 );
            upMI_MM_II += LogOfProduct( tpqMM_1, tpsII, !IsUngapped());
// fprintf(stderr," MI_MM_II: %.4f tpsII=%.4f\n",
// exp(LogOfProduct(tpqMM_1,tpsII,!IsUngapped())/scs),exp(tpsII/scs));

            //-----------------
            //state IM: enter
            leftIM_MI_MM = bestMM;   //subject

            tpqMI = GetLogTransProbsAt( gaps_fst_, P_MI, n-1 );
            leftIM_MI_MM += LogOfProduct( tpqMI, tpsMM_1, !IsUngapped());
// fprintf(stderr," IM_MI_MM: %.4f tpqMI=%.4f\n",
// exp(LogOfProduct(tpqMI,tpsMM_1,!IsUngapped())/scs),exp(tpqMI/scs));

            //state IM: extend
            tpqIM_1 = GetLogTransProbsAt( gaps_fst_, P_IM, n-2 );
            oIM = F[n-1][m-1][stateIM];
            oIM += LogOfProduct( tpqIM_1, tpsMM_1, !IsUngapped());
// fprintf(stderr," MM_IM_MM: %.4f tpqIM_1=%.4f\n",
// exp(LogOfProduct(tpqIM_1,tpsMM_1,!IsUngapped())/scs),exp(tpqIM_1/scs));

            leftIM_II_MM = bestIM;   //subject

            tpqII = GetLogTransProbsAt( gaps_fst_, P_II, n-1 );
            leftIM_II_MM += LogOfProduct( tpqII, tpsMM_1, !IsUngapped());
// fprintf(stderr," IM_II_MM: %.4f tpqII=%.4f\n",
// exp(LogOfProduct(tpqII,tpsMM_1,!IsUngapped())/scs),exp(tpqII/scs));

            //-----------------
            //state DG: enter
            upDG_MD_1 = currentMM;  //query

            tpqMD_1 = GetLogTransProbsAt( gaps_fst_, P_MD, n-2 );
            upDG_MD_1 += LogOfProduct( tpqMD_1, 0.0, !IsUngapped());
// fprintf(stderr," DG_MD: %.4f tpqMD_1=%.4f\n",
// exp(LogOfProduct(tpqMD_1,0.0,!IsUngapped())/scs),exp(tpqMD_1/scs));

            //state DG: extend
            tpqDM_1 = GetLogTransProbsAt( gaps_fst_, P_DM, n-2 );
            oDG = currentDG;
            oDG += LogOfProduct( tpqDM_1, tpsMM_1, !IsUngapped());
// fprintf(stderr," MM_DM_MM: %.4f tpqDM_1=%.4f\n",
// exp(LogOfProduct(tpqDM_1,tpsMM_1,!IsUngapped())/scs),exp(tpqDM_1/scs));

            currentDG = F[n-1][m][stateDG];
            upDG_DD_1 = currentDG;  //query

            tpqDD_1 = GetLogTransProbsAt( gaps_fst_, P_DD, n-2 );
            upDG_DD_1 += LogOfProduct( tpqDD_1, 0.0, !IsUngapped());
// fprintf(stderr," DG_DD: %.4f tpqDD_1=%.4f\n",
// exp(LogOfProduct(tpqDD_1,0.0,!IsUngapped())/scs),exp(tpqDD_1/scs));

            //-----------------
            //state GD: enter
            leftGD_1_MD = bestMM;    //subject

            tpsMD_1 = GetLogTransProbsAt( gaps_sec_, P_MD, m-2 );
            leftGD_1_MD += LogOfProduct( 0.0, tpsMD_1, !IsUngapped());
// fprintf(stderr," GD_MD: %.4f tpsMD_1=%.4f\n",
// exp(LogOfProduct(0.0,tpsMD_1,!IsUngapped())/scs),exp(tpsMD_1/scs));

            //state GD: extend
            tpsDM_1 = GetLogTransProbsAt( gaps_sec_, P_DM, m-2 );
            oGD = F[n-1][m-1][stateGD];
            oGD += LogOfProduct( tpqMM_1, tpsDM_1, !IsUngapped());
// fprintf(stderr," MM_DM_DM: %.4f tpsDM_1=%.4f\n",
// exp(LogOfProduct(tpqMM_1,tpsDM_1,!IsUngapped())/scs),exp(tpsDM_1/scs));

            leftGD_1_DD = bestGD;   //subject

            tpsDD_1 = GetLogTransProbsAt( gaps_sec_, P_DD, m-2 );
            leftGD_1_DD += LogOfProduct( 0.0, tpsDD_1, !IsUngapped());
// fprintf(stderr," GD_DD: %.4f tpsDD_1=%.4f\n",
// exp(LogOfProduct(0.0,tpsDD_1,!IsUngapped())/scs),exp(tpsDD_1/scs));

// // //
// // //
            //process state stateMI
            //.....................
            if( upMI_MM_II < upMI_MM_MI ) {
                maxlog = upMI_MM_MI;
                ptr = dMM; //state (direction) for backtracing
            } else if( upMI_MM_MI < upMI_MM_II ) {
                        maxlog = upMI_MM_II;
                        ptr = dMI; //up direction
                    } else { //upMI_MM_MI == upMI_MM_II
                        maxlog = upMI_MM_MI;
                        ptr = dMM_MI; //diagonal or up
                    }
            upMI_MM_MI -= maxlog;
            upMI_MM_II -= maxlog;
            //save log Z_MI(n,m)
            bestMI = exp(ScBack(upMI_MM_MI,scs)) + exp(ScBack(upMI_MM_II,scs));
            if( 0.0 < bestMI )
                bestMI = maxlog + log( bestMI ) * scs;
            else
                bestMI = LOG_PROB_MIN;
            pointer[n][m][stateMI] = ptr;
            F[n][m][stateMI] = bestMI;

            //process state stateIM
            //.....................
            if( leftIM_II_MM < leftIM_MI_MM ) {
                maxlog = leftIM_MI_MM;
                ptr = dMM; //state (direction) for backtracing
            } else if( leftIM_MI_MM < leftIM_II_MM ) {
                        maxlog = leftIM_II_MM;
                        ptr = dIM; //left direction
                    } else { //leftIM_MI_MM == leftIM_II_MM
                        maxlog = leftIM_MI_MM;
                        ptr = dMM_IM; //diagonal or left
                    }
            leftIM_MI_MM -= maxlog;
            leftIM_II_MM -= maxlog;
            //save log Z_IM(n,m)
            bestIM = exp(ScBack(leftIM_MI_MM,scs)) + exp(ScBack(leftIM_II_MM,scs));
            if( 0.0 < bestIM )
                bestIM = maxlog + log( bestIM ) * scs;
            else
                bestIM = LOG_PROB_MIN;
            pointer[n][m][stateIM] = ptr;
            F[n][m][stateIM] = bestIM;

            //process state stateDG
            //.....................
            if( upDG_DD_1 < upDG_MD_1 ) {
                maxlog = upDG_MD_1;
                ptr = dMM; //backtracing state (direction)
            } else if( upDG_MD_1 < upDG_DD_1 ) {
                        maxlog = upDG_DD_1;
                        ptr = dDG; //up direction
                    } else { //upDG_MD_1 == upDG_DD_1
                        maxlog = upDG_MD_1;
                        ptr = dMM_DG; //diagonal or up
                    }
            upDG_MD_1 -= maxlog;
            upDG_DD_1 -= maxlog;
            //save log Z_DG(n,m)
            bestDG = exp(ScBack(upDG_MD_1,scs)) + exp(ScBack(upDG_DD_1,scs));
            if( 0.0 < bestDG )
                bestDG = maxlog + log( bestDG ) * scs;
            else
                bestDG = LOG_PROB_MIN;
            pointer[n][m][stateDG] = ptr;
            F[n][m][stateDG] = bestDG;

            //process state stateGD
            //.....................
            if( leftGD_1_DD < leftGD_1_MD ) {
                maxlog = leftGD_1_MD;
                ptr = dMM; //backtracing state (direction)
            } else if( leftGD_1_MD < leftGD_1_DD ) {
                        maxlog = leftGD_1_DD;
                        ptr = dGD; //left direction
                    } else { //leftGD_1_MD == leftGD_1_DD
                        maxlog = leftGD_1_MD;
                        ptr = dMM_GD; //diagonal or left
                    }
            leftGD_1_MD -= maxlog;
            leftGD_1_DD -= maxlog;
            //save log Z_GD(n,m)
            bestGD = exp(ScBack(leftGD_1_MD,scs)) + exp(ScBack(leftGD_1_DD,scs));
            if( 0.0 < bestGD )
                bestGD = maxlog + log( bestGD ) * scs;
            else
                bestGD = LOG_PROB_MIN;
            pointer[n][m][stateGD] = ptr;
            F[n][m][stateGD] = bestGD;

            //prepare for processing of state stateMM
            oMM += score;
            oMI += score;
            oIM += score;
            oDG += score;
            oGD += score;

            //process state stateMM
            //.....................
            if( oMI < oMM ) {//1.
                if( oIM < oMM ) {//2.
                    if( oDG < oMM ) {
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM_GD; }
                    }
                    else if( oMM < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMM
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM_DG; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM | dDG | dGD; }
                    }
                }
                else if( oMM < oIM ) {//2.
                    if( oDG < oIM ) {
                        if( oGD < oIM ) { maxlog = oIM; ptr = dIM; }
                        else if( oIM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oIM; ptr = dIM_GD; }
                    }
                    else if( oIM < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oIM
                        if( oGD < oIM ) { maxlog = oIM; ptr = dIM_DG; }
                        else if( oIM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oIM; ptr = dIM | dDG | dGD; }
                    }
                }
                else {//2. oIM == oMM /// ///
                    if( oDG < oMM ) {
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM_IM; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM | dIM | dGD; }
                    }
                    else if( oMM < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMM
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM | dIM | dDG; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM | dIM | dDG | dGD; }
                    }
                }
            }
            else if( oMM < oMI ) {//1.
                if( oIM < oMI ) {//2.
                    if( oDG < oMI ) {
                        if( oGD < oMI ) { maxlog = oMI; ptr = dMI; }
                        else if( oMI < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMI; ptr = dMI_GD; }
                    }
                    else if( oMI < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMI
                        if( oGD < oMI ) { maxlog = oMI; ptr = dMI_DG; }
                        else if( oMI < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMI; ptr = dMI | dDG | dGD; }
                    }
                }
                else if( oMI < oIM ) {//2.
                    if( oDG < oIM ) {
                        if( oGD < oIM ) { maxlog = oIM; ptr = dIM; }
                        else if( oIM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oIM; ptr = dIM_GD; }
                    }
                    else if( oIM < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oIM
                        if( oGD < oIM ) { maxlog = oIM; ptr = dIM_DG; }
                        else if( oIM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oIM; ptr = dIM | dDG | dGD; }
                    }
                }
                else {//2. oIM == oMI /// ///
                    if( oDG < oMI ) {
                        if( oGD < oMI ) { maxlog = oMI; ptr = dMI_IM; }
                        else if( oMI < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMI; ptr = dMI | dIM | dGD; }
                    }
                    else if( oMI < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMI
                        if( oGD < oMI ) { maxlog = oMI; ptr = dMI | dIM | dDG; }
                        else if( oMI < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMI; ptr = dMI | dIM | dDG | dGD; }
                    }
                }
            }
            else {//1. oMI == oMM ///
                if( oIM < oMM ) {//2.
                    if( oDG < oMM ) {
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM_MI; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM | dMI | dGD; }
                    }
                    else if( oMM < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMM
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM | dMI | dDG; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM | dMI | dDG | dGD; }
                    }
                }
                else if( oMM < oIM ) {//2.
                    if( oDG < oIM ) {
                        if( oGD < oIM ) { maxlog = oIM; ptr = dIM; }
                        else if( oIM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oIM; ptr = dIM_GD; }
                    }
                    else if( oIM < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oIM
                        if( oGD < oIM ) { maxlog = oIM; ptr = dIM_DG; }
                        else if( oIM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oIM; ptr = dIM | dDG | dGD; }
                    }
                }
                else {//2. oIM == oMM /// ///
                    if( oDG < oMM ) {
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM | dMI | dIM; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM | dMI | dIM | dGD; }
                    }
                    else if( oMM < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMM
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM | dMI | dIM | dDG; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM | dMI | dIM | dDG | dGD; }
                    }
                }
            }

            if(( GetGlobal())? 0: 1 )
                maxlog = SLC_MAX( maxlog, cllprob );
            oMM -= maxlog;
            oMI -= maxlog;
            oIM -= maxlog;
            oDG -= maxlog;
            oGD -= maxlog;
            //save log Z_MM(n,m)
            bestMM = exp(ScBack(oMM,scs)) + exp(ScBack(oMI,scs)) + exp(ScBack(oIM,scs)) + 
                     exp(ScBack(oDG,scs)) + exp(ScBack(oGD,scs));
            if(( GetGlobal())? 0: 1 )
                bestMM += exp(ScBack(cllprob-maxlog,scs));
            if( 0.0 < bestMM )
                bestMM = maxlog + log( bestMM ) * scs;
            else
                bestMM = LOG_PROB_MIN;
// fprintf(stderr," bestMM=%.6g oMM=%.6g oMI=%.6g oIM=%.6g oDG=%.6g oGD=%.6g maxlog=%g\n",
// exp(bestMM/scs),exp(oMM/scs),exp(oMI/scs),exp(oIM/scs),exp(oDG/scs),exp(oGD/scs),maxlog);

            pointer[n][m][stateMM] = ptr;
            F[n][m][stateMM] = bestMM;
        }
    }
}

// -------------------------------------------------------------------------
// BackwardDP: backward algorithm of dynamic programming to get backward
//  probabilities
//
void MAAlignment::BackwardDP()
{
    const TPAScore  cllprob = GetGlobal()? LOG_PROB_MIN: 0;//log of local probability

    TPAScore    bestMM, currentMM, oMM; //best, current, and operating scores for state stateMM
    TPAScore    bestMI, currentMI, oMI; //best, current, and operating scores for state stateMI
    TPAScore    bestIM, currentIM, oIM; //best, current, and operating scores for state stateIM
    TPAScore    bestDG, currentDG, oDG; //best, current, and operating scores for state stateDG
    TPAScore    bestGD, currentGD, oGD; //best, current, and operating scores for state stateGD

    TPAScore    upMI_MM_MI, upMI_MM_II;     //state MI: score of entering MI state (trans:MM-MI) and extending it (MM-II)
    TPAScore    leftIM_MI_MM, leftIM_II_MM; //state IM: score of entering IM state (trans:MI-MM) and extending it (II-MM)

    TPAScore    upDG_MD_1, upDG_DD_1;       //state DG: score of entering DG state (trans:MD-1) and extending it (DD-1)
    TPAScore    leftGD_1_MD, leftGD_1_DD;   //state GD: score of entering GD state (trans:1-MD) and extending it (1-DD)

    double      tpqMM, tpsMM;   //transition probabilities for query and subject
    double      tpqMI, tpsMI;
    double      tpqIM, tpsIM;
    double      tpqII, tpsII;
    double      tpqMD, tpsMD;
    double      tpqDM, tpsDM;
    double      tpqDD, tpsDD;

    double      maxlog; //maximum value of logarithms

    TPAScore    score;                  //profile score at certain positions
    int         ptr;                    //value of pointer
    int         n, m;
    int         qend = GetRestricted()? GetQueryEnd() + 1: GetQuerySize();
    int         qbeg = GetRestricted()? GetQueryBeg() + 1: 1;
    int         send = GetRestricted()? GetSbjctEnd() + 1: GetSubjectSize();
    int         sbeg = GetRestricted()? GetSbjctBeg() + 1: 1;

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();
    if( !scmatrix )
        throw myruntime_error( "MAAlignment: Unable to align profiles." );

    double          epa = GetExpectPerAlignment();
    double          dyno = 0.0;//dynamic offset
    double          dmul = 0.0;//dynamic multiplier
    double          dgrd = 1.0;
    double          delt = dgrd - gScBase;
    //s=y+x; s^n=ay+bx=a(s-x)+bx=as+cx, a<1,b>1
    if( SLC_SQRT_DBL_MIN < epa ) {
        if( epa < 1.0 ) {
              dyno = gOffset * epa;
              if( 0 < delt )
                  dgrd -= delt * ( 1.0 - epa );
              dmul = gModScMul * ( 1.0 - epa );
        }
        else  dyno = gOffset;
    }

    const double    pmult = scmatrix->GetMultiplier();
    double scs = scmatrix->GetAutoScalingFactor();
    if( gbUSEIMAGE )
        scs = 1.0;


    if( !GetGlobal()) {
        qend--;
        send--;
    }

// // const_cast<GapScheme&>(gaps_fst_).OutputGapScheme(); fprintf(stderr,"\n\n\n");
// // const_cast<GapScheme&>(gaps_sec_).OutputGapScheme();

    //iterate over query positions; 1+1 additional positions at beginning and end
    for( n = qend; qbeg <= n; n-- )
    {
        bestMM      = B_[n][send+1][stateMM];
        bestIM      = B_[n][send+1][stateIM];
        bestGD      = B_[n][send+1][stateGD];

        currentMM   = B_[n+1][send+1][stateMM];
        currentMI   = B_[n+1][send+1][stateMI];
        currentDG   = B_[n+1][send+1][stateDG];

// fprintf(stderr,"\n" );
        //iterate over subject positions
        for( m = send; sbeg <= m; m-- )
        {
            score = 0;
            if(( GetGlobal())? n < qend && m < send: 1 ) {
                if( gbUSEIMAGE ) {
                    if( scmatrix->GetUseModImage() && dmul )
                        score =(scmatrix->GetImageScore(m,n)*dgrd - dyno +
                                scmatrix->GetModImageScore(m,n)*dmul ) * pmult;
                    else
//                         score = scmatrix->GetImageScore( m, n ) * pmult;
                        score =(scmatrix->GetImageScore( m, n ) - dyno) * pmult;
                } else {
                    if( scmatrix->GetUseModImage())
                        score = scmatrix->GetModScore( m, n );
                    else
                        score = scmatrix->GetScore( m, n ); //B_'s offset is one greater
//                         score  = AutocorrScore( scmatrix, m, n ); //autocorrelation score
                }
            }
// fprintf(stderr,"\nBWD: %d %d: score: %.4f\n", n, m, exp(score/scs));

            //-----------------
            //state MM->MM
            tpqMM = GetLogTransProbsAt( gaps_fst_, P_MM, n-1 );//end state at the last pos.
            tpsMM = GetLogTransProbsAt( gaps_sec_, P_MM, m-1 );
            currentMM = B_[n+1][m+1][stateMM];
            oMM = currentMM;
            oMM += score;
            oMM += LogOfProduct( tpqMM, tpsMM, !IsUngapped());
// fprintf(stderr,"BWD: MM_MM_MM: %.4f tpqMM=%.4f tpsMM=%.4f\n", 
// exp(LogOfProduct(tpqMM,tpsMM,!IsUngapped())/scs),exp(tpqMM/scs),exp(tpsMM/scs));

            // MM->MI
            //B_ has EXTRA two rows and columns
            currentMI = B_[n+1][m][stateMI];
            oMI = currentMI; //query
            tpsMI = GetLogTransProbsAt( gaps_sec_, P_MI, m-1 );
            oMI += LogOfProduct( tpqMM, tpsMI, !IsUngapped());
// fprintf(stderr,"BWD: MI_MM_MI: %.4f tpsMI=%.4f\n", 
// exp(LogOfProduct(tpqMM,tpsMI,!IsUngapped())/scs),exp(tpsMI/scs));

            // MM->IM
            oIM = bestIM;   //subject
            tpqMI = GetLogTransProbsAt( gaps_fst_, P_MI, n-1 );
            oIM += LogOfProduct( tpqMI, tpsMM, !IsUngapped());
// fprintf(stderr,"BWD: IM_MI_MM: %.4f tpqMI=%.4f\n", 
// exp(LogOfProduct(tpqMI,tpsMM,!IsUngapped())/scs),exp(tpqMI/scs));

            // MM->DG
            currentDG = B_[n+1][m][stateDG];
            oDG = currentDG;  //query
            tpqMD = GetLogTransProbsAt( gaps_fst_, P_MD, n-1 );
            oDG += LogOfProduct( tpqMD, 0.0, !IsUngapped());
// fprintf(stderr,"BWD: DG_MD: %.4f tpqMD=%.4f\n", 
// exp(LogOfProduct(tpqMD,0.0,!IsUngapped())/scs),exp(tpqMD/scs));

            // MM->GD
            oGD = bestGD;    //subject
            tpsMD = GetLogTransProbsAt( gaps_sec_, P_MD, m-1 );
            oGD += LogOfProduct( 0.0, tpsMD, !IsUngapped());
// fprintf(stderr,"BWD: GD_MD: %.4f tpsMD=%.4f\n", 
// exp(LogOfProduct(0.0,tpsMD,!IsUngapped())/scs),exp(tpsMD/scs));


            //-----------------
            //state MI->MM
            upMI_MM_MI = currentMM; //query
            upMI_MM_MI += score;
            tpsIM = GetLogTransProbsAt( gaps_sec_, P_IM, m-1 );
            upMI_MM_MI += LogOfProduct( tpqMM, tpsIM, !IsUngapped());
// fprintf(stderr,"BWD: MI_MM_MI: %.4f tpsIM=%.4f\n", 
// exp(LogOfProduct(tpqMM,tpsIM,!IsUngapped())/scs),exp(tpsIM/scs));

            // MI->MI
            upMI_MM_II = currentMI; //query
            tpsII = GetLogTransProbsAt( gaps_sec_, P_II, m-1 );
            upMI_MM_II += LogOfProduct( tpqMM, tpsII, !IsUngapped());
// fprintf(stderr,"BWD: MI_MM_II: %.4f tpsII=%.4f\n", 
// exp(LogOfProduct(tpqMM,tpsII,!IsUngapped())/scs),exp(tpsII/scs));

            //-----------------
            //state IM->MM
            leftIM_MI_MM = currentMM;   //subject
            leftIM_MI_MM += score;
            tpqIM = GetLogTransProbsAt( gaps_fst_, P_IM, n-1 );
            leftIM_MI_MM += LogOfProduct( tpqIM, tpsMM, !IsUngapped());
// fprintf(stderr,"BWD: IM_MI_MM: %.4f tpqIM=%.4f\n", 
// exp(LogOfProduct(tpqIM,tpsMM,!IsUngapped())/scs),exp(tpqIM/scs));

            // IM->IM
            leftIM_II_MM = bestIM;   //subject
            tpqII = GetLogTransProbsAt( gaps_fst_, P_II, n-1 );
            leftIM_II_MM += LogOfProduct( tpqII, tpsMM, !IsUngapped());
// fprintf(stderr,"BWD: IM_II_MM: %.4f tpqII=%.4f\n", 
// exp(LogOfProduct(tpqII,tpsMM,!IsUngapped())/scs),exp(tpqII/scs));

            //-----------------
            //state DG->MM
            upDG_MD_1 = currentMM;  //query
            upDG_MD_1 += score;
            tpqDM = GetLogTransProbsAt( gaps_fst_, P_DM, n-1 );
            upDG_MD_1 += LogOfProduct( tpqDM, tpsMM, !IsUngapped());
// fprintf(stderr,"BWD: DG_MD: %.4f tpqDM=%.4f\n", 
// exp(LogOfProduct(tpqDM,tpsMM,!IsUngapped())/scs),exp(tpqDM/scs));

            // DG->DG
            upDG_DD_1 = currentDG;  //query
            tpqDD = GetLogTransProbsAt( gaps_fst_, P_DD, n-1 );
            upDG_DD_1 += LogOfProduct( tpqDD, 0.0, !IsUngapped());
// fprintf(stderr,"BWD: DG_DD: %.4f tpqDD=%.4f\n", 
// exp(LogOfProduct(tpqDD,0.0,!IsUngapped())/scs),exp(tpqDD/scs));

            //-----------------
            //state GD->MM
            leftGD_1_MD = currentMM;    //subject
            leftGD_1_MD += score;
            tpsDM = GetLogTransProbsAt( gaps_sec_, P_DM, m-1 );
            leftGD_1_MD += LogOfProduct( tpqMM, tpsDM, !IsUngapped());
// fprintf(stderr,"BWD: GD_MD: %.4f tpsDM=%.4f\n", 
// exp(LogOfProduct(tpqMM,tpsDM,!IsUngapped())/scs),exp(tpsDM/scs));

            // GD->GD
            leftGD_1_DD = bestGD;   //subject
            tpsDD = GetLogTransProbsAt( gaps_sec_, P_DD, m-1 );
            leftGD_1_DD += LogOfProduct( 0.0, tpsDD, !IsUngapped());
// fprintf(stderr,"BWD: GD_DD: %.4f tpsDD=%.4f\n", 
// exp(LogOfProduct(0.0,tpsDD,!IsUngapped())/scs),exp(tpsDD/scs));

// // //
// // //
            //process state stateMI
            //.....................
            if( upMI_MM_II < upMI_MM_MI ) {
                maxlog = upMI_MM_MI;
                ptr = dMM;
            } else if( upMI_MM_MI < upMI_MM_II ) {
                        maxlog = upMI_MM_II;
                        ptr = dMI; //up direction
                    } else { //upMI_MM_MI == upMI_MM_II
                        maxlog = upMI_MM_MI;
                        ptr = dMM_MI; //diagonal or up
                    }
            upMI_MM_MI -= maxlog;
            upMI_MM_II -= maxlog;
            //save log Z_MI(n,m)
            bestMI = exp(ScBack(upMI_MM_MI,scs)) + exp(ScBack(upMI_MM_II,scs));
            if( 0.0 < bestMI )
                bestMI = maxlog + log( bestMI ) * scs;
            else
                bestMI = LOG_PROB_MIN;

            Bptr_[n][m][stateMI] = ptr;
            B_[n][m][stateMI] = bestMI;

            //process state stateIM
            //.....................
            if( leftIM_II_MM < leftIM_MI_MM ) {
                maxlog = leftIM_MI_MM;
                ptr = dMM;
            } else if( leftIM_MI_MM < leftIM_II_MM ) {
                        maxlog = leftIM_II_MM;
                        ptr = dIM; //left direction
                    } else { //leftIM_MI_MM == leftIM_II_MM
                        maxlog = leftIM_MI_MM;
                        ptr = dMM_IM; //diagonal or left
                    }
            leftIM_MI_MM -= maxlog;
            leftIM_II_MM -= maxlog;
            //save log Z_IM(n,m)
            bestIM = exp(ScBack(leftIM_MI_MM,scs)) + exp(ScBack(leftIM_II_MM,scs));
            if( 0.0 < bestIM )
                bestIM = maxlog + log( bestIM ) * scs;
            else
                bestIM = LOG_PROB_MIN;

            Bptr_[n][m][stateIM] = ptr;
            B_[n][m][stateIM] = bestIM;

            //process state stateDG
            //.....................
            if( upDG_DD_1 < upDG_MD_1 ) {
                maxlog = upDG_MD_1;
                ptr = dMM;
            } else if( upDG_MD_1 < upDG_DD_1 ) {
                        maxlog = upDG_DD_1;
                        ptr = dDG; //up direction
                    } else { //upDG_MD_1 == upDG_DD_1
                        maxlog = upDG_MD_1;
                        ptr = dMM_DG; //diagonal or up
                    }
            upDG_MD_1 -= maxlog;
            upDG_DD_1 -= maxlog;
            //save log Z_DG(n,m)
            bestDG = exp(ScBack(upDG_MD_1,scs)) + exp(ScBack(upDG_DD_1,scs));
            if( 0.0 < bestDG )
                bestDG = maxlog + log( bestDG ) * scs;
            else
                bestDG = LOG_PROB_MIN;

            Bptr_[n][m][stateDG] = ptr;
            B_[n][m][stateDG] = bestDG;

            //process state stateGD
            //.....................
            if( leftGD_1_DD < leftGD_1_MD ) {
                maxlog = leftGD_1_MD;
                ptr = dMM;
            } else if( leftGD_1_MD < leftGD_1_DD ) {
                        maxlog = leftGD_1_DD;
                        ptr = dGD; //left direction
                    } else { //leftGD_1_MD == leftGD_1_DD
                        maxlog = leftGD_1_MD;
                        ptr = dMM_GD; //diagonal or left
                    }
            leftGD_1_MD -= maxlog;
            leftGD_1_DD -= maxlog;
            //save log Z_GD(n,m)
            bestGD = exp(ScBack(leftGD_1_MD,scs)) + exp(ScBack(leftGD_1_DD,scs));
            if( 0.0 < bestGD )
                bestGD = maxlog + log( bestGD ) * scs;
            else
                bestGD = LOG_PROB_MIN;

            Bptr_[n][m][stateGD] = ptr;
            B_[n][m][stateGD] = bestGD;

            //process state stateMM
            //.....................
            if( oMI < oMM ) {//1.
                if( oIM < oMM ) {//2.
                    if( oDG < oMM ) {
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM_GD; }
                    }
                    else if( oMM < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMM
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM_DG; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM | dDG | dGD; }
                    }
                }
                else if( oMM < oIM ) {//2.
                    if( oDG < oIM ) {
                        if( oGD < oIM ) { maxlog = oIM; ptr = dIM; }
                        else if( oIM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oIM; ptr = dIM_GD; }
                    }
                    else if( oIM < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oIM
                        if( oGD < oIM ) { maxlog = oIM; ptr = dIM_DG; }
                        else if( oIM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oIM; ptr = dIM | dDG | dGD; }
                    }
                }
                else {//2. oIM == oMM /// ///
                    if( oDG < oMM ) {
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM_IM; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM | dIM | dGD; }
                    }
                    else if( oMM < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMM
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM | dIM | dDG; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM | dIM | dDG | dGD; }
                    }
                }
            }
            else if( oMM < oMI ) {//1.
                if( oIM < oMI ) {//2.
                    if( oDG < oMI ) {
                        if( oGD < oMI ) { maxlog = oMI; ptr = dMI; }
                        else if( oMI < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMI; ptr = dMI_GD; }
                    }
                    else if( oMI < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMI
                        if( oGD < oMI ) { maxlog = oMI; ptr = dMI_DG; }
                        else if( oMI < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMI; ptr = dMI | dDG | dGD; }
                    }
                }
                else if( oMI < oIM ) {//2.
                    if( oDG < oIM ) {
                        if( oGD < oIM ) { maxlog = oIM; ptr = dIM; }
                        else if( oIM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oIM; ptr = dIM_GD; }
                    }
                    else if( oIM < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oIM
                        if( oGD < oIM ) { maxlog = oIM; ptr = dIM_DG; }
                        else if( oIM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oIM; ptr = dIM | dDG | dGD; }
                    }
                }
                else {//2. oIM == oMI /// ///
                    if( oDG < oMI ) {
                        if( oGD < oMI ) { maxlog = oMI; ptr = dMI_IM; }
                        else if( oMI < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMI; ptr = dMI | dIM | dGD; }
                    }
                    else if( oMI < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMI
                        if( oGD < oMI ) { maxlog = oMI; ptr = dMI | dIM | dDG; }
                        else if( oMI < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMI; ptr = dMI | dIM | dDG | dGD; }
                    }
                }
            }
            else {//1. oMI == oMM ///
                if( oIM < oMM ) {//2.
                    if( oDG < oMM ) {
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM_MI; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM | dMI | dGD; }
                    }
                    else if( oMM < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMM
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM | dMI | dDG; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM | dMI | dDG | dGD; }
                    }
                }
                else if( oMM < oIM ) {//2.
                    if( oDG < oIM ) {
                        if( oGD < oIM ) { maxlog = oIM; ptr = dIM; }
                        else if( oIM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oIM; ptr = dIM_GD; }
                    }
                    else if( oIM < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oIM
                        if( oGD < oIM ) { maxlog = oIM; ptr = dIM_DG; }
                        else if( oIM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oIM; ptr = dIM | dDG | dGD; }
                    }
                }
                else {//2. oIM == oMM /// ///
                    if( oDG < oMM ) {
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM | dMI | dIM; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM | dMI | dIM | dGD; }
                    }
                    else if( oMM < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMM
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM | dMI | dIM | dDG; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM | dMI | dIM | dDG | dGD; }
                    }
                }
            }

            if(( GetGlobal())? 0/*n < qend && m < send*/: 1 )
                maxlog = SLC_MAX( maxlog, cllprob );
            oMM -= maxlog;
            oMI -= maxlog;
            oIM -= maxlog;
            oDG -= maxlog;
            oGD -= maxlog;
            //save log Z_MM(n,m)
            bestMM = exp(ScBack(oMM,scs)) + exp(ScBack(oMI,scs)) + exp(ScBack(oIM,scs)) + 
                     exp(ScBack(oDG,scs)) + exp(ScBack(oGD,scs));
            if(( GetGlobal())? 0/*n < qend && m < send*/: 1 )
                bestMM += exp(ScBack(cllprob-maxlog,scs));

            if( 0.0 < bestMM )
                bestMM = maxlog + log( bestMM ) * scs;
            else
                bestMM = LOG_PROB_MIN;
// fprintf(stderr,"BWD: bestMM=%.6g oMM=%.6g oMI=%.6g oIM=%.6g oDG=%.6g oGD=%.6g\n",
// exp(bestMM/scs),exp(oMM/scs),exp(oMI/scs),exp(oIM/scs),exp(oDG/scs),exp(oGD/scs));

            Bptr_[n][m][stateMM] = ptr;
            B_[n][m][stateMM] = bestMM;
        }
    }
}





// =========================================================================
// ClearFMarginVertScaledAt: set vertical margin at the position of scaled 
//  version of F
//
void MAAlignment::ClearFMarginVertScaledAt( int n, int sbeg )
{
    TPAScore    currentMM;
    TPAScore    currentMI, bestMI;
    TPAScore    currentDG, bestDG;

    TPAScore    upMI_MM_MI, upMI_MM_II;
    TPAScore    upDG_MD_1, upDG_DD_1;

    double      tpqMM_1;//transition probabilities
    double      tpqMD_1;
    double      tpqDD_1;
    double      tpsMI;
    double      tpsII;

    if( n < 1 || GetQuerySize() < n || sbeg < 0 || GetSubjectSize() < sbeg )
        throw myruntime_error("MAAlignment: ClearFMarginVertScaledAt: Memory access error.");

    //-----------------
    //state MI: enter
    tpqMM_1 = gaps_fst_.GetTransProbsAt( P_MM, n-2 );//-1 beg. state

    currentMM = F[n-1][sbeg][stateMM];
    upMI_MM_MI = currentMM; //query

    tpsMI = gaps_sec_.GetTransProbsAt( P_MI, sbeg-1 );
    upMI_MM_MI *= ProbProduct( tpqMM_1, tpsMI, !IsUngapped());

    //state MI: extend
    currentMI = F[n-1][sbeg][stateMI];
    upMI_MM_II = currentMI; //query

    tpsII = gaps_sec_.GetTransProbsAt( P_II, sbeg-1 );
    upMI_MM_II *= ProbProduct( tpqMM_1, tpsII, !IsUngapped());

    //-----------------
    //state DG: enter
    upDG_MD_1 = currentMM;  //query

    tpqMD_1 = gaps_fst_.GetTransProbsAt( P_MD, n-2 );
    upDG_MD_1 *= ProbProduct( tpqMD_1, 1.0, !IsUngapped());

    //state DG: extend
    currentDG = F[n-1][sbeg][stateDG];
    upDG_DD_1 = currentDG;  //query

    tpqDD_1 = gaps_fst_.GetTransProbsAt( P_DD, n-2 );
    upDG_DD_1 *= ProbProduct( tpqDD_1, 1.0, !IsUngapped());

    //process state stateMI
    //.....................
    bestMI = scale_[n] * ( upMI_MM_MI + upMI_MM_II );
    F[n][sbeg][stateMI] = bestMI;

    //process state stateDG
    //.....................
    bestDG = scale_[n] * ( upDG_MD_1 + upDG_DD_1 );
    F[n][sbeg][stateDG] = bestDG;
// fprintf(stderr,"*FWD*  n=%d m=%d   bestMM=%g bestMI=%g bestIM=%g bestDG=%g bestGD=%g   scale=%g\n",
// n,sbeg,F[n][sbeg][stateMM],bestMI,F[n][sbeg][stateIM],bestDG,F[n][sbeg][stateGD],scale_[n]);
}

// -------------------------------------------------------------------------
// ClearFMarginHorzScaledAt: set horizontal margin at the position of scaled 
//  version of F
//
void MAAlignment::ClearFMarginHorzScaledAt( int qbeg, int m )
{
    TPAScore    bestMM;
    TPAScore    bestIM;
    TPAScore    bestGD;

    TPAScore    leftIM_MI_MM, leftIM_II_MM;
    TPAScore    leftGD_1_MD, leftGD_1_DD;

    double      tpsMM_1;//transition probabilities for query and subject
    double      tpsMD_1;
    double      tpsDD_1;
    double      tpqMI;
    double      tpqII;

    if( qbeg < 0 || GetQuerySize() < qbeg || m < 1 || GetSubjectSize() < m )
        throw myruntime_error("MAAlignment: ClearFMarginHorzScaledAt: Memory access error.");

    //-----------------
    //state IM: enter
    tpsMM_1 = gaps_sec_.GetTransProbsAt( P_MM, m-2 );

    bestMM = F[qbeg][m-1][stateMM];
    leftIM_MI_MM = bestMM;   //subject

    tpqMI = gaps_fst_.GetTransProbsAt( P_MI, qbeg-1 );
    leftIM_MI_MM *= ProbProduct( tpqMI, tpsMM_1, !IsUngapped());

    //state IM: extend
    bestIM = F[qbeg][m-1][stateIM];
    leftIM_II_MM = bestIM;   //subject

    tpqII = gaps_fst_.GetTransProbsAt( P_II, qbeg-1 );
    leftIM_II_MM *= ProbProduct( tpqII, tpsMM_1, !IsUngapped());

    //-----------------
    //state GD: enter
    leftGD_1_MD = bestMM;    //subject

    tpsMD_1 = gaps_sec_.GetTransProbsAt( P_MD, m-2 );
    leftGD_1_MD *= ProbProduct( 1.0, tpsMD_1, !IsUngapped());

    //state GD: extend
    bestGD = F[qbeg][m-1][stateGD];
    leftGD_1_DD = bestGD;   //subject

    tpsDD_1 = gaps_sec_.GetTransProbsAt( P_DD, m-2 );
    leftGD_1_DD *= ProbProduct( 1.0, tpsDD_1, !IsUngapped());

    //process state stateIM
    //.....................
    bestIM = ( leftIM_MI_MM + leftIM_II_MM );
    F[qbeg][m][stateIM] = bestIM;

    //process state stateGD
    //.....................
    bestGD = ( leftGD_1_MD + leftGD_1_DD );
    F[qbeg][m][stateGD] = bestGD;
// fprintf(stderr,"*FWD*  n=%d m=%d   bestMM=%g bestMI=%g bestIM=%g bestDG=%g bestGD=%g   scale=%g\n",
// qbeg,m,F[qbeg][m][stateMM],F[qbeg][m][stateMI],bestIM,F[qbeg][m][stateDG],bestGD,scale_[qbeg]);
}

// =========================================================================
// ForwardDPScaled: scaled version of forward algorithm of dynamic 
//  programming to get forward probabilities
//
void MAAlignment::ForwardDPScaled()
{
    TPAScore    cllprob = GetGlobal()? 0.0: 1.0;//local probability
    TPAScore    maxMM;
    TPAScore    bestMM, currentMM, oMM; //best, current, and operating scores for state stateMM
    TPAScore    bestMI, currentMI, oMI; //best, current, and operating scores for state stateMI
    TPAScore    bestIM, currentIM, oIM; //best, current, and operating scores for state stateIM
    TPAScore    bestDG, currentDG, oDG; //best, current, and operating scores for state stateDG
    TPAScore    bestGD, currentGD, oGD; //best, current, and operating scores for state stateGD

    TPAScore    upMI_MM_MI, upMI_MM_II;     //state MI: score of entering MI state (trans:MM-MI) and extending it (MM-II)
    TPAScore    leftIM_MI_MM, leftIM_II_MM; //state IM: score of entering IM state (trans:MI-MM) and extending it (II-MM)

    TPAScore    upDG_MD_1, upDG_DD_1;       //state DG: score of entering DG state (trans:MD-1) and extending it (DD-1)
    TPAScore    leftGD_1_MD, leftGD_1_DD;   //state GD: score of entering GD state (trans:1-MD) and extending it (1-DD)

    double      tpqMM_1, tpsMM_1;   //transition probabilities for query and subject
    double      tpqIM_1, tpsIM_1;
    double      tpqMD_1, tpsMD_1;
    double      tpqDM_1, tpsDM_1;
    double      tpqDD_1, tpsDD_1;
    double      tpqMI, tpsMI;
    double      tpqII, tpsII;

    TPAScore    score;                  //profile score at certain positions
    int         n, m;
    int         qbeg = GetRestricted()? GetQueryBeg() + 1: 1;
    int         qend = GetRestricted()? GetQueryEnd() + 1: GetQuerySize();
    int         sbeg = GetRestricted()? GetSbjctBeg() + 1: 1;
    int         send = GetRestricted()? GetSbjctEnd() + 1: GetSubjectSize();

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();
    if( !scmatrix )
        throw myruntime_error( "MAAlignment: Unable to align profiles." );

    double          epa = GetExpectPerAlignment();
    double          dyno = 0.0;//dynamic offset
    double          dmul = 0.0;//dynamic multiplier
    double          dgrd = 1.0;
    double          delt = dgrd - gScBase;
    //s=y+x; s^n=ay+bx=a(s-x)+bx=as+cx, a<1,b>1
    if( SLC_SQRT_DBL_MIN < epa ) {
        if( epa < 1.0 ) {
              dyno = gOffset * epa;
              if( 0 < delt )
                  dgrd -= delt * ( 1.0 - epa );
              dmul = gModScMul * ( 1.0 - epa );
        }
        else  dyno = gOffset;
    }

    const double    pmult = scmatrix->GetMultiplier();

// // const_cast<GapScheme&>(gaps_fst_).OutputGapScheme(); fprintf(stderr,"\n\n\n");
// // const_cast<GapScheme&>(gaps_sec_).OutputGapScheme();

    powscale_ = 1.0;

    //prepare margins of F
    if( GetGlobal())
        for( m = sbeg; m <= send; m++ )
            ClearFMarginHorzScaledAt( qbeg - 1, m );

    //start full probabilistic (global alignment) calculations
    //iterate over query positions
    for( n = qbeg; n <= qend; n++ )
    {
        maxMM       = 0;
        bestMM      = F[n][sbeg-1][stateMM] *= powscale_;
        bestIM      = F[n][sbeg-1][stateIM] *= powscale_;
        bestGD      = F[n][sbeg-1][stateGD] *= powscale_;

        currentMM   = F[n-1][sbeg-1][stateMM];
        currentMI   = F[n-1][sbeg-1][stateMI];
        currentDG   = F[n-1][sbeg-1][stateDG];

        if( GetGlobal())
            ClearFMarginVertScaledAt( n, sbeg - 1 );

// fprintf(stderr,"\n" );
        //iterate over subject positions
        for( m = sbeg; m <= send; m++ )
        {
            if( scmatrix->GetUseModImage() && dmul )
                score =(scmatrix->GetImageScore(m-1,n-1)*dgrd - dyno +
                        scmatrix->GetModImageScore(m-1,n-1)*dmul ) * pmult;
            else
//                 score = scmatrix->GetImageScore( m-1, n-1 ) * pmult; //score at (m,n)
                score =(scmatrix->GetImageScore( m-1, n-1 ) - dyno) * pmult; //score at (m,n)
//                 score = scmatrix->GetScore( m-1, n-1 ); //scaled score at (m,n)
//                 score  = AutocorrScore( scmatrix, m-1, n-1 ); //TEST:autocorrelation score
            score = exp( score );
// fprintf(stderr,"\n %d %d: score: %.4f\n", n, m, score );

            //-----------------
            //state MI: enter
            tpqMM_1 = gaps_fst_.GetTransProbsAt( P_MM, n-2 );//index value of -1 is ok: beg. state
            tpsMM_1 = gaps_sec_.GetTransProbsAt( P_MM, m-2 );
            oMM = currentMM;
            oMM *= ProbProduct( tpqMM_1, tpsMM_1, !IsUngapped());
// fprintf(stderr," MM_MM_MM: %.4f tpqMM_1=%.4f tpsMM_1=%.4f\n",
// ProbProduct(tpqMM_1,tpsMM_1,!IsUngapped()),tpqMM_1,tpsMM_1);

            //F has EXTRA row and column
            currentMM = F[n-1][m][stateMM];
            upMI_MM_MI = currentMM; //query

            tpsMI = gaps_sec_.GetTransProbsAt( P_MI, m-1 );
            upMI_MM_MI *= ProbProduct( tpqMM_1, tpsMI, !IsUngapped());
// fprintf(stderr," MI_MM_MI: %.4f tpsMI=%.4f\n", ProbProduct(tpqMM_1,tpsMI,!IsUngapped()),tpsMI);

            //state MI: extend
            tpsIM_1 = gaps_sec_.GetTransProbsAt( P_IM, m-2 );
            oMI = currentMI;
            oMI *= ProbProduct( tpqMM_1, tpsIM_1, !IsUngapped());
// fprintf(stderr," MM_MM_IM: %.4f tpsIM_1=%.4f\n", ProbProduct(tpqMM_1,tpsIM_1,!IsUngapped()),tpsIM_1);

            currentMI = F[n-1][m][stateMI];
            upMI_MM_II = currentMI; //query

            tpsII = gaps_sec_.GetTransProbsAt( P_II, m-1 );
            upMI_MM_II *= ProbProduct( tpqMM_1, tpsII, !IsUngapped());
// fprintf(stderr," MI_MM_II: %.4f tpsII=%.4f\n", ProbProduct(tpqMM_1,tpsII,!IsUngapped()),tpsII);

            //-----------------
            //state IM: enter
            leftIM_MI_MM = bestMM;   //subject

            tpqMI = gaps_fst_.GetTransProbsAt( P_MI, n-1 );
            leftIM_MI_MM *= ProbProduct( tpqMI, tpsMM_1, !IsUngapped());
// fprintf(stderr," IM_MI_MM: %.4f tpqMI=%.4f\n", ProbProduct(tpqMI,tpsMM_1,!IsUngapped()),tpqMI);

            //state IM: extend
            tpqIM_1 = gaps_fst_.GetTransProbsAt( P_IM, n-2 );
            oIM = F[n-1][m-1][stateIM];
            oIM *= ProbProduct( tpqIM_1, tpsMM_1, !IsUngapped());
// fprintf(stderr," MM_IM_MM: %.4f tpqIM_1=%.4f\n", ProbProduct(tpqIM_1,tpsMM_1,!IsUngapped()),tpqIM_1);

            leftIM_II_MM = bestIM;   //subject

            tpqII = gaps_fst_.GetTransProbsAt( P_II, n-1 );
            leftIM_II_MM *= ProbProduct( tpqII, tpsMM_1, !IsUngapped());
// fprintf(stderr," IM_II_MM: %.4f tpqII=%.4f\n", ProbProduct(tpqII,tpsMM_1,!IsUngapped()),tpqII);

            //-----------------
            //state DG: enter
            upDG_MD_1 = currentMM;  //query

            tpqMD_1 = gaps_fst_.GetTransProbsAt( P_MD, n-2 );
            upDG_MD_1 *= ProbProduct( tpqMD_1, 1.0, !IsUngapped());
// fprintf(stderr," DG_MD: %.4f tpqMD_1=%.4f\n", ProbProduct(tpqMD_1,1.0,!IsUngapped()),tpqMD_1);

            //state DG: extend
            tpqDM_1 = gaps_fst_.GetTransProbsAt( P_DM, n-2 );
            oDG = currentDG;
            oDG *= ProbProduct( tpqDM_1, tpsMM_1, !IsUngapped());
// fprintf(stderr," MM_DM_MM: %.4f tpqDM_1=%.4f\n", ProbProduct(tpqDM_1,tpsMM_1,!IsUngapped()),tpqDM_1);

            currentDG = F[n-1][m][stateDG];
            upDG_DD_1 = currentDG;  //query

            tpqDD_1 = gaps_fst_.GetTransProbsAt( P_DD, n-2 );
            upDG_DD_1 *= ProbProduct( tpqDD_1, 1.0, !IsUngapped());
// fprintf(stderr," DG_DD: %.4f tpqDD_1=%.4f\n", ProbProduct(tpqDD_1,1.0,!IsUngapped()),tpqDD_1);

            //-----------------
            //state GD: enter
            leftGD_1_MD = bestMM;    //subject

            tpsMD_1 = gaps_sec_.GetTransProbsAt( P_MD, m-2 );
            leftGD_1_MD *= ProbProduct( 1.0, tpsMD_1, !IsUngapped());
// fprintf(stderr," GD_MD: %.4f tpsMD_1=%.4f\n", ProbProduct(1.0,tpsMD_1,!IsUngapped()),tpsMD_1);

            //state GD: extend
            tpsDM_1 = gaps_sec_.GetTransProbsAt( P_DM, m-2 );
            oGD = F[n-1][m-1][stateGD];
            oGD *= ProbProduct( tpqMM_1, tpsDM_1, !IsUngapped());
// fprintf(stderr," MM_DM_DM: %.4f tpsDM_1=%.4f\n", ProbProduct(tpqMM_1,tpsDM_1,!IsUngapped()),tpsDM_1);

            leftGD_1_DD = bestGD;   //subject

            tpsDD_1 = gaps_sec_.GetTransProbsAt( P_DD, m-2 );
            leftGD_1_DD *= ProbProduct( 1.0, tpsDD_1, !IsUngapped());
// fprintf(stderr," GD_DD: %.4f tpsDD_1=%.4f\n", ProbProduct(1.0,tpsDD_1,!IsUngapped()),tpsDD_1);

// // //
// // //
            //process state stateMI
            //.....................
            //scale order is imp.
            bestMI = scale_[n] * ( upMI_MM_MI + upMI_MM_II );
            F[n][m][stateMI] = bestMI;

            //process state stateIM
            //.....................
            bestIM = ( leftIM_MI_MM + leftIM_II_MM );
            F[n][m][stateIM] = bestIM;

            //process state stateDG
            //.....................
            bestDG = scale_[n] * ( upDG_MD_1 + upDG_DD_1 );
            F[n][m][stateDG] = bestDG;

            //process state stateGD
            //.....................
            bestGD = ( leftGD_1_MD + leftGD_1_DD );
            F[n][m][stateGD] = bestGD;

            //process state stateMM
            //.....................
            //scale order is imp.
            bestMM = scale_[n] * ( cllprob + score * ( oMM + oMI + oIM + oDG + oGD ));
//             bestMM = score * scale_[n] * ( cllprob + oMM + oMI + oIM + oDG + oGD );
// fprintf(stderr," bestMM=%.6g oMM=%.6g oMI=%.6g oIM=%.6g oDG=%.6g oGD=%.6g\n",bestMM,oMM,oMI,oIM,oDG,oGD);
            F[n][m][stateMM] = bestMM;
            if( maxMM < bestMM )
                maxMM = bestMM;
// fprintf(stderr,"*FWD*  n=%d m=%d   bestMM=%g bestMI=%g bestIM=%g bestDG=%g bestGD=%g   score=%g scale=%g\n",
// n,m,bestMM,bestMI,bestIM,bestDG,bestGD,score,scale_[n]);
        }

        cllprob *= scale_[n];
        scale_[n+1] = 1.0;
        if( maxMM )
            scale_[n+1] /= maxMM;
        powscale_ *= scale_[n+1];
    }
}

// -------------------------------------------------------------------------
// BackwardDPScaled: scaled version of backward algorithm of dynamic 
//  programming to get backward probabilities
//
void MAAlignment::BackwardDPScaled()
{
    TPAScore    cllprob = GetGlobal()? 0.0: 1.0;//local probability

    TPAScore    bestMM, currentMM, oMM; //best, current, and operating scores for state stateMM
    TPAScore    bestMI, currentMI, oMI; //best, current, and operating scores for state stateMI
    TPAScore    bestIM, currentIM, oIM; //best, current, and operating scores for state stateIM
    TPAScore    bestDG, currentDG, oDG; //best, current, and operating scores for state stateDG
    TPAScore    bestGD, currentGD, oGD; //best, current, and operating scores for state stateGD

    TPAScore    upMI_MM_MI, upMI_MM_II;     //state MI: score of entering MI state (trans:MM-MI) and extending it (MM-II)
    TPAScore    leftIM_MI_MM, leftIM_II_MM; //state IM: score of entering IM state (trans:MI-MM) and extending it (II-MM)

    TPAScore    upDG_MD_1, upDG_DD_1;       //state DG: score of entering DG state (trans:MD-1) and extending it (DD-1)
    TPAScore    leftGD_1_MD, leftGD_1_DD;   //state GD: score of entering GD state (trans:1-MD) and extending it (1-DD)

    double      tpqMM, tpsMM;   //transition probabilities for query and subject
    double      tpqMI, tpsMI;
    double      tpqIM, tpsIM;
    double      tpqII, tpsII;
    double      tpqMD, tpsMD;
    double      tpqDM, tpsDM;
    double      tpqDD, tpsDD;

    TPAScore    score;                  //profile score at certain positions
    int         n, m;
    int         qend = GetRestricted()? GetQueryEnd() + 1: GetQuerySize();
    int         qbeg = GetRestricted()? GetQueryBeg() + 1: 1;
    int         send = GetRestricted()? GetSbjctEnd() + 1: GetSubjectSize();
    int         sbeg = GetRestricted()? GetSbjctBeg() + 1: 1;

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();
    if( !scmatrix )
        throw myruntime_error( "MAAlignment: Unable to align profiles." );

    double          epa = GetExpectPerAlignment();
    double          dyno = 0.0;//dynamic offset
    double          dmul = 0.0;//dynamic multiplier
    double          dgrd = 1.0;
    double          delt = dgrd - gScBase;
    //s=y+x; s^n=ay+bx=a(s-x)+bx=as+cx, a<1,b>1
    if( SLC_SQRT_DBL_MIN < epa ) {
        if( epa < 1.0 ) {
              dyno = gOffset * epa;
              if( 0 < delt )
                  dgrd -= delt * ( 1.0 - epa );
              dmul = gModScMul * ( 1.0 - epa );
        }
        else  dyno = gOffset;
    }

    const double    pmult = scmatrix->GetMultiplier();


    powscale_ = scale_[qend+1];
    if( GetGlobal())
        powscale_ = 1.0;

    if( !GetGlobal()) {
        cllprob *= scale_[qend+1];
        qend--;
        send--;
    }

    for( m = send+1; sbeg <= m; m-- ) {
        B_[qend+1][m][stateMM] *= powscale_;
        B_[qend+1][m][stateMI] *= powscale_;
        B_[qend+1][m][stateIM] *= powscale_;
        B_[qend+1][m][stateDG] *= powscale_;
        B_[qend+1][m][stateGD] *= powscale_;
    }

//  fprintf(stderr," powscale=%g\n",powscale_);

// // const_cast<GapScheme&>(gaps_fst_).OutputGapScheme(); fprintf(stderr,"\n\n\n");
// // const_cast<GapScheme&>(gaps_sec_).OutputGapScheme();

    //iterate over query positions; 1+1 additional positions at beginning and end
    for( n = qend; qbeg <= n; n-- )
    {
        powscale_ *= scale_[n+1];
        cllprob *= scale_[n+1];

        bestMM      = B_[n][send+1][stateMM] *= powscale_;
        bestIM      = B_[n][send+1][stateIM] *= powscale_;
        bestGD      = B_[n][send+1][stateGD] *= powscale_;

        currentMM   = B_[n+1][send+1][stateMM];
        currentMI   = B_[n+1][send+1][stateMI];
        currentDG   = B_[n+1][send+1][stateDG];

// fprintf(stderr,"\n" );
        //iterate over subject positions
        for( m = send; sbeg <= m; m-- )
        {
            score = 1;
            if(( GetGlobal())? n < qend && m < send: 1 ) {
//             if( n < qend && m < send ) {
                if( scmatrix->GetUseModImage() && dmul )
                    score =(scmatrix->GetImageScore(m,n)*dgrd - dyno +
                            scmatrix->GetModImageScore(m,n)*dmul ) * pmult;
                else
//                     score = scmatrix->GetImageScore( m, n ) * pmult; //B_'s offset is one greater
                    score =(scmatrix->GetImageScore( m, n ) - dyno ) * pmult; //B_'s offset is one greater
//                     score = scmatrix->GetScore( m, n ); //B_'s offset is one greater
//                     score  = AutocorrScore( scmatrix, m, n ); //autocorrelation score
                score = exp( score );
            }
// fprintf(stderr,"\nBWD: %d %d: score: %.4f\n", n, m, score );

            //-----------------
            //state MM->MM
            oMM = 0;
            tpqMM = gaps_fst_.GetTransProbsAt( P_MM, n-1 );//end state at the last pos.
            tpsMM = gaps_sec_.GetTransProbsAt( P_MM, m-1 );
            currentMM = B_[n+1][m+1][stateMM];
            if( score ) {
                oMM = currentMM;
                oMM *= score;
                oMM *= ProbProduct( tpqMM, tpsMM, !IsUngapped());
            }
// fprintf(stderr,"BWD: MM_MM_MM: %.4f tpqMM=%.4f tpsMM=%.4f\n", ProbProduct(tpqMM,tpsMM,!IsUngapped()),tpqMM,tpsMM);

            // MM->MI
            //B_ has EXTRA two rows and columns
            currentMI = B_[n+1][m][stateMI];
            oMI = currentMI; //query
            tpsMI = gaps_sec_.GetTransProbsAt( P_MI, m-1 );
            oMI *= ProbProduct( tpqMM, tpsMI, !IsUngapped());
// fprintf(stderr,"BWD: MI_MM_MI: %.4f tpsMI=%.4f\n", ProbProduct(tpqMM,tpsMI,!IsUngapped()),tpsMI);

            // MM->IM
            oIM = bestIM;   //subject
            tpqMI = gaps_fst_.GetTransProbsAt( P_MI, n-1 );
            oIM *= ProbProduct( tpqMI, tpsMM, !IsUngapped());
// fprintf(stderr,"BWD: IM_MI_MM: %.4f tpqMI=%.4f\n", ProbProduct(tpqMI,tpsMM,!IsUngapped()),tpqMI);

            // MM->DG
            currentDG = B_[n+1][m][stateDG];
            oDG = currentDG;  //query
            tpqMD = gaps_fst_.GetTransProbsAt( P_MD, n-1 );
            oDG *= ProbProduct( tpqMD, 1.0, !IsUngapped());
// fprintf(stderr,"BWD: DG_MD: %.4f tpqMD=%.4f\n", ProbProduct(tpqMD,1.0,!IsUngapped()),tpqMD);

            // MM->GD
            oGD = bestGD;    //subject
            tpsMD = gaps_sec_.GetTransProbsAt( P_MD, m-1 );
            oGD *= ProbProduct( 1.0, tpsMD, !IsUngapped());
// fprintf(stderr,"BWD: GD_MD: %.4f tpsMD=%.4f\n", ProbProduct(1.0,tpsMD,!IsUngapped()),tpsMD);


            //-----------------
            //state MI->MM
            upMI_MM_MI = 0;
            tpsIM = gaps_sec_.GetTransProbsAt( P_IM, m-1 );
            if( score ) {
                upMI_MM_MI = currentMM; //query
                upMI_MM_MI *= score;
                upMI_MM_MI *= ProbProduct( tpqMM, tpsIM, !IsUngapped());
            }
// fprintf(stderr,"BWD: MI_MM_MI: %.4f tpsIM=%.4f\n", ProbProduct(tpqMM,tpsIM,!IsUngapped()),tpsIM);

            // MI->MI
            upMI_MM_II = currentMI; //query
            tpsII = gaps_sec_.GetTransProbsAt( P_II, m-1 );
            upMI_MM_II *= ProbProduct( tpqMM, tpsII, !IsUngapped());
// fprintf(stderr,"BWD: MI_MM_II: %.4f tpsII=%.4f\n", ProbProduct(tpqMM,tpsII,!IsUngapped()),tpsII);

            //-----------------
            //state IM->MM
            leftIM_MI_MM = 0;
            tpqIM = gaps_fst_.GetTransProbsAt( P_IM, n-1 );
            if( score ) {
                leftIM_MI_MM = currentMM;   //subject
                leftIM_MI_MM *= score;
                leftIM_MI_MM *= ProbProduct( tpqIM, tpsMM, !IsUngapped());
            }
// fprintf(stderr,"BWD: IM_MI_MM: %.4f tpqIM=%.4f\n", ProbProduct(tpqIM,tpsMM,!IsUngapped()),tpqIM);

            // IM->IM
            leftIM_II_MM = bestIM;   //subject
            tpqII = gaps_fst_.GetTransProbsAt( P_II, n-1 );
            leftIM_II_MM *= ProbProduct( tpqII, tpsMM, !IsUngapped());
// fprintf(stderr,"BWD: IM_II_MM: %.4f tpqII=%.4f\n", ProbProduct(tpqII,tpsMM,!IsUngapped()),tpqII);

            //-----------------
            //state DG->MM
            upDG_MD_1 = 0;
            tpqDM = gaps_fst_.GetTransProbsAt( P_DM, n-1 );
            if( score ) {
                upDG_MD_1 = currentMM;  //query
                upDG_MD_1 *= score;
                upDG_MD_1 *= ProbProduct( tpqDM, tpsMM, !IsUngapped());
            }
// fprintf(stderr,"BWD: DG_MD: %.4f tpqDM=%.4f\n", ProbProduct(tpqDM,tpsMM,!IsUngapped()),tpqDM);

            // DG->DG
            upDG_DD_1 = currentDG;  //query
            tpqDD = gaps_fst_.GetTransProbsAt( P_DD, n-1 );
            upDG_DD_1 *= ProbProduct( tpqDD, 1.0, !IsUngapped());
// fprintf(stderr,"BWD: DG_DD: %.4f tpqDD=%.4f\n", ProbProduct(tpqDD,1.0,!IsUngapped()),tpqDD);

            //-----------------
            //state GD->MM
            leftGD_1_MD = 0;
            tpsDM = gaps_sec_.GetTransProbsAt( P_DM, m-1 );
            if( score ) {
                leftGD_1_MD = currentMM;    //subject
                leftGD_1_MD *= score;
                leftGD_1_MD *= ProbProduct( tpqMM, tpsDM, !IsUngapped());
            }
// fprintf(stderr,"BWD: GD_MD: %.4f tpsDM=%.4f\n", ProbProduct(tpqMM,tpsDM,!IsUngapped()),tpsDM);

            // GD->GD
            leftGD_1_DD = bestGD;   //subject
            tpsDD = gaps_sec_.GetTransProbsAt( P_DD, m-1 );
            leftGD_1_DD *= ProbProduct( 1.0, tpsDD, !IsUngapped());
// fprintf(stderr,"BWD: GD_DD: %.4f tpsDD=%.4f\n", ProbProduct(1.0,tpsDD,!IsUngapped()),tpsDD);

// // //
// // //
            //process state stateMI
            //.....................
            bestMI = scale_[n+1] * ( upMI_MM_MI + upMI_MM_II );
            B_[n][m][stateMI] = bestMI;

            //process state stateIM
            //.....................
            bestIM = ( scale_[n+1] * leftIM_MI_MM + leftIM_II_MM );
            B_[n][m][stateIM] = bestIM;

            //process state stateDG
            //.....................
            bestDG = scale_[n+1] * ( upDG_MD_1 + upDG_DD_1 );
            B_[n][m][stateDG] = bestDG;

            //process state stateGD
            //.....................
            bestGD = ( scale_[n+1] * leftGD_1_MD + leftGD_1_DD );
            B_[n][m][stateGD] = bestGD;

            //process state stateMM
            //.....................
            bestMM = cllprob + scale_[n+1]*(oMM + oMI + oDG) + (oIM + oGD);
            B_[n][m][stateMM] = bestMM;
// fprintf(stderr,"BWD: bestMM=%.6g oMM=%.6g oMI=%.6g oIM=%.6g oDG=%.6g oGD=%.6g\n",bestMM,oMM,oMI,oIM,oDG,oGD);
// fprintf(stderr,"*BWD*  n=%d m=%d   bestMM=%g bestMI=%g bestIM=%g bestDG=%g bestGD=%g   score=%g scale=%g\n",
// n,m,bestMM,bestMI,bestIM,bestDG,bestGD,score,scale_[n+1]);
        }
    }
}





// =========================================================================
// AccuracyDP: dynamic programming to fill up accuracy DP matrices
//
void MAAlignment::AccuracyDPObs()
{
    double      areg = GetAReg();
    TPAScore    bestMM, currentMM, oMM; //best, current, and operating scores for state stateMM
    TPAScore    bestMI, currentMI, oMI; //best, current, and operating scores for state stateMI
    TPAScore    bestIM, currentIM, oIM; //best, current, and operating scores for state stateIM
    TPAScore    bestDG, currentDG, oDG; //best, current, and operating scores for state stateDG
    TPAScore    bestGD, currentGD, oGD; //best, current, and operating scores for state stateGD

    TPAScore    upMI_MM_MI, upMI_MM_II;     //state MI: score of entering MI state (trans:MM-MI) and extending it (MM-II)
    TPAScore    leftIM_MI_MM, leftIM_II_MM; //state IM: score of entering IM state (trans:MI-MM) and extending it (II-MM)

    TPAScore    upDG_MD_1, upDG_DD_1;       //state DG: score of entering DG state (trans:MD-1) and extending it (DD-1)
    TPAScore    leftGD_1_MD, leftGD_1_DD;   //state GD: score of entering GD state (trans:1-MD) and extending it (1-DD)

    double      tpqMM_1, tpsMM_1;   //transition probabilities for query and subject
    double      tpqIM_1, tpsIM_1;
    double      tpqMD_1, tpsMD_1;
    double      tpqDM_1, tpsDM_1;
    double      tpqDD_1, tpsDD_1;
    double      tpqMI, tpsMI;
    double      tpqII, tpsII;

    TPAScore    accMM, accMI, accIM, accDG, accGD, accavg;
    bool        gppd = !IsUngapped();

    double      lFE = GetFE();//full probability of all alignments
    double      score;
    double      maxlog;//maximum value of logarithms
    int         ptr;//pointer value
    int         n, m;
    int         qbeg = GetRestricted()? GetQueryBeg() + 1: 1;
    int         qend = GetRestricted()? GetQueryEnd() + 1: GetQuerySize();
    int         sbeg = GetRestricted()? GetSbjctBeg() + 1: 1;
    int         send = GetRestricted()? GetSbjctEnd() + 1: GetSubjectSize();

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();
    if( !scmatrix )
        throw myruntime_error( "MAAlignment: Unable to align profiles." );

    double scs = scmatrix->GetAutoScalingFactor();
    if( gbUSEIMAGE )
        scs = 1.0;

#ifdef SCALEDFWDBWD
    if( lFE <= 0.0 )
        return;
#endif
#ifdef MAATESTPRINT
    //average sum of all state values per ALIGNMENT position should be 1!
    double sum = 0;
    double msum[send+1];
    memset( msum, 0, ( send + 1 )* sizeof( double ));
#endif

    //start maximum accuracy calculations
    //iterate over query positions
    for( n = qbeg; n <= qend; n++ )
    {
        bestMM      = A_[n][sbeg-1][stateMM];
        bestIM      = A_[n][sbeg-1][stateIM];
        bestGD      = A_[n][sbeg-1][stateGD];

        currentMM   = A_[n-1][sbeg-1][stateMM];
        currentMI   = A_[n-1][sbeg-1][stateMI];
        currentDG   = A_[n-1][sbeg-1][stateDG];

        //iterate over subject positions
        for( m = sbeg; m <= send; m++ )
        {
            if( scmatrix->GetUseModImage())
                score = scmatrix->GetModScore( m-1, n-1 );
            else
                score = scmatrix->GetScore( m-1, n-1 );

            //get transition probabilities
            tpqMM_1 = gaps_fst_.GetLogTransProbsAt( P_MM, n-2 );
            tpsMM_1 = gaps_sec_.GetLogTransProbsAt( P_MM, m-2 );
            tpsIM_1 = gaps_sec_.GetLogTransProbsAt( P_IM, m-2 );
            tpqIM_1 = gaps_fst_.GetLogTransProbsAt( P_IM, n-2 );
            tpqDM_1 = gaps_fst_.GetLogTransProbsAt( P_DM, n-2 );
            tpsDM_1 = gaps_sec_.GetLogTransProbsAt( P_DM, m-2 );

            tpsMI = gaps_sec_.GetLogTransProbsAt( P_MI, m-1 );
            tpsII = gaps_sec_.GetLogTransProbsAt( P_II, m-1 );

            tpqMI = gaps_fst_.GetLogTransProbsAt( P_MI, n-1 );
            tpqII = gaps_fst_.GetLogTransProbsAt( P_II, n-1 );

            tpqMD_1 = gaps_fst_.GetLogTransProbsAt( P_MD, n-2 );
            tpqDD_1 = gaps_fst_.GetLogTransProbsAt( P_DD, n-2 );

            tpsMD_1 = gaps_sec_.GetLogTransProbsAt( P_MD, m-2 );
            tpsDD_1 = gaps_sec_.GetLogTransProbsAt( P_DD, m-2 );

            //calculate posterior probabilities
#ifdef SCALEDFWDBWD
            accMM = F[n][m][stateMM] * B_[n][m][stateMM] / lFE;// - areg;
            accMI = F[n][m][stateMI] * B_[n][m][stateMI] / lFE;// - areg;
            accIM = F[n][m][stateIM] * B_[n][m][stateIM] / lFE;// - areg;
            accDG = F[n][m][stateDG] * B_[n][m][stateDG] / lFE;// - areg;
            accGD = F[n][m][stateGD] * B_[n][m][stateGD] / lFE;// - areg;
            accavg = (accMM+accMI+accIM+accDG+accGD)*.2;
            accMM -= areg;
            accavg = pow(accavg,.6);
            accMI -= accavg + areg;
            accIM -= accavg + areg;
            accDG -= accavg + areg;
            accGD -= accavg + areg;
            //similarly behaved: multiplication (summation in log space) of 
            //independent positional probabilities with a small weight for 
            //match probabilities (accMM/(from .25 to .5));
            //in that case, local alignment should be changed to semi-global and
            //MakeAlignmentPath accordingly modified.
#else
            accMM = exp(ScBack( F[n][m][stateMM] + B_[n][m][stateMM] - lFE, scs )) - areg;
            accMI = exp(ScBack( F[n][m][stateMI] + B_[n][m][stateMI] - lFE, scs )) - areg;
            accIM = exp(ScBack( F[n][m][stateIM] + B_[n][m][stateIM] - lFE, scs )) - areg;
            accDG = exp(ScBack( F[n][m][stateDG] + B_[n][m][stateDG] - lFE, scs )) - areg;
            accGD = exp(ScBack( F[n][m][stateGD] + B_[n][m][stateGD] - lFE, scs )) - areg;
#endif
#ifdef MAATESTPRINT
//             if(n==43||n==67||n==83)
//                 fprintf(stderr,"   m=%d accMM=%g accMI=%g accIM=%g accDG=%g accGD=%g\n",
//                         m, accMM, accMI, accIM, accDG, accGD );
            sum += accMM + accMI + accIM + accDG + accGD;
            msum[m] += accMM + accMI + accIM + accDG + accGD;
#endif

            //state MM
            oMM = currentMM + accMM;
            oMI = currentMI + accMM;
            oIM = A_[n-1][m-1][stateIM] + accMM;
            oDG = currentDG + accMM;
            oGD = A_[n-1][m-1][stateGD] + accMM;

            //-----------------
            //state MI: enter
            currentMM = A_[n-1][m][stateMM];
            upMI_MM_MI = currentMM; //query
            upMI_MM_MI += accMI;

            //state MI: extend
            currentMI = A_[n-1][m][stateMI];
            upMI_MM_II = currentMI; //query
            upMI_MM_II += accMI;

            //-----------------
            //state IM: enter
            leftIM_MI_MM = bestMM;   //subject
            leftIM_MI_MM += accIM;

            //state IM: extend
            leftIM_II_MM = bestIM;   //subject
            leftIM_II_MM += accIM;

            //-----------------
            //state DG: enter
            upDG_MD_1 = currentMM;  //query
            upDG_MD_1 += accDG;

            //state DG: extend
            currentDG = A_[n-1][m][stateDG];
            upDG_DD_1 = currentDG;  //query
            upDG_DD_1 += accDG;

            //-----------------
            //state GD: enter
            leftGD_1_MD = bestMM;    //subject
            leftGD_1_MD += accGD;

            //state GD: extend
            leftGD_1_DD = bestGD;   //subject
            leftGD_1_DD += accGD;

// // //
// // //
            //process state stateMI
            //.....................
            if( upMI_MM_II < upMI_MM_MI ) {
                maxlog = upMI_MM_MI;
                ptr = dMM; //state (direction) for backtracing
            } else if( upMI_MM_MI < upMI_MM_II ) {
                        maxlog = upMI_MM_II;
                        ptr = dMI; //up direction
                    } else { //upMI_MM_MI == upMI_MM_II
                        maxlog = upMI_MM_MI;
                        ptr = dMM_MI; //diagonal or up
                    }
            //save log A_MI(n,m)
            bestMI = maxlog;
            if( 0 < bestMI ) {
                Aptr_[n][m][stateMI] = ptr;
                A_[n][m][stateMI] = bestMI;
            } else {
                Aptr_[n][m][stateMI] = ptr = dNo;
                A_[n][m][stateMI] = bestMI = 0;
            }

            if( ptr & dMM )      S_[n][m][stateMI] = S_[n-1][m][stateMM] + LogOfProduct(tpqMM_1,tpsMI,gppd);
            else if( ptr & dMI ) S_[n][m][stateMI] = S_[n-1][m][stateMI] + LogOfProduct(tpqMM_1,tpsII,gppd);

            //process state stateIM
            //.....................
            if( leftIM_II_MM < leftIM_MI_MM ) {
                maxlog = leftIM_MI_MM;
                ptr = dMM; //state (direction) for backtracing
            } else if( leftIM_MI_MM < leftIM_II_MM ) {
                        maxlog = leftIM_II_MM;
                        ptr = dIM; //left direction
                    } else { //leftIM_MI_MM == leftIM_II_MM
                        maxlog = leftIM_MI_MM;
                        ptr = dMM_IM; //diagonal or left
                    }
            //save log A_IM(n,m)
            bestIM = maxlog;
            if( 0 < bestIM ) {
                Aptr_[n][m][stateIM] = ptr;
                A_[n][m][stateIM] = bestIM;
            } else {
                Aptr_[n][m][stateIM] = ptr = dNo;
                A_[n][m][stateIM] = bestIM = 0;
            }

            if( ptr & dMM )      S_[n][m][stateIM] = S_[n][m-1][stateMM] + LogOfProduct(tpqMI,tpsMM_1,gppd);
            else if( ptr & dIM ) S_[n][m][stateIM] = S_[n][m-1][stateIM] + LogOfProduct(tpqII,tpsMM_1,gppd);

            //process state stateDG
            //.....................
            if( upDG_DD_1 < upDG_MD_1 ) {
                maxlog = upDG_MD_1;
                ptr = dMM; //backtracing state (direction)
            } else if( upDG_MD_1 < upDG_DD_1 ) {
                        maxlog = upDG_DD_1;
                        ptr = dDG; //up direction
                    } else { //upDG_MD_1 == upDG_DD_1
                        maxlog = upDG_MD_1;
                        ptr = dMM_DG; //diagonal or up
                    }
            //save log A_DG(n,m)
            bestDG = maxlog;
            if( 0 < bestDG ) {
                Aptr_[n][m][stateDG] = ptr;
                A_[n][m][stateDG] = bestDG;
            } else {
                Aptr_[n][m][stateDG] = ptr = dNo;
                A_[n][m][stateDG] = bestDG = 0;
            }

            if( ptr & dMM )      S_[n][m][stateDG] = S_[n-1][m][stateMM] + LogOfProduct(tpqMD_1,0.0,gppd);
            else if( ptr & dDG ) S_[n][m][stateDG] = S_[n-1][m][stateDG] + LogOfProduct(tpqDD_1,0.0,gppd);

            //process state stateGD
            //.....................
            if( leftGD_1_DD < leftGD_1_MD ) {
                maxlog = leftGD_1_MD;
                ptr = dMM; //backtracing state (direction)
            } else if( leftGD_1_MD < leftGD_1_DD ) {
                        maxlog = leftGD_1_DD;
                        ptr = dGD; //left direction
                    } else { //leftGD_1_MD == leftGD_1_DD
                        maxlog = leftGD_1_MD;
                        ptr = dMM_GD; //diagonal or left
                    }
            //save log A_GD(n,m)
            bestGD = maxlog;
            if( 0 < bestGD ) {
                Aptr_[n][m][stateGD] = ptr;
                A_[n][m][stateGD] = bestGD;
            } else {
                Aptr_[n][m][stateGD] = ptr = dNo;
                A_[n][m][stateGD] = bestGD = 0;
            }

            if( ptr & dMM )      S_[n][m][stateGD] = S_[n][m-1][stateMM] + LogOfProduct(0.0,tpsMD_1,gppd);
            else if( ptr & dGD ) S_[n][m][stateGD] = S_[n][m-1][stateGD] + LogOfProduct(0.0,tpsDD_1,gppd);

            //process state stateMM
            //.....................
            if( oMI < oMM ) {//1.
                if( oIM < oMM ) {//2.
                    if( oDG < oMM ) {
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM_GD; }
                    }
                    else if( oMM < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMM
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM_DG; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM | dDG | dGD; }
                    }
                }
                else if( oMM < oIM ) {//2.
                    if( oDG < oIM ) {
                        if( oGD < oIM ) { maxlog = oIM; ptr = dIM; }
                        else if( oIM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oIM; ptr = dIM_GD; }
                    }
                    else if( oIM < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oIM
                        if( oGD < oIM ) { maxlog = oIM; ptr = dIM_DG; }
                        else if( oIM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oIM; ptr = dIM | dDG | dGD; }
                    }
                }
                else {//2. oIM == oMM /// ///
                    if( oDG < oMM ) {
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM_IM; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM | dIM | dGD; }
                    }
                    else if( oMM < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMM
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM | dIM | dDG; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM | dIM | dDG | dGD; }
                    }
                }
            }
            else if( oMM < oMI ) {//1.
                if( oIM < oMI ) {//2.
                    if( oDG < oMI ) {
                        if( oGD < oMI ) { maxlog = oMI; ptr = dMI; }
                        else if( oMI < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMI; ptr = dMI_GD; }
                    }
                    else if( oMI < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMI
                        if( oGD < oMI ) { maxlog = oMI; ptr = dMI_DG; }
                        else if( oMI < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMI; ptr = dMI | dDG | dGD; }
                    }
                }
                else if( oMI < oIM ) {//2.
                    if( oDG < oIM ) {
                        if( oGD < oIM ) { maxlog = oIM; ptr = dIM; }
                        else if( oIM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oIM; ptr = dIM_GD; }
                    }
                    else if( oIM < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oIM
                        if( oGD < oIM ) { maxlog = oIM; ptr = dIM_DG; }
                        else if( oIM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oIM; ptr = dIM | dDG | dGD; }
                    }
                }
                else {//2. oIM == oMI /// ///
                    if( oDG < oMI ) {
                        if( oGD < oMI ) { maxlog = oMI; ptr = dMI_IM; }
                        else if( oMI < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMI; ptr = dMI | dIM | dGD; }
                    }
                    else if( oMI < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMI
                        if( oGD < oMI ) { maxlog = oMI; ptr = dMI | dIM | dDG; }
                        else if( oMI < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMI; ptr = dMI | dIM | dDG | dGD; }
                    }
                }
            }
            else {//1. oMI == oMM ///
                if( oIM < oMM ) {//2.
                    if( oDG < oMM ) {
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM_MI; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM | dMI | dGD; }
                    }
                    else if( oMM < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMM
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM | dMI | dDG; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM | dMI | dDG | dGD; }
                    }
                }
                else if( oMM < oIM ) {//2.
                    if( oDG < oIM ) {
                        if( oGD < oIM ) { maxlog = oIM; ptr = dIM; }
                        else if( oIM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oIM; ptr = dIM_GD; }
                    }
                    else if( oIM < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oIM
                        if( oGD < oIM ) { maxlog = oIM; ptr = dIM_DG; }
                        else if( oIM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oIM; ptr = dIM | dDG | dGD; }
                    }
                }
                else {//2. oIM == oMM /// ///
                    if( oDG < oMM ) {
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM | dMI | dIM; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM | dMI | dIM | dGD; }
                    }
                    else if( oMM < oDG ) {
                        if( oGD < oDG ) { maxlog = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMM
                        if( oGD < oMM ) { maxlog = oMM; ptr = dMM | dMI | dIM | dDG; }
                        else if( oMM < oGD ) { maxlog = oGD; ptr = dGD; }
                        else { maxlog = oMM; ptr = dMM | dMI | dIM | dDG | dGD; }
                    }
                }
            }

            //save log A_MM(n,m)
            bestMM = maxlog;
            if( 0 < bestMM ) {
                Aptr_[n][m][stateMM] = ptr;
                A_[n][m][stateMM] = bestMM;
            } else {
                Aptr_[n][m][stateMM] = ptr = dNo;
                A_[n][m][stateMM] = bestMM = 0;
            }

            if( ptr != dNo ) {
              if( ptr & dMM )      S_[n][m][stateMM] = S_[n-1][m-1][stateMM] + LogOfProduct(tpqMM_1,tpsMM_1,gppd);
              else if( ptr & dMI ) S_[n][m][stateMM] = S_[n-1][m-1][stateMI] + LogOfProduct(tpqMM_1,tpsIM_1,gppd);
              else if( ptr & dIM ) S_[n][m][stateMM] = S_[n-1][m-1][stateIM] + LogOfProduct(tpqIM_1,tpsMM_1,gppd);
              else if( ptr & dDG ) S_[n][m][stateMM] = S_[n-1][m-1][stateDG] + LogOfProduct(tpqDM_1,tpsMM_1,gppd);
              else if( ptr & dGD ) S_[n][m][stateMM] = S_[n-1][m-1][stateGD] + LogOfProduct(tpqMM_1,tpsDM_1,gppd);
              S_[n][m][stateMM] += score;
            }
        }
#ifdef MAATESTPRINT
        fprintf( stderr," n=%d sum=%g\n", n, sum ); sum = 0;
#endif
    }
#ifdef MAATESTPRINT
    fprintf( stderr, "\n\n" );
    for( m = sbeg; m <= send; m++ ) fprintf( stderr," m=%d msum=%g\n", m, msum[m]);
#endif
}

// -------------------------------------------------------------------------
// AccuracyDP: dynamic programming to fill up accuracy DP matrices
//
void MAAlignment::AccuracyDP()
{
    double      areg = GetAReg();//0--extend the alignment maximally
    double      aext = GetARegExt();//0--no gap penalties
    TPAScore    bestMM, currentMM, oMM, oMI, oIM;

    double      tpqMM_1, tpsMM_1;   //transition probabilities for query and subject
    double      tpqIM_1, tpsIM_1;
    double      tpqMD_1, tpsMD_1;
    double      tpqDM_1, tpsDM_1;
    double      tpqDD_1, tpsDD_1;
    double      tpqMI, tpsMI;
    double      tpqII, tpsII;

    TPAScore    accMM;//, accMI, accIM, accDG, accGD;
    bool        gppd = !IsUngapped();

    double      lFE = GetFE();//full probability of all alignments
    double      score;
    double      maxlog;//maximum value of logarithms
    int         ptr, pptr;//pointer valueS
    int         margin;
    int         n, m, mb, me;
    int         qbeg = GetRestricted()? GetQueryBeg() + 1: 1;
    int         qend = GetRestricted()? GetQueryEnd() + 1: GetQuerySize();
    int         sbeg = GetRestricted()? GetSbjctBeg() + 1: 1;
    int         send = GetRestricted()? GetSbjctEnd() + 1: GetSubjectSize();

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();
    if( !scmatrix )
        throw myruntime_error( "MAAlignment: Unable to align profiles." );

    double          epa = GetExpectPerAlignment();
    double          dyno = 0.0;//dynamic offset
    double          dmul = 0.0;//dynamic multiplier
    double          dgrd = 1.0;
    double          delt = dgrd - gScBase;
    //s=y+x; s^n=ay+bx=a(s-x)+bx=as+cx, a<1,b>1
    if( SLC_SQRT_DBL_MIN < epa ) {
        if( epa < 1.0 ) {
              dyno = gOffset * epa;
              if( 0 < delt )
                  dgrd -= delt * ( 1.0 - epa );
              dmul = gModScMul * ( 1.0 - epa );
        }
        else  dyno = gOffset;
    }
    if( scmatrix->GetAutoScaling()) {
        if( dyno )
            dyno *= scmatrix->GetAutoScalingFactor();
        if( dgrd )
            dgrd *= scmatrix->GetAutoScalingFactor();
        if( dmul )
            dmul *= scmatrix->GetAutoScalingFactor();
    }

    double scs = scmatrix->GetAutoScalingFactor();
    if( gbUSEIMAGE )
        scs = 1.0;

#ifdef SCALEDFWDBWD
    if( lFE <= 0.0 )
        return;
#endif

    margin = ( int )rint( gdMARGINPRC * ( double )SLC_MIN( GetQuerySize(), GetSubjectSize()));
    if( gnMAXMARGIN < margin )
        margin = gnMAXMARGIN;

    //start maximum accuracy calculations
    //iterate over query positions
    for( n = qbeg; n <= qend; n++ )
    {
        bestMM      = A_[n][sbeg-1][stateMM];

        currentMM   = A_[n-1][sbeg-1][stateMM];

        //set boundaries of DP matrix triangles formed by both margins
        mb = 1; me = GetSubjectSize();
        if( GetQuerySize() < n + margin )
            mb += margin + n - GetQuerySize();
        if( n <= margin )
            me -= margin - n + 1;
        mb = SLC_MAX( mb, sbeg );
        me = SLC_MIN( me, send );

        //iterate over subject positions
        for( m = mb; m <= me; m++ )
        {
            if( scmatrix->GetUseModImage() && dmul )
                score = scmatrix->GetScore(m-1,n-1)*dgrd - dyno +
                        scmatrix->GetModScore(m-1,n-1)*dmul;
            else
//                 score = scmatrix->GetScore( m-1, n-1 );
                score = scmatrix->GetScore( m-1, n-1 ) - dyno;

            //get transition probabilities
            tpqMM_1 = gaps_fst_.GetLogTransProbsAt( P_MM, n-2 );
            tpsMM_1 = gaps_sec_.GetLogTransProbsAt( P_MM, m-2 );
            tpsIM_1 = gaps_sec_.GetLogTransProbsAt( P_IM, m-2 );
            tpqIM_1 = gaps_fst_.GetLogTransProbsAt( P_IM, n-2 );
            tpqDM_1 = gaps_fst_.GetLogTransProbsAt( P_DM, n-2 );
            tpsDM_1 = gaps_sec_.GetLogTransProbsAt( P_DM, m-2 );

            tpsMI = gaps_sec_.GetLogTransProbsAt( P_MI, m-1 );
            tpsII = gaps_sec_.GetLogTransProbsAt( P_II, m-1 );

            tpqMI = gaps_fst_.GetLogTransProbsAt( P_MI, n-1 );
            tpqII = gaps_fst_.GetLogTransProbsAt( P_II, n-1 );

            tpqMD_1 = gaps_fst_.GetLogTransProbsAt( P_MD, n-2 );
            tpqDD_1 = gaps_fst_.GetLogTransProbsAt( P_DD, n-2 );

            tpsMD_1 = gaps_sec_.GetLogTransProbsAt( P_MD, m-2 );
            tpsDD_1 = gaps_sec_.GetLogTransProbsAt( P_DD, m-2 );

            //calculate posterior probabilities
#ifdef SCALEDFWDBWD
            accMM = F[n][m][stateMM] * B_[n][m][stateMM] / lFE;
#else
            accMM = exp(ScBack( F[n][m][stateMM] + B_[n][m][stateMM] - lFE, scs ));
#endif
            //state MM
            oMM = currentMM + accMM - areg;
// fprintf(stderr,"*ACCURACY*  n=%d m=%d  oMM=%g accMM=%g\n",n,m,oMM,accMM);
// fprintf(stderr,"*ACCURACY*  n=%d m=%d  oMM=%g accMM=%g   F=%g B=%g FE=%g\n",
// n,m,oMM,accMM,F[n][m][stateMM],B_[n][m][stateMM],lFE);

            //-----------------
            //state UP
            currentMM = A_[n-1][m][stateMM];
            oMI = currentMM - aext;//query

            //-----------------
            //state IM
            oIM = bestMM - aext;//subject

// // //
// // //
            //process state stateMM
            //.....................
            if( oMI < oMM ) {
                if( oIM < oMM ) { maxlog = oMM; ptr = dMM; }
                else if( oMM < oIM ) { maxlog = oIM; ptr = dIM; }
                else { maxlog = oMM; ptr = dMM_IM; }
            }
            else if( oMM < oMI ) {
                if( oIM < oMI ) { maxlog = oMI; ptr = dMI; }
                else if( oMI < oIM ) { maxlog = oIM; ptr = dIM; }
                else { maxlog = oMI; ptr = dMI_IM; }
            }
            else {// oMI == oMM
                if( oIM < oMM ) { maxlog = oMM; ptr = dMM_MI; }
                else if( oMM < oIM ) { maxlog = oIM; ptr = dIM; }
                else { maxlog = oMM; ptr = dMM | dMI | dIM; }
            }

            //save A_MM(n,m)
            bestMM = maxlog;
            if( 0 < bestMM ) {
                Aptr_[n][m][stateMM] = ptr;
                A_[n][m][stateMM] = bestMM;
            } else {
                Aptr_[n][m][stateMM] = ptr = dNo;
                A_[n][m][stateMM] = bestMM = 0;
            }

            maxlog = 0.0;
            if( ptr & dMM ) {
                pptr = Aptr_[n-1][m-1][stateMM];
                if( pptr & dMM || pptr == dNo )
                    maxlog = LogOfProduct(tpqMM_1,tpsMM_1,gppd);
                else if( pptr & dMI )
                    maxlog = SLC_MAX( LogOfProduct(tpqMM_1,tpsIM_1,gppd), LogOfProduct(tpqDM_1,tpsMM_1,gppd));
                else if( pptr & dIM )
                    maxlog = SLC_MAX( LogOfProduct(tpqIM_1,tpsMM_1,gppd), LogOfProduct(tpqMM_1,tpsDM_1,gppd));
                S_[n][m][stateMM] = S_[n-1][m-1][stateMM] + maxlog;
                S_[n][m][stateMM] += score;
            }
            else if( ptr & dMI ) {
                pptr = Aptr_[n-1][m][stateMM];
                if( pptr & dMM || pptr == dNo )
                    maxlog = SLC_MAX( LogOfProduct(tpqMM_1,tpsMI,gppd), LogOfProduct(tpqMD_1,0.0,gppd));
                else if( pptr & dMI )
                    maxlog = SLC_MAX( LogOfProduct(tpqMM_1,tpsII,gppd), LogOfProduct(tpqDD_1,0.0,gppd));
                else //if( pptr & dIM )
                    maxlog = LOG_PROB_MIN;
//                     throw myruntime_error("MAAlignment: AccuracyDP: Invalid transition: IM->MI.");
                S_[n][m][stateMM] = S_[n-1][m][stateMM] + maxlog;
            }
            else if( ptr & dIM ) {
                pptr = Aptr_[n][m-1][stateMM];
                if( pptr & dMM || pptr == dNo )
                    maxlog = SLC_MAX( LogOfProduct(tpqMI,tpsMM_1,gppd), LogOfProduct(0.0,tpsMD_1,gppd));
                else if( pptr & dIM )
                    maxlog = SLC_MAX( LogOfProduct(tpqII,tpsMM_1,gppd), LogOfProduct(0.0,tpsDD_1,gppd));
                else //if( pptr & dMI )
                    maxlog = LOG_PROB_MIN;
//                     throw myruntime_error("MAAlignment: AccuracyDP: Invalid transition: MI->IM.");
                S_[n][m][stateMM] = S_[n][m-1][stateMM] + maxlog;
            }
        }
    }
}

// -------------------------------------------------------------------------
// AccuracyDP: dynamic programming to fill up accuracy DP matrices
//
void MAAlignment::AccuracyDP2()
{
    double      areg = GetAReg();
    double      aext = GetARegExt();
    TPAScore    bestMM, currentMM, oMM;
    TPAScore    bestMI, currentMI, oMI;
    TPAScore    bestIM, currentIM, oIM;

    TPAScore    upMI_MM_MI, upMI_MM_II;
    TPAScore    leftIM_MI_MM, leftIM_II_MM;

    double      tpqMM_1, tpsMM_1;
    double      tpqIM_1, tpsIM_1;
    double      tpqMD_1, tpsMD_1;
    double      tpqDM_1, tpsDM_1;
    double      tpqDD_1, tpsDD_1;
    double      tpqMI, tpsMI;
    double      tpqII, tpsII;

    TPAScore    accMM;//, accMI, accIM, accDG, accGD;
    bool        gppd = !IsUngapped();

    double      lFE = GetFE();//full probability of all alignments
    double      score;
    double      maxlog;//maximum value of logarithms
    int         ptr, pptr;//pointer valueS
    int         margin;
    int         n, m, mb, me;
    int         qbeg = GetRestricted()? GetQueryBeg() + 1: 1;
    int         qend = GetRestricted()? GetQueryEnd() + 1: GetQuerySize();
    int         sbeg = GetRestricted()? GetSbjctBeg() + 1: 1;
    int         send = GetRestricted()? GetSbjctEnd() + 1: GetSubjectSize();

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();
    if( !scmatrix )
        throw myruntime_error( "MAAlignment: Unable to align profiles." );

    double          epa = GetExpectPerAlignment();
    double          dyno = 0.0;//dynamic offset
    double          dmul = 0.0;//dynamic multiplier
    double          dgrd = 1.0;
    double          delt = dgrd - gScBase;
    //s=y+x; s^n=ay+bx=a(s-x)+bx=as+cx, a<1,b>1
    if( SLC_SQRT_DBL_MIN < epa ) {
        if( epa < 1.0 ) {
              dyno = gOffset * epa;
              if( 0 < delt )
                  dgrd -= delt * ( 1.0 - epa );
              dmul = gModScMul * ( 1.0 - epa );
        }
        else  dyno = gOffset;
    }
    if( scmatrix->GetAutoScaling()) {
        if( dyno )
            dyno *= scmatrix->GetAutoScalingFactor();
        if( dgrd )
            dgrd *= scmatrix->GetAutoScalingFactor();
        if( dmul )
            dmul *= scmatrix->GetAutoScalingFactor();
    }

    double scs = scmatrix->GetAutoScalingFactor();
    if( gbUSEIMAGE )
        scs = 1.0;

#ifdef SCALEDFWDBWD
    if( lFE <= 0.0 )
        return;
#endif

    margin = ( int )rint( gdMARGINPRC * ( double )SLC_MIN( GetQuerySize(), GetSubjectSize()));
    if( gnMAXMARGIN < margin )
        margin = gnMAXMARGIN;

    //start maximum accuracy calculations
    //iterate over query positions
    for( n = qbeg; n <= qend; n++ )
    {
        bestMM      = A_[n][sbeg-1][stateMM];

        currentMM   = A_[n-1][sbeg-1][stateMM];

        //set boundaries of DP matrix triangles formed by both margins
        mb = 1; me = GetSubjectSize();
        if( GetQuerySize() < n + margin )
            mb += margin + n - GetQuerySize();
        if( n <= margin )
            me -= margin - n + 1;
        mb = SLC_MAX( mb, sbeg );
        me = SLC_MIN( me, send );

        //iterate over subject positions
        for( m = mb; m <= me; m++ )
        {
            if( scmatrix->GetUseModImage() && dmul )
                score = scmatrix->GetScore(m-1,n-1)*dgrd - dyno +
                        scmatrix->GetModScore(m-1,n-1)*dmul;
            else
//                 score = scmatrix->GetScore( m-1, n-1 );
                score = scmatrix->GetScore( m-1, n-1 ) - dyno;

            //get transition probabilities
            tpqMM_1 = gaps_fst_.GetLogTransProbsAt( P_MM, n-2 );
            tpsMM_1 = gaps_sec_.GetLogTransProbsAt( P_MM, m-2 );
            tpsIM_1 = gaps_sec_.GetLogTransProbsAt( P_IM, m-2 );
            tpqIM_1 = gaps_fst_.GetLogTransProbsAt( P_IM, n-2 );
            tpqDM_1 = gaps_fst_.GetLogTransProbsAt( P_DM, n-2 );
            tpsDM_1 = gaps_sec_.GetLogTransProbsAt( P_DM, m-2 );

            tpsMI = gaps_sec_.GetLogTransProbsAt( P_MI, m-1 );
            tpsII = gaps_sec_.GetLogTransProbsAt( P_II, m-1 );

            tpqMI = gaps_fst_.GetLogTransProbsAt( P_MI, n-1 );
            tpqII = gaps_fst_.GetLogTransProbsAt( P_II, n-1 );

            tpqMD_1 = gaps_fst_.GetLogTransProbsAt( P_MD, n-2 );
            tpqDD_1 = gaps_fst_.GetLogTransProbsAt( P_DD, n-2 );

            tpsMD_1 = gaps_sec_.GetLogTransProbsAt( P_MD, m-2 );
            tpsDD_1 = gaps_sec_.GetLogTransProbsAt( P_DD, m-2 );

            //calculate posterior probabilities
#ifdef SCALEDFWDBWD
            accMM = F[n][m][stateMM] * B_[n][m][stateMM] / lFE;
#else
            accMM = exp(ScBack( F[n][m][stateMM] + B_[n][m][stateMM] - lFE, scs ));
#endif

            //-----------------
            //state MM
            oMM = currentMM + accMM - areg;
            oMI = currentMI + accMM - areg;
            oIM = A_[n-1][m-1][stateIM] + accMM - areg;

            //-----------------
            //state UP: enter
            currentMM = A_[n-1][m][stateMM];
            upMI_MM_MI = currentMM - aext;//query

            //state MI: extend
            currentMI = A_[n-1][m][stateMI];
            upMI_MM_II = currentMI - aext * 0.8;

            //-----------------
            //state IM: enter
            leftIM_MI_MM = bestMM - aext;//subject

            //state IM: extend
            leftIM_II_MM = bestIM - aext * 0.8;

// // //
// // //
            //process state stateMI
            //.....................
            if( upMI_MM_II < upMI_MM_MI ) {
                maxlog = upMI_MM_MI;
                ptr = dMM; //state for backtracing
            } else if( upMI_MM_MI < upMI_MM_II ) {
                        maxlog = upMI_MM_II;
                        ptr = dMI; //up direction
                    } else { //upMI_MM_MI == upMI_MM_II
                        maxlog = upMI_MM_MI;
                        ptr = dMM_MI; //diagonal or up
                    }
            //save log A_MI(n,m)
            bestMI = maxlog;
            if( 0 < bestMI ) {
                Aptr_[n][m][stateMI] = ptr;
                A_[n][m][stateMI] = bestMI;
            } else {
                Aptr_[n][m][stateMI] = ptr = dNo;
                A_[n][m][stateMI] = bestMI = 0;
            }


            if( ptr & dMM ) {
                maxlog = SLC_MAX( LogOfProduct(tpqMM_1,tpsMI,gppd), LogOfProduct(tpqMD_1,0.0,gppd));
                S_[n][m][stateMI] = S_[n-1][m][stateMM] + maxlog;
            }
            else if( ptr & dMI ) {
                maxlog = SLC_MAX( LogOfProduct(tpqMM_1,tpsII,gppd), LogOfProduct(tpqDD_1,0.0,gppd));
                S_[n][m][stateMI] = S_[n-1][m][stateMI] + maxlog;
            }

            //process state stateIM
            //.....................
            if( leftIM_II_MM < leftIM_MI_MM ) {
                maxlog = leftIM_MI_MM;
                ptr = dMM; //state for backtracing
            } else if( leftIM_MI_MM < leftIM_II_MM ) {
                        maxlog = leftIM_II_MM;
                        ptr = dIM; //left direction
                    } else { //leftIM_MI_MM == leftIM_II_MM
                        maxlog = leftIM_MI_MM;
                        ptr = dMM_IM; //diagonal or left
                    }
            //save log A_IM(n,m)
            bestIM = maxlog;
            if( 0 < bestIM ) {
                Aptr_[n][m][stateIM] = ptr;
                A_[n][m][stateIM] = bestIM;
            } else {
                Aptr_[n][m][stateIM] = ptr = dNo;
                A_[n][m][stateIM] = bestIM = 0;
            }

            if( ptr & dMM ) {
                maxlog = SLC_MAX( LogOfProduct(tpqMI,tpsMM_1,gppd), LogOfProduct(0.0,tpsMD_1,gppd));
                S_[n][m][stateIM] = S_[n][m-1][stateMM] + maxlog;
            }
            else if( ptr & dIM ) {
                maxlog = SLC_MAX( LogOfProduct(tpqII,tpsMM_1,gppd), LogOfProduct(0.0,tpsDD_1,gppd));
                S_[n][m][stateIM] = S_[n][m-1][stateIM] + maxlog;
            }

            //process state stateMM
            //.....................
            if( oMI < oMM ) {
                if( oIM < oMM ) { maxlog = oMM; ptr = dMM; }
                else if( oMM < oIM ) { maxlog = oIM; ptr = dIM; }
                else { maxlog = oMM; ptr = dMM_IM; }
            }
            else if( oMM < oMI ) {
                if( oIM < oMI ) { maxlog = oMI; ptr = dMI; }
                else if( oMI < oIM ) { maxlog = oIM; ptr = dIM; }
                else { maxlog = oMI; ptr = dMI_IM; }
            }
            else {// oMI == oMM
                if( oIM < oMM ) { maxlog = oMM; ptr = dMM_MI; }
                else if( oMM < oIM ) { maxlog = oIM; ptr = dIM; }
                else { maxlog = oMM; ptr = dMM | dMI | dIM; }
            }

            bestMM = maxlog;
            if( 0 < bestMM ) {
                Aptr_[n][m][stateMM] = ptr;
                A_[n][m][stateMM] = bestMM;
            } else {
                Aptr_[n][m][stateMM] = ptr = dNo;
                A_[n][m][stateMM] = bestMM = 0;
            }

            if( ptr != dNo ) {
                if( ptr & dMM )
                    S_[n][m][stateMM] = S_[n-1][m-1][stateMM] + LogOfProduct(tpqMM_1,tpsMM_1,gppd);
                else if( ptr & dMI ) {
                    maxlog = SLC_MAX( LogOfProduct(tpqMM_1,tpsIM_1,gppd), LogOfProduct(tpqDM_1,tpsMM_1,gppd));
                    S_[n][m][stateMM] = S_[n-1][m-1][stateMI] + maxlog;
                }
                else if( ptr & dIM ) {
                    maxlog = SLC_MAX( LogOfProduct(tpqIM_1,tpsMM_1,gppd), LogOfProduct(tpqMM_1,tpsDM_1,gppd));
                    S_[n][m][stateMM] = S_[n-1][m-1][stateIM] + maxlog;
                }
                S_[n][m][stateMM] += score;
            }
        }
    }
}

// -------------------------------------------------------------------------
// MakeAlignmentPath: make alignment path between two profiles given 
//  filled up dynamic programming matrices
//
void MAAlignment::MakeAlignmentPathObs()
{
    const bool  lSEARCH = true;
    TPAScore    lprob;
    TPAScore    score;
    int         row = 0;
    int         col = 0;
    int         prow, pcol;
    State       nextstate;//last state backtraced
    State       state;//state of back-tracing
    State       vstate;
    int         step, idns, psts, mtched, valid;
    int         qbeg = GetQueryBeg();
    int         sbeg = GetSbjctBeg();
    int         qend = GetQueryEnd();
    int         send = GetSbjctEnd();
    int         qsize = GetQuerySize();
    int         ssize = GetSubjectSize();
    int n, m;

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();
    if( !scmatrix )
        throw myruntime_error( "MAAlignment: MakeAlignmentPath: Null score matrix." );

    if( !A_ || !S_ )
        throw myruntime_error("MAAlignment: MakeAlignmentPath: Unable to draw alignment path.");

    if( !Aptr_ )
        throw myruntime_error("MAAlignment: MakeAlignmentPath: Unable to make alignment path.");

    if( qend < 0 || send < 0 || qsize < qend || ssize < send ||
        qend < qbeg || send < sbeg )
        throw myruntime_error("MAAlignment: MakeAlignmentPath: Invalid start-end positions.");

    nextstate = state = stateMM;
    prow = GetRestricted()? qend+1: qsize;
    pcol = GetRestricted()? send+1: ssize;
    lprob = 0;
    if( lSEARCH ) {
        //find nearest maximum value
        for( n = 1; n <= prow; n++ )
            for( m = 1; m <= pcol; m++ ) {
                if( lprob < A_[n][m][stateMM] ) {
                    lprob = A_[n][m][stateMM];
                    row = n;
                    col = m;
                }
            }
        if( !GetRestricted())
            qbeg = sbeg = 0;//beg. at (0,0)
    }
    else {
        row = prow;
        col = pcol;
    }
    score = idns = psts = mtched = 0;
    valid = -1;
    step = 0;
    nextstate = stateMM;

    while( qbeg <= row && sbeg <= col )
    {
        state = nextstate;
        if( noStates <= state )
            break;
        lprob = A_[row][col][state];
        prow = row;
        pcol = col;
        nextstate = GetState( Aptr_[row][col][state] );
        switch( nextstate ) {
            case stateMM:   row--; col--; 
                            if( qbeg <= row && sbeg <= col ) {
                                if( logo_fst_[row] == logo_sec_[col])
                                    idns++;
                                if( 0 < scmatrix->GetScore(col,row))
                                    psts++;
                                mtched++;
                            }
                            break;
            case stateMI:   row--; break;
            case stateIM:   col--; break;
            case stateDG:   row--; break;
            case stateGD:   col--; break;
            default: break;
        };

        if( state == stateMM ) {
            if( step <= 0 )
                score = S_[prow][pcol][state];
            vstate = state;
            valid = step;
        }
        if( 0 <= valid ) {
            path[step][nQuery] = prow;
            path[step][nSbjct] = pcol;
            step++;
        }
    }

    step = valid;
    if( 0 <= step ) {
        row = path[step][nQuery];
        col = path[step][nSbjct];

        if( row < 0 || col < 0 || noStates <= vstate )
            throw myruntime_error("MAAlignment: Invalid final values of path.");

        nextstate = GetState( Aptr_[row][col][vstate]);
        switch( vstate ) {
            case stateMM:   row--; col--; break;
            case stateMI:   row--; break;
            case stateIM:   col--; break;
            case stateDG:   row--; break;
            case stateGD:   col--; break;
            default: break;
        };
        if( nextstate < noStates && 0 <= row && 0 <= col )
            score -= S_[row][col][nextstate];
    }
    SetAlnScore( score );
    SetAlnSteps( SLC_MAX( 0, step ));
    SetNoMatched( mtched );
    SetNoIdentities( idns );
    SetNoPositives( psts );
}

// -------------------------------------------------------------------------
// MakeAlignmentPath: make alignment path between two profiles given 
//  filled up dynamic programming matrices
//
void MAAlignment::MakeAlignmentPath()
{
    const bool  lSEARCH = true;
    TPAScore    lprob;
    TPAScore    score;
    int         row = 0;
    int         col = 0;
    int         prow, pcol;
    State       nextstate;//last state backtraced
    State       state;//state of back-tracing
    State       vstate;
    int         step, idns, psts, mtched, valid;
    int         qbeg = GetQueryBeg();
    int         sbeg = GetSbjctBeg();
    int         qend = GetQueryEnd();
    int         send = GetSbjctEnd();
    int         qsize = GetQuerySize();
    int         ssize = GetSubjectSize();
    int n, m, s;

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();
    if( !scmatrix )
        throw myruntime_error( "MAAlignment: MakeAlignmentPath: Null score matrix." );

    if( !A_ || !S_ )
        throw myruntime_error("MAAlignment: MakeAlignmentPath: Unable to draw alignment path.");

    if( !Aptr_ )
        throw myruntime_error("MAAlignment: MakeAlignmentPath: Unable to make alignment path.");

    if( qend < 0 || send < 0 || qsize < qend || ssize < send ||
        qend < qbeg || send < sbeg )
        throw myruntime_error("MAAlignment: MakeAlignmentPath: Invalid start-end positions.");

    prow = GetRestricted()? qend+1: qsize;
    pcol = GetRestricted()? send+1: ssize;
    lprob = 0;
    if( lSEARCH ) {
        //find nearest maximum value
        for( n = 1; n <= prow; n++ )
            for( m = 1; m <= pcol; m++ ) {
                nextstate = GetState( Aptr_[n][m][stateMM]);
                if( lprob < A_[n][m][stateMM] && nextstate == stateMM ) {
                    lprob = A_[n][m][stateMM];
                    row = n;
                    col = m;
                }
            }
        if( !GetRestricted())
            qbeg = sbeg = 0;//beg. at (0,0)
    }
    else {
        row = prow;
        col = pcol;
    }
    score = 0;
    valid = -1;
    step = idns = psts = mtched = 0;
    nextstate = stateMM;
    //traceback through path and find a segment with 
    //match states of >50% at both ends
    while( qbeg <= row && sbeg <= col )
    {
        state = nextstate;
        if( noStates <= state )
            break;
        lprob = A_[row][col][stateMM/*state*/];
        prow = row;
        pcol = col;
        nextstate = GetState( Aptr_[row][col][stateMM/*state*/] );
        switch( nextstate ) {
            case stateMM:   row--; col--; 
                            if( qbeg <= row && sbeg <= col ) {
                                if( logo_fst_[row] == logo_sec_[col])
                                    idns++;
                                if( 0 < scmatrix->GetScore(col,row))
                                    psts++;
                                mtched++;
                            }
                            break;
            case stateMI:   row--; break;
            case stateIM:   col--; break;
            case stateDG:   row--; break;
            case stateGD:   col--; break;
            default: break;
        };

        if( state == stateMM ) {
            if( step <= 0 )
                score = S_[prow][pcol][stateMM/*state*/];
            vstate = state;
            valid = step;
        }
        if( 0 <= valid ) {
            path[step][nQuery] = prow;
            path[step][nSbjct] = pcol;
            step++;
        }
    }

    step = valid;
    if( 0 <= step ) {
        row = path[step][nQuery];
        col = path[step][nSbjct];

        if( row < 0 || col < 0 || noStates <= vstate )
            throw myruntime_error("MAAlignment: Invalid final values of path.");

        nextstate = GetState( Aptr_[row][col][stateMM/*vstate*/]);
        switch( vstate ) {
            case stateMM:   row--; col--; break;
            case stateMI:   row--; break;
            case stateIM:   col--; break;
            case stateDG:   row--; break;
            case stateGD:   col--; break;
            default: break;
        };
        if( nextstate < noStates && 0 <= row && 0 <= col )
            score -= S_[row][col][stateMM/*nextstate*/];
    }
    SetAlnScore( score );
    SetAlnSteps( SLC_MAX( 0, step ));
    SetNoMatched( mtched );
    SetNoIdentities( idns );
    SetNoPositives( psts );
}

// -------------------------------------------------------------------------
// MakeAlignmentPath2: make alignment path between two profiles from 
//  non-penalized probabilities
//
void MAAlignment::MakeAlignmentPath2()
{
    double      areg = GetAReg();//prob. threshold to regulate alignment length
    double      aext = GetARegExt();
    TPAScore    lprob, segcprb, maxsegcprb;//segment's cumulative probability
    TPAScore    score, begscore, maxscore;
    int         segbeg, segend;//segment's start and end positions
    int         maxsegbeg, maxsegend;
    int         row = 0;
    int         col = 0;
    int         prow, pcol;
    int         step, idns, psts, mtched;
    int         maxidns, maxpsts, maxmtched;
    int         (*locpath)[noInds] = NULL;//auxiliary alignment path
    mystring    errstr;

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();
    if( !scmatrix )
        throw myruntime_error( "MAAlignment: MakeAlignmentPath2: Null score matrix." );

    if( !A_ || !S_ )
        throw myruntime_error("MAAlignment: MakeAlignmentPath2: Unable to draw alignment path.");

    locpath = (int(*)[noInds])malloc(((querySize_+1)+(subjectSize_+1))*noInds*sizeof(int));
    if( !locpath )
        throw myruntime_error( "MAAlignment: MakeAlignmentPath2: Not enough memory." );
    memset( locpath, 0, (GetQuerySize()+GetSubjectSize())*noInds*sizeof(int));

    try {
        MakeAlignmentPath();
        //
        if( GetAReg() <= 0.)
            return;
        if( GetAlnSteps() < 1 )
            return;
        //
        idns = psts = mtched = 0;
        maxidns = maxpsts = maxmtched = 0;
        segcprb = 0;
        maxsegcprb = 0;
        score = begscore = maxscore = 0;
        segbeg = segend = -1;
        maxsegbeg = maxsegend = -1;
        for( step = GetAlnSteps()-1; 0 <= step; step-- ) {
            if( step == GetAlnSteps()-1 ||
               (path[step+1][nQuery] != path[step][nQuery] && 
                path[step+1][nSbjct] != path[step][nSbjct])) {
                row = path[step][nQuery];
                col = path[step][nSbjct];
                prow = row-1;
                pcol = col-1;
                if( prow < 0 || pcol < 0 )
                    continue;
                lprob = A_[row][col][stateMM] - A_[prow][pcol][stateMM];
                segcprb += lprob - areg;
                if( segcprb <= 0 ) {
                    score = begscore = 0;
                    segcprb = 0;
                    segbeg = segend = -1;
                    idns = psts = mtched = 0;
                }
                else {
                    if( segend < 0 ) {
                        segend = step;
                        begscore = S_[prow][pcol][stateMM];
                    }
                    segbeg = step;
                    score = S_[row][col][stateMM];
                    //stat.variables
                    if( logo_fst_[prow] == logo_sec_[pcol])
                        idns++;
                    if( 0 < scmatrix->GetScore(pcol,prow))
                        psts++;
                    mtched++; 
                    if( maxsegcprb < segcprb ) {
                        maxsegcprb = segcprb;
                        maxsegbeg = segbeg;
                        maxsegend = segend;
                        maxscore = score - begscore;
                        maxidns = idns;
                        maxpsts = psts;
                        maxmtched = mtched;
                    }
                }
            }
        }
        //save and rewrite path
        if( 0 <= maxsegbeg ) {
            for( step = 0; maxsegbeg <= maxsegend; step++, maxsegbeg++ ) {
                locpath[step][nQuery] = path[maxsegbeg][nQuery];
                locpath[step][nSbjct] = path[maxsegbeg][nSbjct];
            }
        }
        SetAlnSteps( SLC_MAX(0,step));
        for( step = 0; step < GetAlnSteps(); step++ ) {
            row = path[step][nQuery] = locpath[step][nQuery];
            col = path[step][nSbjct] = locpath[step][nSbjct];
        }
        SetAlnScore( maxscore );
        SetNoMatched( maxmtched );
        SetNoIdentities( maxidns );
        SetNoPositives( maxpsts );
    } catch( myexception const& ex ) {
        errstr = ex.what();
    }

    free(locpath);
    locpath = NULL;

    if( !errstr.empty())
        throw myruntime_error( errstr );
}

