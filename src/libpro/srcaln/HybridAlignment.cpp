/***************************************************************************
 *   Copyright (C) 2012 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "data.h"
#include "mystring.h"
#include "myexcept.h"
#include "ext/psl.h"
#include "HybridAlignment.h"


////////////////////////////////////////////////////////////////////////////
// CLASS HybridAlignment
//
// constructor
//
HybridAlignment::HybridAlignment(
        const FrequencyMatrix& freq_fst, const LogOddsMatrix& logo_fst, const GapScheme& gaps_fst,
        const FrequencyMatrix& freq_sec, const LogOddsMatrix& logo_sec, const GapScheme& gaps_sec,
        const AbstractScoreMatrix*  usc_system,
        bool ungapped
    )
:   ProfileAlignment( freq_fst, logo_fst, gaps_fst, 
                      freq_sec, logo_sec, gaps_sec, usc_system, ungapped )
{
}

// -------------------------------------------------------------------------
// default constructor is not allowed
//
HybridAlignment::HybridAlignment()
:   ProfileAlignment()
{
}

// -------------------------------------------------------------------------
// destructor: deallocate memory
//
HybridAlignment::~HybridAlignment()
{
}

// -------------------------------------------------------------------------
// ClearF: clear dynamic programming structure
//
void HybridAlignment::ClearF()
{
    int n, m;
    if( F == NULL )
        return;

    for( n = 0; n < GetQuerySize() + 1; n++ ) {
        if( F[n] == NULL )
            throw myruntime_error("HybridAlignment: Null vector of DPS.");
        for( m = 0; m < GetSubjectSize() + 1; m++ ) {
            if( n == 0 || m == 0 )
                F[n][m][stateMM] = 0;
            else
                F[n][m][stateMM] = LOG_PROB_MIN;
            F[n][m][stateMI] = LOG_PROB_MIN;
            F[n][m][stateIM] = LOG_PROB_MIN;
            F[n][m][stateDG] = LOG_PROB_MIN;
            F[n][m][stateGD] = LOG_PROB_MIN;
        }
    }
}

// -------------------------------------------------------------------------
// Run: prepare and start alignment of profiles
//
void HybridAlignment::Run()
{
    Clear();

    AlignProfiles();                //1.
    MakeAlignmentPath();            //2.
    PostProcess();
    SetFinalScore( GetAlnScore());
    ComputeStatistics();            //3.
// fprintf( stderr, "%12.4g\n", GetRawExpectation());
}

// -------------------------------------------------------------------------
// PostProcess: post process alignment
//
void HybridAlignment::PostProcess()
{
}

// -------------------------------------------------------------------------
// AlignProfiles: dynamic programming to align a pair of profiles;
//  semi-probabilistic approach to alignment
//
void HybridAlignment::AlignProfiles()
{
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
    double      tpqMI_1, tpsMI_1;

    double      deltasMM;   //delta score at state MM
    double      deltasMI;   //delta score at state MI
    double      deltasIM;   //delta score at state IM
    double      deltasDG;   //delta score at state DG
    double      deltasGD;   //delta score at state GD
    double      maxlog; //maximum value of logarithms
    double      logt1, logt2, logt3, logt4, logt5;

    TPAScore    score;                  //profile score at certain positions
    int         ptr;                    //value of pointer
    int         n, m;

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();

    if( !scmatrix )
        throw myruntime_error( "HybridAlignment: Unable to align profiles." );


// // const_cast<GapScheme&>(gaps_fst_).OutputGapScheme(); fprintf(stderr,"\n\n\n");
// // const_cast<GapScheme&>(gaps_sec_).OutputGapScheme();

    //iterate over query positions
    for( n = 1; n < GetQuerySize() + 1; n++ )
    {
        bestMM      = F[n][0][stateMM];
        bestIM      = F[n][0][stateIM];
        bestGD      = F[n][0][stateGD];

        currentMM   = F[n-1][0][stateMM];
        currentMI   = F[n-1][0][stateMI];
        currentDG   = F[n-1][0][stateDG];

// // fprintf(stderr,"\n" );
        //iterate over subject positions
        for( m = 1; m < GetSubjectSize() + 1; m++ )
        {
            score = scmatrix->GetImageScore( m-1, n-1 ); //profile score at positions m and n
//             score  = AutocorrScore( scmatrix, m-1, n-1 ); //autocorrelation score
// // fprintf(stderr,"\n %d %d: score: %.4f ", n, m, score/32.0 );

            //-----------------
            //state MI: enter
            tpqMM_1 = gaps_fst_.GetLogTransProbsAt( P_MM, n-2 );//index value of -1 is ok: beg. state
            tpsMM_1 = gaps_sec_.GetLogTransProbsAt( P_MM, m-2 );
            oMM = currentMM;
            oMM += LogOfProduct( tpqMM_1, tpsMM_1, !IsUngapped());
// // fprintf(stderr," MM_MM_MM: %.4f ", LogOfProduct( tpqMM_1, tpsMM_1, !IsUngapped())/32.0);

            //F has EXTRA row and column
            currentMM = F[n-1][m][stateMM];
            upMI_MM_MI = currentMM; //query

            tpsMI_1 = gaps_sec_.GetLogTransProbsAt( P_MI, m-2 );
            tpsMI = gaps_sec_.GetLogTransProbsAt( P_MI, m-1 );
            upMI_MM_MI += LogOfProduct( tpqMM_1, tpsMI, !IsUngapped());
// // fprintf(stderr," MI_MM_MI: %.4f ", LogOfProduct( tpqMM_1, tpsMI, !IsUngapped())/32.0);

            //state MI: extend
            tpsIM_1 = gaps_sec_.GetLogTransProbsAt( P_IM, m-2 );
            oMI = currentMI;
            oMI += LogOfProduct( tpqMM_1, tpsIM_1, !IsUngapped());
// fprintf(stderr," MM_MM_IM: %.4f ", LogOfProduct( tpqMM_1, tpsIM_1, !IsUngapped())/32.0);

            currentMI = F[n-1][m][stateMI];
            upMI_MM_II = currentMI; //query

            tpsII = gaps_sec_.GetLogTransProbsAt( P_II, m-1 );
            upMI_MM_II += LogOfProduct( tpqMM_1, tpsII, !IsUngapped());
// fprintf(stderr," MI_MM_II: %.4f ", LogOfProduct( tpqMM_1, tpsII, !IsUngapped())/32.0);

            //-----------------
            //state IM: enter
            leftIM_MI_MM = bestMM;   //subject

            tpqMI_1 = gaps_fst_.GetLogTransProbsAt( P_MI, n-2 );
            tpqMI = gaps_fst_.GetLogTransProbsAt( P_MI, n-1 );
            leftIM_MI_MM += LogOfProduct( tpqMI, tpsMM_1, !IsUngapped());
// // fprintf(stderr," IM_MI_MM: %.4f ", LogOfProduct( tpqMI, tpsMM_1, !IsUngapped())/32.0);

            //state IM: extend
            tpqIM_1 = gaps_fst_.GetLogTransProbsAt( P_IM, n-2 );
            oIM = F[n-1][m-1][stateIM];
            oIM += LogOfProduct( tpqIM_1, tpsMM_1, !IsUngapped());
// fprintf(stderr," MM_IM_MM: %.4f ", LogOfProduct( tpqIM_1, tpsMM_1, !IsUngapped())/32.0);

            leftIM_II_MM = bestIM;   //subject

            tpqII = gaps_fst_.GetLogTransProbsAt( P_II, n-1 );
            leftIM_II_MM += LogOfProduct( tpqII, tpsMM_1, !IsUngapped());
// fprintf(stderr," IM_II_MM: %.4f ", LogOfProduct( tpqII, tpsMM_1, !IsUngapped())/32.0);

            //-----------------
            //state DG: enter
            upDG_MD_1 = currentMM;  //query

            tpqMD_1 = gaps_fst_.GetLogTransProbsAt( P_MD, n-2 );
            upDG_MD_1 += LogOfProduct( tpqMD_1, 0.0, !IsUngapped());
// // fprintf(stderr," DG_MD: %.4f ", LogOfProduct( tpqMD_1, 0.0, !IsUngapped())/32.0);

            //state DG: extend
            tpqDM_1 = gaps_fst_.GetLogTransProbsAt( P_DM, n-2 );
            oDG = currentDG;
            oDG += LogOfProduct( tpqDM_1, tpsMM_1, !IsUngapped());
// fprintf(stderr," MM_DM_MM: %.4f ", LogOfProduct( tpqDM_1, tpsMM_1, !IsUngapped())/32.0);

            currentDG = F[n-1][m][stateDG];
            upDG_DD_1 = currentDG;  //query

            tpqDD_1 = gaps_fst_.GetLogTransProbsAt( P_DD, n-2 );
            upDG_DD_1 += LogOfProduct( tpqDD_1, 0.0, !IsUngapped());
// fprintf(stderr," DG_DD: %.4f ", LogOfProduct( tpqDD_1, 0.0, !IsUngapped())/32.0);

            //-----------------
            //state GD: enter
            leftGD_1_MD = bestMM;    //subject

            tpsMD_1 = gaps_sec_.GetLogTransProbsAt( P_MD, m-2 );
            leftGD_1_MD += LogOfProduct( 0.0, tpsMD_1, !IsUngapped());
// // fprintf(stderr," GD_MD: %.4f ", LogOfProduct( 0.0, tpsMD_1, !IsUngapped())/32.0);

            //state GD: extend
            tpsDM_1 = gaps_sec_.GetLogTransProbsAt( P_DM, m-2 );
            oGD = F[n-1][m-1][stateGD];
            oGD += LogOfProduct( tpqMM_1, tpsDM_1, !IsUngapped());
// fprintf(stderr," MM_DM_DM: %.4f ", LogOfProduct( tpqMM_1, tpsDM_1, !IsUngapped())/32.0);

            leftGD_1_DD = bestGD;   //subject

            tpsDD_1 = gaps_sec_.GetLogTransProbsAt( P_DD, m-2 );
            leftGD_1_DD += LogOfProduct( 0.0, tpsDD_1, !IsUngapped());
// fprintf(stderr," GD_DD: %.4f ", LogOfProduct( 0.0, tpsDD_1, !IsUngapped())/32.0);

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
            bestMI = exp( upMI_MM_MI ) + exp( upMI_MM_II );
            if( 0.0 < bestMI )
                bestMI = maxlog + log( bestMI );
            else
                bestMI = LOG_PROB_MIN;
            if( 0 < bestMI ) {
                pointer[n][m][stateMI] = ptr;
            } else {
                ///bestMI = 0;
                pointer[n][m][stateMI] = dNo;
            }
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
            bestIM = exp( leftIM_MI_MM ) + exp( leftIM_II_MM );
            if( 0.0 < bestIM )
                bestIM = maxlog + log( bestIM );
            else
                bestIM = LOG_PROB_MIN;
            if( 0 < bestIM ) {
                pointer[n][m][stateIM] = ptr;
            } else {
                ///bestIM = 0;
                pointer[n][m][stateIM] = dNo;
            }
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
            bestDG = exp( upDG_MD_1 ) + exp( upDG_DD_1 );
            if( 0.0 < bestDG )
                bestDG = maxlog + log( bestDG );
            else
                bestDG = LOG_PROB_MIN;
            if( 0 < bestDG ) {
                pointer[n][m][stateDG] = ptr;
            } else {
                ///bestDG = 0;
                pointer[n][m][stateDG] = dNo;
            }
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
            bestGD = exp( leftGD_1_MD ) + exp( leftGD_1_DD );
            if( 0.0 < bestGD )
                bestGD = maxlog + log( bestGD );
            else
                bestGD = LOG_PROB_MIN;
            if( 0 < bestGD ) {
                pointer[n][m][stateGD] = ptr;
            } else {
                ///bestGD = 0;
                pointer[n][m][stateGD] = dNo;
            }
            F[n][m][stateGD] = bestGD;

            //prepare for processing of state stateMM
            //calculate and adjust scores for each state
            logt1 = -tpqMM_1 - tpsMM_1;
            logt2 = tpqMD_1 + logt1;
            logt3 = tpsMD_1 + logt1;
            logt4 = tpqMI_1 - tpqMM_1;
            logt5 = tpsMI_1 - tpsMM_1;
            maxlog = SLC_MAX5( logt1, logt2, logt3, logt4, logt5 );
            logt1 -= maxlog;
            logt2 -= maxlog;
            logt3 -= maxlog;
            logt4 -= maxlog;
            logt5 -= maxlog;
            deltasMM = exp(logt1) - exp(logt2) - exp(logt3) - exp(logt4) - exp(logt5);
            if( 0.0 < deltasMM )
                deltasMM = maxlog + log( deltasMM );
            else
                deltasMM = LOG_PROB_MIN;

            deltasMI = -tpqMM_1 - tpsIM_1;
            //make sure argument of logarithm is positive
            if( tpqMM_1 < 0.0 || SLC_LOG_DBL_MIN < tpsIM_1 )
                deltasMI += log( 1.0 - exp( tpqMM_1 ) + exp( tpqMM_1 + tpsIM_1 ));

            deltasIM = -tpqIM_1 - tpsMM_1;
            //make sure argument of logarithm is positive
            if( tpsMM_1 < 0.0 || SLC_LOG_DBL_MIN < tpqIM_1 )
                deltasIM += log( 1.0 - exp( tpsMM_1 ) + exp( tpqIM_1 + tpsMM_1 ));;

            deltasDG = -tpsMM_1;
            deltasGD = -tpqMM_1;

            oMM += deltasMM + score;
            oMI += deltasMI + score;
            oIM += deltasIM + score;
            oDG += deltasDG + score;
            oGD += deltasGD + score;

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
            if( maxlog < 0.0 ) {
                maxlog = 0.0;
//                 ptr = dNo;
            }
            oMM -= maxlog;
            oMI -= maxlog;
            oIM -= maxlog;
            oDG -= maxlog;
            oGD -= maxlog;
            //save log Z_MM(n,m)
            bestMM = exp(-maxlog) + exp(oMM) + exp(oMI) + exp(oIM) + exp(oDG) + exp(oGD);
            if( 0.0 < bestMM )
                bestMM = maxlog + log( bestMM );
            else
                bestMM = LOG_PROB_MIN;

            if( 0 < bestMM ) {
                pointer[n][m][stateMM] = ptr;
            } else {
                ///bestMM = 0;
                pointer[n][m][stateMM] = dNo;
            }
            F[n][m][stateMM] = bestMM;
       }
    }
}

// -------------------------------------------------------------------------
// GetMax: return the first maximum entry found in array of values
//
void HybridAlignment::GetMax( TPAScore* zss, int noes, TPAScore* maxv, int* indx )
{
    int n;
    if( noes <= 0 || maxv == NULL )
        return;
    *maxv = zss[0];
    if( indx )
        *indx = 0;
    for( n = 1; n < noes; n++ ) {
        if( *maxv < zss[n]) {
            *maxv = zss[n];
            if( indx )
                *indx = n;
        }
    }
}

// -------------------------------------------------------------------------
// MakeAlignmentPath: make alignment path between two profiles given 
//  filled up dynamic programming matrices
//
void HybridAlignment::MakeAlignmentPath()
{
    TPAScore    score, maxscore, maxlog;//score, maximum scores
    TPAScore    zss[noStates];          //DPM values for each state
    int         row = 0;                //row index of maximum score
    int         col = 0;                //column index of maximum score
    int         prow, pcol;
    State       nextstate = stateMM;    //next state of back-tracing
    State       state = nextstate;      //state of back-tracing
    int         step;                   //index for variable path
    int n, m;

    if( !F )
        throw myruntime_error("HybridAlignment: Null DP matrices.");

    if( !path )
        throw myruntime_error("HybridAlignment: Null path structure.");

    maxscore = score = maxlog = 0;
    //find the most distant maximum value 
    for( n = GetQuerySize(); n; n-- )
        for( m = GetSubjectSize(); m; m-- ) {
            zss[stateMM] = F[n][m][stateMM];
            zss[stateMI] = F[n][m][stateMI];
            zss[stateIM] = F[n][m][stateIM];
            zss[stateDG] = F[n][m][stateDG];
            zss[stateGD] = F[n][m][stateGD];
            GetMax( zss, noStates, &maxlog, ( int* )&state );
            if( noStates <= state )
                throw myruntime_error("HybridAlignment: MakeAlignmentPath failed.");

//             if( state != stateMM )
//                 continue;
            zss[stateMM] -= maxlog;
            zss[stateMI] -= maxlog;
            zss[stateIM] -= maxlog;
            zss[stateDG] -= maxlog;
            zss[stateGD] -= maxlog;
            //save log Z_MM(n,m)
            score = exp(zss[stateMM]) + exp(zss[stateMI]) + 
                    exp(zss[stateIM]) + exp(zss[stateDG]) + exp(zss[stateGD]);
            if( 0.0 < score )
                score = maxlog + log( score );
            else
                score = LOG_PROB_MIN;

            if( maxscore < score ) {
                nextstate = state;
                maxscore = score;
                row = n;
                col = m;
            }
        }

    if( maxscore <= 0 )
        return;

    step = 0;
    while( 0 < row && 0 < col ) {
        state = nextstate;
        prow = row;
        pcol = col;
        zss[stateMM] = F[row][col][stateMM];
        zss[stateMI] = F[row][col][stateMI];
        zss[stateIM] = F[row][col][stateIM];
        zss[stateDG] = F[row][col][stateDG];
        zss[stateGD] = F[row][col][stateGD];
        GetMax( zss, noStates, &maxlog, ( int* )&nextstate );
        if( noStates <= nextstate )
            throw myruntime_error("HybridAlignment: MakeAlignmentPath: Invalid state value.");
        switch( state ) {
            case stateMM:   row--; col--; break;
            case stateMI:   row--; break;
            case stateIM:   col--; break;
            case stateDG:   row--; break;
            case stateGD:   col--; break;
            default:        break;
        };
        if( row < 0 || col < 0 )
            throw myruntime_error("HybridAlignment: MakeAlignmentPath: Invalid indices.");
        if( 0 < step || state == stateMM ) {
            path[step][nQuery] = prow;
            path[step][nSbjct] = pcol;
            step++;
        }
        if( nextstate == stateMM && zss[stateMM] <= 0.0 )
            break;
    }

    SetAlnScore( maxscore );
    SetAlnSteps( step - 1 ); //substract: beginning points to cell of score 0
}
