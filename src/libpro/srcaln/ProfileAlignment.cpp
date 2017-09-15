/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
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
#include "libpro/srcpro/SEGProfile.h"
#include "libpro/srcaln/MAAlignment.h"
#include "ProfileAlignment.h"

//global file variables
int     gnMAXMARGIN = 60;   //maximum margin length at both ends of profile
double  gdMARGINPRC = 0.34; //percentage of profile length, set as margin length

double  ProfileAlignment::s_information_thrsh = 0.0;    //information content threshold
double  ProfileAlignment::s_segdistance_thrsh = 0.0;    //SEG distance threshold

FrequencyMatrix dummyfreq;
LogOddsMatrix   dummylogo;
GapScheme       dummygaps;

////////////////////////////////////////////////////////////////////////////
// CLASS ProfileAlignment
//
// Constructor
//  Parameters are frequency, log-odds matrices, and gap cost schemes
//
ProfileAlignment::ProfileAlignment(
        const FrequencyMatrix& freq_fst, const LogOddsMatrix& logo_fst, const GapScheme& gaps_fst,
        const FrequencyMatrix& freq_sec, const LogOddsMatrix& logo_sec, const GapScheme& gaps_sec,
        const AbstractScoreMatrix*  usc_system,
        bool ungapped
    )
:
    model( StatModel::PC_DIST_HOM_0_2 ),

    F( NULL ),
    pointer( NULL ),
    path( NULL ),

    querySize_( 0 ),
    subjectSize_( 0 ),

    freq_fst_( freq_fst ),
    logo_fst_( logo_fst ),
    gaps_fst_( gaps_fst ),

    freq_sec_( freq_sec ),
    logo_sec_( logo_sec ),
    gaps_sec_( gaps_sec ),

    proprobs_( NULL ),
    scoreSystem( usc_system ),

    bitscore( -1.0 ),
    ref_expectation( -1.0 ),
    expectperaln_( -1.0 ),
    raw_expect_( -1.0 ),
    expectation( -1.0 ),
    p_value( -1.0 ),

    finscore_( 0.0 ),
    alnscore_( 0.0 ),
    alnsteps_( 0 ),
    mtcsteps_( 0 ),
    idnsteps_( 0 ),
    pststeps_( 0 ),

    ungapped_( ungapped )
{
    if( !logo_fst_.IsCompatible( freq_fst_ ) || !logo_fst_.IsCompatible( gaps_fst_ ) ||
        !logo_sec_.IsCompatible( freq_sec_ ) || !logo_sec_.IsCompatible( gaps_sec_ ) )
            throw myruntime_error( "Profile matrices are incompatible." );

    querySize_ = freq_fst_.GetColumns();
    subjectSize_ = freq_sec_.GetColumns();

    Initialize();
}

// -------------------------------------------------------------------------
// default constructor is invalid for initialization
//
ProfileAlignment::ProfileAlignment()
:   model( StatModel::PC_DIST_HOM_0_2 ),

    F( 0 ),
    pointer( 0 ),
    path( 0 ),

    querySize_( 0 ),
    subjectSize_( 0 ),

    freq_fst_( dummyfreq ),
    logo_fst_( dummylogo ),
    gaps_fst_( dummygaps ),

    freq_sec_( dummyfreq ),
    logo_sec_( dummylogo ),
    gaps_sec_( dummygaps ),

    proprobs_( NULL ),
    scoreSystem( NULL ),

    bitscore( -1.0 ),
    ref_expectation( -1.0 ),
    expectperaln_( -1.0 ),
    raw_expect_( -1.0 ),
    expectation( -1.0 ),
    p_value( 1.0 ),

    finscore_( 0.0 ),
    alnscore_( 0.0 ),
    alnsteps_( 0 ),
    mtcsteps_( 0 ),
    idnsteps_( 0 ),
    pststeps_( 0 )
{
    throw myruntime_error( "ProfileAlignment: Default construction is not allowed." );
}

// -------------------------------------------------------------------------
// Destructor
//
ProfileAlignment::~ProfileAlignment()
{
    Destroy();
}

// -------------------------------------------------------------------------
// Destroy: deallocate memory used by class structures
//
void ProfileAlignment::Destroy()
{
    int n;

    if( path ) {
        free( path );
        path = NULL;
    }
    if( F ) {
        for( n = 0; n < querySize_ + 1; n++ )
            if( F[n])
                free( F[n]);
        free( F );
        F = NULL;
    }
    if( pointer ) {
        for( n = 0; n < querySize_ + 1; n++ )
            if( pointer[n])
                free( pointer[n]);
        free( pointer );
        pointer = NULL;
    }
}

// -------------------------------------------------------------------------
// Initialize: create and initialize related structures and variables
//
void ProfileAlignment::Initialize()
{
    int n;

    if(  querySize_ < 1 || MAXCOLUMNS < querySize_ ||
         subjectSize_ < 1 || MAXCOLUMNS < subjectSize_ )
            throw myruntime_error( "Invalid profile size." );

    Destroy();

    //reserve with extra size of 1...
    F       = (TPAScore(**)[noStates] )malloc( sizeof( TPAScore* )*( querySize_ + 1 ));
    pointer = (     int(**)[noStates] )malloc( sizeof( int*      )*( querySize_ + 1 ));
    path    = (     int( *)[noInds] )malloc( sizeof( int )*(( querySize_+1 ) + ( subjectSize_+1 )) * noInds );

    if( !F || !pointer || !path )
        throw myruntime_error( "ProfileAlignment: Not enough memory." );

    for( n = 0; n < querySize_ + 1; n++ )
    {
        F[n]        = (TPAScore(*)[noStates] )malloc( sizeof( TPAScore )*( subjectSize_ + 1 ) * noStates );
        pointer[n]  = ( int    (*)[noStates] )malloc( sizeof( int      )*( subjectSize_ + 1 ) * noStates );

        if( !F[n] || !pointer[n] )
            throw myruntime_error( "ProfileAlignment: Not enough memory." );
    }
}


// -------------------------------------------------------------------------
// ClearPath: clear alignment path structure
//
void ProfileAlignment::ClearPath()
{
    if( path )
        memset( path, 0, sizeof( int )*( GetQuerySize() + GetSubjectSize()) * noInds );
}

// -------------------------------------------------------------------------
// ClearF: clear dynamic programming structure
//
void ProfileAlignment::ClearF()
{
    int n;
    if( F )
        for( n = 0; n < GetQuerySize() + 1; n++ )
            if( F[n] )
                memset( F[n], 0,  sizeof( TPAScore )*( GetSubjectSize() + 1 ) * noStates );
}

// -------------------------------------------------------------------------
// ClearPointer: clear back-tracing pointer structure
//
void ProfileAlignment::ClearPointer()
{
    int n;
    if( pointer )
        for( n = 0; n < GetQuerySize() + 1; n++ )
            if( pointer[n] )
                memset( pointer[n], 0, sizeof( int )*( GetSubjectSize() + 1 ) * noStates );
}

// -------------------------------------------------------------------------
// Clear: clear information structures
//
void ProfileAlignment::Clear()
{
    ClearPath();
    ClearF();
    ClearPointer();

    SetBitScore( -1.0 );
//     SetReferenceExpectation( -1.0 );
//     SetExpectPerAlignment( -1.0 );
//     SetRawExpectation( -1.0 );
//     SetExpectation( -1.0 );
//     SetPvalue( 1.0 );

    finscore_ = 0.0;
    SetAlnScore( 0.0 );
    SetAlnSteps( 0 );
    SetNoMatched( 0 );
}

// -------------------------------------------------------------------------
// ComputeStatistics: compute e-value and p-value
// -------------------------------------------------------------------------

void ProfileAlignment::ComputeStatistics()
{
    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();
    if( scmatrix == NULL )
        throw myruntime_error("ProfileAlignment: ComputeStatistics: Null score matrix.");

    double  ref_expect = -1.0;
    double  pure_expect = -1.0;
    double  pair_expect = -1.0;
    double  locbit_score = -1.0;
    double  expectation = scmatrix->ComputeExpectation( GetScore(), &ref_expect, &pure_expect, &pair_expect, &locbit_score );

    SetBitScore( locbit_score );
    SetRawExpectation( pure_expect );
    SetExpectPerAlignment( pair_expect );

    if( expectation < 0.0 ) {
        SetReferenceExpectation( ref_expect );
        return;
    }

    SetReferenceExpectation( expectation ); //set this value equal to the available original value
    SetExpectation( expectation );
    SetPvalue( ComputePvalue( GetExpectation()));
}

// -------------------------------------------------------------------------
// Run: steps to produce alignment between profiles
// -------------------------------------------------------------------------

void ProfileAlignment::Run( size_t nopros )
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
// MARealign: realign by initiating maximum accuracy calculations
//
void ProfileAlignment::MARealign()
{
    int         nosteps, mtched, idns, psts, s;
    int         qbeg, sbeg;
    int         qend, send;
    int         qsize = GetQuerySize();
    int         ssize = GetSubjectSize();
    double      epa;//expect per pair alignment

    nosteps = GetAlnSteps();
    mtched = GetNoMatched();
    idns = GetNoIdentities();
    psts = GetNoPositives();
    if( nosteps < 1 )
        //may be cases when scores are barely positive but 
        //transitions probabilities outweigh them resulting in 
        //negative final scores; so no. steps may be 0
        return;
    if( qsize + ssize < nosteps )
        throw myruntime_error("ProfileAlignment: MARealign: Invalid no. steps.");

    //determine alignment boundaries of this alignment
    qbeg = path[nosteps-1][nQuery]-1;
    sbeg = path[nosteps-1][nSbjct]-1;
    qend = path[0][nQuery]-1;
    send = path[0][nSbjct]-1;

    if( qend < 0 || send < 0 || qend < qbeg || send < sbeg ||
        qsize <= qend || ssize <= send )
        throw myruntime_error("ProfileAlignment: MARealign: Invalid alignment boundaries.");

    epa = GetExpectPerAlignment();
    if( SLC_SQRT_DBL_MIN < epa && epa < SLC_SQRT_DBL_MAX )
        epa *= exp( 26.0 );

    MAAlignment maa( freq_fst_, logo_fst_, gaps_fst_, 
                     freq_sec_, logo_sec_, gaps_sec_,
                     qbeg, qend, sbeg, send,
                     GetScoreSystem());
    maa.SetBitScore( GetBitScore());
    maa.SetRawExpectation( GetRawExpectation());
    maa.SetExpectPerAlignment( epa );
    maa.SetReferenceExpectation( GetReferenceExpectation()); //set this value equal to the available original value
    maa.SetExpectation( GetExpectation());
    maa.Run();

    //import MAA path data
    nosteps = maa.GetAlnSteps();
    mtched = maa.GetNoMatched();
    idns = maa.GetNoIdentities();
    psts = maa.GetNoPositives();

    if( nosteps <= 0 ) {
        SetAlnScore( 0 );
        SetFinalScore( 0 );
        SetAlnSteps( 0 );
        SetNoMatched( 0 );
        SetNoIdentities( 0 );
        SetNoPositives( 0 );
        return;
    }

    if( qsize + ssize < nosteps )
        throw myruntime_error("ProfileAlignment: MARealign: Invalid no. MAA steps.");

    for( s = 0; s < nosteps; s++ ) {
        path[s][nQuery] = maa.path[s][nQuery];
        path[s][nSbjct] = maa.path[s][nSbjct];
    }

    SetAlnSteps( nosteps );
    SetNoMatched( mtched );
//     SetNoIdentities( idns );
//     SetNoPositives( psts );
//     SetAlnScore( maa.GetAlnScore());
//     PostProcess();
//     SetFinalScore( GetAlnScore());
    ComputeStatistics();
}

// -------------------------------------------------------------------------
// AdjustScore: adjust score and recompute its statistical significance
//
void ProfileAlignment::AdjustScore( double value )
{
    finscore_ = value;
    ComputeStatistics();
}

// -------------------------------------------------------------------------
// PostProcess: post process alignment; change score and statistical 
//     significance
//
void ProfileAlignment::PostProcess()
{
    if( GetScoreMatrix() == NULL )
        return;
    if( GetScoreMatrix()->GetScoreAdjRelEnt()) RelEntropyAdjustment();
    if( GetScoreMatrix()->GetScoreAdjCorrel()) CorrelatedScoreAdjustment();
}

// -------------------------------------------------------------------------
// AutocorrScoreCmb: calculate combined autocorrelation value of diagonal 
//  scores given window size
//
double ProfileAlignment::AutocorrScoreCmb( const AbstractScoreMatrix* scmatrix, int sbjctpos, int querypos )
{
    const int   winsz = 4;
    double  acorr, value;
    int w;
    value = 0.0;
    for( w = 1; w <= winsz; w++ ) {
        acorr = AutocorrScore( scmatrix, sbjctpos, querypos, 0, 0, GetSubjectSize(), GetQuerySize(), w );
        value += acorr;
    }
    return value / ( double )winsz;
}

// -------------------------------------------------------------------------
// Autocorrelation: compute sum of autocorrelation of diagonal scores given
//     window size
//
double ProfileAlignment::AutocorrScore(
    const AbstractScoreMatrix* scmatrix, int sbjctpos, int querypos )
{
    int winsz = 4;
    return AutocorrScore( scmatrix, sbjctpos, querypos, 0, 0, GetSubjectSize(), GetQuerySize(), winsz );
}

// Autocorrelation: overloaded
//

double ProfileAlignment::AutocorrScore(
    const AbstractScoreMatrix* scmatrix,
    int sbjctpos, int querypos,
    int sbjctmin, int querymin,
    int sbjctlen, int querylen,
    int winsz )
{
    if( scmatrix == NULL ||
        sbjctpos <  sbjctmin || querypos <  querymin ||
        sbjctmin + sbjctlen <= sbjctpos || querymin + querylen <= querypos )
        throw myruntime_error( mystring( "ProfileAlignment: Memory access error." ));

    int         midwin = winsz >> 1;
    int         winpar = winsz ^ 1;

    double      ro = 0.0;       //sum of autocorrelations
    int         sign = 1;       //1 -- '+', -1 -- '-'
    size_t      no_pairs = 0;   //number of pairs, which is winsz! / ((winsz-2)! 2! )
    double      score1, score2;

    if( winsz <= 1 || sbjctlen < winsz || querylen < winsz ) {
        score1 = scmatrix->GetScore( sbjctpos, querypos );
        return score1;
    }

    int         spos = sbjctpos - midwin;
    int         qpos = querypos - midwin;

    //shifts have to be performed synchronously
    if( sbjctlen < spos + winsz ) { qpos -= spos + winsz - sbjctlen; spos = sbjctlen - winsz; }
    if( querylen < qpos + winsz ) { spos -= qpos + winsz - querylen; qpos = querylen - winsz; }

    if( spos < sbjctmin ) { qpos += sbjctmin - spos; spos = sbjctmin; }
    if( qpos < querymin ) { spos += querymin - qpos; qpos = querymin; }


    for( int i = 0; i < winsz - 1 && spos + i < sbjctmin + sbjctlen && qpos + i < querymin + querylen; i++ )
        for( int j = i + 1; j < winsz && spos + j < sbjctmin + sbjctlen && qpos + j < querymin + querylen; j++ ) {
            score1 = scmatrix->GetScore( spos + i, qpos + i );
            score2 = scmatrix->GetScore( spos + j, qpos + j );
            if( score1 == 0.0 || score2 == 0.0 )
                continue;
            ro += ( score1 < 0.0 && score2 < 0.0 )? -score1 * score2: score1 * score2;
            no_pairs++;
        }

    if( !no_pairs ) {
        score1 = scmatrix->GetScore( sbjctpos, querypos );
        return score1;
    }

    if( ro == 0.0 )
        return ro;

    if( ro < 0.0 )
        sign = -1;

    ro = sqrt( fabs( ro / ( double )no_pairs ));

    if( sign < 0 )
        return -ro;
    return ro;
}

// -------------------------------------------------------------------------
// Product: product of two probabilities
//
double ProfileAlignment::ProbProduct( double term1, double term2, bool apply )
{
    double val = 0.0;//min prob.
    if( !apply )
        return val;
    if( val < term1 && val < term2 )
        val = term1 * term2;
    return val;
}

// -------------------------------------------------------------------------
// LogOfProduct: log of product
//
double ProfileAlignment::LogOfProduct( double term1, double term2, bool apply )
{
    double val = SCORE_MIN;//low negative value
    if( !apply )
        return val;
    if( val < term1 && val < term2 )
        val = term1 + term2;
    return val;
}

// -------------------------------------------------------------------------
// AlignProfiles: dynamic programming to align a pair of profiles;
//  optimizational approach to alignment
//
void ProfileAlignment::AlignProfiles()
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

    int         extimeleft = 0;         //left gap extension time 
    int         extimeup = 0;           //'up' gap extension time
    int         margin;

    TPAScore    score;                  //profile score at certain positions
    int         ptr;                    //value of pointer
    int         n, m, mb, me;

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();

    if( !scmatrix )
        throw myruntime_error( "ProfileAlignment: Unable to align profiles." );

    double      gapadj = scmatrix->GetDeletionCoefficient();    //deletion coefficient
    double      deladj;                 //deletion probability
    double      insadj;                 //insertion probability
    double      tprobq, tprobs; //logs of transition probabilities of query and subject
    double      logtprod;       //log of product of transition probabilities

    margin = ( int )rint( gdMARGINPRC * ( double )SLC_MIN( GetQuerySize(), GetSubjectSize()));
    if( gnMAXMARGIN < margin )
        margin = gnMAXMARGIN;

// // const_cast<GapScheme&>(gaps_fst_).OutputGapScheme(); fprintf(stderr,"\n\n\n");
// // const_cast<GapScheme&>(gaps_sec_).OutputGapScheme();

    //iterate over query positions
    for( n = 1; n <= GetQuerySize(); n++ )
    {
        bestMM      = F[n][0][stateMM];
        bestIM      = F[n][0][stateIM];
        bestGD      = F[n][0][stateGD];

        currentMM   = F[n-1][0][stateMM];
        currentMI   = F[n-1][0][stateMI];
        currentDG   = F[n-1][0][stateDG];

// // fprintf(stderr,"\n" );

        //set boundaries of DP matrix triangles formed by both margins
        mb = 1; me = GetSubjectSize();
        if( GetQuerySize() < n + margin )
            mb += margin + n - GetQuerySize();
        if( n <= margin )
            me -= margin - n + 1;

        //iterate over subject positions
        for( m = mb; m <= me; m++ )
        {
            score  = scmatrix->GetScore( m-1, n-1 ); //profile score at positions m and n
//             score  = AutocorrScore( scmatrix, m-1, n-1 ); //TEST:autocorrelation score
// // fprintf(stderr,"\n %d %d: score: %.4f ", n, m, score/32.0 );
//  static double ww[]={0.25,0.5,0.25};
//  static double ww[]={0.08464703,0.09766335,0.1126812,0.1300084,0.15,0.1300084,0.1126812,0.09766335,0.08464703};
//  static double ww[]={0.04312742,0.06691911,0.1038357,0.1611178,0.25,0.1611178,0.1038357,0.06691911,0.04312742};
//  static double ww[]={0.009838449,0.02558575,0.06653798,0.1730378,0.45,0.1730378,0.06653798,0.02558575,0.009838449};
//  const int csz = sizeof(ww)/sizeof(double) >> 1;
//if( csz<m && m<GetSubjectSize()-csz && csz<n && n<GetQuerySize()-csz){
//  score = 0.0;
//  for(int f=-csz;f<=csz;f++)
//    score += ww[csz+f] * scmatrix->GetScore( m-1+f, n-1+f );
//}

            //-----------------
            //state MI: enter
            tprobq = gaps_fst_.GetLogTransProbsAt( P_MM, n-2 );//index value of -1 is ok: beg. state
            tprobs = gaps_sec_.GetLogTransProbsAt( P_MM, m-2 );
            oMM = currentMM;
            oMM += LogOfProduct( tprobq, tprobs, !IsUngapped());
// // fprintf(stderr," MM_MM_MM: %.4f ", LogOfProduct( tprobq, tprobs, !IsUngapped())/32.0);

            currentMM = F[n-1][m][stateMM];
            upMI_MM_MI = currentMM; //query

            tprobq = gaps_fst_.GetLogTransProbsAt( P_MM, n-2 );//index value of -1 is ok: beg. state
            tprobs = gaps_sec_.GetLogTransProbsAt( P_MI, m-1 );
            upMI_MM_MI += LogOfProduct( tprobq, tprobs, !IsUngapped());
// // fprintf(stderr," MI_MM_MI: %.4f ", LogOfProduct( tprobq, tprobs, !IsUngapped())/32.0);

            //state MI: extend
            tprobq = gaps_fst_.GetLogTransProbsAt( P_MM, n-2 );
            tprobs = gaps_sec_.GetLogTransProbsAt( P_IM, m-2 );
            oMI = currentMI;
            oMI += LogOfProduct( tprobq, tprobs, !IsUngapped());
// fprintf(stderr," MM_MM_IM: %.4f ", LogOfProduct( tprobq, tprobs, !IsUngapped())/32.0);

            currentMI = F[n-1][m][stateMI];
            upMI_MM_II = currentMI; //query

            tprobq = gaps_fst_.GetLogTransProbsAt( P_MM, n-2 );
            tprobs = gaps_sec_.GetLogTransProbsAt( P_II, m-1 );
            upMI_MM_II += LogOfProduct( tprobq, tprobs, !IsUngapped());
// fprintf(stderr," MI_MM_II: %.4f ", LogOfProduct( tprobq, tprobs, !IsUngapped())/32.0);

            //-----------------
            //state IM: enter
            leftIM_MI_MM = bestMM;   //subject

            tprobq = gaps_fst_.GetLogTransProbsAt( P_MI, n-1 );
            tprobs = gaps_sec_.GetLogTransProbsAt( P_MM, m-2 );//index value of -1 is ok: beg. state
            leftIM_MI_MM += LogOfProduct( tprobq, tprobs, !IsUngapped());
// // fprintf(stderr," IM_MI_MM: %.4f ", LogOfProduct( tprobq, tprobs, !IsUngapped())/32.0);

            //state IM: extend
            tprobq = gaps_fst_.GetLogTransProbsAt( P_IM, n-2 );
            tprobs = gaps_sec_.GetLogTransProbsAt( P_MM, m-2 );
            oIM = F[n-1][m-1][stateIM];
            oIM += LogOfProduct( tprobq, tprobs, !IsUngapped());
// fprintf(stderr," MM_IM_MM: %.4f ", LogOfProduct( tprobq, tprobs, !IsUngapped())/32.0);

            leftIM_II_MM = bestIM;   //subject

            tprobq = gaps_fst_.GetLogTransProbsAt( P_II, n-1 );
            tprobs = gaps_sec_.GetLogTransProbsAt( P_MM, m-2 );
            leftIM_II_MM += LogOfProduct( tprobq, tprobs, !IsUngapped());
// fprintf(stderr," IM_II_MM: %.4f ", LogOfProduct( tprobq, tprobs, !IsUngapped())/32.0);

            //-----------------
            //state DG: enter
            upDG_MD_1 = currentMM;  //query

            tprobq = gaps_fst_.GetLogTransProbsAt( P_MD, n-2 );
            upDG_MD_1 += LogOfProduct( tprobq, 0.0, !IsUngapped());
// // fprintf(stderr," DG_MD: %.4f ", LogOfProduct( tprobq, 0.0, !IsUngapped())/32.0);

            //state DG: extend
            tprobq = gaps_fst_.GetLogTransProbsAt( P_DM, n-2 );
            tprobs = gaps_sec_.GetLogTransProbsAt( P_MM, m-2 );
            oDG = currentDG;
            oDG += LogOfProduct( tprobq, tprobs, !IsUngapped());
// fprintf(stderr," MM_DM_MM: %.4f ", LogOfProduct( tprobq, tprobs, !IsUngapped())/32.0);

            currentDG = F[n-1][m][stateDG];
            upDG_DD_1 = currentDG;  //query

            tprobq = gaps_fst_.GetLogTransProbsAt( P_DD, n-2 );
            upDG_DD_1 += LogOfProduct( tprobq, 0.0, !IsUngapped());
// fprintf(stderr," DG_DD: %.4f ", LogOfProduct( tprobq, 0.0, !IsUngapped())/32.0);

            //-----------------
            //state GD: enter
            leftGD_1_MD = bestMM;    //subject

            tprobs = gaps_sec_.GetLogTransProbsAt( P_MD, m-2 );
            leftGD_1_MD += LogOfProduct( 0.0, tprobs, !IsUngapped());
// // fprintf(stderr," GD_MD: %.4f ", LogOfProduct( 0.0, tprobs, !IsUngapped())/32.0);

            //state GD: extend
            tprobq = gaps_fst_.GetLogTransProbsAt( P_MM, n-2 );
            tprobs = gaps_sec_.GetLogTransProbsAt( P_DM, m-2 );
            oGD = F[n-1][m-1][stateGD];
            oGD += LogOfProduct( tprobq, tprobs, !IsUngapped());
// fprintf(stderr," MM_DM_DM: %.4f ", LogOfProduct( tprobq, tprobs, !IsUngapped())/32.0);

            leftGD_1_DD = bestGD;   //subject

            tprobs = gaps_sec_.GetLogTransProbsAt( P_DD, m-2 );
            leftGD_1_DD += LogOfProduct( 0.0, tprobs, !IsUngapped());
// fprintf(stderr," GD_DD: %.4f ", LogOfProduct( 0.0, tprobs, !IsUngapped())/32.0);

// // //
// // //
            //process state stateMI
            //.....................
            if( upMI_MM_II < upMI_MM_MI ) {
                extimeup = 0;
                bestMI = upMI_MM_MI;
                ptr = dMM; //state (direction) for backtracing
            } else if( upMI_MM_MI < upMI_MM_II ) {
                        extimeup ++;
                        bestMI = upMI_MM_II;
                        ptr = dMI; //up direction
                    } else { //upMI_MM_MI == upMI_MM_II
                        extimeup = 0;
                        bestMI = upMI_MM_MI;
                        ptr = dMM_MI; //diagonal or up
                    }
            if( 0 < bestMI ) {
                F[n][m][stateMI] = bestMI;
                pointer[n][m][stateMI] = ptr;
            } else {
                bestMI = 0;
                pointer[n][m][stateMI] = dNo;
            }

            //process state stateIM
            //.....................
            if( leftIM_II_MM < leftIM_MI_MM ) {
                extimeleft = 0;
                bestIM = leftIM_MI_MM;
                ptr = dMM; //state (direction) for backtracing
            } else if( leftIM_MI_MM < leftIM_II_MM ) {
                        extimeleft ++;
                        bestIM = leftIM_II_MM;
                        ptr = dIM; //left direction
                    } else { //leftIM_MI_MM == leftIM_II_MM
                        extimeleft = 0;
                        bestIM = leftIM_MI_MM;
                        ptr = dMM_IM; //diagonal or left
                    }
            if( 0 < bestIM ) {
                F[n][m][stateIM] = bestIM;
                pointer[n][m][stateIM] = ptr;
            } else {
                bestIM = 0;
                pointer[n][m][stateIM] = dNo;
            }

            //process state stateDG
            //.....................
            if( upDG_DD_1 < upDG_MD_1 ) {
                extimeup = 0;
                bestDG = upDG_MD_1;
                ptr = dMM; //backtracing state (direction)
            } else if( upDG_MD_1 < upDG_DD_1 ) {
                        extimeup ++;
                        bestDG = upDG_DD_1;
                        ptr = dDG; //up direction
                    } else { //upDG_MD_1 == upDG_DD_1
                        extimeup = 0;
                        bestDG = upDG_MD_1;
                        ptr = dMM_DG; //diagonal or up
                    }
            if( 0 < bestDG ) {
                F[n][m][stateDG] = bestDG;
                pointer[n][m][stateDG] = ptr;
            } else {
                bestDG = 0;
                pointer[n][m][stateDG] = dNo;
            }

            //process state stateGD
            //.....................
            if( leftGD_1_DD < leftGD_1_MD ) {
                extimeleft = 0;
                bestGD = leftGD_1_MD;
                ptr = dMM; //backtracing state (direction)
            } else if( leftGD_1_MD < leftGD_1_DD ) {
                        extimeleft ++;
                        bestGD = leftGD_1_DD;
                        ptr = dGD; //left direction
                    } else { //leftGD_1_MD == leftGD_1_DD
                        extimeleft = 0;
                        bestGD = leftGD_1_MD;
                        ptr = dMM_GD; //diagonal or left
                    }
            if( 0 < bestGD ) {
                F[n][m][stateGD] = bestGD;
                pointer[n][m][stateGD] = ptr;
            } else {
                bestGD = 0;
                pointer[n][m][stateGD] = dNo;
            }

            //process state stateMM
            //.....................
            if( oMI < oMM ) {//1.
                if( oIM < oMM ) {//2.
                    if( oDG < oMM ) {
                        if( oGD < oMM ) { bestMM = oMM; ptr = dMM; }
                        else if( oMM < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oMM; ptr = dMM_GD; }
                    }
                    else if( oMM < oDG ) {
                        if( oGD < oDG ) { bestMM = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMM
                        if( oGD < oMM ) { bestMM = oMM; ptr = dMM_DG; }
                        else if( oMM < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oMM; ptr = dMM | dDG | dGD; }
                    }
                }
                else if( oMM < oIM ) {//2.
                    if( oDG < oIM ) {
                        if( oGD < oIM ) { bestMM = oIM; ptr = dIM; }
                        else if( oIM < oGD ) { bestMM = oGD; ptr = dGD; }
                        else {bestMM = oIM; ptr = dIM_GD; }
                    }
                    else if( oIM < oDG ) {
                        if( oGD < oDG ) { bestMM = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oIM
                        if( oGD < oIM ) { bestMM = oIM; ptr = dIM_DG; }
                        else if( oIM < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oIM; ptr = dIM | dDG | dGD; }
                    }
                }
                else {//2. oIM == oMM /// ///
                    if( oDG < oMM ) {
                        if( oGD < oMM ) { bestMM = oMM; ptr = dMM_IM; }
                        else if( oMM < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oMM; ptr = dMM | dIM | dGD; }
                    }
                    else if( oMM < oDG ) {
                        if( oGD < oDG ) { bestMM = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMM
                        if( oGD < oMM ) { bestMM = oMM; ptr = dMM | dIM | dDG; }
                        else if( oMM < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oMM; ptr = dMM | dIM | dDG | dGD; }
                    }
                }
            }
            else if( oMM < oMI ) {//1.
                if( oIM < oMI ) {//2.
                    if( oDG < oMI ) {
                        if( oGD < oMI ) { bestMM = oMI; ptr = dMI; }
                        else if( oMI < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oMI; ptr = dMI_GD; }
                    }
                    else if( oMI < oDG ) {
                        if( oGD < oDG ) { bestMM = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMI
                        if( oGD < oMI ) { bestMM = oMI; ptr = dMI_DG; }
                        else if( oMI < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oMI; ptr = dMI | dDG | dGD; }
                    }
                }
                else if( oMI < oIM ) {//2.
                    if( oDG < oIM ) {
                        if( oGD < oIM ) { bestMM = oIM; ptr = dIM; }
                        else if( oIM < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oIM; ptr = dIM_GD; }
                    }
                    else if( oIM < oDG ) {
                        if( oGD < oDG ) { bestMM = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oIM
                        if( oGD < oIM ) { bestMM = oIM; ptr = dIM_DG; }
                        else if( oIM < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oIM; ptr = dIM | dDG | dGD; }
                    }
                }
                else {//2. oIM == oMI /// ///
                    if( oDG < oMI ) {
                        if( oGD < oMI ) { bestMM = oMI; ptr = dMI_IM; }
                        else if( oMI < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oMI; ptr = dMI | dIM | dGD; }
                    }
                    else if( oMI < oDG ) {
                        if( oGD < oDG ) { bestMM = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMI
                        if( oGD < oMI ) { bestMM = oMI; ptr = dMI | dIM | dDG; }
                        else if( oMI < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oMI; ptr = dMI | dIM | dDG | dGD; }
                    }
                }
            }
            else {//1. oMI == oMM ///
                if( oIM < oMM ) {//2.
                    if( oDG < oMM ) {
                        if( oGD < oMM ) { bestMM = oMM; ptr = dMM_MI; }
                        else if( oMM < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oMM; ptr = dMM | dMI | dGD; }
                    }
                    else if( oMM < oDG ) {
                        if( oGD < oDG ) { bestMM = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMM
                        if( oGD < oMM ) { bestMM = oMM; ptr = dMM | dMI | dDG; }
                        else if( oMM < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oMM; ptr = dMM | dMI | dDG | dGD; }
                    }
                }
                else if( oMM < oIM ) {//2.
                    if( oDG < oIM ) {
                        if( oGD < oIM ) { bestMM = oIM; ptr = dIM; }
                        else if( oIM < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oIM; ptr = dIM_GD; }
                    }
                    else if( oIM < oDG ) {
                        if( oGD < oDG ) { bestMM = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oIM
                        if( oGD < oIM ) { bestMM = oIM; ptr = dIM_DG; }
                        else if( oIM < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oIM; ptr = dIM | dDG | dGD; }
                    }
                }
                else {//2. oIM == oMM /// ///
                    if( oDG < oMM ) {
                        if( oGD < oMM ) { bestMM = oMM; ptr = dMM | dMI | dIM; }
                        else if( oMM < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oMM; ptr = dMM | dMI | dIM | dGD; }
                    }
                    else if( oMM < oDG ) {
                        if( oGD < oDG ) { bestMM = oDG; ptr = dDG; }
                        else if( oDG < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oDG; ptr = dDG_GD; }
                    }
                    else {// oDG == oMM
                        if( oGD < oMM ) { bestMM = oMM; ptr = dMM | dMI | dIM | dDG; }
                        else if( oMM < oGD ) { bestMM = oGD; ptr = dGD; }
                        else { bestMM = oMM; ptr = dMM | dMI | dIM | dDG | dGD; }
                    }
                }
            }

            bestMM += score;
            if( 0 < bestMM ) {
                F[n][m][stateMM] = bestMM;
                pointer[n][m][stateMM] = ptr;
            } else {
                bestMM = 0;
                pointer[n][m][stateMM] = dNo;
            }
        }
    }
}

// -------------------------------------------------------------------------
// MakeAlignmentPath: make alignment path between two profiles given 
//  filled up dynamic programming matrices
//
void ProfileAlignment::MakeAlignmentPath()
{
    TPAScore    score = 0;              //maximum score from dynamic programming matrix
    int         row = 0;                //row index of maximum score
    int         col = 0;                //column index of maximum score
    State       laststate = stateMM;    //last state of back-tracing
    State       state = laststate;      //state of back-tracing
    int         step = 0;               //path index
    int         mtched = 0;
    int         identities = 0;
    int         positives = 0;
    int         margin;
    int n, m, mb, me;

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();

    if( !scmatrix )
        throw myruntime_error( "ProfileAlignment: No scoring system." );

    if( !F )
        throw myruntime_error( "ProfileAlignment: Unable to derive alignment path." );

    if( !path )
        throw myruntime_error( mystring( "ProfileAlignment: Unable to make alignment path." ));

    margin = ( int )rint( gdMARGINPRC * ( double )SLC_MIN( GetQuerySize(), GetSubjectSize()));
    if( gnMAXMARGIN < margin )
        margin = gnMAXMARGIN;

    //find the maximum value 
    for( n = 1; n <= GetQuerySize(); n++ ) {
        //set boundaries of DP matrix triangles formed by both margins
        mb = 1; me = GetSubjectSize();
        if( GetQuerySize() < n + margin )
            mb += margin + n - GetQuerySize();
        if( n <= margin )
            me -= margin - n + 1;
        for( m = mb; m <= me; m++ )
            if( score < F[n][m][stateMM]) {//ending with gap is not allowed
                score = F[n][m][stateMM];
                row = n;
                col = m;
            }
    }
    if( score <= 0 )
        return;

    path[step][nQuery] = row;
    path[step][nSbjct] = col;
    step++;

    while( 0 < row && 0 < col ) {
        state = laststate;
        laststate = GetState( pointer[row][col][state] );
        if( laststate == noStates )
            break;  //end of path
        switch( state ) {
            case stateMM:   row--; col--; 
                            if( logo_fst_[row] == logo_sec_[col])
                                identities++;
                            if( 0 < scmatrix->GetScore(col,row))
                                positives++;
                            mtched++; 
                            break;
            case stateMI:   row--; break;
            case stateIM:   col--; break;
            case stateDG:   row--; break;
            case stateGD:   col--; break;
            default:        break;
        };
        path[step][nQuery] = row;
        path[step][nSbjct] = col;

        step++;
    }

    SetAlnScore( score );
    SetAlnSteps( step - 1 ); //substract: beginning points to cell of score 0
    SetNoMatched( mtched );
    SetNoIdentities( identities );
    SetNoPositives( positives );
}

// -------------------------------------------------------------------------
// CorrelatedScoreAdjustment: adjust alignment scores by quadratically 
//     superposing adjacent scores of aligned match positions
//
void ProfileAlignment::CorrelatedScoreAdjustment()
{
    if( GetScoreMatrix() == NULL )
        return;

    const int       nopps = 4;
    const bool      doscale = GetScoreMatrix()->GetAutoScaling();
    const double    scale = GetScoreMatrix()->GetAutoScalingFactor();

//     static double ww[]={0.08464703,0.09766335,0.1126812,0.1300084,0.15,0.1300084,0.1126812,0.09766335,0.08464703};
//     static double ww[]={0.04312742,0.06691911,0.1038357,0.1611178,0.25,0.1611178,0.1038357,0.06691911,0.04312742};
//     static double ww[]={0.009838449,0.02558575,0.06653798,0.1730378,0.45,0.1730378,0.06653798,0.02558575,0.009838449};
//     static double ww[]={0.07637643,0.1986235,0.45,0.1986235,0.07637643};
//     static double ww[]={0.25,0.5,0.25};
    int qin = -1, sin = -1, len = 0;
    int step = 0;//path index
    int n;
    TPAScore sc, prev[nopps];
    TPAScore mmsc = 0, crsc = 0;
    TPAScore score = GetAlnScore();

    for( n = 0; n < nopps; n++ )
        prev[n] = 0;

    for( step = GetAlnSteps() - 1; 0 <= step; step-- ) {
        //if previous index coincides with current index
        if(  step < GetAlnSteps() - 1  &&  path[step][nQuery] == path[step+1][nQuery] )
            continue; // gap in query
        qin = path[step][nQuery] - 1;//-1 for DPM indices are one greater...

        if(  step < GetAlnSteps() - 1  &&  path[step][nSbjct] == path[step+1][nSbjct] )
            continue; // gap in subject
        sin = path[step][nSbjct] - 1;

        if( qin < 0 || sin < 0 )
            continue;

        sc = GetScoreMatrix()->GetImageScore( sin, qin );
//         mmsc += sc;

        for( n = 0; n < nopps; n++ )
            crsc += sc * prev[n];
// fprintf(stderr,"sc=%g crsc=%g  p=%g,%g,%g,%g\n",sc,crsc,prev[0],prev[1],prev[2],prev[3]);

        for( n = nopps-1; 0 < n; n-- )
            prev[n] = prev[n-1];
        prev[0] = sc;
//         for( n = 0; n < nopps; n++ )
//             crsc += ww[n] * prev[n];
        len++;
    }

    if( crsc ) {
        crsc /= ( double )( nopps >> 1 )/** len*/;
        crsc = ( 0.0 <= crsc )? sqrt( crsc ): -sqrt( -crsc );
        crsc -= mmsc;
        if( doscale )
            crsc *= scale;
        score += crsc;
        if( score < 0.0 )
            score = 0.0;
        SetAlnScore( score );
    }
}

// -------------------------------------------------------------------------
// RelEntropyAdjustment: adjust alignment score using mutual similarity
//     measure after alignment is produced. For this, all ungapped alignment
//     segments within original alignment are processed seperately to obtain
//     partial adjustment sums. For an ungapped pair segment a relative
//     entropy is computed using probabilities of run within the segment.
//     The obtained relative entropy is divided by the ungapped
//     alignment length to obtain the mean distance per position.
//
void ProfileAlignment::RelEntropyAdjustment()
{
    if( GetQuerySize() <= 0 || GetSubjectSize() <= 0 )
        return;

    if( !path )
        return;

    int     qin = -1;       //index for query
    int     sin = -1;       //index for subject
    double  colscore = 0.0; //column pair score

    int     step;
    int     alnlength = 0;
    bool    gapcur = false;

    const int       no_cols = 100;

    FrequencyMatrix loc_fst_freq;   loc_fst_freq.Reserve( no_cols );
    FrequencyMatrix loc_sec_freq;   loc_sec_freq.Reserve( no_cols );
    LogOddsMatrix   loc_fst_pssm;   loc_fst_pssm.Reserve( no_cols );
    LogOddsMatrix   loc_sec_pssm;   loc_sec_pssm.Reserve( no_cols );

    try {
        for( step = GetAlnSteps() - 1; 0 <= step; step-- )
        {
            //if previous index coincides with current index
            if(  step < GetAlnSteps() - 1  &&  path[step][nQuery] == path[step+1][nQuery] )
                    gapcur = true; // gap in query
                    //-1 for DPM indices are one greater...
            else    qin = path[step][nQuery] - 1;

            if(  step < GetAlnSteps() - 1  &&  path[step][nSbjct] == path[step+1][nSbjct] )
                    gapcur = true; // gap in subject
            else    sin = path[step][nSbjct] - 1;

            if( gapcur ) {
                gapcur = false;
                if( 0 < alnlength ) {
//                     ProcessUngappedSegment( loc_fst_freq, loc_fst_pssm, loc_sec_freq, loc_sec_pssm );

                    loc_fst_freq.Clear(); loc_fst_pssm.Clear();
                    loc_sec_freq.Clear(); loc_sec_pssm.Clear();
                    alnlength = 0;
                }
                continue;
            }

            if( qin < 0 || sin < 0 )
                continue;

            colscore = GetScoreMatrix()->GetScore( sin, qin );

            if( logo_fst_.GetInformationAt( qin ) < GetScoreMatrix()->GetInformationThreshold() ||
                logo_sec_.GetInformationAt( sin ) < GetScoreMatrix()->GetInformationThreshold())
                    SubtractScore( colscore );
//             else
//                 if( 0.0 < colscore && colscore < 0.5 * PP_SCALE_CONSTANT )
//                     SubtractScore( colscore );

//             CopyFrequencyColumn( freq_fst_, loc_fst_freq, qin );
//             CopyFrequencyColumn( freq_sec_, loc_sec_freq, sin );

//             CopyProfileColumn( logo_fst_, loc_fst_pssm, qin );
//             CopyProfileColumn( logo_sec_, loc_sec_pssm, sin );

            alnlength++;
        }

        if( 0 < alnlength ) {
//             ProcessUngappedSegment( loc_fst_freq, loc_fst_pssm, loc_sec_freq, loc_sec_pssm );
        }

// fprintf( stderr, "----*\n\n" );

    } catch( myexception const& ex ) {
        error( ex.what());
        return;
    }
}

// -------------------------------------------------------------------------
// CopyFrequencyColumn: Copy frequency vector column from one matrix to a
//     new one
// -------------------------------------------------------------------------

void ProfileAlignment::CopyFrequencyColumn( const FrequencyMatrix& freq_from, FrequencyMatrix& freq_to, int ind ) const
{
    if( freq_from.GetColumns() <= ind || ind < 0 )
        throw myruntime_error( mystring( "ProfileAlignment: Memory access error." ));

    char    residue = freq_from.GetResidueAt( ind );

    if( residue == GAP )
        //cannot be gap
        throw myruntime_error( mystring( "ProfileAlignment: Gap found in ungapped segment." ));

    //obtain column attributes
    const double  ( *weights )[NUMALPH] = freq_from.GetVectorAt( ind );

    freq_to.Push( *weights, residue );
}

// -------------------------------------------------------------------------
// CopyProfileColumn: Copy profile column from one pssm to a new one
// -------------------------------------------------------------------------

void ProfileAlignment::CopyProfileColumn( const LogOddsMatrix& pssm_from, LogOddsMatrix& pssm_to, int ind ) const
{
    if( pssm_from.GetColumns() <= ind || ind < 0 )
        throw myruntime_error( mystring( "ProfileAlignment: Memory access error." ));

    char    residue = pssm_from.GetResidueAt( ind );

    if( residue == GAP )
        //cannot be gap
        throw myruntime_error( mystring( "ProfileAlignment: Gap found in ungapped segment." ));

    //obtain PSSM column attributes
    const double  ( *scores )[NUMALPH] = pssm_from.GetVectorAt( ind );
    double           freqweight = pssm_from.GetFrequencyWeightAt( ind );
    double           information = pssm_from.GetInformationAt( ind );
    double           expMIDs[PS_NSTATES];

    expMIDs[PS_M] = pssm_from.GetMIDExpNoObservationsAt( ind, PS_M );
    expMIDs[PS_I] = pssm_from.GetMIDExpNoObservationsAt( ind, PS_I );
    expMIDs[PS_D] = pssm_from.GetMIDExpNoObservationsAt( ind, PS_D );

    pssm_to.Push( *scores, residue, freqweight, information, expMIDs );
}

// -------------------------------------------------------------------------
// ProcessUngappedSegment: Process one segment of ungapped alignment by
//     computing its relative entropy;
// return absolute value of difference of entropies of the aligned segments
// -------------------------------------------------------------------------

void ProfileAlignment::ProcessUngappedSegment(
    const FrequencyMatrix& f1, const LogOddsMatrix& l1,
    const FrequencyMatrix& f2, const LogOddsMatrix& l2 )
{
    if( !l1.IsCompatible( f1 ) || !l2.IsCompatible( f2 ) || f1.GetColumns() != f2.GetColumns())
        throw myruntime_error( mystring( "ProfileAlignment: Ungapped segments are incompatible." ));

    SEGProfile  segdiff( f1, f2 );

//     SEGProfile  segpro1( f1, l1, f1.GetColumns());  segpro1.SetDistance( GetSEGdistanceThreshold());
//     SEGProfile  segpro2( f2, l2, f2.GetColumns());  segpro2.SetDistance( GetSEGdistanceThreshold());
// 
//     double      ent1 = segpro1.Entropy();
//     double      ent2 = segpro2.Entropy();

//     SubtractScore( fabs( ent1 - ent2 ));

//     double      logP1 = segpro1.LogProbability();
//     double      logP2 = segpro2.LogProbability();


// int n;
// fprintf( stderr, "%.6g\n", fabs( ent1 - ent2 ) /* * f1.GetColumns() */);
// 
// for( n = 0; n < f1.GetColumns(); n++ )
//     fprintf( stderr, "%c", DehashCode( f1.GetResidueAt( n )));
// fprintf( stderr, "\n" );
// 
// for( n = 0; n < f1.GetColumns(); n++ )
//     fprintf( stderr, "%c", DehashCode( f2.GetResidueAt( n )));
// fprintf( stderr, "\n\n" );

}

// -------------------------------------------------------------------------
// GetAlnLength: return alignment length after alignment has been produced
// -------------------------------------------------------------------------

int ProfileAlignment::GetAlnLength() const
{
    if( GetQuerySize() <= 0 || GetSubjectSize() <= 0 )
        return -1;

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();

    if( !path || !scmatrix )
        return -1;

    int     qin = -1;   //index of query
    int     sin = -1;   //index of subject

    int     step;
    int     alnlength = 0;
    bool    gapcur = false;

    for( step = GetAlnSteps() - 1; 0 <= step; step-- )
    {
        if( step < GetAlnSteps() - 1  &&  path[step][nQuery] == path[step+1][nQuery] )
              gapcur = true; // gap in query
              //-1 for DPM indices are one greater...
        else  qin = path[step][nQuery] - 1;

        if( step < GetAlnSteps() - 1  &&  path[step][nSbjct] == path[step+1][nSbjct] )
              gapcur = true; // gap in subject
        else  sin = path[step][nSbjct] - 1;

        if( gapcur ) {
            alnlength++;
            gapcur = false;
            continue;
        }

        if( 0 <= qin && 0 <= sin )
            if( scmatrix->GetScore( sin, qin ) != 0 )
                alnlength++;
    }

    return alnlength;
}

// -------------------------------------------------------------------------
// Output: output alignment and statistical information
// -------------------------------------------------------------------------

void ProfileAlignment::Output( const char* filename )
{
    FILE*   fp = stdout;

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw myruntime_error( "Failed to open file for writing." );

    Print( fp );

    if( fp != stdout )
        fclose( fp );
}

// -------------------------------------------------------------------------
// GetMaxReqSizeForAlignment: predict amount of memory to contain alignment 
//     information
// -------------------------------------------------------------------------

size_t ProfileAlignment::GetMaxReqSizeForAlignment() const
{
    static size_t   nl_sz = 2;                  //size of newline
    static size_t   no_lines = 4;               //number of lines per alignment fragment
    static size_t   no_infol = 4;               //number of info lines coming along with alignment (scores, etc.)
    static size_t   num_space = 6;              //space for numbers
    static size_t   indent = OUTPUTINDENT;      //indentation
    static size_t   textwidth = OUTPUTWIDTH;    //output width
    static size_t   max_foot_lines = 6;         //maximum number of lines comprising footer

    if( indent < num_space )
        indent = num_space;

    size_t  max_line_width = textwidth + 3 * indent + 2 * nl_sz; //maximum line width

    if( max_line_width < 80 )
        max_line_width = 80;

    size_t  max_footer_size = max_foot_lines * max_line_width;

    if( !textwidth )
        return max_footer_size;

    size_t  alnlen = GetAlnSteps();             //alignment length
    size_t  alnfrags = alnlen / textwidth + 1;  //number of alignment fragments

    size_t  aln_size = ( no_infol + alnfrags * no_lines ) * max_line_width;

    return aln_size + max_footer_size;
}

// -------------------------------------------------------------------------
// Print: print alignment information to string stream supposed to have 
//     enough space allocated before
//

void ProfileAlignment::Print( char* sp, bool showpars )
{
    if( !sp )
        return;
    *sp = 0; //to ensure printing to the end of stream
    Print( &string_print, sp, showpars );
}

// -------------------------------------------------------------------------
// Print: print alignment information to file
//

void ProfileAlignment::Print( FILE* fp, bool showpars )
{
    Print( &file_print, fp, showpars );
}

// -------------------------------------------------------------------------
// Print: print alignment information to stream pointed by vpn which can be
//     either file or string
// -------------------------------------------------------------------------

void ProfileAlignment::Print( TPrintFunction print_func, void* vpn, bool showpars )
{
    if( vpn == NULL )
        return;

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();

    if( GetQuerySize() <= 0 || GetSubjectSize() <= 0 )
        throw myruntime_error( "ProfileAlignment: Unable to output alignment." );

    if( !path )
        throw myruntime_error( "ProfileAlignment: No tracing path generated." );

    if( !scmatrix )
        throw myruntime_error( "ProfileAlignment: No score matrix." );

    char    qaa;    //query amino acid
    char    saa;    //subject amino acid

    int     qin;    //index for query 
    int     sin;    //index for subject

    int     step;
    int     identities = 0, positives = 0, gaps = 0; //useful statistics

    mystring    query; query.reserve( KBYTE );
    mystring    sbjct; sbjct.reserve( KBYTE );
    mystring    match; match.reserve( KBYTE );

    for( step = GetAlnSteps() - 1; 0 <= step; step-- )
    {
        qaa = 0; saa = 0;
        if( step < GetAlnSteps()-1  &&  path[step][nQuery] == path[step+1][nQuery] )
              query.append( '-' );
              //-1 for DPM indices are one greater...
        else  query.append( qaa = DehashCode( logo_fst_[ qin=path[step][nQuery]-1 ] ));

        if( step < GetAlnSteps()-1  &&  path[step][nSbjct] == path[step+1][nSbjct] )
              sbjct.append( '-' );
        else  sbjct.append( saa = DehashCode( logo_sec_[ sin=path[step][nSbjct]-1 ] ));

        if( qaa && saa && qaa == saa ) {
            match.append( qaa );
            identities++;
        } else if( qaa && saa && 0 < scmatrix->GetScore( sin, qin )) {
                    match.append( '+' );
                    positives++;
                } else {
                    match.append( ' ' );
                    if( !qaa || !saa )
                        gaps++;
                }
    }

    print_func( vpn, "\n" );
    logo_sec_.PrintDescription( print_func, vpn );

    print_func( vpn, "  Length/enos: Query = %d/%.1f, Sbjct = %d/%.1f\n", 
                logo_fst_.GetColumns(), logo_fst_.GetEffNoSequences(), 
                logo_sec_.GetColumns(), logo_sec_.GetEffNoSequences());
    print_func( vpn, "\n" );

    if( GetBitScore() < 0.0 )
        print_func( vpn, " Score = %.2f, ", GetScore());
    else
        print_func( vpn, " Score = %.2f (%.1f bits), ", GetScore(), GetBitScore());
    if( GetExpectation() < 0.0 )
//         print_func( vpn, " Expect = n/a, P-value = n/a\n" );
        print_func( vpn, " Expect = %.2g, P-value = %.2g\n", GetReferenceExpectation(), ComputePvalue( GetReferenceExpectation()));
    else
        print_func( vpn, " Expect = %.2g, P-value = %.2g\n", GetExpectation(), GetPvalue());

    if( identities )
        print_func( vpn, " Identities = %d/%d (%d%%)", identities, GetAlnSteps(), int( identities * 100 / GetAlnSteps() ));
    if( positives ) {
        if( identities )
            print_func( vpn, "," );
        print_func( vpn, " Positives = %d/%d (%d%%)", positives, GetAlnSteps(), int( positives * 100 / GetAlnSteps() ));
    }
    if( gaps ) {
        if( identities || positives )
            print_func( vpn, "," );
        print_func( vpn, " Gaps = %d/%d (%d%%)", gaps, GetAlnSteps(), int( gaps * 100 / GetAlnSteps() ));
    }


    print_func( vpn, "\n\n" );

    const int   width = OUTPUTWIDTH;
    int         begp, endp, ind;

    for( int n = GetAlnSteps() - 1; 0 <= n; n -= width )
    {
        begp = path[n][nQuery];
        endp = path[ ind = ( n - width + 1 > 0 )? n - width + 1: 0 ][ nQuery ];
        if( endp == begp && ind != n )
            endp++;
        if( n < GetAlnSteps() - 1 && begp == path[n+1][nQuery] )
            begp++;

        print_func( vpn, "Query: %5d %s %-5d\n",
                begp,
                ( query.substr( GetAlnSteps() - n - 1, width )).c_str(),
                endp );

        print_func( vpn, "%12c %s\n", 32, ( match.substr( GetAlnSteps() - n - 1, width )).c_str());

        begp = path[n][nSbjct];
        endp = path[ ind = ( n - width + 1 > 0 )? n - width + 1: 0 ][ nSbjct ];
        if( endp == begp && ind != n )
            endp++;
        if( n < GetAlnSteps() - 1 && begp == path[n+1][nSbjct] )
            begp++;

        print_func( vpn, "Sbjct: %5d %s %-5d\n\n",
                begp,
                ( sbjct.substr( GetAlnSteps() - n - 1, width )).c_str(),
                endp );
    }

    if( showpars )
        scmatrix->PrintParameterTable( print_func, vpn );
}

// -------------------------------------------------------------------------
// OutputScoringMatrix: output profile scoring system
// -------------------------------------------------------------------------

void ProfileAlignment::OutputScoringMatrix( const char* filename )
{
    const AbstractScoreMatrix* scmatrix = GetScoreMatrix();

    FILE*   fp = stdout;

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw myruntime_error( "Failed to open file for writing." );

    const_cast<AbstractScoreMatrix*>( scmatrix )->PrintScoringMatrix( fp );

    if( fp != stdout )
        fclose( fp );
}
