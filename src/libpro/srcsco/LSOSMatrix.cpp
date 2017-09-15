/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "mystring.h"
#include "myexcept.h"

#include "rc.h"
#include "data.h"
#include "ext/pslcodes.h"
#include "libpro/srcpro/SSSScores.h"
#include "libpro/srcpro/CVS2Scores.h"
#include "libpro/tfopt/TargetFreqOptimizerH.h"
#include "LSOSMatrix.h"


FrequencyMatrix smp_lsodummyfreq;
LogOddsMatrix   smp_lsodummylogo;

const bool      gcbAVGBPRB = true;//false;//averaged backgrounds
const double    g_scoreX = 0.0; //-1.0;

// -------------------------------------------------------------------------
// constructor: frequency matrices, log-odds matrices for the first and
//   second profiles are given with the parameters
//
LSOSMatrix::LSOSMatrix(
    const FrequencyMatrix&  freq_fst, const LogOddsMatrix& logo_fst,
    const FrequencyMatrix&  freq_sec, const LogOddsMatrix& logo_sec,
    double          infrm_threshold,
    int             thick_number,
    double          thick_percnt,
    double          mask_percnt,
    Configuration   config[NoSchemes],
    TBehaviour      s_behaviour,
    TScaling        a_scaling,
    TMask           c_masking )
:
    AbstractScoreMatrix( ProfSpecLSO, config, s_behaviour, a_scaling, c_masking ),

    freq_fst_( freq_fst ),
    logo_fst_( logo_fst ),

    freq_sec_( freq_sec ),
    logo_sec_( logo_sec ),

    hdpbase_( NULL ),
    hdpctbase_( NULL ),
    qryhdpbuf_( NULL ),

    thickness_number( thick_number ),
    thickness_percnt( thick_percnt ),
    maskscale_percnt( mask_percnt ),

    scores_( freq_fst_.GetColumns() * freq_sec_.GetColumns()),
    rprobs_( freq_fst_.GetColumns()),
    cprobs_( freq_sec_.GetColumns()),

    pqry_( logo_fst.GetTrgFreqs()),
    psub_( NULL ),

    cvs2scores_( NULL )
{
    SetUseModImage( true );
    PrivateInit();
    SetInformationThreshold( infrm_threshold );
    SetScoreAdjCorrel( true );
    SetSupportOptimFreq( false );
}

// -------------------------------------------------------------------------
// constructor: frequency matrices, log-odds matrices, HDPbaase given
//
LSOSMatrix::LSOSMatrix(
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
    TMask           c_masking )
:
    AbstractScoreMatrix( ProfSpecLSO, config, s_behaviour, a_scaling, c_masking ),

    freq_fst_( freq_fst ),
    logo_fst_( logo_fst ),

    freq_sec_( freq_sec ),
    logo_sec_( logo_sec ),

    hdpbase_( hdpbase ),
    hdpctbase_( hdpctbase ),
    qryhdpbuf_( NULL ),

    thickness_number( thick_number ),
    thickness_percnt( thick_percnt ),
    maskscale_percnt( mask_percnt ),

    scores_( freq_fst_.GetColumns() * freq_sec_.GetColumns()),
    rprobs_( freq_fst_.GetColumns()),
    cprobs_( freq_sec_.GetColumns()),

    pqry_( logo_fst.GetTrgFreqs()),
    psub_( NULL ),

    cvs2scores_( NULL )
{
    SetUseModImage( true );
    PrivateInit();
    SetInformationThreshold( infrm_threshold );
    SetScoreAdjCorrel( true );
    SetSupportOptimFreq( false );
}

// -------------------------------------------------------------------------
// constructor: overloaded
//
LSOSMatrix::LSOSMatrix(
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
    TType           type )
:
    AbstractScoreMatrix( type, config, s_behaviour, a_scaling, c_masking ),

    freq_fst_( freq_fst ),
    logo_fst_( logo_fst ),

    freq_sec_( freq_sec ),
    logo_sec_( logo_sec ),

    hdpbase_( NULL ),
    hdpctbase_( NULL ),
    qryhdpbuf_( NULL ),

    thickness_number( thick_number ),
    thickness_percnt( thick_percnt ),
    maskscale_percnt( mask_percnt ),

    scores_( freq_fst_.GetColumns() * freq_sec_.GetColumns()),
    rprobs_( freq_fst_.GetColumns()),
    cprobs_( freq_sec_.GetColumns()),

    pqry_( logo_fst.GetTrgFreqs()),
    psub_( NULL ),

    cvs2scores_( NULL )
{
    SetUseModImage( true );
    PrivateInit();
    SetInformationThreshold( infrm_threshold );
}

// -------------------------------------------------------------------------
// default constructor is invalid for initialization
//
LSOSMatrix::LSOSMatrix()
:
    AbstractScoreMatrix(),

    freq_fst_( smp_lsodummyfreq ),
    logo_fst_( smp_lsodummylogo ),

    freq_sec_( smp_lsodummyfreq ),
    logo_sec_( smp_lsodummylogo ),

    hdpbase_( NULL ),
    hdpctbase_( NULL ),
    qryhdpbuf_( NULL ),

    thickness_number( 0 ),
    thickness_percnt( 0.0 ),
    maskscale_percnt( 0.0 ),

    scores_( 1 ),
    rprobs_( 1 ),
    cprobs_( 1 ),

    pqry_( NULL ),
    psub_( NULL ),

    cvs2scores_( NULL )
{
    memset( mixbp_, 0, NUMAA * sizeof( double ));
}

// -------------------------------------------------------------------------
// destructor: deallocation of resources
//
LSOSMatrix::~LSOSMatrix()
{
    if( psub_ ) { free( psub_ ); psub_ = NULL; }
    if( qryhdpbuf_ ) { free( qryhdpbuf_ ); qryhdpbuf_ = NULL; }
    DestroyCVS2scores();
}

// -------------------------------------------------------------------------
// PrivateInit: initialization method private for this class
//
void LSOSMatrix::PrivateInit()
{
    const double (*tfrqs)[NUMAA];
    double bprob;
    int size;
    int p, a;
    const double    ensqry = logo_fst_.GetEffNoSequences();
    const double    enssub = logo_sec_.GetEffNoSequences();

    if( !logo_fst_.IsCompatible( freq_fst_ ) || !logo_sec_.IsCompatible( freq_sec_ ))
        throw myruntime_error("LSOSMatrix: Profile matrices are incompatible.");

    if( freq_fst_.GetColumns() < 1 || freq_fst_.GetColumns() > MAXCOLUMNS ||
        freq_sec_.GetColumns() < 1 || freq_sec_.GetColumns() > MAXCOLUMNS )
        throw myruntime_error("LSOSMatrix: Wrong profile matrices.");

    if( GetMaskscalePercents() < 0.0 )
        throw myruntime_error("LSOSMatrix: Negative scaling factor of masked positions.");

    Init( freq_fst_.GetColumns(), freq_sec_.GetColumns());

    if( GetHDPBase()) {
        size = 0;
        if( GetHDPBase()->GetMenu())
            size = GetHDPBase()->GetMenu()->GetSize();
        if( 0 < size ) {
            qryhdpbuf_ = ( double* )malloc( size * sizeof( double ));
            if( qryhdpbuf_ == NULL )
                throw myruntime_error("LSOSMatrix: Not enough memory.");
       }
    }

    for( int r = 0; r < NUMAA; r++ ) {
        SetMixBPAt( r, (logo_fst_.GetBackProbsAt(r)+logo_sec_.GetBackProbsAt(r)) * 0.5 );
        //SetMixBPAt( r, (ensqry*logo_fst_.GetBackProbsAt(r)+enssub*logo_sec_.GetBackProbsAt(r)) / (ensqry + enssub));
        //SetMixBPAt( r, logo_fst_.GetBackProbsAt(r) * logo_sec_.GetBackProbsAt(r) / GetMixBPAt(r));//sim
        //SetMixBPAt( r, logo_fst_.GetBackProbsAt(r) * logo_sec_.GetBackProbsAt(r) / LOSCORES.PROBABility(r));//nope
    }

    psub_ = (double(*)[NUMAA])malloc( logo_sec_.GetColumns()*NUMAA*sizeof( double ));
    if( psub_ == NULL )
        throw myruntime_error("LSOSMatrix: Not enough memory.");

    for( p = 0; p < logo_sec_.GetColumns(); p++ ) {
        tfrqs = GetSbjctLogo().GetTrgFreqsAt(p);
        for( a = 0; a < NUMAA; a++ ) {
            bprob = gcbAVGBPRB? GetMixBPAt(a): LOSCORES.PROBABility(a);
            if( bprob <= 0.0 )
                throw myruntime_error("LSOSMatrix: Invalid back. probabilities.");
            psub_[p][a] = (*tfrqs)[a] / bprob;
        }
    }

    //initialize cvs2scores table
    if( 1 < CVS.AVGLEN && CVS.AVGCWGT < 1.0 && 0.0 < CVS.AVGCWGT ) {
        DestroyCVS2scores();
        cvs2scores_ = ( double** )malloc( GetSubjectSize() * sizeof(double*));
        if( !cvs2scores_ )
            throw myruntime_error("LSOSMatrix: Not enough memory.");
        for( int m = 0; m < GetSubjectSize(); m++ ) {
            cvs2scores_[m] = ( double* )malloc( GetQuerySize() * sizeof(double));
            if( !cvs2scores_[m])
                throw myruntime_error("LSOSMatrix: Not enough memory.");
            memset( cvs2scores_[m], 0, GetQuerySize() * sizeof( double ));
        }
    }

#ifdef USEPROFBACKPROBS
    LOSCORES.StoreProbabilities_1( logo_fst_.GetBackProbs());
    LOSCORES.StoreProbabilities_2( logo_sec_.GetBackProbs());
#endif
}

// -------------------------------------------------------------------------
// ComputeHDPStatScore: compute HDP score within statistical HDP framework
//
double LSOSMatrix::ComputeHDPStatScore( int m, int n, bool sson, double* modscore ) const
{
    static const double minscore = -3.0;
    static const double maxscore = 3.0;
    const HDPbase*  hdpbase = GetHDPBase();
    const Menu*     menu = hdpbase? hdpbase->GetMenu(): NULL;
    const int       nodshs = menu? menu->GetSize(): 0;

    if( modscore )
        *modscore = 0.0;
    if( hdpbase == NULL || menu == NULL || nodshs < 2 )
        return 0.0;
    if( qryhdpbuf_ == NULL )
        return 0.0;

    const double    bppqry = logo_fst_.GetBckPPProbAt(n);
    const double*   ppsqry = logo_fst_.GetPPProbsAt(n);
    const int*      nxsqry = logo_fst_.GetPPPIndsAt(n);
    size_t          npsqry = logo_fst_.GetNoPPProbsAt(n);
    double          prprqr = 0.0;
    //
    const double    bppsub = logo_sec_.GetBckPPProbAt(m);
    const double*   ppssub = logo_sec_.GetPPProbsAt(m);
    const int*      nxssub = logo_sec_.GetPPPIndsAt(m);
    size_t          npssub = logo_sec_.GetNoPPProbsAt(m);
    double          prprsb = 0.0;
    //
    int     nosupcls = hdpbase->GetNoSupClusters();

    if( bppqry <= 0.0 || bppsub <= 0.0 || npsqry < 1 || npssub < 1 )
        return 0.0;

    double  score = 0.0;
    int sx;
    int i;

    memset( qryhdpbuf_, 0, nodshs * sizeof( double ));

    prprqr = ppsqry[npsqry-1];//prior probability is last
    prprsb = ppssub[npssub-1];//prior probability is last

    for( i = 0; i < npsqry; i++ ) {
        if( nxsqry[i] < 0 )
            continue;
        if( nxsqry[i] < nodshs )
            qryhdpbuf_[nxsqry[i]] = ppsqry[i];
    }

    if( prprqr && prprsb )
        score += prprqr * prprsb * menu->GetPriorProbNewDish();
    for( i = 0; i < npssub; i++ ) {
        sx = nxssub[i];
        if( sx < 0 )
            continue;
        if( sx < nodshs && qryhdpbuf_[sx] && ppssub[i])
            score += qryhdpbuf_[sx] * ppssub[i] * menu->GetPriorProbAt(sx);
    }

    if( score )
        score /= bppqry * bppsub;
    else
        return 0.0;

    if( !isfinite(score)) {
        warning("HDP score NaN.");
        return 0.0;
    }
    if( score < 0.0 ) {
        warning("Negative HDP odds ratio.");
        return 0.0;
    }
    else if( score == 0.0 )
        score = minscore;
    else
        score = log(score);

//     score = SLC_MAX( minscore, score );
    score = SLC_MAX( 0.0, score );
    score = SLC_MIN( maxscore, score );
    if( modscore )
        *modscore = score;
    return score;
}

// -------------------------------------------------------------------------
// ComputeHDPDrvdScore: get HDP score derived from distribution along 
//  structure alignments
//
double LSOSMatrix::ComputeHDPDrvdScore( int m, int n, bool sson, double* modscore ) const
{
    static const double minscore = -5.0;
    static const double maxscore = 5.0;
//     static const double scd_ensthr = 3.0;
    const HDPbase*  hdpbase = GetHDPBase();
    const HDPscores* hdpscores = hdpbase? hdpbase->GetScores(): NULL;
    const Menu*     menu = hdpbase? hdpbase->GetMenu(): NULL;
    const int       nodshs = menu? menu->GetSize(): 0;

    if( modscore )
        *modscore = 0.0;
    if( hdpbase == NULL || hdpscores == NULL || 
        menu == NULL || nodshs < 2 )
        return 0.0;

    static const double offset = 0;//0.05;
    static const double sssdowns = 0.45;//0.95;//downscale for SS states
    //////
    const double*   ppsqry = logo_fst_.GetPPProbsAt(n);
    const int*      nxsqry = logo_fst_.GetPPPIndsAt(n);
    size_t          npsqry = logo_fst_.GetNoPPProbsAt(n);
//     const double    ensqry = logo_fst_.GetMIDExpNoObservationsAt(n,PS_M);
    const double    ensqry = logo_fst_.GetEffNoSequences();
    double          prbqry = 0.0;//cluster membership probability for query
    int             dishqr = 0;
    //
    const double*   ppssub = logo_sec_.GetPPProbsAt(m);
    const int*      nxssub = logo_sec_.GetPPPIndsAt(m);
    size_t          npssub = logo_sec_.GetNoPPProbsAt(m);
//     const double    enssub = logo_sec_.GetMIDExpNoObservationsAt(m,PS_M);
    const double    enssub = logo_sec_.GetEffNoSequences();
    double          prbsub = 0.0;//cluster membership probability for sbjct
    int             dishsb = 0;
    double          score = 0.0;
    double          tmpmul;
    //
    //double          xxicov =  0.6793478;//cov=5.0
    //double          xyicov =  0.5706522;//cov=-4.2
    double          zzicov =  0.47;
    double          zwicov =  0.47;
    double          xxicov =  0.37;//0.47
    double          xyicov =  1.17;//0.77

//     if( ensqry < scd_ensthr || enssub < scd_ensthr )
//         return 0.0;

    if( npsqry < 1 || npssub < 1 )
        return 0.0;

    prbqry = ppsqry[0];
    prbsub = ppssub[0];

    dishqr = nxsqry[0];
    dishsb = nxssub[0];

    if( dishqr < 0 || dishsb < 0 )
        return 0.0;

    if( nodshs <= dishqr || nodshs <= dishsb )
        throw myruntime_error("LSOSMatrix: ComputeHDPDrvdScore: Invalid dish index.");

    score = hdpscores->GetScore( dishqr, dishsb, prbqry, prbsub, ensqry, enssub );

    if( hdpscores->ScoreNA( score ))
        //score NA
        return 0.0;

    if( sson ) {
        if( modscore )
            *modscore = score;
        score *= sssdowns;
        if( modscore && *modscore ) {
            *modscore -= offset;
            *modscore *= sssdowns;
        }
    }

    if( score ) {
//         tmpmul = exp( -0.5*2.e-5*(pow(ensqry-21.,4.)+pow(enssub-21.,4.)));
///        tmpmul = 0.4+0.6*exp( -0.5*3.e-5*(pow(ensqry-23.,4.)+pow(enssub-23.,4.)));
        tmpmul = 0.4+0.6*exp( -0.5*3.e-5*(
            zzicov*pow(ensqry-23.,4.) + zzicov*pow(enssub-23.,4.) + 
                2.*zwicov*SQUARE(ensqry-23.)*SQUARE(enssub-23.)));
        score *= tmpmul;
        if( modscore && *modscore )
            *modscore *= tmpmul;

        //tmpmul = exp( -0.5*15.*(pow(prbqry-1.,2.)+pow(prbsub-1.,2.)));
        tmpmul = exp( -0.5*15.*(
            xxicov*(prbqry-1.)*(prbqry-1.) + xxicov*(prbsub-1.)*(prbsub-1.) + 
                2.*xyicov*(prbqry-1.)*(prbsub-1.)));
        score *= tmpmul;
        if( modscore && *modscore )
            *modscore *= tmpmul;
        //
    }
    return score;
}

// -------------------------------------------------------------------------
// ComputeHDPScore: compute HDP score given positions of query and subject
//
double LSOSMatrix::ComputeHDPScore( int m, int n, bool sson, double* modscore ) const
{
    const HDPbase*  hdpbase = GetHDPBase();
    if( modscore )
        *modscore = 0.0;
    if( hdpbase == NULL )
        return 0.0;

    double  scwght = hdpbase->GetAdjWeight();
    double  score;

    if( scwght <= 0.0 )
        return 0.0;

    if( hdpbase->GetScores())
        score = ComputeHDPDrvdScore( m, n, sson, modscore );
    else
        score = ComputeHDPStatScore( m, n, sson, modscore );

    if( score )
        score *= scwght;
    if( modscore && *modscore )
        *modscore *= scwght;

    return score;
}

// -------------------------------------------------------------------------
// ComputeHDPctDrvdScore: get HDP ctx score derived from structure 
//  alignments
//
double LSOSMatrix::ComputeHDPctDrvdScore( int m, int n, bool sson, double* modscore ) const
{
    const HDPbase*  hdpctbase = GetHDPctBase();
    const HDPscores* hdpctscores = hdpctbase? hdpctbase->GetScores(): NULL;
    const Menu*     menu = hdpctbase? hdpctbase->GetMenu(): NULL;
    const int       nodshs = menu? menu->GetSize(): 0;

    if( modscore )
        *modscore = 0.0;
    if( hdpctbase == NULL || hdpctscores == NULL ||
        menu == NULL || nodshs < 2 )
        return 0.0;

    static const double offset = 0;//0.05;
    static const double sssdowns = 0.45;//0.95;//downscale for SS states
    //////
    const double*   ppsqry = logo_fst_.GetctPPProbsAt(n);
    const int*      nxsqry = logo_fst_.GetctPPPIndsAt(n);
    size_t          npsqry = logo_fst_.GetNoctPPProbsAt(n);
//     const double    ensqry = logo_fst_.GetMIDExpNoObservationsAt(n,PS_M);
    const double    ensqry = logo_fst_.GetEffNoSequences();
    double          prbqry = 0.0;//cluster membership probability for query
    int             dishqr = 0;
    //
    const double*   ppssub = logo_sec_.GetctPPProbsAt(m);
    const int*      nxssub = logo_sec_.GetctPPPIndsAt(m);
    size_t          npssub = logo_sec_.GetNoctPPProbsAt(m);
//     const double    enssub = logo_sec_.GetMIDExpNoObservationsAt(m,PS_M);
    const double    enssub = logo_sec_.GetEffNoSequences();
    double          prbsub = 0.0;//cluster membership probability for sbjct
    int             dishsb = 0;
    double          score = 0.0;
    double          tmpmul;
    //
    //double          xxicov =  0.6793478;//cov=5.0
    //double          xyicov =  0.5706522;//cov=-4.2
    double          zzicov =  0.47;
    double          zwicov =  0.47;
    double          xxicov =  0.37;//0.47
    double          xyicov =  1.17;//0.77


    if( npsqry < 1 || npssub < 1 )
        return 0.0;

    prbqry = ppsqry[0];
    prbsub = ppssub[0];

    dishqr = nxsqry[0];
    dishsb = nxssub[0];

    if( dishqr < 0 || dishsb < 0 )
        return 0.0;

    if( nodshs <= dishqr || nodshs <= dishsb )
        throw myruntime_error("LSOSMatrix: ComputeHDPctDrvdScore: Invalid dish index.");

    score = hdpctscores->GetScore( dishqr, dishsb, prbqry, prbsub, ensqry, enssub );

    if( hdpctscores->ScoreNA( score ))
        //score NA
        return 0.0;

    if( sson ) {
        if( modscore )
            *modscore = score;
        score *= sssdowns;
        if( modscore && *modscore ) {
            *modscore -= offset;
            *modscore *= sssdowns;
        }
    }

    if( score ) {
//         tmpmul = exp( -0.5*2.e-5*(pow(ensqry-21.,4.)+pow(enssub-21.,4.)));
///        tmpmul = 0.4+0.6*exp( -0.5*3.e-5*(pow(ensqry-23.,4.)+pow(enssub-23.,4.)));
        tmpmul = 0.4+0.6*exp( -0.5*3.e-5*(
            zzicov*pow(ensqry-23.,4.) + zzicov*pow(enssub-23.,4.) +
                2.*zwicov*SQUARE(ensqry-23.)*SQUARE(enssub-23.)));
        score *= tmpmul;
        if( modscore && *modscore )
            *modscore *= tmpmul;

        //tmpmul = exp( -0.5*15.*(pow(prbqry-1.,2.)+pow(prbsub-1.,2.)));
        tmpmul = exp( -0.5*15.*(
            xxicov*(prbqry-1.)*(prbqry-1.) + xxicov*(prbsub-1.)*(prbsub-1.) +
                2.*xyicov*(prbqry-1.)*(prbsub-1.)));
        score *= tmpmul;
        if( modscore && *modscore )
            *modscore *= tmpmul;
        //
    }
    return score;
}

// -------------------------------------------------------------------------
// ComputeHDPctScore: compute HDP ctx score given query and sbjct positions 
//
double LSOSMatrix::ComputeHDPctScore( int m, int n, bool sson, double* modscore ) const
{
    const HDPbase*  hdpctbase = GetHDPctBase();
    if( modscore )
        *modscore = 0.0;
    if( hdpctbase == NULL )
        return 0.0;

    double  scwght = hdpctbase->GetAdjWeight();
    double  score = 0.0;

    if( scwght <= 0.0 )
        return 0.0;

    if( hdpctbase->GetScores())
        score = ComputeHDPctDrvdScore( m, n, sson, modscore );

    if( score )
        score *= scwght;
    if( modscore && *modscore )
        *modscore *= scwght;

    return score;
}

// -------------------------------------------------------------------------
// ComputeSSSScore: compute SS state score given positions of query and 
//  subject
//
double LSOSMatrix::ComputeSSSScore( int m, int n, double* modscore ) const
{
    double  scwght = SSSSCORES.GetSSSWeight();
    if( modscore )
        *modscore = 0.0;
    if( scwght <= 0.0 )
        return 0.0;
    if( !logo_fst_.GetSSSSet() || !logo_sec_.GetSSSSet())
        return 0.0;
    if( logo_fst_.GetSSSP3() != logo_sec_.GetSSSP3()) {
        warning("Not using inconsistent SS predictions in profiles.");
        return 0.0;
    }

    static const double  offset = 0.05;//0.4;0.34
    ////
    const double    ensqry = logo_fst_.GetEffNoSequences();
    char            sssqry = logo_fst_.GetSSStateAt(n);
    double          sspqry = logo_fst_.GetSSStateProbAt(n);
    ////
    const double    enssub = logo_sec_.GetEffNoSequences();
    char            ssssub = logo_sec_.GetSSStateAt(m);
    double          sspsub = logo_sec_.GetSSStateProbAt(m);
    double          sc, score = 0.0;
    int t1, t2;
    int t10 = 0, t1n = SS_NSTATES-1;
    int t20 = 0, t2n = SS_NSTATES-1;
    ////
//     double          zzicov =  0.47;
//     double          zwicov =  0.47;

    if( SS_NSTATES <= sssqry || SS_NSTATES <= ssssub )
        throw myruntime_error("LSOSMatrix: ComputeSSSScore: Invalid SS class.");

    if( logo_fst_.GetSSSP3())
    {
        for( t1 = t10; t1 <= t1n; t1++ )
            for( t2 = t20; t2 <= t2n; t2++ ) {
                sssqry = t1;
                ssssub = t2;
                sspqry = logo_fst_.GetSSStateProbAt(n, t1 );
                sspsub = logo_sec_.GetSSStateProbAt(m, t2 );
                sc = SSSSCORES.GetScore( sssqry, ssssub, sspqry, sspsub, ensqry, enssub );
                if( SSSSCORES.ScoreNA( sc ))
                    return 0.0;
                score += sspqry * sspsub * sc;
            }
    }
    else
        score = SSSSCORES.GetScore( sssqry, ssssub, sspqry, sspsub, ensqry, enssub );

    if( SSSSCORES.ScoreNA( score ))
        //score NA
        return 0.0;

//     score *= 0.4+0.6*exp( -0.5*3.e-5*(
//         zzicov*pow(ensqry-16.,4.) + zzicov*pow(enssub-16.,4.) + 
//             2.*zwicov*SQUARE(ensqry-16.)*SQUARE(enssub-16.)));


    if( sssqry == ssssub  &&  sssqry != SS_C )
        score += (sspqry+.1)*(sspsub+.1) / 9.;//90.;

    if( modscore )
        *modscore = score - offset;

    score -= offset;

    if( score )
        score *= scwght;
    if( modscore && *modscore )
        *modscore *= scwght;

    return score;
}

// -------------------------------------------------------------------------
// ComputeSSSScore: compute SS state score given positions of query and 
//  subject
//
double LSOSMatrix::RegulateHDPSSSS( int m, int n, double* modscore ) const
{
    double  scwght = 
#if defined( USEiHDPSSSSCORES )
        iHDPSSSSCORES.GetiHDPSSSWeight();
#else
        HDPSSSSCORES.GetHDPSSSWeight();
#endif
    if( modscore )
        *modscore = 0.0;
    if( scwght <= 0.0 )
        return 0.0;
    if( !logo_fst_.GetSSSSet() || !logo_sec_.GetSSSSet())
        return 0.0;

    const double*   ppsqry = logo_fst_.GetPPProbsAt(n);
    const int*      nxsqry = logo_fst_.GetPPPIndsAt(n);
    size_t          npsqry = logo_fst_.GetNoPPProbsAt(n);
    double          prbqry = 0.0;//cluster membership probability for query
    int             dishqr = 0;
    //
    const char      sssqry = logo_fst_.GetSSStateAt(n);
    const double    sspqry = logo_fst_.GetSSStateProbAt(n);
    const double    ensqry = logo_fst_.GetEffNoSequences();
    ////////
    const double*   ppssub = logo_sec_.GetPPProbsAt(m);
    const int*      nxssub = logo_sec_.GetPPPIndsAt(m);
    size_t          npssub = logo_sec_.GetNoPPProbsAt(m);
    double          prbsub = 0.0;//cluster membership probability for sbjct
    int             dishsb = 0;
    //
    const char      ssssub = logo_sec_.GetSSStateAt(m);
    const double    sspsub = logo_sec_.GetSSStateProbAt(m);
    const double    enssub = logo_sec_.GetEffNoSequences();
    double          score = 0.0, tmpval;
    ////////
    double          zzicov =  0.47;
    double          zwicov =  0.47;
    ////////

    if( SS_NSTATES <= sssqry || SS_NSTATES <= ssssub )
        throw myruntime_error("LSOSMatrix: RegulateHDPSSSS: Invalid SS class.");

    if( npsqry < 1 || npssub < 1 )
        return 0.0;

    prbqry = ppsqry[0];
    prbsub = ppssub[0];

    dishqr = nxsqry[0];
    dishsb = nxssub[0];

    //if priors
    if( dishqr < 0 || dishsb < 0 )
        return 0.0;

//     if( nodshs <= dishqr || nodshs <= dishsb )
//         throw myruntime_error("LSOSMatrix: RegulateHDPSSSS: Invalid dish index.");

#if defined( USEiHDPSSSSCORES )
    score = iHDPSSSSCORES.GetScore( dishqr, sssqry, prbqry, sspqry, ensqry );
    if( iHDPSSSSCORES.ScoreNA( score ))
        return 0.0;
    tmpval = iHDPSSSSCORES.GetScore( dishsb, ssssub, prbsub, sspsub, enssub );
    if( iHDPSSSSCORES.ScoreNA( tmpval ))
        return 0.0;
    score += tmpval;
#else
    score = HDPSSSSCORES.GetScore( 
        dishqr, dishsb, 
        sssqry, ssssub, 
        prbqry, prbsub, 
        sspqry, sspsub, 
        ensqry, enssub 
    );
    if( HDPSSSSCORES.ScoreNA( score ))
        //score NA
        return 0.0;
#endif

//     score *= 0.4+0.6*exp( -0.5*3.e-5*(
//         zzicov*pow(ensqry-23.,4.) + zzicov*pow(enssub-23.,4.) +
//             2.*zwicov*SQUARE(ensqry-23.)*SQUARE(enssub-23.)));

    if( score )
        score *= scwght;
    if( modscore )
        *modscore = score;

    return score;
}

// -------------------------------------------------------------------------
// ComputeCVS2ScoreHelper: compute log odds score of normal vectors
//
double LSOSMatrix::ComputeCVS2ScoreHelper( int m, int n, bool sson ) const
{
    const char      lppqry = logo_fst_.GetCtxVecLpprobAt(n);
    const double    nm2qry = logo_fst_.GetCtxVecNorm2At(n);
    const double*   vecqry = logo_fst_.GetCtxVecAt(n);
    const int       vszqry = logo_fst_.GetCtxVecSize();
    ////
    const char      lppsub = logo_sec_.GetCtxVecLpprobAt(m);
    const double    nm2sub = logo_sec_.GetCtxVecNorm2At(m);
    const double*   vecsub = logo_sec_.GetCtxVecAt(m);
    const int       vszsub = logo_sec_.GetCtxVecSize();
    double          score = 0.0;
    double          dot = 0.0, det = 0.0;
    int             i;
    ////

    if( vecqry == NULL || vecsub == NULL )
        throw myruntime_error("LSOSMatrix: ComputeCVS2ScoreHelper: Null vectors.");
    if( vszqry != vszsub || vszqry != CVS.DIM )
        throw myruntime_error("LSOSMatrix: ComputeCVS2ScoreHelper: Invalid vector sizes.");

    for( i = 0; i < vszqry; i++ )
        dot += vecqry[i] * vecsub[i];

    dot += CVS.loKAPPA0;//add off-diagonal entry of C_2

    //calculate det(C_2+H'H), where C_2 shifted centering matrix, H represents observations
    det = (nm2qry+1.+CVS.loKAPPA0)*(nm2sub+1.+CVS.loKAPPA0) - dot*dot;
    if( det <= 0.0 )
        throw myruntime_error("LSOSMatrix: ComputeCVS2ScoreHelper: Non-positive determinant.");

    //log odds score of normal vectors
    score = CVS.CTERM + CVS.PowerNU0*log(det) - lppqry - lppsub;

    return score;
}

// -------------------------------------------------------------------------
// ComputeCVS2Score: translate log odds score of normal vectors to alignment 
//  score
//
double LSOSMatrix::ComputeCVS2Score( int m, int n, bool sson, double* modscore ) const
{
    double  scwght = CVS2SCORES.GetCVSWeight();
    double  score = 0.0;

    if( modscore )
        *modscore = 0.0;
    if( scwght <= 0.0 )
        return 0.0;
    if( !logo_fst_.GetCtxVecSet() || !logo_sec_.GetCtxVecSet())
        return 0.0;

    static const double  offset = 0.5;///0.3;//0.0;
    const double    ensqry = logo_fst_.GetEffNoSequences();
    const double    enssub = logo_sec_.GetEffNoSequences();
//     double          zzicov =  0.47;
//     double          zwicov =  0.47;
    ////
    double  s;
    int i, j, k;

    if( 1 < CVS.AVGLEN && CVS.AVGCWGT < 1.0 && 0.0 < CVS.AVGCWGT ) {
        if( cvs2scores_ == NULL )
            throw myruntime_error("LSOSMatrix: ComputeCVS2Score: Memory access error.");
        //fill the table on the fly while calculating log odds scores
        for( i=m-CVS.hAVGLEN, j=n-CVS.hAVGLEN, k=0; k < CVS.AVGLEN; k++, i++, j++ ) {
            if( i < 0 || j < 0 )
                continue;
            if( GetSubjectSize() <= i || GetQuerySize() <= j )
                break;
            if( cvs2scores_[i][j])
                s = cvs2scores_[i][j];
            else {
                s = ComputeCVS2ScoreHelper( i, j, sson );
                cvs2scores_[i][j] = s;
            }
            score += s * CVS.AVGCOEFFS.GetCoefficientAt(k);
        }
    }
    else
        score = ComputeCVS2ScoreHelper( m, n, sson );

    //translate score
    score = CVS2SCORES.GetScore( score, ensqry, enssub );

    if( !score )
        return 0.0;

//     score *= 0.4+0.6*exp( -0.5*3.e-5*(
//         zzicov*pow(ensqry-23.,4.) + zzicov*pow(enssub-23.,4.) + 
//             2.*zwicov*SQUARE(ensqry-23.)*SQUARE(enssub-23.)));

    if( modscore )
        *modscore = score - offset;

//     score -= offset;

    if( score )
        score *= scwght;
    if( modscore && *modscore )
        *modscore *= scwght;

    return score;
}

// -------------------------------------------------------------------------
// ComputeScore: compute score given positions of query and subject
//
double LSOSMatrix::ComputeScore( int m, int n, bool final, double* modscore ) const
{
    double  hsc, ssssc, tmpv = 0.0;
    double  scwght = SSSSCORES.GetSSSWeight();
    bool    usingss = logo_fst_.GetSSSSet() && logo_sec_.GetSSSSet() &&( 0.0 < scwght );
    double  score = 0.0;
    double  offset = 0.0003;
    double  scscale = 0.95;//score scale
    int a;

    if( final ) {
    }

    for( a = 0; a < NUMAA; a++ )
        score += pqry_[n][a]*psub_[m][a];

    score = log(score) - offset;

    if( usingss ) {
        //score *= scscale;
        tmpv = 0.0;
        ssssc = ComputeSSSScore( m, n, GetUseModImage()? &tmpv: NULL );
        if( ssssc )
            score += ssssc;
        if( modscore && tmpv )
            *modscore += tmpv;
        tmpv = 0.0;
        hsc = RegulateHDPSSSS( m, n, NULL );
        if( hsc )
            score += hsc;
        if( modscore && tmpv )
            *modscore += tmpv;
    }
    {   tmpv = 0.0;
        hsc = ComputeHDPScore( m, n, usingss, NULL );
        if( hsc )
            score += hsc;
        if( modscore && tmpv )
            *modscore += tmpv;
    }
    {   tmpv = 0.0;
        hsc = ComputeHDPctScore( m, n, usingss, NULL );
        if( hsc )
            score += hsc;
        if( modscore && tmpv )
            *modscore += tmpv;
    }
    {   tmpv = 0.0;
        hsc = ComputeCVS2Score( m, n, usingss, GetUseModImage()? &tmpv: NULL );
        if( hsc )
            score += hsc;
        if( modscore && tmpv )
            *modscore += tmpv;
    }
    if( score )
        score *= GetMultiplier();
    if( modscore )
        *modscore *= GetMultiplier();
    return score;
}

// -------------------------------------------------------------------------
// ComputeProfileScoringMatrix: compute profile-profile scoring matrix
//
void LSOSMatrix::ComputeProfileScoringMatrix( bool final )
{
    if( GetQuerySize() <= 0 || GetSubjectSize() <= 0 )
        throw myruntime_error( mystring( "LSOSMatrix: Unable to compute matrix of scores." ));

    if( !IsValid())
        throw myruntime_error( mystring( "LSOSMatrix: Failed to compute matrix of scores." ));

    int     m, n;
    int     no_elems = 0;
    int     negatives = 0;
    TMask   currentmask = Unmasked;
    double  score = 0.0;
    double  modscore = 0.0;


    //fill matrix with values
    for( m = 0; m < GetSubjectSize(); m++ ) {
        for( n = 0; n < GetQuerySize(); n++ )
        {
            score = 0.0;
            modscore = 0.0;
            currentmask = Unmasked;

            if( logo_fst_.GetInformationAt( n ) < GetInformationThreshold() ||
                logo_sec_.GetInformationAt( m ) < GetInformationThreshold() ||
                ! ThicknessConstraintsMet( m, n ))
            {
                //mask the score possibly excluding it from the ongoing computational statistics
                SetMasked( currentmask = GetMaskingApproach(), m, n );
            }

            if( 1 || currentmask != MaskToIgnore )
            {
                score = ComputeScore( m, n, final, GetUseModImage()? &modscore: NULL );

                if( currentmask != MaskToIgnore ) {
                    no_elems++;
                    if( score < 0.0 ) negatives++;    //increase counter for each negative score found
                }
                if( currentmask == MaskToConsider || currentmask == MaskToIgnore ) {
                    score *= GetMaskscalePercents();
                    modscore *= GetMaskscalePercents();
                }
            }
            else {
                //penalize score at the positions;
                score = g_scoreX;
                modscore = g_scoreX;
            }

            SetScore( score, m, n );
            if( GetUseModImage())
                SetModScore( modscore, m, n );
        }
    }


    if( negatives == no_elems ) {
        SetAllNegatives( true );
        return;
    }
}

// -------------------------------------------------------------------------
// ComputeExpectation: compute e-value given calculated statistical params
//
double LSOSMatrix::ComputeExpectation(
        double score,
        double* ref_expect,
        double* pure_expect,
        double* pair_expect,
        double* bitscore ) const
{
    double  expect = -1.0;
    return AbstractScoreMatrix::ComputeExpectation( score, ref_expect, pure_expect, pair_expect, bitscore );
}

// -------------------------------------------------------------------------
// ComputePositionalScoreProbs: compute probabilities of scores at each
//     query and subject position
//
void LSOSMatrix::ComputePositionalScoreProbs( AttributableScores* PATTR_SCORES )
{
    if( PATTR_SCORES == NULL )
        throw myruntime_error( mystring( "LSOSMatrix: Unable to compute probabilities: Wrong argument." ));

    double  sum = 0.0;          //initial sum
    double  consval = 0.0;      //probability conservation value
    double  exp2sum = 0.0;      //e to sum
    double  normterm = 0.0;     //normalizing term
    double* querynorm = NULL;   //array of normalizing terms
    TScore  loc_score, sc;
    int     m, n;

    if( 0 < GetQuerySize()) {
        querynorm = ( double* )malloc( sizeof( double ) * GetQuerySize());
        if( querynorm == NULL )
            throw myruntime_error( mystring( "LSOSMatrix: ComputePositionalScoreProbs: Not enough memory." ));
        for( n = 0; n < GetQuerySize(); n++ )
            querynorm[n] = 0.0;
    }
    try{
        for( m = 0; m < GetSubjectSize(); m++ ) {

            normterm = 0.0;
            for( n = 0; n < GetQuerySize(); n++ )
            {
                sum = 0.0;

                if( PATTR_SCORES->GetMaskedToIgnore( m, n ))
                    continue;

                loc_score = ( TScore )rint( PATTR_SCORES->GetScore( m, n ));

                if( loc_score <= SCORE_MIN )
                    continue;

                exp2sum = 1.0;
                PATTR_SCORES->IncQueryInfProbOf( exp2sum, loc_score, n );
                PATTR_SCORES->IncSbjctInfProbOf( exp2sum, loc_score, m );
                normterm += exp2sum;
                querynorm[n] += exp2sum;
            }
            if( normterm ) {
                consval = 0.0;
                PATTR_SCORES->SetSbjctInfProbNorm( normterm, m );
                //normalize probabilities for query positions
                for( sc = PATTR_SCORES->GetSbjctMinScoreAt( m ); sc <= PATTR_SCORES->GetSbjctMaxScoreAt( m ); sc++ )
                    consval += PATTR_SCORES->DivSbjctInfProbOf( normterm, sc, m );
                if( consval < 0.999 || consval > 1.001 )
                    warning( "LSOSMatrix: ComputePositionalScoreProbs: Probabilities not conserved." );
            }
        }

        for( n = 0; n < GetQuerySize(); n++ )
        {   //normalize probabilities
            if( querynorm[n] ) {
                consval = 0.0;
                PATTR_SCORES->SetQueryInfProbNorm( querynorm[n], n );
                for( sc = PATTR_SCORES->GetQueryMinScoreAt( n ); sc <= PATTR_SCORES->GetQueryMaxScoreAt( n ); sc++ )
                    consval += PATTR_SCORES->DivQueryInfProbOf( querynorm[n], sc, n );
                if( consval < 0.999 || consval > 1.001 )
                    warning( "LSOSMatrix: ComputePositionalScoreProbs: Probabilities not conserved." );
            }
        }
    } catch( myexception const& ex ) {
        if( querynorm ) free( querynorm );
        throw myruntime_error( ex.what(), ex.eclass());
    }

    if( querynorm ) free( querynorm );
}

// -------------------------------------------------------------------------
// ComputeScoreProbabilities: compute score probabilities;
//  score probabilities are all equally expected according to the LSOS model
//
void LSOSMatrix::ComputeScoreProbabilities( AttributableScores* PATTR_SCORES )
{
    if( PATTR_SCORES == NULL )
        throw myruntime_error( mystring( "LSOSMatrix: Unable to compute probabilities: Wrong argument." ));

    const double    accuracy = 1.0e-4;

    int     n, m, nn, mm, r;
    char    strbuf[KBYTE];
    double  sum = 0.0;      //initial sum
    double  proquery = 0.0; //non-normalized probability of query position
    double  prosbjct = 0.0; //non-normalized probability of subject position
    double  exp2sum = 0.0;  //e to sum
    double  normprob = 0.0; //normalizing term
    double  normterm = 0.0; //normalizing term
    double  avgscore = 0.0;
    double  consv;
    TScore  loc_score  = 0;

    scores_.Clear();
    rprobs_.Clear();
    cprobs_.Clear();


    for( n = 0, nn = 0; n < GetQuerySize(); n++ ) {
        if(0) {
            sum = 0.0;
            for( r = 0; r < NUMAA; r++ ) {
                if( LOSCORES.PROBABility( r ) <= 0.0 )
                    continue;
                sum += freq_fst_( n, r ) * LOSCORES.LogPROBABILITY_1( r );
            }
            rprobs_.AddValueAt( nn++, exp( sum ));
        } else
            rprobs_.SetValueAt( nn++, 1.0 );
    }

    for( m = 0, mm = 0; m < GetSubjectSize(); m++ ) {
        if(0) {
            sum = 0.0;
            for( r = 0; r < NUMAA; r++ ) {
                if( LOSCORES.PROBABility( r ) <= 0.0 )
                    continue;
                sum += freq_sec_( m, r ) * LOSCORES.LogPROBABILITY_2( r );
            }
            cprobs_.AddValueAt( mm++, exp( sum ));
        } else
            cprobs_.SetValueAt( mm++, 1.0 );
    }
    
    avgscore = normterm = 0.0;
    for( n = 0, nn = 0; n < GetQuerySize(); n++ )
    {
        proquery = rprobs_.GetValueAt( nn++ );

        for( m = 0, mm = 0; m < GetSubjectSize(); m++ )
        {
            prosbjct = cprobs_.GetValueAt( mm++ );

            loc_score = ( TScore )rint( PATTR_SCORES->GetScore( m, n ));

            exp2sum = 1.;//proquery * prosbjct;//same

            if( loc_score <= SCORE_MIN )
                    scores_.Push( g_scoreX );
            else    scores_.Push( loc_score );

            if( loc_score <= SCORE_MIN )
                continue;

            if( PATTR_SCORES->GetMaskedToIgnore( m, n ))
                continue;

            PATTR_SCORES->IncProbabilityOf( exp2sum, loc_score );
            normterm += exp2sum;
        }
    }
    //normalize probabilities
    if( normterm )
        for( int sc = PATTR_SCORES->GetMinScore(); sc <= PATTR_SCORES->GetMaxScore(); sc++ ) {
            PATTR_SCORES->DivideProbabilityOf( normterm, sc );
            avgscore += sc * PATTR_SCORES->GetProbabilityOf( sc );
        }

    normprob = rprobs_.Sum();
    if( 0.0 < normprob ) {
        rprobs_.MultiplyBy( 1.0 / normprob );
        consv = rprobs_.Sum();
        if( consv < 1.0 - accuracy || consv > 1.0 + accuracy ) {
            sprintf( strbuf, "Background probabilities not conserved: %f", consv );
            throw myruntime_error( strbuf );
        }
    }
    normprob = cprobs_.Sum();
    if( 0.0 < normprob ) {
        cprobs_.MultiplyBy( 1.0 / normprob );
        consv = cprobs_.Sum();
        if( consv < 1.0 - accuracy || consv > 1.0 + accuracy ) {
            sprintf( strbuf, "Background probabilities not conserved: %f", consv );
            throw myruntime_error( strbuf );
        }
    }

    PATTR_SCORES->SetExpectedScore( avgscore );
}

// -------------------------------------------------------------------------
// OptimizeTargetFrequencies: Optimize target frequencies
//
void LSOSMatrix::OptimizeTargetFrequencies()
{
    double  lambda;

    if( GetAllNegatives()) {
        //warning( "Unable to optimize target frequencies." );
        return;
    }

    ComputeStatisticalParameters();

    if( scores_.GetSize() < 1 || rprobs_.GetSize() < 1 || cprobs_.GetSize() < 1 )
        return;

    if( 0.0 <= GetExpectedScore() || GetLambda() <= 0.0 ) {
        warning( "Unable to optimize target frequencies." );
        return;
    }

    lambda = GetLambda();
    if( GetAutoScaling())
        lambda = GetScaledScores()->GetLambda();

    scores_.MultiplyBy( lambda );


    TargetFreqOptimizerH    optimizer( scores_, rprobs_, cprobs_ );
    int                     status;

    //no need to set lambda if scores were previously multiplied by
//     optimizer.SetLambda( GetLambda());
    optimizer.SetConstrainedH( GetRefH());
    optimizer.SetNMResidTol( TARGET_FREQ_OPT_TOLERANCE );
    optimizer.SetMaxNoNMIterations( TARGET_FREQ_OPT_MAXITERATIONS );

    status = optimizer.Optimize();

    if( status != PSL_SUCCESS /*&& status != PSL_MAXITERATS*/ ) {
#ifdef TFOPTESTPRINT
        warning( TranslatePSLError( status ));
#else
        warning( "Target frequencies not optimized." );
#endif
        return;
    }

    if( status != PSL_SUCCESS )
        optimizer.Normalize();

    if( optimizer.Negatives()) {
        warning( "Invalid target frequencies after optimization." );
        return;
    }

    IncorporateTargetFrequencies( optimizer.GetOptTargetFreqs());
}

// -------------------------------------------------------------------------
// IncorporateTargetFrequencies: Calculate new scores based on optimized
//     target frequencies
//
void LSOSMatrix::IncorporateTargetFrequencies( const Pslvector& optfreqs )
{
    if( rprobs_.GetSize() < 1 || cprobs_.GetSize() < 1 ||
        optfreqs.GetSize() != rprobs_.GetSize() * cprobs_.GetSize())
        return;

    const double    accuracy = 1.0e-6;
    int     n, m, nn, mm, ind;
    double  tfval = 0.0;    //target frequency value
    double  proquery = 0.0; //background probability of query position
    double  prosbjct = 0.0; //background probability of subject position
    double  loc_score = 0.0;
    double  score = 0.0;
    double  consv;

    consv = optfreqs.Sum();
    if( consv < 1.0 - accuracy || consv > 1.0 + accuracy ) {
        warning( "Target frequencies after optimization not conserved." );
        return;
    }

    ind = 0;
    for( n = 0, nn = 0; n < GetQuerySize(); n++ )
    {
        proquery = rprobs_.GetValueAt( nn++ );

        for( m = 0, mm = 0; m < GetSubjectSize(); m++ )
        {
            prosbjct = cprobs_.GetValueAt( mm++ );

            tfval = optfreqs.GetValueAt( ind++ );

            if( proquery <= 0.0 || prosbjct <= 0.0 )
                throw myruntime_error( "Invalid probabilities after optimizing of target frequencies." );

            score = log( tfval / ( proquery * prosbjct ));

            loc_score = GetScore( m, n );

            if( loc_score <= SCORE_MIN )
                continue;

//             if( GetMaskedToIgnore( m, n ))
//                 continue;

            SetScore( score, m, n );
        }
    }
}



// =========================================================================
// PRINT ROUTINES
//
// PrintParameterTable: print table of computed statistical parameter
//     values to stream
//
void LSOSMatrix::PrintParameterTable( TPrintFunction print_func, void* vpn ) const
{
    if( vpn == NULL )
        return;

    if( 0.0 <= GetExpectedScore()) {
        print_func( vpn, "Expected score per position is non-negative, %.4f!\n\n", GetExpectedScore());
        return;
    }

    char    _sk[BUF_MAX], _slambda[BUF_MAX];
    char    _sdk[BUF_MAX], _sdlambda[BUF_MAX];
    if( 0.0 < GetK())   sprintf( _sk, "%6.4f", GetK());
        else            sprintf( _sk, "n/a" );
    if( 0.0 < GetLambda())  sprintf( _slambda, "%6.4f", GetLambda());
        else                sprintf( _slambda, "n/a" );
    if( 0.0 < GetDerivedGappedK())  sprintf( _sdk, "%6.4f", GetDerivedGappedK());
        else                        sprintf( _sdk, "n/a" );
    if( 0.0 < GetDerivedGappedLambda()) sprintf( _sdlambda, "%6.4f", GetDerivedGappedLambda());
        else                            sprintf( _sdlambda, "n/a" );

    print_func( vpn, "%-25s  %-6s   %-6s\n", " ", "K", "Lambda" );
//    print_func( vpn, "%-25s  %6.4f   %6.4f\n", "Reference ungapped,", GetRefK(),            GetRefLambda());
//    print_func( vpn, "%-25s  %6.4f   %6.4f\n", "Reference gapped,",   GetExpGappedK(),      GetExpGappedLambda());
    print_func( vpn, "%-25s  %6s   %6s\n", "Computed  ungapped,", _sk, _slambda );
    print_func( vpn, "%-25s  %6s   %6s\n", "Estimated gapped,", _sdk, _sdlambda );

    print_func( vpn, "Entropy, %6.4f; Expected, %6.4f; Min/Max, %d/%-d\n\n",
                GetEntropy(), GetExpectedScore(), GetMinScore(), GetMaxScore());
}

// PrintFinal: nothing to print in the final
//
void LSOSMatrix::PrintFinal( TPrintFunction print_func, void* vpn ) const
{
    if( vpn == NULL )
        return;

    PrintReferenceParameterTable( print_func, vpn );
}

// -------------------------------------------------------------------------
// PrintScoringMatrix: output computed scoring system
//
void LSOSMatrix::PrintScoringMatrix( FILE* fp )
{
    int l = 0;
    int m, n;
    const double* const queryentrops = GetQueryInfContent();
    const double* const sbjctentrops = GetSbjctInfContent();

    if( fp == NULL )
        return;

    fprintf( fp,"%12c Position-specific profile scoring matrix\n", 32 );

    fprintf( fp, "%9c", 32 );
    for( m = 0; m < GetSubjectSize(); m++ )
        fprintf( fp, " %4d", m + 1 );

    fprintf( fp, "\n%9c", 32 );
    for( m = 0; m < GetSubjectSize(); m++ )
        fprintf( fp, " %4c", DehashCode( logo_sec_[m] ));

    for( n = 0; n < GetQuerySize(); n++ ) {

        fprintf( fp, "\n%5d %c   ", ++l, DehashCode( logo_fst_[n] ));

        for( m = 0; m < GetSubjectSize(); m++ )
            fprintf( fp, "%4d ", ( int )rint( GetImageScore/*GetScoreBase*/( m, n ) ));

        if( queryentrops )
            fprintf( fp, " %4.2f", queryentrops[n] );
    }

    if( sbjctentrops ) {
        fprintf( fp, "\n\n          " );
        for( m = 0; m < GetSubjectSize(); m++ )
            fprintf( fp, "%4.2f ", sbjctentrops[m] );
    }

    fprintf( fp, "\n" );
}

