/***************************************************************************
 *   Copyright (C) 2012 by Mindaugas Margelevicius                         *
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
#include "libpro/tfopt/TargetFreqOptimizerH.h"
#include "HDP0ScoreMatrix.h"


FrequencyMatrix smp_hdpdummyfreq_;
LogOddsMatrix   smp_hdpdummylogo_;

const double    g_scoreX = 0.0; //-1.0;

//if true, use matrix variate t-distribution
bool        HDP0ScoreMatrix::mtxvart_ = true;
//if true, background probabilities are multivariate t-dist.
bool        HDP0ScoreMatrix::hdpbprobs_ = false;
//default basin size
const int   HDP0ScoreMatrix::s_defbasize = 50;

//local file constants
double HDPSM_kappa0 = 1.0;//10.0;
double HDPSM_Qdegfr = 0.0;//21.0;


// -------------------------------------------------------------------------
// constructor: frequency matrices, log-odds matrices for the first and
//   second profiles are given with the parameters
//
HDP0ScoreMatrix::HDP0ScoreMatrix(
    const HDPbase*  hdpparent,
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
    AbstractScoreMatrix( HDPProfileSpecific, config, s_behaviour, a_scaling, c_masking ),

    hdpparent_( hdpparent ),
    freq_fst_( freq_fst ),
    logo_fst_( logo_fst ),

    freq_sec_( freq_sec ),
    logo_sec_( logo_sec ),

    thickness_number( thick_number ),
    thickness_percnt( thick_percnt ),
    maskscale_percnt( mask_percnt ),

    scores_( freq_fst_.GetColumns() * freq_sec_.GetColumns()),
    rprobs_( freq_fst_.GetColumns()),
    cprobs_( freq_sec_.GetColumns()),

    hdpbase_( NULL ),
    quebasin_( NULL ),
    subbasin_( NULL )
{
    PrivateInit();
    SetInformationThreshold( infrm_threshold );
    SetScoreAdjCorrel( true );
    SetSupportOptimFreq( false );
}

// -------------------------------------------------------------------------
// constructor: overloaded
//
HDP0ScoreMatrix::HDP0ScoreMatrix(
    const HDPbase*  hdpparent,
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

    hdpparent_( hdpparent ),
    freq_fst_( freq_fst ),
    logo_fst_( logo_fst ),

    freq_sec_( freq_sec ),
    logo_sec_( logo_sec ),

    thickness_number( thick_number ),
    thickness_percnt( thick_percnt ),
    maskscale_percnt( mask_percnt ),

    scores_( freq_fst_.GetColumns() * freq_sec_.GetColumns()),
    rprobs_( freq_fst_.GetColumns()),
    cprobs_( freq_sec_.GetColumns()),

    hdpbase_( NULL ),
    quebasin_( NULL ),
    subbasin_( NULL )
{
    PrivateInit();
    SetInformationThreshold( infrm_threshold );
}

// -------------------------------------------------------------------------
// default constructor is invalid for initialization
//
HDP0ScoreMatrix::HDP0ScoreMatrix()
:
    AbstractScoreMatrix(),

    hdpparent_( NULL ),
    freq_fst_( smp_hdpdummyfreq_ ),
    logo_fst_( smp_hdpdummylogo_ ),

    freq_sec_( smp_hdpdummyfreq_ ),
    logo_sec_( smp_hdpdummylogo_ ),

    thickness_number( 0 ),
    thickness_percnt( 0.0 ),
    maskscale_percnt( 0.0 ),

    scores_( 1 ),
    rprobs_( 1 ),
    cprobs_( 1 ),

    hdpbase_( NULL ),
    quebasin_( NULL ),
    subbasin_( NULL )
{
}

// -------------------------------------------------------------------------
// destructor: deallocation of resources
//
HDP0ScoreMatrix::~HDP0ScoreMatrix()
{
    DestroyHDPbase();
    DestroyQueryBasin();
    DestroySbjctBasin();
}

// -------------------------------------------------------------------------
// PrivateInit: initialization method private for this class
//
void HDP0ScoreMatrix::PrivateInit()
{
    if( !logo_fst_.IsCompatible( freq_fst_ ) || !logo_sec_.IsCompatible( freq_sec_ ))
            throw myruntime_error( mystring( "HDP0ScoreMatrix: Profile matrices are incompatible." ));

    if(  freq_fst_.GetColumns() < 1 || freq_fst_.GetColumns() > MAXCOLUMNS ||
         freq_sec_.GetColumns() < 1 || freq_sec_.GetColumns() > MAXCOLUMNS )
            throw myruntime_error( mystring( "HDP0ScoreMatrix: Wrong profile matrices." ));

    if( GetMaskscalePercents() < 0.0 )
            throw myruntime_error( mystring( "HDP0ScoreMatrix: Negative scaling factor of masked positions." ));

    Init( freq_fst_.GetColumns(), freq_sec_.GetColumns());

#ifdef USEPROFBACKPROBS
    LOSCORES.StoreProbabilities_1( logo_fst_.GetBackProbs());
    LOSCORES.StoreProbabilities_2( logo_sec_.GetBackProbs());
#endif
}


// -------------------------------------------------------------------------
// InitHDPVecobs: initialize hdp-related structures; 
//  overloaded for making alignments symmetric
//
void HDP0ScoreMatrix::InitHDPVecobs( const LogOddsMatrix& fst, const LogOddsMatrix& sec,
                              Pslvector& fstvec, Pslvector& secvec,
                              Basin& fstbas, Basin& secbas, HDPbase& basHDP )
{
    //reference lambda
    const double    lambda = LOSCORES.StatisParam( Ungapped, Lambda );
    const double    tol = 0.1;//error tolerance for target probabilities
    const int       dim = NUMAA - 1;
    const int       upperdof = 20;//upper limit for degrees of freedom
    const int       quesize = fst.GetColumns();
    const int       subsize = sec.GetColumns();
    const double    queryenos = fst.GetEffNoSequences();
    const double    sbjctenos = sec.GetEffNoSequences();
    const double    querydenm = 19.0;//queryenos /( double )upperdof;
    double  tprob;//target probability from profile
    double  pscor;//score from profile
    double  lp;//log-probability
    Dish*   dsh;
    Pslvector*  ps, *bps;
    bool    added = false;//noise added
    int propsize;//subject proposal size
    int bnpos, npos, kpos, dpos;
    int n, m, a, r;

    if( quesize < 1 || subsize < 1 )
        return;

    //one position is for background prob. vector
    basHDP.ReserveMenu( quesize + subsize + 1 );
    basHDP.ReserveBasin( quesize + subsize + 1 );
    basHDP.GetBasin()->SetDestroy( false );

    for( n = 0; n < quesize; n++ ) {
        ps = new Pslvector( NUMAA+0 );
        if( ps == NULL )
            throw myruntime_error("HDP0ScoreMatrix: InitHDPVec: Not enough memory.");
        for( a = 0; a < NUMAA; a++ ) {
            pscor = fst( n, a );
            tprob = exp( pscor * lambda ) * LOSCORES.PROBABility( a );
            ps->SetValueAt( a, tprob );
        }
//         ps->SetValueAt( NUMAA, 1.e-3 );
        HDPbase::LogitNormal2Normal( ps, tol );
        //reduce dimensionality
        ps->DecDim();
        bnpos = fstbas.NewValue( ps );
//         npos = GetHDPbase()->AddToBasin( ps );
    }

    for( m = 0; m <= subsize; m++ ) {
        dsh = new Dish( 1 );
        ps = new Pslvector( NUMAA+0 );
        if( ps == NULL || dsh == NULL )
            throw myruntime_error("HDP0ScoreMatrix: InitHDPVec: Not enough memory.");
//         dsh->SetBasin( GetHDPbase()->GetBasin());
        dsh->SetBasin( &secbas );
        for( a = 0; a < NUMAA; a++ ) {
            tprob = LOSCORES.PROBABility( a );
            //the last position is for background probability vector
            if( m < subsize ) {
                pscor = sec( m, a );
                tprob *= exp( pscor * lambda );
            }
            ps->SetValueAt( a, tprob );
        }
//         ps->SetValueAt( NUMAA, 1.e-3 );
        HDPbase::LogitNormal2Normal( ps, tol );
        //reduce dimensionality
        ps->DecDim();
        bnpos = secbas.NewValue( ps );
        npos = basHDP.AddToBasin( ps );
        kpos = basHDP.GetMenu()->NewDish( dsh );
//         dpos = dsh->NewVectorNInd( npos );
        if( m < subsize )
//             propsize = ( int )SLC_MIN(
//                 SLC_MAX( rint( sec.GetMIDExpNoObservationsAt(m,PS_M)/ querydenm), 1.0), 21 );
            propsize = 1;
        else
//             propsize = ( int )SLC_MIN( SLC_MAX( rint( sbjctenos / querydenm ), 1.0 ), 21 );
            propsize = 1;
        dpos = dsh->NewVectorNInd( bnpos );
        dsh->SetReadSize( propsize );
    }

    //{{prepare HDPbase
    basHDP.SetDegFAdjustment( HDPSM_Qdegfr );
    basHDP.GetMenu()->SetKappa0( HDPSM_kappa0 );//IMP.!
    basHDP.GetMenu()->SetNu0( HDPSM_kappa0 );
    basHDP.GetMenu()->SetDim( ps->GetSize());
    try {
        if( 0.0 < HDPSM_kappa0 )
            basHDP.CalcPriorParams();
        //reduce weight of prior parameters after calculating them
        basHDP.GetMenu()->SetKappa0( HDPSM_kappa0 );
        basHDP.GetMenu()->SetNu0( HDPSM_kappa0 );
        basHDP.GetMenu()->SetDim( ps->GetSize());
        //calculate parameters for each dish
        basHDP.RecalcMenuParams();
    } catch( myexception const& ex ) {
        //calculate prior parameters based on self profile vectors;
        //when based on both sets of vectors, noise does not imply
        //high fluctuations (even with amplification of kappa0), and
        //scores tend to be lower
//         GetHDPbase()->GetBasin()->Clear();
        for( n = 0; n < quesize; n++ ) {
            ps = fstbas.GetValueAt( n );
            if( ps == NULL )
                throw myruntime_error("HDP0ScoreMatrix: InitHDPVec: Memory access error.");
            npos = GetHDPbase()->AddToBasin( ps );
        }
        //add noise when vectors are linearly dependent
        added = basHDP.CalcPriorParams( true, 0.0001/*std. fact.*/);
        //scores are reduced when prior params. are calculated from self vectors
        if( 0&&added ) {
            //if kappa0 is small, scores in cases may increase substantially
            HDPSM_kappa0 = SLC_MAX( 0.0, HDPSM_kappa0 - 10.0 );
            HDPSM_Qdegfr = SLC_MAX( 0.0, HDPSM_Qdegfr - 10.0 );
        }
        basHDP.GetMenu()->SetKappa0( HDPSM_kappa0 );
        basHDP.GetMenu()->SetNu0( HDPSM_kappa0 );
        //recalculate dish parameters
        basHDP.RecalcMenuParams();
    }
    //}}

    //precalculate log-probabilities for query positions
    fstvec.Reserve( quesize );
    for( n = 0; n < quesize; n++ ) {
        ps = fstbas.GetValueAt( n );
        if( ps == NULL )
            throw myruntime_error("HDP0ScoreMatrix: InitHDPVec: Memory access error.");
//         basHDP.SetDegFAdjustment(
//             SLC_MIN( SLC_MAX( rint( fst.GetMIDExpNoObservationsAt( n, PS_M )/ querydenm ), 1.0), 21 ));
        basHDP.SetDegFAdjustment( HDPSM_Qdegfr );
//         basHDP.PriorProbVec( ps, &lp );//P(Q)-- prior probability
        basHDP.ProbVecOfDish( ps, subsize, &lp );//P(Q|B), where B is background vector
        fstvec.SetValueAt( n, lp );
    }

    if( GetUseHDPBProbs()) {
        //precalculate log-probabilities for subject positions
        secvec.Reserve( subsize );
        for( m = 0; m < subsize; m++ ) {
            ps = secbas.GetValueAt( m );
            if( ps == NULL )
                throw myruntime_error("HDP0ScoreMatrix: InitHDPVec: Memory access error.");
//             basHDP.SetDegFAdjustment(
//                 SLC_MIN( SLC_MAX( rint( sec.GetMIDExpNoObservationsAt( m, PS_M )/ querydenm ), 1.0), 21));
            basHDP.SetDegFAdjustment( 0.0 );
            basHDP.PriorProbVec( ps, &lp );
            secvec.SetValueAt( m, lp );
        }
    }
}

// -------------------------------------------------------------------------
// InitHDPVec: initialize hdp-related structures; 
//  overloaded for making alignments symmetric
//
void HDP0ScoreMatrix::InitHDPVec( const LogOddsMatrix& fst, const LogOddsMatrix& sec,
                              Pslvector& fstvec, Pslvector& secvec,
                              Basin& fstbas, Basin& secbas, HDPbase& basHDP )
{
    const HDPbase*  parentHDP = GetHDPparent();
    const double    lambda = LOSCORES.StatisParam( Ungapped, Lambda );//reference lambda
    const double    tol = 0.1;//error tolerance for target probabilities
    const int       dim = NUMAA - 1;
    const int       upperdof = 20;//upper limit for degrees of freedom
    const int       quesize = fst.GetColumns();
    const int       subsize = sec.GetColumns();
    const double    queryenos = fst.GetEffNoSequences();
    const double    sbjctenos = sec.GetEffNoSequences();
    const double    querydenm = 19.0;//queryenos /( double )upperdof;
    double  tprob;//target probability from profile
    double  pscor;//score from profile
    double  lp;//log-probability
    Dish*   dsh;
    Pslvector*  ps, *bps;
    bool    added = false;//noise added
    int propsize;//subject proposal size
    int bnpos, npos, kpos, dpos;
    int n, m, a, r;

    if( parentHDP == NULL )
        throw myruntime_error("HDP0ScoreMatrix: InitHDPVec: Null parent HDP object.");
    if( parentHDP->GetMenu() == NULL )
        throw myruntime_error("HDP0ScoreMatrix: InitHDPVec: Null parent HDP menu.");

    if( quesize < 1 || subsize < 1 )
        return;

    //one position is for background prob. vector
    basHDP.ReserveMenu( quesize + subsize + 1 );
    basHDP.ReserveBasin( quesize + subsize + 1 );
    basHDP.GetBasin()->SetDestroy( false );

    for( n = 0; n < quesize; n++ ) {
        ps = new Pslvector( NUMAA+0 );
        if( ps == NULL )
            throw myruntime_error("HDP0ScoreMatrix: InitHDPVec: Not enough memory.");
        for( a = 0; a < NUMAA; a++ ) {
            pscor = fst( n, a );
            tprob = exp( pscor * lambda ) * LOSCORES.PROBABility( a );
            ps->SetValueAt( a, tprob );
        }
//         ps->SetValueAt( NUMAA, 1.e-3 );
        HDPbase::LogitNormal2Normal( ps, tol );
        //reduce dimensionality
        ps->DecDim();
        bnpos = fstbas.NewValue( ps );
    }

    for( m = 0; m <= subsize; m++ ) {
        dsh = new Dish( 1 );
        ps = new Pslvector( NUMAA+0 );
        if( ps == NULL || dsh == NULL )
            throw myruntime_error("HDP0ScoreMatrix: InitHDPVec: Not enough memory.");
        dsh->SetBasin( &secbas );
        for( a = 0; a < NUMAA; a++ ) {
            tprob = LOSCORES.PROBABility( a );
            //the last position is for background probability vector
            if( m < subsize ) {
                pscor = sec( m, a );
                tprob *= exp( pscor * lambda );
            }
            ps->SetValueAt( a, tprob );
        }
//         ps->SetValueAt( NUMAA, 1.e-3 );
        HDPbase::LogitNormal2Normal( ps, tol );
        //reduce dimensionality
        ps->DecDim();
        bnpos = secbas.NewValue( ps );
        npos = basHDP.AddToBasin( ps );
        kpos = basHDP.GetMenu()->NewDish( dsh );
        propsize = 1;
        dpos = dsh->NewVectorNInd( bnpos );
        dsh->SetReadSize( propsize );
    }

    //{{prepare HDPbase
    basHDP.GetMenu()->SoftCopyPriors( *parentHDP->GetMenu());
    basHDP.SetDegFAdjustment( HDPSM_Qdegfr );
    //reduce weight of prior parameters after calculating them
    basHDP.GetMenu()->SetKappa0( HDPSM_kappa0 );//IMP.!
    basHDP.GetMenu()->SetNu0( HDPSM_kappa0 );
    basHDP.CalcPriorProbFact();
    //calculate parameters for each dish
    basHDP.RecalcMenuParams();
    //to leverage probability order
    for( n = 0; n < basHDP.GetMenu()->GetSize(); n++ )
        if( basHDP.GetMenu()->GetDishAt(n))
            basHDP.GetMenu()->GetDishAt(n)->SetReadSize(0);
    //}}

    //precalculate log-probabilities for query positions
    fstvec.Reserve( quesize );
    for( n = 0; n < quesize; n++ ) {
        ps = fstbas.GetValueAt( n );
        if( ps == NULL )
            throw myruntime_error("HDP0ScoreMatrix: InitHDPVec: Memory access error.");
        basHDP.SetDegFAdjustment( HDPSM_Qdegfr );
        basHDP.PriorProbVec( ps, &lp );//P(Q)-- prior probability
//         basHDP.ProbVecOfDish( ps, subsize, &lp );//P(Q|B), where B is background vector
        fstvec.SetValueAt( n, lp );
    }

    if( GetUseHDPBProbs()) {
        //precalculate log-probabilities for subject positions
        secvec.Reserve( subsize );
        for( m = 0; m < subsize; m++ ) {
            ps = secbas.GetValueAt( m );
            if( ps == NULL )
                throw myruntime_error("HDP0ScoreMatrix: InitHDPVec: Memory access error.");
            basHDP.SetDegFAdjustment( HDPSM_Qdegfr );
            basHDP.PriorProbVec( ps, &lp );
            secvec.SetValueAt( m, lp );
        }
    }
}

// -------------------------------------------------------------------------
// InitHDPVec: initialize hdp-related structures
//
void HDP0ScoreMatrix::InitHDPVec()
{
    const int       quesize = GetQuerySize();
    const int       subsize = GetSubjectSize();
    const double    queryenos = GetQueryLogo().GetEffNoSequences();
    const double    sbjctenos = GetSbjctLogo().GetEffNoSequences();

    if( quesize < 1 || subsize < 1 )
        return;

    InitHDPbase();
    InitQueryBasin();
    InitSbjctBasin();

    if( sbjctenos <= queryenos )
        InitHDPVec( GetQueryLogo(), GetSbjctLogo(),
                 GetQueValues(), GetSubValues(),
                 *GetQueryBasin(), *GetSbjctBasin(), *GetHDPbase());
    else
        InitHDPVec( GetSbjctLogo(), GetQueryLogo(),
                 GetSubValues(), GetQueValues(),
                 *GetSbjctBasin(), *GetQueryBasin(), *GetHDPbase());
}





// =========================================================================
// InitHDPMtx: initialize hdp-related structures
//
void HDP0ScoreMatrix::InitHDPMtx()
{
    const int       quesize = GetQuerySize();
    const int       subsize = GetSubjectSize();
    const double    queryenos = GetQueryLogo().GetEffNoSequences();
    const double    sbjctenos = GetSbjctLogo().GetEffNoSequences();

    if( quesize < 1 || subsize < 1 )
        return;

    InitHDPbase();
    InitQueryBasin();
    InitSbjctBasin();

    const LogOddsMatrix& fst = GetQueryLogo();
    const LogOddsMatrix& sec = GetSbjctLogo();
    Pslvector&  fstvec = GetQueValues();
    Pslvector&  secvec = GetSubValues();
    Basin&      fstbas = *GetQueryBasin();
    Basin&      secbas = *GetSbjctBasin();
    HDPbase&    basHDP = *GetHDPbase();

    const HDPbase*  parentHDP = GetHDPparent();
    const double    lambda = LOSCORES.StatisParam( Ungapped, Lambda );//reference lambda
    const double    tol = 0.1;//error tolerance for target probabilities
    const int       dim = NUMAA - 1;
//     const int       upperdof = 20;//upper limit for degrees of freedom
//     const int       quesize = fst.GetColumns();
//     const int       subsize = sec.GetColumns();
//     const double    queryenos = fst.GetEffNoSequences();
//     const double    sbjctenos = sec.GetEffNoSequences();
//     const double    querydenm = 19.0;//queryenos /( double )upperdof;
    double  tprob;//target probability from profile
    double  pscor;//score from profile
    double  lp;//log-probability
//     Dish*   dsh;
    Pslvector*  ps, *bps;
//     bool    added = false;//noise added
//     int propsize;//subject proposal size
    int bnpos;//, npos, kpos, dpos;
    int n, m, a;

    if( parentHDP == NULL )
        throw myruntime_error("HDP0ScoreMatrix: InitHDPMtx: Null parent HDP object.");
    if( parentHDP->GetMenu() == NULL )
        throw myruntime_error("HDP0ScoreMatrix: InitHDPMtx: Null parent HDP menu.");

    if( quesize < 1 || subsize < 1 )
        return;

    //one position is for background prob. vector
    basHDP.ReserveMenu( quesize + subsize + 1 );
    basHDP.ReserveBasin( quesize + subsize + 1 );
    basHDP.GetBasin()->SetDestroy( false );

    for( n = 0; n < quesize; n++ ) {
        ps = new Pslvector( NUMAA+0 );
        if( ps == NULL )
            throw myruntime_error("HDP0ScoreMatrix: InitHDPMtx: Not enough memory.");
        for( a = 0; a < NUMAA; a++ ) {
            pscor = fst( n, a );
            tprob = exp( pscor * lambda ) * LOSCORES.PROBABility( a );
            ps->SetValueAt( a, tprob );
        }
//         ps->SetValueAt( NUMAA, 1.e-3 );
        HDPbase::LogitNormal2Normal( ps, tol );
        //reduce dimensionality
        ps->DecDim();
        bnpos = fstbas.NewValue( ps );
    }

    for( m = 0; m <= subsize; m++ ) {
//         dsh = new Dish( 1 );
        ps = new Pslvector( NUMAA+0 );
        if( ps == NULL /*|| dsh == NULL */)
            throw myruntime_error("HDP0ScoreMatrix: InitHDPMtx: Not enough memory.");
//         dsh->SetBasin( &secbas );
        for( a = 0; a < NUMAA; a++ ) {
            tprob = LOSCORES.PROBABility( a );
            //the last position is for background probability vector
            if( m < subsize ) {
                pscor = sec( m, a );
                tprob *= exp( pscor * lambda );
            }
            ps->SetValueAt( a, tprob );
        }
//         ps->SetValueAt( NUMAA, 1.e-3 );
        HDPbase::LogitNormal2Normal( ps, tol );
        //reduce dimensionality
        ps->DecDim();
        bnpos = secbas.NewValue( ps );
//         npos = basHDP.AddToBasin( ps );
//         kpos = basHDP.GetMenu()->NewDish( dsh );
//         propsize = 1;
//         dpos = dsh->NewVectorNInd( bnpos );
//         dsh->SetReadSize( propsize );
    }

    //{{prepare HDPbase
    basHDP.GetMenu()->SoftCopyPriors( *parentHDP->GetMenu());
    basHDP.SetDegFAdjustment( HDPSM_Qdegfr );
    //reduce weight of prior parameters after calculating them
    basHDP.GetMenu()->SetKappa0( HDPSM_kappa0 );//IMP.!
    basHDP.GetMenu()->SetNu0( HDPSM_kappa0 );
    basHDP.CalcPriorProbFact();
//     //calculate parameters for each dish
//     basHDP.RecalcMenuParams();
//     //to leverage probability order
//     for( n = 0; n < basHDP.GetMenu()->GetSize(); n++ )
//         if( basHDP.GetMenu()->GetDishAt(n))
//             basHDP.GetMenu()->GetDishAt(n)->SetReadSize(0);
//     //}}

    //precalculate log-probabilities for query positions
    fstvec.Reserve( quesize );
    for( n = 0; n < quesize; n++ ) {
        ps = fstbas.GetValueAt( n );
        if( ps == NULL )
            throw myruntime_error("HDP0ScoreMatrix: InitHDPMtx: Memory access error.");
        basHDP.SetDegFAdjustment( HDPSM_Qdegfr );
        basHDP.PriorProbVec( ps, &lp );//P(Q)-- prior probability
//         basHDP.ProbVecOfDish( ps, subsize, &lp );//P(Q|B), where B is background vector
        fstvec.SetValueAt( n, lp );
    }

    //precalculate log-probabilities for subject positions
    secvec.Reserve( subsize );
    for( m = 0; m < subsize; m++ ) {
        ps = secbas.GetValueAt( m );
        if( ps == NULL )
            throw myruntime_error("HDP0ScoreMatrix: InitHDPMtx: Memory access error.");
        basHDP.SetDegFAdjustment( HDPSM_Qdegfr );
        basHDP.PriorProbVec( ps, &lp );//P(S)-- prior probability
        secvec.SetValueAt( m, lp );
    }
}





// =========================================================================
// InitHDP: initialize hdp-related structures
//
void HDP0ScoreMatrix::InitHDP()
{
    if( GetUseMtxVarT())
        InitHDPMtx();
    else
        InitHDPVec();
}





// =========================================================================
// ComputeScoreVec: compute score at a query-subject position;
//  overloaded for making alignments symmetric
//
double HDP0ScoreMatrix::ComputeScoreVec( int m, int n, bool final,
                                     const LogOddsMatrix& sec, const LogOddsMatrix& fst, 
                                     Pslvector& secvec, Pslvector& fstvec,
                                     Basin& secbas, Basin& fstbas, HDPbase& basHDP )
{
    const double    queryenos = fst.GetEffNoSequences();
    const double    sbjctenos = sec.GetEffNoSequences();
    const double    querydenm = 19.0;//queryenos /( double )upperdof;
    double  score = 0.0;
    double  scadj = 0.0;
    double  qlp, slp;//log probabilities of query, subject
    double  qslp, qblp;//log probabilities of query and subject, and query and background
    int dim = NUMAA - 1;
    int d;
    Pslvector*  qv, *sv;
    Pslmatrix   qsv( dim, 2 );

    qv = fstbas.GetValueAt( n );
    if( qv == NULL )
        throw myruntime_error("HDP0ScoreMatrix: ComputeScoreVec: Memory access error.");

//     basHDP.SetDegFAdjustment(
//         SLC_MIN( SLC_MAX( rint( fst.GetMIDExpNoObservationsAt( n, PS_M )/ querydenm ), 1.0), 21 ));
    basHDP.SetDegFAdjustment( HDPSM_Qdegfr );

    basHDP.ProbVecOfDish( qv, m, &qslp );//P(Q|S)
    qblp = qlp = fstvec.GetValueAt( n );//P(Q|B)

    if( final ) {
    }

    score = qslp - qblp;
    //score range depends on prior distribution; if subject is non-informative,
    //scores may be lifted several times (probabilities are much smaller)
    score = SLC_MIN( 10.0, score );

    return ( score - scadj ) * GetMultiplier();
}

// -------------------------------------------------------------------------
// ComputeScoreVec: compute score at a query-subject position
//
double HDP0ScoreMatrix::ComputeScoreVec( int m, int n, bool final )
{
    const int       qsize = GetQuerySize();
    const int       ssize = GetSubjectSize();
    const double    queryenos = GetQueryLogo().GetEffNoSequences();
    const double    sbjctenos = GetSbjctLogo().GetEffNoSequences();

    if( GetHDPbase() == NULL || GetQueryBasin() == NULL || GetSbjctBasin() == NULL )
        throw myruntime_error("Null HDP structures.");

    if( n < 0 || GetQueryBasin()->GetSize() <= n ||
        m < 0 || GetSbjctBasin()->GetSize() <= m )
        throw myruntime_error("HDP0ScoreMatrix: ComputeScoreVec: Invalid positions.");

    if( sbjctenos <= queryenos ) {
        if( GetHDPbase()->GetMenu()->GetSize() <= m )
            throw myruntime_error("HDP0ScoreMatrix: ComputeScoreVec: Invalid sbjct position.");
        if( GetQueValues().GetSize() <= n )
            throw myruntime_error("HDP0ScoreMatrix: ComputeScoreVec: Invalid query position.");
        return ComputeScoreVec( m, n, final,
                  GetSbjctLogo(), GetQueryLogo(),
                  GetSubValues(), GetQueValues(),
                  *GetSbjctBasin(), *GetQueryBasin(), *GetHDPbase());
    } else {
        if( GetHDPbase()->GetMenu()->GetSize() <= n )
            throw myruntime_error("HDP0ScoreMatrix: ComputeScoreVec: Invalid query position.");
        if( GetSubValues().GetSize() <= m )
            throw myruntime_error("HDP0ScoreMatrix: ComputeScoreVec: Invalid sbjct position.");
        return ComputeScoreVec( n, m, final,
                  GetQueryLogo(), GetSbjctLogo(),
                  GetQueValues(), GetSubValues(),
                  *GetQueryBasin(), *GetSbjctBasin(), *GetHDPbase());
    }
}





// =========================================================================
// ComputeScoreMtx: compute score at a query-subject position
//
double HDP0ScoreMatrix::ComputeScoreMtx( int m, int n, bool final )
{
    const LogOddsMatrix& sec = GetSbjctLogo();
    const LogOddsMatrix& fst = GetQueryLogo();
    Pslvector&  secvec = GetSubValues();
    Pslvector&  fstvec = GetQueValues();
    Basin&      secbas = *GetSbjctBasin();
    Basin&      fstbas = *GetQueryBasin();
    HDPbase&    basHDP = *GetHDPbase();

    const double    queryenos = fst.GetEffNoSequences();
    const double    sbjctenos = sec.GetEffNoSequences();
    double  score = 0.0;
    double  scadj = 0.0;
    double  qlp, slp;//log probabilities of query, subject
    double  qslp;//log probabilities of query and subject
    Pslvector*  qv, *sv;
    int d;

    sv = secbas.GetValueAt( m );
    qv = fstbas.GetValueAt( n );
    if( sv == NULL || qv == NULL )
        throw myruntime_error("HDP0ScoreMatrix: ComputeScoreMtx: Memory access error.");

    const int dim = sv->GetSize();
    Pslmatrix   qsv( dim, 2 );

    for( d = 0; d < dim; d++ )  qsv.SetValueAt( d, 0, sv->GetValueAt( d ));
    for( d = 0; d < dim; d++ )  qsv.SetValueAt( d, 1, qv->GetValueAt( d ));

    basHDP.SetDegFAdjustment( HDPSM_Qdegfr );

    basHDP.PriorProbMtx( &qsv, &qslp );//P(S,Q)
    slp = secvec.GetValueAt( m );//P(S)
    qlp = fstvec.GetValueAt( n );//P(Q)

    if( final ) {
    }

    score = qslp - qlp - slp;
    return ( score - scadj ) * GetMultiplier();
}





// =========================================================================
// ComputeScore: compute score at a query-subject position
//
double HDP0ScoreMatrix::ComputeScore( int m, int n, bool final )
{
    if( GetUseMtxVarT())
        return ComputeScoreMtx( m, n, final );
    else
        return ComputeScoreVec( m, n, final );
}





// =========================================================================
// ComputeProfileScoringMatrix: compute scoring matrix of a pair of profiles
//
void HDP0ScoreMatrix::ComputeProfileScoringMatrix( bool final )
{
    if( GetQuerySize() <= 0 || GetSubjectSize() <= 0 )
        throw myruntime_error( mystring( "HDP0ScoreMatrix: Unable to compute matrix of scores." ));

    if( !IsValid())
        throw myruntime_error( mystring( "HDP0ScoreMatrix: Not valid score matrix." ));

    int     m, n;
    int     no_elems = 0;
    int     negatives = 0;
    TMask   currentmask = Unmasked;
    double  score = 0.0;

    if( !final )
        InitHDP();

    //fill matrix with values
    for( m = 0; m < GetSubjectSize(); m++ ) {
        for( n = 0; n < GetQuerySize(); n++ )
        {
            //give constant penalties at the positions of X
//             if( GetQueryFreq().GetResidueAt( n ) == X ||
//                 GetSbjctFreq().GetResidueAt( m ) == X )
//             {
//                 SetScore( g_scoreX, m, n );
//                 continue;
//             }

            score = 0.0;
            currentmask = Unmasked;

//             if( logo_fst_.GetInformationAt( n ) < GetInformationThreshold() ||
//                 logo_sec_.GetInformationAt( m ) < GetInformationThreshold() ||
//                 ! ThicknessConstraintsMet( m, n ))
//             {
//                 //mask the score possibly excluding it from the ongoing computational statistics
//                 SetMasked( currentmask = GetMaskingApproach(), m, n );
//             }

            score = ComputeScore( m, n, final );

            if( currentmask != MaskToIgnore ) {
                no_elems++;
                if( score < 0.0 ) negatives++;
            }
            if( currentmask == MaskToConsider || currentmask == MaskToIgnore )
                score *= GetMaskscalePercents();

            SetScore( score, m, n );
        }
    }

    if( negatives == no_elems ) {
        SetAllNegatives( true );
        return;
    }

    if( GetAutoScaling())
            GetScaledScores()->CheckforExpectedScore();
    else    GetCorresScores()->CheckforExpectedScore();
}

// -------------------------------------------------------------------------
// ReduceScoresBy: reduce matrix scores by given value
//
bool HDP0ScoreMatrix::ReduceScores( AttributableScores* PATTR_SCORES )
{
return false;
    const double    adjconst = 0.1;
    double          expected = PATTR_SCORES->GetExpectedScore();
    double          ascale = PATTR_SCORES->GetAutoScalingFactor();
    double          delta;
    if( ascale )
        expected /= ascale;
    //when expected score is just below 0, lambda can be extremely small with
    //no potential to scale scores
    delta = expected + expected + adjconst;
    ReduceScoresBy( delta );
    return true;
}

// -------------------------------------------------------------------------
// ReduceScoresBy: reduce matrix scores by given value
//
void HDP0ScoreMatrix::ReduceScoresBy( double delta )
{
    int     m, n;
    int     no_elems = 0;
    int     negatives = 0;
    TMask   currentmask = Unmasked;
    double  score = 0.0;

    for( m = 0; m < GetSubjectSize(); m++ ) {
        for( n = 0; n < GetQuerySize(); n++ )
        {
            score = 0.0;
            currentmask = Unmasked;

            score = GetImageScore( m, n );
            if( SCORE_MIN < score )
                score -= delta * GetMultiplier();

            if( currentmask != MaskToIgnore ) {
                no_elems++;
                if( score < 0.0 ) negatives++;
            }
            SetScore( score, m, n );
        }
    }

    if( negatives == no_elems ) {
        SetAllNegatives( true );
        return;
    }
}

// -------------------------------------------------------------------------
// ComputeExpectation: compute e-value given previously computed
//     statistical parameters
//
double HDP0ScoreMatrix::ComputeExpectation(
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
void HDP0ScoreMatrix::ComputePositionalScoreProbs( AttributableScores* PATTR_SCORES )
{
    throw myruntime_error( "HDP0ScoreMatrix: ComputePositionalScoreProbs is not to be called." );
}

// -------------------------------------------------------------------------
// ComputeScoreProbabilities: compute probabilities to observe scores at
//     each position (i,j)
//
void HDP0ScoreMatrix::ComputeScoreProbabilities( AttributableScores* PATTR_SCORES )
{
    if( GetFPScaling())
        ComputeScoreProbabilitiesFPI( PATTR_SCORES );
    else
        ComputeScoreProbabilitiesII( PATTR_SCORES );
}

// -------------------------------------------------------------------------
// ComputeScoreProbabilitiesII: computes probabilities to observe scores at
//     each position (i,j)
//
void HDP0ScoreMatrix::ComputeScoreProbabilitiesII( AttributableScores* PATTR_SCORES )
{
    if( PATTR_SCORES == NULL )
        throw myruntime_error( mystring( "HDP0ScoreMatrix: Unable to compute probabilities: Wrong argument." ));

    const double    accuracy = 1.0e-4;

    int     n, m, nn, mm;
    char    strbuf[KBYTE];
    double  val = 0.0;      //initial sum
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


    if( GetUseHDPBProbs()) {
        normterm = SLC_MAX( GetQueValues().Max(), GetSubValues().Max());

        if( GetSubjectSize() != GetSubValues().GetSize())
            throw myruntime_error( 
                "HDP0ScoreMatrix: ComputeScoreProbabilities: Invalid sbjct size." );
        if( GetQuerySize() != GetQueValues().GetSize())
            throw myruntime_error( 
                "HDP0ScoreMatrix: ComputeScoreProbabilities: Invalid query size." );
    }
    for( n = 0, nn = 0; n < GetQuerySize(); n++ ) {
//         if( freq_fst_[n] == X )
//             continue;
        if( GetUseHDPBProbs()) {
            val = GetQueValues().GetValueAt( n );
            rprobs_.SetValueAt( nn++, exp( val - normterm ));
        } else 
            rprobs_.SetValueAt( nn++, 1.0 );
    }

    for( m = 0, mm = 0; m < GetSubjectSize(); m++ ) {
//         if( freq_sec_[m] == X )
//             continue;
        if( GetUseHDPBProbs()) {
            val = GetSubValues().GetValueAt( m );
            cprobs_.SetValueAt( mm++, exp( val - normterm ));
        } else 
            cprobs_.SetValueAt( mm++, 1.0 );
    }

    avgscore = normterm = 0.0;
    for( n = 0, nn = 0; n < GetQuerySize(); n++ )
    {
//         if( freq_fst_[n] == X )
//             continue;
        proquery = rprobs_.GetValueAt( nn++ );

        for( m = 0, mm = 0; m < GetSubjectSize(); m++ ) {
//             if( freq_sec_[m] == X )
//                 continue;
            prosbjct = cprobs_.GetValueAt( mm++ );

            loc_score = ( TScore )rint( PATTR_SCORES->GetScore( m, n ));

            exp2sum = proquery * prosbjct;

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
// ComputeScoreProbabilitiesFPI: computes probabilities to observe scores at
//     each position (i,j)
//
void HDP0ScoreMatrix::ComputeScoreProbabilitiesFPI( AttributableScores* PATTR_SCORES )
{
    if( PATTR_SCORES == NULL )
        throw myruntime_error( mystring( "HDP0ScoreMatrix: Unable to compute probabilities: Wrong argument." ));

    const double    accuracy = 1.0e-4;

    int     ind;
    int     n, m, nn, mm, r;
    char    strbuf[KBYTE];
    double  val = 0.0;      //initial sum
    double  proquery = 0.0; //non-normalized probability of query position
    double  prosbjct = 0.0; //non-normalized probability of subject position
    double  exp2sum = 0.0;  //e to sum
    double  normprob = 0.0; //normalizing term
    double  normterm = 0.0; //normalizing term
    double  avgscore = 0.0;
    double  consv;
    double  loc_score  = 0;

    scores_.Clear();
    rprobs_.Clear();
    cprobs_.Clear();


    if( GetUseHDPBProbs()) {
        normterm = SLC_MAX( GetQueValues().Max(), GetSubValues().Max());

        if( GetSubjectSize() != GetSubValues().GetSize())
            throw myruntime_error( 
                "HDP0ScoreMatrix: ComputeScoreProbabilities: Invalid sbjct size." );
        if( GetQuerySize() != GetQueValues().GetSize())
            throw myruntime_error( 
                "HDP0ScoreMatrix: ComputeScoreProbabilities: Invalid query size." );
    }
    for( n = 0, nn = 0; n < GetQuerySize(); n++ ) {
//         if( freq_fst_[n] == X )
//             continue;
        if( GetUseHDPBProbs()) {
            val = GetQueValues().GetValueAt( n );
            rprobs_.SetValueAt( nn++, exp( val - normterm ));
        } else 
            rprobs_.SetValueAt( nn++, 1.0 );
    }

    for( m = 0, mm = 0; m < GetSubjectSize(); m++ ) {
//         if( freq_sec_[m] == X )
//             continue;
        if( GetUseHDPBProbs()) {
            val = GetSubValues().GetValueAt( m );
            cprobs_.SetValueAt( mm++, exp( val - normterm ));
        } else 
            cprobs_.SetValueAt( mm++, 1.0 );
    }


    normterm = 0.0;
    for( n = 0, nn = 0; n < GetQuerySize(); n++ )
    {
//         if( freq_fst_[n] == X )
//             continue;

        proquery = rprobs_.GetValueAt( nn++ );

        for( m = 0, mm = 0; m < GetSubjectSize(); m++ ) {
//             if( freq_sec_[m] == X )
//                 continue;

            prosbjct = cprobs_.GetValueAt( mm++ );

            loc_score = PATTR_SCORES->GetScore( m, n );

            exp2sum = proquery * prosbjct;

            if( loc_score <= SCORE_MIN )
                    scores_.Push( g_scoreX );
            else    scores_.Push( loc_score );

            ind = n * GetSubjectSize() + m;
            if( loc_score <= SCORE_MIN ||
                PATTR_SCORES->GetMaskedToIgnore( m, n )) {
                PATTR_SCORES->IncProbabilityOf( 0.0, ( TScore )ind );
                continue;
            }
            PATTR_SCORES->IncProbabilityOf( exp2sum, ( TScore )ind );
            avgscore += exp2sum * loc_score;
            normterm += exp2sum;
        }
    }

    if( normterm ) {
        avgscore /= normterm;
        for( ind = PATTR_SCORES->GetMinScore(); ind <= PATTR_SCORES->GetMaxScore(); ind++ ) {
            PATTR_SCORES->DivideProbabilityOf( normterm, ( TScore )ind );
        }
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





// =========================================================================
// OptimizeTargetFrequencies: set target frequencies so that they become 
//  valid probabilities
//
void HDP0ScoreMatrix::OptimizeTargetFrequencies()
{
    double      lambda;
    mystring    errstr;
    int         m, n;
    double      score;

    if( GetAllNegatives())
        return;

    ComputeStatisticalParameters();

    if( 0.0 <= GetExpectedScore() || GetLambda() <= 0.0 ) {
        warning( "Unable to adjust target frequencies." );
        return;
    }

    lambda = GetLambda();
    if( GetAutoScaling())
        lambda = GetScaledScores()->GetLambda();

    for( m = 0; m < GetSubjectSize(); m++ ) {
        for( n = 0; n < GetQuerySize(); n++ )
        {
            score = GetImageScore( m, n );
            if( score <= SCORE_MIN )
                continue;
            if( score )
                SetScore( lambda * score, m, n );
        }
    }
}

// =========================================================================
// OptimizeTargetFrequencies: Optimize target frequencies so that marginal 
//  probabilities be compatible with null model probabilities
//
void HDP0ScoreMatrix::OptimizeTargetFrequenciesObs()
{
    const double    c_tol = 1.e-6;
    const int   c_maxiters = 200;
    double      lambda;
    mystring    errstr;

    if( GetAllNegatives()) {
        warning( "Unable to optimize target frequencies." );
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
//     optimizer.SetConstrainedH( GetRefH());
    optimizer.SetNMResidTol( c_tol );
    optimizer.SetMaxNoNMIterations( c_maxiters );

    status = optimizer.Optimize();

    if( status != PSL_SUCCESS /*&& status != PSL_MAXITERATS*/ ) {
        errstr = "Target frequencies not optimized: ";
        errstr += TranslatePSLError( status );
        warning( errstr.c_str());
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
// IncorporateTargetFrequencies: calculate new scores based on optimized
//     target frequencies
//
void HDP0ScoreMatrix::IncorporateTargetFrequencies( const Pslvector& optfreqs )
{
    if( rprobs_.GetSize() < 1 || cprobs_.GetSize() < 1 ||
        optfreqs.GetSize() != rprobs_.GetSize() * cprobs_.GetSize()) {
        warning( "Scores not updated: dimensionality error." );
        return;
    }

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
        warning( "New target frequencies not conserved." );
        return;
    }

    ind = 0;
    for( n = 0, nn = 0; n < GetQuerySize(); n++ )
    {
        if( freq_fst_[n] == X )
            continue;

        proquery = rprobs_.GetValueAt( nn++ );

        for( m = 0, mm = 0; m < GetSubjectSize(); m++ ) {
            if( freq_sec_[m] == X )
                continue;

            prosbjct = cprobs_.GetValueAt( mm++ );

            tfval = optfreqs.GetValueAt( ind++ );

            if( proquery <= 0.0 || prosbjct <= 0.0 )
                throw myruntime_error(
                    "IncorporateTargetFrequencies: Invalid null model probabilities." );

            score = log( tfval /( proquery * prosbjct ));

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
//     values to a stream
//

void HDP0ScoreMatrix::PrintParameterTable( TPrintFunction print_func, void* vpn ) const
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

void HDP0ScoreMatrix::PrintFinal( TPrintFunction print_func, void* vpn ) const
{
    if( vpn == NULL )
        return;

    PrintReferenceParameterTable( print_func, vpn );
}

// -------------------------------------------------------------------------
// PrintScoringMatrix: output the computed scoring system
//
void HDP0ScoreMatrix::PrintScoringMatrix( FILE* fp )
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

