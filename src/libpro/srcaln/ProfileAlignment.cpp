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
#include "ext/gammainc.h"
#include "ext/betainc.h"
#include "libpro/srcpro/SEGProfile.h"
#include "libpro/srcaln/AlignmentProbe.h"
#include "libpro/srcaln/MAAlignment.h"
#include "libpro/srcaln/ProfileAlignment.h"

//global file variables
int     gnMAXMARGIN = 60;   //maximum margin length at both ends of profile
double  gdMARGINPRC = 0.34; //percentage of profile length, set as margin length

double  ProfileAlignment::s_information_thrsh = 0.0;    //information content threshold
double  ProfileAlignment::s_segdistance_thrsh = 0.0;    //SEG distance threshold

FrequencyMatrix dummyfreq;
LogOddsMatrix   dummylogo;
GapScheme       dummygaps;

//values of pseudo hyperparameters for model SS18
double  gSS18hypers_1[ProfileAlignment::noSS18hypers] = {0.1, 12.,  1., 4.,  0.05, 1.,   0.1,  0.35,  0.65};
double  gSS18hypers_2[ProfileAlignment::noSS18hypers] = {0.6, 16.,  1.,12.,  0.3, -1.7,  0.1,  0.45,  0.65};

//comments:
//scaleE=.1*muE+.84;//linnear fit
//scaleR=.17*muR+.09;//linnear fit
//scaleR=.05*muR;//works similarly as .1*(exp(1.*scaleR)-1.) for mss18_1

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
    modelSS18_( 0 ),
    hypersSS18_( NULL ),
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

    SetModelSS18( mss18_1 );

    querySize_ = freq_fst_.GetColumns();
    subjectSize_ = freq_sec_.GetColumns();

    Initialize();
}

// -------------------------------------------------------------------------
// default constructor is invalid for initialization
//
ProfileAlignment::ProfileAlignment()
:   modelSS18_( 0 ),
    hypersSS18_( NULL ),

    model( StatModel::PC_DIST_HOM_0_2 ),

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
// SetModelSS18: set model for statistical significance estimation, SS18
//
void ProfileAlignment::SetModelSS18( int mod )
{
    if( noSS18models < mod )
        throw myruntime_error( "ProfileAlignment: SetModelSS18: Invalid SS18 model." );
    modelSS18_ = mod;
    switch( modelSS18_ ) {
        case mss18_1: hypersSS18_ = gSS18hypers_1; break;
        case mss18_2: hypersSS18_ = gSS18hypers_2; break;
    }
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

// *************************************************************************
// PARAMETERS FOR the ESTIMATION of STATISTICAL SIGNIFICANCE 
// -------------------------------------------------------------------------
// *** ********************************** ***
// COMBINED PREDICTION of PARAMETERS      ***
// *** ********************************** ***
// *** MODEL: genpro-fit2-rnd3-fr13-r0.05 ***
// *** ********************************** ***
// *** ONE LAYER: TRAINING UP TO n09      ***
const double WEIGHTS_MUSCALE[] = {//rows -- inputs (0th -- biases), columns -- units
 -0.2521174,  -0.09077898,  0.4872214, -0.2937255,
 -0.2945740,   0.10689298, -0.4693488,  1.1206345,
 -1.6767546, -51.75890070, -1.1090034, -0.7146433,
 -3.1166613,   6.21859310, -9.2319566, -0.8462229,
 47.3575084,   2.59408759, -1.0808783, -0.7300542,

  0.9882636, -0.9937430,
 -0.1616453,  0.6158327,
  0.1249525, -0.4921280,
  0.2313532, -3.1451028,
 -3.8387104,  1.3241476
};

// SEPARATE PREDICTION for each PARAMETER ***
enum TGEVDpar {
    nWGT_MU,
    nWGT_SCALE,
    nWGT
};

// *** *************************************************** ***
// *** RE MODEL: genpro-fit2-rnd3-fr9-r0.03-RE-altseed     ***
// *** *************************************************** ***
// *** ONE LAYER                                           ***
// *** *************************************************** ***
const double WEIGHTS_MU[] = {//rows -- inputs (0th -- biases), columns -- units
     3.1179029,     2.8802990,    -2.3533863,    -0.1763256,
     0.2732032,    -0.0845794,    -0.1646748,     0.3270687,
    -1.7865388,     0.7335855,     2.6470928,     0.5918245,
    -2.7370650,    -0.4155596,     1.6489798,     0.4900585,
    -2.4113938,    -3.5943765,    -0.7481717,     0.0253580,

     2.0547640,
     1.4051396,
    -1.7874640,
     2.0048461,
     0.6947985,
};
const double WEIGHTS_SCALE[] = {//rows -- inputs (0th -- biases), columns -- units
    -0.7879687,    -4.9900681,    -0.8238729,    -4.0875038,
    -1.2976974,    -0.1579375,    -0.4557760,     0.3549362,
    -1.5414541,     1.1291076,    -0.7614433,     1.2263836,
    -2.2530756,     4.1729244,     2.1426744,     6.7480938,
     3.6836031,     1.1692273,     2.8791774,     0.5062199,

     1.6914424,
    -1.1273527,
     3.0416709,
     1.6852105,
    -1.5032083,
};/**/
// *** K reference distribution                        ***
// *** ONE LAYER                                       ***
const double WEIGHTS_K[] = {//rows -- inputs (0th -- biases), columns -- units
    -1.7721228,     2.4816519,    -8.6518419,    -0.4259275,
     9.0964014,    -6.8232719,     9.5204826,    -6.0912179,
     0.1275992,    -0.0631605,     0.5367407,     0.0549147,
     4.7687464,   -16.5030831,     0.5409633,     5.6129098,
    -0.2514888,     0.3401646,     0.3956280,     0.5282228,

     5.7578105,
    -7.1751879,
    -3.9589345,
    -4.3828591,
    -2.1614059,
};
// *** regularized                                     ***
/*const double WEIGHTS_K[] = {//rows -- inputs (0th -- biases), columns -- units
    -2.8845544,    -1.5211136,    -0.1359915,     0.0383737,
     1.8756024,     2.0403454,     0.4074908,    -0.7979058,
     0.1621456,     0.0739473,     0.0175347,    -0.0838067,
     3.4892843,    -3.5095957,     0.1623601,    -0.2353806,
     0.2068843,     0.1709714,    -0.0249227,     0.0665686,

    -0.0146893,
    -4.2199248,
    -3.0887579,
    -0.4580971,
     0.8015521,
};*/
// *** LAMBDA distributed as gamma distribution        ***
// *** ONE LAYER                                       ***
// *** regularized                                     ***
const double WEIGHTS_SHAPE[] = {//rows -- inputs (0th -- biases), columns -- units
    -2.9018623,    -0.9684774,     1.2608146,     1.5975289,
    -1.9899775,     0.2388696,    -2.3885081,    -1.4792473,
     1.8689068,    -2.3583094,     0.2683977,    -2.4600280,
    -2.2387877,     2.0502527,    -4.5687282,    -1.0528608,
     1.6114651,     1.2506443,     1.3080572,     2.8845918,

     1.4792966,
     2.8305984,
    -2.2532871,
    -1.3726210,
     1.8922552,
};
// *** regularized                                     ***
const double WEIGHTS_RATE[] = {//rows -- inputs (0th -- biases), columns -- units
    -0.4075754,     3.2312478,     1.2221065,    -0.2575275,
     0.6070305,     3.0940547,    -1.5283114,     2.1143662,
    -2.8166513,    -1.5014809,    -2.5348326,     1.0289814,
     1.4609513,     3.6779055,    -0.3311446,     1.7141028,
     1.1267811,    -2.3566856,     3.1777027,    -1.2797282,

     1.2463668,
    -1.7424291,
    -2.6175948,
     1.7659680,
     1.7149519,
};
// *** Identities distributed as neg. binomial distr.  ***
// *** ONE LAYER                                       ***
// *** unregularized                                   ***
const double WEIGHTS_IDNS_SHAPE[] = {//rows -- inputs (0th -- biases), columns -- units
    -8.4636515,    16.9797299,    -9.1250755,     0.3564449,
    -4.2519708,     1.6559860,    -0.2791682,     1.4218709,
     7.7580154,   -43.6828921,    -5.3104691,     0.1543695,
     2.0993262,     0.2204461,    -0.1575909,     0.5074642,
     6.3921483,    10.1670416,    23.7761161,    -1.9629601,

     1.5204094,
    -0.2345974,
     0.4274742,
    -0.7237101,
    -0.2228414,
};
// *** unregularized                                   ***
const double WEIGHTS_IDNS_PROB[] = {//rows -- inputs (0th -- biases), columns -- units
     4.8073506,    -8.6539151,    -2.8968828,     8.0184503,
    -2.1758273,    -1.5768951,     1.2869805,    -4.4719919,
     4.0526030,     5.7376477,     7.4444679,    -5.3911772,
     1.6753047,     0.5279203,    -0.8111873,     2.4874509,
    -7.9253355,     6.1001819,    -6.2181541,    -3.9124716,

    -0.0137954,
     0.7790411,
     0.6318541,
    -0.8145518,
    -0.4963973,
};
// *** Positives distributed as neg. binomial distr.   ***
// *** ONE LAYER                                       ***
// *** unregularized                                   ***
const double WEIGHTS_PSTS_SHAPE[] = {//rows -- inputs (0th -- biases), columns -- units
    -2.2536918,     6.5498621,   -17.7679277,     8.9469180,
     0.2918674,     0.4063074,    -0.6413057,    -0.1579744,
     2.2383033,    -6.6635263,    39.9775323,     3.1968203,
    -2.9466323,     6.3040974,    -0.0988215,     0.7525636,
     1.7852864,    -4.2965767,    -6.6078174,   -21.7463742,

     1.3728666,
     0.3601257,
     0.1532177,
    -0.4205272,
     0.8366193,
};
// *** unregularized                                   ***
const double WEIGHTS_PSTS_PROB[] = {//rows -- inputs (0th -- biases), columns -- units
    10.7701016,     5.2423821,    -6.0872594,     2.2645375,
    -0.3471798,    -0.9742062,     2.1039062,    -4.3148497,
    -8.6466264,   -16.8524608,    -6.1815851,    -2.1717738,
     0.5593822,     0.1871712,    -1.1179758,     3.9230042,
    -6.3451519,     6.1641695,    18.0493775,     1.5164023,

     0.5530306,
    -0.7280781,
    -0.9750970,
     0.9860431,
     0.5354816,
};
// *** EVD of scores grouped by RED                    ***
// *** Exclusively-RED-based                           ***
// *** ONE LAYER                                       ***
const double WEIGHTS_R_MU[] = {//rows -- inputs (0th -- biases), columns -- units
    -3.9242536,     4.4161349,     1.7215862,     0.1724395,
     0.5252747,    -1.0654045,    -5.4591037,     4.3266981,
     4.6103076,     0.1210349,    -0.5459632,     1.5179339,
    -1.2452664,    -3.0690073,     1.6370711,    -4.8518849,

     2.1171617,
     0.4062795,
    -1.6235214,
     0.6316156,
     0.5419556
};
// *** regularized ***
const double WEIGHTS_R_SCALE[] = {//rows -- inputs (0th -- biases), columns -- units
    -0.0026744,    -1.8972487,    -0.2211457,    -1.3897139,
    -0.0002680,     3.3353436,    -1.5897148,     0.5816330,
     0.0004945,    -0.5830618,    -0.4223673,    -0.5388343,
     0.0044131,    -1.1410319,     5.1129024,     2.3233640,

     0.1723871,
    -0.0018545,
     1.7041458,
     2.3123308,
    -1.0130396
};

// *** EVD of scores grouped by Lmb                    ***
// *** Exclusively-Lmb-based                           ***
// *** ONE LAYER                                       ***
const double WEIGHTS_L_MU[] = {//rows -- inputs (0th -- biases), columns -- units
    -0.0605839,    -0.0173019,     0.0600411,    -5.9555713,
     0.0533123,    -0.0051053,    -0.0527767,    -3.1437676,
    -1.0046299,     0.6215832,     1.0036328,     1.8754493,
     0.9682557,    -0.5038284,    -0.9667071,     3.4776044,

     2.8086869,
    -0.5250554,
     0.3085814,
     0.5244239,
     2.4166220,
};
// *** regularized ***
const double WEIGHTS_L_SCALE[] = {//rows -- inputs (0th -- biases), columns -- units
    -4.4935289,     1.4710020,    -1.4710009,     4.8137710,
   -11.5416831,     0.7746900,    -0.7746948,     1.0000337,
     0.4601384,     0.5248285,    -0.5248293,     2.0837171,
     4.1707612,    -2.2621661,     2.2621667,    -5.3104535,

     2.5137293,
     4.2316344,
    -0.8620832,
     0.8620837,
     3.1901578,
};/**/

// *** EVD of scores grouped by Idns                   ***
// *** Exclusively-Idns-based                          ***
// *** ONE LAYER                                       ***
const double WEIGHTS_I_MU[] = {//rows -- inputs (0th -- biases), columns -- units
     0.0505477,    -0.9116769,    -5.3143876,     1.8500664,
    -0.9814986,     0.9831398,    -0.0997495,    -2.0215871,
    -2.3893879,     0.4754569,     1.7375035,     1.5275100,
     4.4739267,    -0.4448464,     2.4874111,    -1.8889199,

     2.7438446,
    -0.5897648,
    -0.5601991,
     2.1557387,
    -0.3400550,
};
// *** not regularized ***
const double WEIGHTS_I_SCALE[] = {//rows -- inputs (0th -- biases), columns -- units
    -3.1160686,     5.3274713,    -0.9449234,     3.0163482,
    -2.4690218,    -2.3375725,     1.3030284,    -0.0944874,
     2.1655387,    -1.4349413,     3.4371325,    -1.8755715,
     4.8969044,     3.8531494,    -7.4919642,    -2.3814889,

     1.9843227,
    -0.3898042,
    -1.6549375,
    -0.6424486,
    -0.8907515,
};

// =========================================================================

#define NN_UNIT_ACT_SIGM 0
#define NN_UNIT_ACT_TANH 1

// my_NN_evdpar: interface to general prediction function
int my_NN_evdpar( int inputs, int layers, int hunits, double* input, int inmax,
                  const double* cmbwgts, const double* par1wgts, const double* par2wgts,
                  int activn,
                  double* par1, double* par2 );


// my_NN_Kref: Predicting K
int my_NN_Kref( double E1, int L1, double E2, int L2, double* K )
{
    int tint;
    double tdbl;
    if( E1 < E2 ) {
        //change places of sequences: NN trained asymmetrically
        tdbl = E2; E2 = E1; E1 = tdbl; 
        tint = L2; L2 = L1; L1 = tint;
    }
    else if( E1 == E2 && L1 < L2 ) {
        tint = L2; L2 = L1; L1 = tint;
    }
    int       inputs = 4;
    const int layers = 1;//2;
    const int hunits = 4;//8;
    const int inmax = 1+4;//SLC_MAX( inputs, hunits );
    //double    input[inmax] = { 1.0, 0.05*E1, 0.001*((double)L1), 0.05*E2, 0.001*((double)L2) };
    double    input[inmax] = { 1.0, 0.05*E1, log((double)L1)/SLC_LN1K, 0.05*E2, log((double)L2)/SLC_LN1K };
    return my_NN_evdpar( inputs, layers, hunits, input, inmax, 
                         NULL, WEIGHTS_K, NULL, NN_UNIT_ACT_SIGM, K, NULL );
}

// my_NN_evdpar_E: Predicting EVD's mu and scale based on profile attributes, 
//  effnos and length 
int my_NN_evdpar_E( double E1, int L1, double E2, int L2, double* mu, double* scale )
{
    int tint;
    double tdbl;
    if( E1 < E2 ) {
        //change places of sequences: NN trained asymmetrically
        tdbl = E2; E2 = E1; E1 = tdbl; 
        tint = L2; L2 = L1; L1 = tint;
    }
    else if( E1 == E2 && L1 < L2 ) {
        tint = L2; L2 = L1; L1 = tint;
    }
    int       inputs = 4;
    const int layers = 1;//2;
    const int hunits = 4;//8;
    const int inmax = 1+4;//SLC_MAX( inputs, hunits );
    //double    input[inmax] = { 1.0, 0.05*E1, 0.001*((double)L1), 0.05*E2, 0.001*((double)L2) };
    double    input[inmax] = { 1.0, 0.05*E1, log((double)L1)/SLC_LN1K, 0.05*E2, log((double)L2)/SLC_LN1K };
    return my_NN_evdpar( inputs, layers, hunits, input, inmax, 
                         WEIGHTS_MUSCALE, WEIGHTS_MU, WEIGHTS_SCALE, NN_UNIT_ACT_TANH, mu, scale );
}

// my_NN_gammadpar_Lambda: Predicting the parameters of a Gamma distribution governing the 
//  distribution of Lambda 
int my_NN_gammadpar_Lambda( double E1, int L1, double E2, int L2, double* shape, double* rate )
{
    int tint;
    double tdbl;
    if( E1 < E2 ) {
        //change places of sequences: NN trained asymmetrically
        tdbl = E2; E2 = E1; E1 = tdbl; 
        tint = L2; L2 = L1; L1 = tint;
    }
    else if( E1 == E2 && L1 < L2 ) {
        tint = L2; L2 = L1; L1 = tint;
    }
    int       inputs = 4;
    const int layers = 1;//2;
    const int hunits = 4;//8;
    const int inmax = 1+4;//SLC_MAX( inputs, hunits );
    //double    input[inmax] = { 1.0, 0.05*E1, 0.001*((double)L1), 0.05*E2, 0.001*((double)L2) };
    double    input[inmax] = { 1.0, 0.05*E1, log((double)L1)/SLC_LN1K, 0.05*E2, log((double)L2)/SLC_LN1K };
    return my_NN_evdpar( inputs, layers, hunits, input, inmax, 
                         NULL, WEIGHTS_SHAPE, WEIGHTS_RATE, NN_UNIT_ACT_TANH, shape, rate );
}

// my_NN_nbdpar_Idns: Predicting the parameters of a negative binomial distribution governing the 
//  distribution of #Identities
int my_NN_nbdpar_Idns( double E1, int L1, double E2, int L2, double* shape, double* prob )
{
    int tint;
    double tdbl;
    if( E1 < E2 ) {
        //change places of sequences: NN trained asymmetrically
        tdbl = E2; E2 = E1; E1 = tdbl; 
        tint = L2; L2 = L1; L1 = tint;
    }
    else if( E1 == E2 && L1 < L2 ) {
        tint = L2; L2 = L1; L1 = tint;
    }
    int       inputs = 4;
    const int layers = 1;//2;
    const int hunits = 4;//8;
    const int inmax = 1+4;//SLC_MAX( inputs, hunits );
    //double    input[inmax] = { 1.0, 0.05*E1, 0.001*((double)L1), 0.05*E2, 0.001*((double)L2) };
    double    input[inmax] = { 1.0, 0.05*E1, log((double)L1)/SLC_LN1K, 0.05*E2, log((double)L2)/SLC_LN1K };
    return my_NN_evdpar( inputs, layers, hunits, input, inmax, 
                         NULL, WEIGHTS_IDNS_SHAPE, WEIGHTS_IDNS_PROB, NN_UNIT_ACT_TANH, shape, prob );
}

// my_NN_nbdpar_Psts: Predicting the parameters of a negative binomial distribution governing the 
//  distribution of #Positives
int my_NN_nbdpar_Psts( double E1, int L1, double E2, int L2, double* shape, double* prob )
{
    int tint;
    double tdbl;
    if( E1 < E2 ) {
        //change places of sequences: NN trained asymmetrically
        tdbl = E2; E2 = E1; E1 = tdbl; 
        tint = L2; L2 = L1; L1 = tint;
    }
    else if( E1 == E2 && L1 < L2 ) {
        tint = L2; L2 = L1; L1 = tint;
    }
    int       inputs = 4;
    const int layers = 1;//2;
    const int hunits = 4;//8;
    const int inmax = 1+4;//SLC_MAX( inputs, hunits );
    //double    input[inmax] = { 1.0, 0.05*E1, 0.001*((double)L1), 0.05*E2, 0.001*((double)L2) };
    double    input[inmax] = { 1.0, 0.05*E1, log((double)L1)/SLC_LN1K, 0.05*E2, log((double)L2)/SLC_LN1K };
    return my_NN_evdpar( inputs, layers, hunits, input, inmax, 
                         NULL, WEIGHTS_PSTS_SHAPE, WEIGHTS_PSTS_PROB, NN_UNIT_ACT_TANH, shape, prob );
}

// my_NN_evdpar_R: Predicting the EVD's parameters for scores grouped according to RED 
int my_NN_evdpar_R( double R, int L1, int L2, double* mu, double* scale )
{
    int tint;
    if( 10. < R )
        R = 10.;
    if( L1 < L2 ) {
        //NN trained asymmetrically
        tint = L2; L2 = L1; L1 = tint;
    }
    int       inputs = 3;
    const int layers = 1;//2;
    const int hunits = 4;//8;
    const int inmax = 1+4;//SLC_MAX( inputs, hunits );
    double    input[inmax] = { 1.0, 0.1*R, log((double)L1)/SLC_LN1K, log((double)L2)/SLC_LN1K };
    return my_NN_evdpar( inputs, layers, hunits, input, inmax, 
                         NULL, WEIGHTS_R_MU, WEIGHTS_R_SCALE, NN_UNIT_ACT_TANH, mu, scale );
}

// my_NN_evdpar_L: Predicting the EVD's parameters for scores grouped according to the 
//  value of Lambda
int my_NN_evdpar_L( double Lmb, int L1, int L2, double* mu, double* scale )
{
    int tint;
    Lmb = round(Lmb*10.)/10.;
    if( Lmb <= 0. )
        Lmb = 1.;
    if( 10. < Lmb )
        Lmb = 10.;
    if( L1 < L2 ) {
        //NN trained asymmetrically
        tint = L2; L2 = L1; L1 = tint;
    }
    int       inputs = 3;
    const int layers = 1;//2;
    const int hunits = 4;//8;
    const int inmax = 1+4;//SLC_MAX( inputs, hunits );
    double    input[inmax] = { 1.0, 0.2*Lmb, log((double)L1)/SLC_LN1K, log((double)L2)/SLC_LN1K };
    return my_NN_evdpar( inputs, layers, hunits, input, inmax, 
                         NULL, WEIGHTS_L_MU, WEIGHTS_L_SCALE, NN_UNIT_ACT_TANH, mu, scale );
}
// my_NN_evdpar_I: Predicting the EVD's parameters for scores grouped according to 
//  #Identities
int my_NN_evdpar_I( double idns, int L1, int L2, double* mu, double* scale )
{
    int tint;
    if( idns < 0. )
        idns = 0.;
    if( 50. < idns )
        idns = 50.;
    if( L1 < L2 ) {
        //NN trained asymmetrically
        tint = L2; L2 = L1; L1 = tint;
    }
    int       inputs = 3;
    const int layers = 1;//2;
    const int hunits = 4;//8;
    const int inmax = 1+4;//SLC_MAX( inputs, hunits );
    double    input[inmax] = { 1.0, 0.2*idns, log((double)L1)/SLC_LN1K, log((double)L2)/SLC_LN1K };
    return my_NN_evdpar( inputs, layers, hunits, input, inmax, 
                         NULL, WEIGHTS_I_MU, WEIGHTS_I_SCALE, NN_UNIT_ACT_TANH, mu, scale );
}

// my_NN_evdpar_RE: Predicting the EVD's parameters for scores grouped according to 
//  both profile attributes and RED
int my_NN_evdpar_RE( double R, double E1, int L1, double E2, int L2, double* mu, double* scale )
{
    int tint;
    double tdbl;
    if( E1 < E2 ) {
        //change places of sequences: NN trained asymmetrically
        tdbl = E2; E2 = E1; E1 = tdbl; 
        tint = L2; L2 = L1; L1 = tint;
    }
    else if( E1 == E2 && L1 < L2 ) {
        tint = L2; L2 = L1; L1 = tint;
    }
    int       inputs = 5;
    const int layers = 2;//1;//2;
    const int hunits = 8;//4;//8;
    const int inmax = 1+8;//5;//8;//SLC_MAX( inputs, hunits );
    double    input[inmax] = { 1.0, 0.1*R, 
    0.05*E1, log((double)L1)/SLC_LN1K, 0.05*E2, log((double)L2)/SLC_LN1K };
    return my_NN_evdpar( inputs, layers, hunits, input, inmax, 
                         WEIGHTS_MUSCALE, WEIGHTS_MU, WEIGHTS_SCALE, NN_UNIT_ACT_SIGM, mu, scale );
}

// .........................................................................
// my_NN_evdpar: implementation of the general prediction function
int my_NN_evdpar( int inputs, int layers, int hunits, double* input, int inmax,
                  const double* cmbwgts, const double* par1wgts, const double* par2wgts,
                  int activn,
                  double* par1, double* par2 )
{
    if( par1 == NULL && par2 == NULL ) {
        error("my_NN_evdpar: Null artguments.");
        throw myruntime_error("my_NN_evdpar: Null artguments.");
    }
    if( input == NULL || inmax != 1+SLC_MAX(inputs,hunits)) {
        error("my_NN_evdpar: Invalid inputs.");
        throw myruntime_error("my_NN_evdpar: Invalid inputs.");
    }

    const bool linout = false;//linear output
    const bool comb = par1 && par2;//combined output of mu and scale
    const double* weights = comb? cmbwgts: ( par1? par1wgts: par2wgts );
    const int outputs = comb? 2: 1;//number of outputs
    double midinp[hunits];
    int inputs1 = inputs+1;
    int lpsd, ulim;
    int l, u, i, j;

    if( comb && hunits < 2 ) {
        error("my_NN_evdpar: Invalid NN architecture.");
        throw myruntime_error("my_NN_evdpar: Invalid NN architecture.");
    }
    if( !weights ) {
        error("my_NN_evdpar: Null weights.");
        throw myruntime_error("my_NN_evdpar: Null weights.");
    }

    ulim = hunits;
    for( l = 0, lpsd = 0; l < layers+1; l++ ) {
        if( layers <= l ) ulim = outputs;
        memset( midinp, 0, hunits*sizeof(double));
        for( u = 0; u < ulim; u++ )
            for( i = 0; i < inputs1; i++ )
                midinp[u] += input[i] * weights[lpsd+u+ulim*i];
        if( !linout || l < layers )
            for( u = 0; u < ulim; u++ ) {
                input[u+1] = 1./(1.+exp(-midinp[u]));
                if( activn == NN_UNIT_ACT_SIGM );
                else if( activn == NN_UNIT_ACT_TANH )
                    input[u+1] += input[u+1] - 1.; //tanh(midinp[u]*.5); .5, steepness
                else {
                    error("my_NN_evdpar: Unknown activation function.");
                    throw myruntime_error("my_NN_evdpar: Unknown activation function.");
                }
            }
        lpsd += hunits * inputs1;
        inputs = ulim;
        inputs1 = inputs+1;
    }
    if( linout )
        for( u = 0; u < ulim; u++ )
            input[u+1] = midinp[u];

    if( par1 )
        *par1 = input[1];
    if( par2 )
        *par2 = input[outputs];
    return 0;
}

// =========================================================================

// -------------------------------------------------------------------------
// ComputeStatisticsHelper: helper method to compute statistical parameters
//
void ProfileAlignment::ComputeStatisticsHelper( 
    double qnos, int qlen, double snos, int slen,
    double red, double lmbd, int idns,
    double* mu, double* scale )
{
    double  muE = 0.0, muR = 0.0;
    double  scaleE = 0.0, scaleR = 0.0;
    double  alpha = 0.0, ibeta = 0.0;
    char    strbuf[KBYTE];
    const double fctmu = 30.;//output scale for mu
    const double minmu = 3.;//minimum value for mu
    const double fctscale = 5.;//output scale for EVD's scale
    const double maxscale = 5.;//maximum value for scale
    const double minscale = .1;//minimum value for scale

    if( hypersSS18_ == NULL )
        throw myruntime_error("ProfileAlignment: ComputeStatisticsHelper: Model params: Memory access error.");
    if( mu == NULL || scale == NULL )
        throw myruntime_error("ProfileAlignment: ComputeStatisticsHelper: Memory access error.");

    my_NN_evdpar_E( qnos, qlen, snos, slen, &muE, NULL );
    my_NN_evdpar_L( lmbd, qlen, slen, &muR, NULL );
    muE *= fctmu; muR *= fctmu;//NN predictions are scaled
    if( muE < minmu ) muE = minmu;
    if( muR < minmu ) muR = minmu;

    scaleE = hypersSS18_[ss18scaleE_S] * muE + hypersSS18_[ss18scaleE_I];
    if( scaleE < minscale ) scaleE = minscale;
    if( maxscale < scaleE ) scaleE = maxscale;

    my_NN_evdpar_L( lmbd, qlen, slen, NULL, &scaleR );
    scaleR *= fctscale;//NN predictions are scaled
    if( scaleR < minscale ) scaleR = minscale;
    scaleR = hypersSS18_[ss18scaleR_F] * ( exp(scaleR) - 1. );
    if( scaleR < minscale ) scaleR = minscale;
    if( maxscale < scaleR ) scaleR = maxscale;

    muE = hypersSS18_[ss18muE_S] * muE + hypersSS18_[ss18muE_I]; 
    muR = hypersSS18_[ss18muR_S] * muR + hypersSS18_[ss18muR_I];

//     if( muE < 3.0 ) muE = 3.0;
//     if( muR < 3.0 ) muR = 3.0;

    alpha = hypersSS18_[ss18alpha];
    //alpha += ( 1. - alpha ) * 0.1 * red;
    alpha -= alpha * 0.1 * SLC_MIN(10.,1./lmbd);
    //conditional mean estimate of mu
    *mu = alpha * muE + ( 1. - alpha ) * muR;

    ibeta = hypersSS18_[ss18ibeta];
    //ibeta += ( 1. - ibeta ) * 0.1 * red;
    ibeta += ( 1. - ibeta ) * 0.1 * SLC_MIN(10.,1./lmbd);
    //conditional mean estimate of scale
    *scale = ( 1. - ibeta ) * scaleE + ibeta * scaleR;

    if( *mu <= 0.0 ) {
        sprintf(strbuf, "ProfileAlignment: ComputeStatisticsHelper: Invalid location: %g.", *mu );
        throw myruntime_error( strbuf );
    }
    if( *scale <= 0.0 ) {
        sprintf(strbuf, "ProfileAlignment: ComputeStatisticsHelper: Invalid scale: %g.", *scale );
        throw myruntime_error( strbuf );
    }
}

// -------------------------------------------------------------------------
// ComputeKrefStatisticsHelper: helper method to compute parameter Kref
//
void ProfileAlignment::ComputeKrefStatisticsHelper( 
    double qnos, int qlen, double snos, int slen,
    double red, double K,
    double* Kref )
{
    char   strbuf[KBYTE];

    if( Kref == NULL )
        throw myruntime_error("ProfileAlignment: ComputeKrefStatisticsHelper: Memory access error.");

    my_NN_Kref( qnos, qlen, snos, slen, Kref );

    if( *Kref < 0.0  || 1.0 < *Kref ) {
        sprintf(strbuf, "ProfileAlignment: ComputeKrefStatisticsHelper: Invalid Kref: %g.", *Kref );
        throw myruntime_error( strbuf );
    }
}

// -------------------------------------------------------------------------
// ComputeLambdaStatisticsHelper: helper method to compute statistical 
//  parameters of gamma distribution representing the distribution of lambda
//
void ProfileAlignment::ComputeLambdaStatisticsHelper( 
    double qnos, int qlen, double snos, int slen, 
    double red, double lmbd, 
    double* shape, double* rate )
{
    double shapeE = 0.0;
    double rateE = 0.0;
    char   strbuf[KBYTE];

    if( shape == NULL || rate == NULL )
        throw myruntime_error("ProfileAlignment: ComputeLambdaStatisticsHelper: Memory access error.");

    my_NN_gammadpar_Lambda( qnos, qlen, snos, slen, &shapeE, NULL );
    if( 9. < shapeE ) {
        sprintf(strbuf, "ProfileAlignment: ComputeLambdaStatisticsHelper: "
                        "Invalid shape prediction: %g.", shapeE );
        throw myruntime_error( strbuf );
    }
    shapeE = exp(9.*shapeE);
    shapeE*= 0.05;
    *shape = shapeE;

    my_NN_gammadpar_Lambda( qnos, qlen, snos, slen, NULL, &rateE );
    if( 9. < rateE ) {
        sprintf(strbuf, "ProfileAlignment: ComputeLambdaStatisticsHelper: "
                        "Invalid rate prediction: %g.", rateE );
        throw myruntime_error( strbuf );
    }
    rateE = exp(9.*rateE);
    rateE*= 0.05;
    *rate = rateE;

    if( *shape <= 0.0 ) {
        sprintf(strbuf, "ProfileAlignment: ComputeLambdaStatisticsHelper: "
                        "Invalid shape parameter: %g.", *shape );
        throw myruntime_error( strbuf );
    }
    if( *rate <= 0.0 ) {
        sprintf(strbuf, "ProfileAlignment: ComputeLambdaStatisticsHelper: "
                        "Invalid rate parameter: %g.", *rate );
        throw myruntime_error( strbuf );
    }
}

// -------------------------------------------------------------------------
// ComputeIdnsStatisticsHelper: helper method to compute statistical 
//  parameters of neg. binomial distribution that represents the 
//  distribution of #identities
//
void ProfileAlignment::ComputeIdnsStatisticsHelper( 
    double qnos, int qlen, double snos, int slen, 
    double red, double lmbd, 
    double* shape, double* prob )
{
    double shapeE = 0.0;
    double probE = 0.0;
    char   strbuf[KBYTE];

    if( shape == NULL || prob == NULL )
        throw myruntime_error("ProfileAlignment: ComputeIdnsStatisticsHelper: Memory access error.");

    my_NN_nbdpar_Idns( qnos, qlen, snos, slen, &shapeE, NULL );
    if( 1. < shapeE || shapeE < 0. ) {
        sprintf(strbuf, "ProfileAlignment: ComputeIdnsStatisticsHelper: "
                        "Invalid shape prediction: %g.", shapeE );
        throw myruntime_error( strbuf );
    }
    shapeE = 15.*shapeE;
    //shapeE *= 10.;
    *shape = shapeE;

    my_NN_nbdpar_Idns( qnos, qlen, snos, slen, NULL, &probE );
    if( probE < .01 )
        probE = .01;
    if( .9 < probE )
        probE = .9;
    if( 1. < probE || probE < 0. ) {
        sprintf(strbuf, "ProfileAlignment: ComputeIdnsStatisticsHelper: "
                        "Invalid prob prediction: %g.", probE );
        throw myruntime_error( strbuf );
    }
    *prob = probE;

    if( *shape <= 0.0 ) {
        sprintf(strbuf, "ProfileAlignment: ComputeIdnsStatisticsHelper: "
                        "Invalid shape parameter: %g.", *shape );
        throw myruntime_error( strbuf );
    }
    if( 1. < *prob || *prob <= 0.0 ) {
        sprintf(strbuf, "ProfileAlignment: ComputeIdnsStatisticsHelper: "
                        "Invalid prob parameter: %g.", *prob );
        throw myruntime_error( strbuf );
    }
}

// -------------------------------------------------------------------------
// ComputeIdnsStatisticsHelper: helper method to compute statistical 
//  parameters of neg. binomial distribution that represents the 
//  distribution of #positives
//
void ProfileAlignment::ComputePstsStatisticsHelper( 
    double qnos, int qlen, double snos, int slen, 
    double red, double lmbd, 
    double* shape, double* prob )
{
    double shapeE = 0.0;
    double probE = 0.0;
    char   strbuf[KBYTE];

    if( shape == NULL || prob == NULL )
        throw myruntime_error("ProfileAlignment: ComputePstsStatisticsHelper: Memory access error.");

    my_NN_nbdpar_Psts( qnos, qlen, snos, slen, &shapeE, NULL );
    if( 1. < shapeE || shapeE < 0. ) {
        sprintf(strbuf, "ProfileAlignment: ComputePstsStatisticsHelper: "
                        "Invalid shape prediction: %g.", shapeE );
        throw myruntime_error( strbuf );
    }
    shapeE = 15.*shapeE;
    //shapeE *= 10.;
    *shape = shapeE;

    my_NN_nbdpar_Psts( qnos, qlen, snos, slen, NULL, &probE );
    if( probE < .01 )
        probE = .01;
    if( .999 < probE )
        probE = .999;
    if( 1. < probE || probE < 0. ) {
        sprintf(strbuf, "ProfileAlignment: ComputePstsStatisticsHelper: "
                        "Invalid prob prediction: %g.", probE );
        throw myruntime_error( strbuf );
    }
    *prob = probE;

    if( *shape <= 0.0 ) {
        sprintf(strbuf, "ProfileAlignment: ComputePstsStatisticsHelper: "
                        "Invalid shape parameter: %g.", *shape );
        throw myruntime_error( strbuf );
    }
    if( 1. < *prob || *prob <= 0.0 ) {
        sprintf(strbuf, "ProfileAlignment: ComputePstsStatisticsHelper: "
                        "Invalid prob parameter: %g.", *prob );
        throw myruntime_error( strbuf );
    }
}

// -------------------------------------------------------------------------
// ComputeStatisticsHelper2: helper method to compute e-value
//
void ProfileAlignment::ComputeExpectHelper2( 
    double qnos, int qlen, double snos, int slen, 
    double red, 
    double score, double* expect )
{
    mystring preamb = "ProfileAlignment: ComputeExpectHelper2: ";
    const AbstractScoreMatrix* scmatrix = GetScoreMatrix();
    if( scmatrix == NULL )
        throw myruntime_error( preamb + "Null score matrix.");

    const bool   bDEPENDENTKREF = false;//calculate Kref in the profile context
    //
    static const double mn_act = (double)scmatrix->GetSearchSpace();
    static const double inopairs = 0.0001;//inverse to the number of pairs used in experiments

    //norm: psts/pow(l1*l2,1.5)
    const double psts_factor = 100000.;//factor of psts in the numerator
    const double psts_exp = 1.5;//exponent in the denominator when normalizing psts 
    const double psts_nbd_shape_exp = 0.1;
    const double psts_nbd_prob_exp = 0.757213;
    //Fisher's generalized method
    const double dfo2 = 1.346206;//Fisher's test: degrees of freedom * 0.5
    const double c = 1.485656;//Fisher's test: scale factor for the test

    int    idns = GetNoIdentities();
    int    psts = GetNoPositives();
    double lmbd = scmatrix->GetLambda();//calculated ungapped
    double lmbdref = scmatrix->GetRefLambda();//ungapped reference
    double K = scmatrix->GetK();//calculated ungapped
    double Kref = scmatrix->GetRefK();//ungapped reference

    double mu_exp = 0.0;//empirical mu
    double scale_exp = 0.0;//empirical scale
    double scale_der;//scale derived
    double logK_der, logK_exp;//K derived and empirical
    double logKmn_der, logexpect;//K derived * mn and log expect
    double gamma_shape_exp = 0.0;//empirical shape of lambda's gammad
    double gamma_rate_exp = 0.0;//empirical rate of lambda's gammad
    double nbd_shape_exp = 0.0;//empirical shape of idnts/psts' nbd
    double nbd_prob_exp = 0.0;//empirical prob of idnts/psts' nbd
    double fct = 0.;
    double err;
    int    code;
    char   strbuf[KBYTE];

    if( expect == NULL )
        throw myruntime_error( preamb + "Null arguments.");

    ComputeStatisticsHelper( qnos, qlen, snos, slen, red, lmbd, idns, &mu_exp, &scale_exp );

    if( mu_exp <= 0.0 || scale_exp <= 0.0 ) {
        sprintf(strbuf, "%sInvalid evd parameters: %g %g.", preamb.c_str(), mu_exp, scale_exp );
        throw myruntime_error( strbuf );
    }

    if( bDEPENDENTKREF ) {
        //Kref in the profile context: slightly better but in the limit of error
        ComputeKrefStatisticsHelper( qnos, qlen, snos, slen, red, K, &Kref );
        if( Kref <= 0.0 || 1.0 <= Kref ) {
            sprintf(strbuf, "%sIgnoring reference-K parameter: %g.", preamb.c_str(), Kref );
            warning( strbuf );
        }
    }

    scale_der = scale_exp;
    if( 0&& 0. < lmbd && 0. < lmbdref ) {
        //unused...
        fct = lmbdref/lmbd;
        fct = SLC_MAX(1.,SLC_MIN(1.1,fct));
        scale_der = fct*scale_exp;
        //scale_der = lmbdref*scale_exp/lmbd;//s_d=s*s_e/s_u; s/s_u=l_u/l
    }

    logKmn_der = mu_exp / scale_exp;
    if( K > Kref && 0. < K && K < 1.0 && 0. < Kref && Kref < 1.0 )
        logKmn_der += log(K/Kref);//K_d=K*K_e/K_u

    logexpect = logKmn_der - score / scale_der;

    if( SLC_LOG_DBL_MAX < logexpect )
        *expect = SLC_DBL_MAX;
    else if( logexpect < SLC_LOG_DBL_MIN )
        *expect = 0.;
    else
        *expect = exp(logexpect);

    if( *expect <= 0.0 )
        *expect = SLC_DBL_MIN;

    //{{combining dependent significance
    int mod = GetModelSS18();//model SS18
    double pairexp = *expect;//expect calculated based on profile attributes and compositional similarity
    double pairpaval;//alignment p-value
    double pairpcval,pairpcvali;//dependent p-value //NOTE:
    double pairproduct;//p-value product
    double pairpuval;//unified p-value
    double pairexpunif;//unified pair expect
    double expunif;//unified expect

    if( mod == mss18_2 ) {
        //this follows from the assumption in prediction that 
        //pairexp=expect*qlen*slen/(qlen*slen*#pairs_in_experiments)=expect/#pairs_in_experiments;
        //note, however, an inverse dependence of mu on the prediction in this case
        pairexp *= inopairs;
    }

    if( pairexp < 0.01 )
        pairpaval = pairexp;
    else
        pairpaval = 1. - exp(-pairexp);

    pairpcval = 1.;

//     //{{LAMBDA
//     ComputeLambdaStatisticsHelper( qnos, qlen, snos, slen, red, lmbd, &gamma_shape_exp, &gamma_rate_exp );
//     if( gamma_shape_exp <= 0.0 || gamma_rate_exp <= 0.0 ) {
//         sprintf(strbuf, "%sInvalid gammad parameters: %g %g.", 
//                 preamb.c_str(), gamma_shape_exp, gamma_rate_exp );
//         throw myruntime_error( strbuf );
//     }
//     if( lmbd <= 0. )
//         lmbd = 0.1;
//     if(( code = psl_gammainc_P_e( gamma_shape_exp, gamma_rate_exp*(lmbd), &pairpcval, &err )) != PSL_OK ) {
//         throw myruntime_error( preamb + "lower inc. gamma failed: " + TranslatePSLError( code ));
//     }
//     //}}

    //{{PSTS normalized, general case
    psts = (int)rint( (double)psts*psts_factor / pow((double)(qlen*slen),psts_exp) );
    if(( code = psl_betainc_e( psts+1, psts_nbd_shape_exp, psts_nbd_prob_exp, &pairpcval, &err )) != PSL_OK ) {
          throw myruntime_error( preamb + "regularized inc. beta failed: " + TranslatePSLError( code ));
    }
    //}}
    //{{Fisher's generalized method for arbitrary k
    //s,psts
    if(( code = psl_gammainc_Q_e( dfo2, (-log(pairpaval)-log(pairpcval))/c, &pairpuval, &err )) != PSL_OK ) {
        throw myruntime_error( preamb + "upper inc. gamma failed: " + TranslatePSLError( code ));
    }
    //}}

    if( pairpuval < 0.01 )
        pairexpunif = pairpuval;
    else
        pairexpunif = -log(1.-pairpuval);

    expunif = pairexpunif * mn_act / ( qlen*slen );
    if( expunif <= 0.0 )
        expunif = SLC_DBL_MIN;
    *expect = expunif;
}

// -------------------------------------------------------------------------
// ComputeStatisticsHelper: helper method to compute e-value
//
void ProfileAlignment::ComputeExpectHelper( 
    double qnos, int qlen, double snos, int slen, 
    double red, double lmbd, int idns, 
    double score, double* expect )
{
    double mu = 0.0;
    double scale = 0.0;
    double explambda, expKmn;//lambda and Kmn
    char   strbuf[KBYTE];

    if( expect == NULL )
        throw myruntime_error("ProfileAlignment: ComputeExpectHelper: Null arguments.");

    ComputeStatisticsHelper( qnos, qlen, snos, slen, red, lmbd, idns, &mu, &scale );

    if( mu <= 0.0 || scale <= 0.0 ) {
        sprintf(strbuf, "ProfileAlignment: ComputeExpectHelper: "
                        "Invalid evd parameters: %g %g.", mu, scale );
        throw myruntime_error( strbuf );
    }

    *expect = exp((mu-score)/scale);
}

// -------------------------------------------------------------------------
// ComputeMixtureExpect: calculate e-value of the mixture
//
void ProfileAlignment::ComputeMixtureExpect( 
    double qnos, int qlen, double snos, int slen, 
    double red, double lmbd, 
    double score, double* expect )
{
    const DBProfileProbs* proprobs = GetProProbs();
    const bool probsready = proprobs && proprobs->ProbsReady();
    int    idns = GetNoIdentities();
    double prob, sumprob;
    double exprt;//partial expect
    int e1, e2, l1, l2;//level indices
    int n;

    if( expect == NULL )
        throw myruntime_error("ProfileAlignment: ComputeMixtureExpect: Null arguments.");

    *expect = sumprob = 0.;
    n = 0;

    if( !probsready ) {
        ComputeExpectHelper( qnos, qlen, snos, slen, red, lmbd, idns, score, expect );
    }
    else {
        proprobs->GetLevIndices( qnos, qlen, &e1, &l1 );
        //qnos = proprobs->GetEffLvs().GetValueAt(e1);
        //qlen = proprobs->GetLenLvs().GetValueAt(l1);
        for( e2 = 0; e2 <= e1; e2++ ) 
        {
            snos = proprobs->GetEffLvs().GetValueAt(e2);
            for( l2 = 0; l2 < proprobs->GetLenLvs().GetSize(); l2++ ) 
            {
                if( e1 <= e2 && l1 < l2 )
                  break;
                slen = proprobs->GetLenLvs().GetValueAt(l2);

                prob = proprobs->GetProbAt( e1, l1, e2, l2 ); n++;
                if( prob <= 0. )
                    continue;
                ComputeExpectHelper( qnos, qlen, snos, slen, red, lmbd, idns, score, &exprt );
                *expect += prob * exprt;
                sumprob += prob;
            }
        }
        snos = qnos;
        slen = qlen;
        l1++;
        if( proprobs->GetLenLvs().GetSize() <= l1 ) {
            e1++; l1 = 0;
        }
        for( --e2, --l2; e1 < proprobs->GetEffLvs().GetSize(); e1++, l1 = 0 ) 
        {
            qnos = proprobs->GetEffLvs().GetValueAt(e1);
            for( ; l1 < proprobs->GetLenLvs().GetSize(); l1++ ) 
            {
                qlen = proprobs->GetLenLvs().GetValueAt(l1);

                prob = proprobs->GetProbAt( e1, l1, e2, l2 ); n++;
                if( prob <= 0. )
                    continue;
                ComputeExpectHelper( qnos, qlen, snos, slen, red, lmbd, idns, score, &exprt );
                *expect += prob * exprt;
                sumprob += prob;
            }
        }

//       for( e1 = 0; e1 < proprobs->GetEffLvs().GetSize(); e1++ ) 
//       {
//         qnos = proprobs->GetEffLvs().GetValueAt(e1);
//         for( l1 = 0; l1 < proprobs->GetLenLvs().GetSize(); l1++ ) 
//         {
//           qlen = proprobs->GetLenLvs().GetValueAt(l1);
//           for( e2 = 0; e2 <= e1; e2++ ) 
//           {
//             snos = proprobs->GetEffLvs().GetValueAt(e2);
//             for( l2 = 0; l2 < proprobs->GetLenLvs().GetSize(); l2++ ) 
//             {
//               if( e1 <= e2 && l1 < l2 )
//                 break;
//               slen = proprobs->GetLenLvs().GetValueAt(l2);
// 
//               prob = proprobs->GetProbAt(n++);
//               if( prob <= 0. )
//                   continue;
//               ComputeExpectHelper( qnos, qlen, snos, slen, red, lmbd, idns, score, &exprt );
//               *expect += prob * exprt;
//             }
//           }
//         }
//       }

    }

    if( probsready && n != proprobs->GetCardinality())
        throw myruntime_error("ProfileAlignment: ComputeStatisticsHelper: Some probabilities unprocessed.");

    if( 0. < *expect && sumprob )
        *expect /= sumprob;
}

// -------------------------------------------------------------------------
// ComputeStatisticsSS18: calculate e-value and related measures using 
//  model SS18
//
void ProfileAlignment::ComputeStatisticsSS18()
{
    const AbstractScoreMatrix* scmatrix = GetScoreMatrix();
    if( scmatrix == NULL )
        throw myruntime_error("ProfileAlignment: ComputeStatisticsHelper: Null score matrix.");

    const int effnoress = NUMAA;
    double lmbd = scmatrix->GetLambda();
    double lmbdref = scmatrix->GetRefLambda();
    double K = scmatrix->GetK();
    double Kref = scmatrix->GetRefK();
    int qlen = GetQuerySize();
    int slen = GetSubjectSize();
    int pspace = qlen * slen;
    double qnos = logo_fst_.GetEffNoSequences();
    double snos = logo_sec_.GetEffNoSequences();
    double minos = SLC_MIN(qnos,snos);
    double score = GetScore();
    double expect, pair_expect, bitscore;//e-value and bit score
    double rescale = 10.;
    double pwscale = 80.;
    double re = 0.0;
    double red, val;
    double p1, p2;
    const double (*qtp)[NUMAA], (*stp)[NUMAA];
    double qap[NUMAA], sap[NUMAA];
    char   strbuf[KBYTE];
    int qin, sin;
    int n, c, r;

    red = 0.0;
    if( 0 ) {
        //uncomment when in use
        for( r = 0; r < effnoress; r++ ) {
            p1 = logo_fst_.GetPostProbsAt(r);
            p2 = logo_sec_.GetPostProbsAt(r);
            if( p1 <= 0. || p2 <= 0. )
                throw myruntime_error("ProfileAlignment: ComputeStatisticsHelper: "
                                      "Invalid profile posterior probabilities.");
            red += p1*log(p1/p2);
            red += p2*log(p2/p1);
        }
        red = rescale * exp(-pwscale*0.5*red);
    }
    if( 0 && 0 < GetAlnSteps()) {
        //not in use
        for( qnos = 0.0, c = 0, n = path[GetAlnSteps()-1][nQuery]-1; n < path[0][nQuery]; n++, c++ )
            qnos += logo_fst_.GetMIDExpNoObservationsAt( n, PS_M );
        qnos /= double(c);
        for( snos = 0.0, c = 0, n = path[GetAlnSteps()-1][nSbjct]-1; n < path[0][nSbjct]; n++, c++ )
            snos += logo_sec_.GetMIDExpNoObservationsAt( n, PS_M );
        snos /= double(c);
        minos = SLC_MIN(qnos,snos);
    }

    //{{Helper methods:
    //ComputeMixtureExpect( qnos, qlen, snos, slen, red, lmbd, score, &expect );
    ComputeExpectHelper2( qnos, qlen, snos, slen, red, score, &expect );
    //}}

    if( expect < 0. ) {
        sprintf(strbuf, "ProfileAlignment: ComputeStatisticsHelper: Invalid expect: %g.", expect );
        throw myruntime_error( strbuf );
    }
    if( expect <= 0. )
        bitscore = SLC_DBL_MAX;
    else
        bitscore = (log((double)scmatrix->GetSearchSpace())-log(expect))/LN2;


    if( lmbd <= 0.0 || K < 0.0 )
        pair_expect = Kref * pspace * exp(-lmbdref*score);
    else
        pair_expect = K * pspace * exp(-lmbd*score);


    SetBitScore( bitscore );
    SetRawExpectation( expect );//NOTE:changed
    SetExpectPerAlignment( pair_expect );

    SetReferenceExpectation( expect );
    SetExpectation( expect );
    SetPvalue( ComputePvalue( GetExpectation()));
}



// -------------------------------------------------------------------------
// ComputeStatisticsObs: compute e-value and p-value
//
void ProfileAlignment::ComputeStatisticsObs()
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
// ComputeStatistics: statistical significance
//
void ProfileAlignment::ComputeStatistics()
{
    int mod = GetModelSS18();
    if( mod == mss18_unused )
        ComputeStatisticsObs();
    else
        ComputeStatisticsSS18();
    return;
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
//     ComputeStatisticsTEST( nopros );
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
