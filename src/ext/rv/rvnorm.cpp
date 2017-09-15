/***************************************************************************
 *   Copyright (C) 2010 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rvnorm.h"

// =========================================================================
// RVNorm
//
RVNorm::RVNorm( Rng& rng, TRVN met )
:   RVar( rng ),
    met_( met ),
    std_( 1.0 ),
    mean_( 0.0 )
{
}

RVNorm::~RVNorm()
{
}

// -------------------------------------------------------------------------
// GenPolar: generate random value of normally distributed variable;
//  Polar form of the Box-Mueller method (GNU SL)
//  (Knuth TAOCP v2 p.122)
//
int RVNorm::GenPolar( double* rv1, double* rv2 )
{
    const int maxn = 100;
    double  r1, r2;
    double  x, y, R, f;
    bool    in = false;
    int n;

    if( rv1 == NULL )
        return PSL_ERR_ADDRESS;

    //1.27 times on average (Knuth TAOCP v2)
    for( n = 0; !in && n < maxn; n++ ) {
        //get x,y in uniform square: (-1,1)
        r1 = GetRng().GetDouble();
        r2 = GetRng().GetDouble();
        x = -1.0 + r1 + r1;
        y = -1.0 + r2 + r2;
        R = x * x + y * y;
        //check if in the unit circle
        in = R <= 1.0 && R != 0.0;
    }

    if( !in )
        r1 = 0.01;
    else
        r1 = log( R )/ R;
    //Box-Muller transform
    f = sqrt ( -r1 - r1 );
    *rv1 = x * f;
    if( rv2 )
        *rv2 = y * f;
    if( !in )
        return PSL_MAXITERATS;
    return PSL_OK;
}

// -------------------------------------------------------------------------
// GenRatio: generate random value of normally distributed variable;
//  Ratio method (Kinderman-Monahan), Leva implementation (adopted from GSL)
//  Leva JL. (1992) ACM Trans Math Software 18, 449-455
//  (Knuth TAOCP v2 p.130)
//
int RVNorm::GenRatio( double* rv )
{
    const int maxn = 100;
    double u, v, x, y, Q;
    const double s = 0.449871; //Leva constants
    const double t = -0.386595;
    const double a = 0.19600;
    const double b = 0.25472;
    const double r1 = 0.27597;
    const double r2 = 0.27846;
    int n = 0;

    if( rv == NULL )
        return PSL_ERR_ADDRESS;

    //1.369 times on average (GSL)
    do {
        if( maxn <= n++ )
            break;
        //generate point P=(u,v) uniform in rectangle
        //enclosing K&M region v^2 <= - 4 u^2 log(u)

        //u in (0,1] to avoid singularity at u=0
        u = 1.0 - GetRng().GetDouble1();
        if( !u )
            return PSL_ERR_DOMAIN;

        //v is in asymmetric interval [-0.5,0.5);  however v=-0.5 is
        //rejected in the last part of the while clause;
        //resulting normal variate is strictly symmetric about 0
        //(provided that v is symmetric once v = -0.5 is excluded)
        v = GetRng().GetDouble() - 0.5;

        //constant 1.7156 > sqrt(8/e) (for accuracy);
        //but not by too much (for efficiency)
        v *= 1.7156;

        //compute Leva's quadratic form Q
        x = u - s;
        y = fabs ( v ) - t;
        Q = x*x + y*( a*y - b*x );

        //accept P if Q < r1; /Leva/
        //reject P if Q > r2; /Leva/
        //accept if v^2 <= -4 u^2 log(u); /K&M/
        //this final test is executed 0.012 times on average (GSL)
    } while( r1 <= Q && ( r2 < Q || -4.0*u*u*log(u) < v*v ));

    //save slope
    if( rv )
        *rv = v / u;
    if( maxn < n )
        return PSL_MAXITERATS;
    return PSL_OK;
}

// =========================================================================
// DATA FOR ZIGGURAT METHOD BY MARSAGLIA & TSANG; GSL implementation
//
//position of right-most step
static const double grvn_ZM_R = 3.44428647676;//R
static const double grvn_1oZM_R = 0.29033589591;//1/R

//tabulated values for the heigt of the Ziggurat levels
static const double grvn_ytab[128] = {
  1, 0.963598623011, 0.936280813353, 0.913041104253,
  0.892278506696, 0.873239356919, 0.855496407634, 0.838778928349,
  0.822902083699, 0.807732738234, 0.793171045519, 0.779139726505,
  0.765577436082, 0.752434456248, 0.739669787677, 0.727249120285,
  0.715143377413, 0.703327646455, 0.691780377035, 0.68048276891,
  0.669418297233, 0.65857233912, 0.647931876189, 0.637485254896,
  0.62722199145, 0.617132611532, 0.607208517467, 0.597441877296,
  0.587825531465, 0.578352913803, 0.569017984198, 0.559815170911,
  0.550739320877, 0.541785656682, 0.532949739145, 0.524227434628,
  0.515614886373, 0.507108489253, 0.498704867478, 0.490400854812,
  0.482193476986, 0.47407993601, 0.466057596125, 0.458123971214,
  0.450276713467, 0.442513603171, 0.434832539473, 0.427231532022,
  0.419708693379, 0.41226223212, 0.404890446548, 0.397591718955,
  0.390364510382, 0.383207355816, 0.376118859788, 0.369097692334,
  0.362142585282, 0.355252328834, 0.348425768415, 0.341661801776,
  0.334959376311, 0.328317486588, 0.321735172063, 0.31521151497,
  0.308745638367, 0.302336704338, 0.29598391232, 0.289686497571,
  0.283443729739, 0.27725491156, 0.271119377649, 0.265036493387,
  0.259005653912, 0.253026283183, 0.247097833139, 0.241219782932,
  0.235391638239, 0.229612930649, 0.223883217122, 0.218202079518,
  0.212569124201, 0.206983981709, 0.201446306496, 0.195955776745,
  0.190512094256, 0.185114984406, 0.179764196185, 0.174459502324,
  0.169200699492, 0.1639876086, 0.158820075195, 0.153697969964,
  0.148621189348, 0.143589656295, 0.138603321143, 0.133662162669,
  0.128766189309, 0.123915440582, 0.119109988745, 0.114349940703,
  0.10963544023, 0.104966670533, 0.100343857232, 0.0957672718266,
  0.0912372357329, 0.0867541250127, 0.082318375932, 0.0779304915295,
  0.0735910494266, 0.0693007111742, 0.065060233529, 0.0608704821745,
  0.056732448584, 0.05264727098, 0.0486162607163, 0.0446409359769,
  0.0407230655415, 0.0368647267386, 0.0330683839378, 0.0293369977411,
  0.0256741818288, 0.0220844372634, 0.0185735200577, 0.0151490552854,
  0.0118216532614, 0.00860719483079, 0.00553245272614, 0.00265435214565
};

//tabulated values for 2^24 times x[i]/x[i+1],
//used to accept for U*x[i+1]<=x[i] without any floating point operations
static const unsigned long  grvn_ktab[128] = {
  0, 12590644, 14272653, 14988939,
  15384584, 15635009, 15807561, 15933577,
  16029594, 16105155, 16166147, 16216399,
  16258508, 16294295, 16325078, 16351831,
  16375291, 16396026, 16414479, 16431002,
  16445880, 16459343, 16471578, 16482744,
  16492970, 16502368, 16511031, 16519039,
  16526459, 16533352, 16539769, 16545755,
  16551348, 16556584, 16561493, 16566101,
  16570433, 16574511, 16578353, 16581977,
  16585398, 16588629, 16591685, 16594575,
  16597311, 16599901, 16602354, 16604679,
  16606881, 16608968, 16610945, 16612818,
  16614592, 16616272, 16617861, 16619363,
  16620782, 16622121, 16623383, 16624570,
  16625685, 16626730, 16627708, 16628619,
  16629465, 16630248, 16630969, 16631628,
  16632228, 16632768, 16633248, 16633671,
  16634034, 16634340, 16634586, 16634774,
  16634903, 16634972, 16634980, 16634926,
  16634810, 16634628, 16634381, 16634066,
  16633680, 16633222, 16632688, 16632075,
  16631380, 16630598, 16629726, 16628757,
  16627686, 16626507, 16625212, 16623794,
  16622243, 16620548, 16618698, 16616679,
  16614476, 16612071, 16609444, 16606571,
  16603425, 16599973, 16596178, 16591995,
  16587369, 16582237, 16576520, 16570120,
  16562917, 16554758, 16545450, 16534739,
  16522287, 16507638, 16490152, 16468907,
  16442518, 16408804, 16364095, 16301683,
  16207738, 16047994, 15704248, 15472926
};

//tabulated values of 2^{-24}*x[i]
static const double grvn_wtab[128] = {
  1.62318314817e-08, 2.16291505214e-08, 2.54246305087e-08, 2.84579525938e-08,
  3.10340022482e-08, 3.33011726243e-08, 3.53439060345e-08, 3.72152672658e-08,
  3.8950989572e-08, 4.05763964764e-08, 4.21101548915e-08, 4.35664624904e-08,
  4.49563968336e-08, 4.62887864029e-08, 4.75707945735e-08, 4.88083237257e-08,
  5.00063025384e-08, 5.11688950428e-08, 5.22996558616e-08, 5.34016475624e-08,
  5.44775307871e-08, 5.55296344581e-08, 5.65600111659e-08, 5.75704813695e-08,
  5.85626690412e-08, 5.95380306862e-08, 6.04978791776e-08, 6.14434034901e-08,
  6.23756851626e-08, 6.32957121259e-08, 6.42043903937e-08, 6.51025540077e-08,
  6.59909735447e-08, 6.68703634341e-08, 6.77413882848e-08, 6.8604668381e-08,
  6.94607844804e-08, 7.03102820203e-08, 7.11536748229e-08, 7.1991448372e-08,
  7.2824062723e-08, 7.36519550992e-08, 7.44755422158e-08, 7.52952223703e-08,
  7.61113773308e-08, 7.69243740467e-08, 7.77345662086e-08, 7.85422956743e-08,
  7.93478937793e-08, 8.01516825471e-08, 8.09539758128e-08, 8.17550802699e-08,
  8.25552964535e-08, 8.33549196661e-08, 8.41542408569e-08, 8.49535474601e-08,
  8.57531242006e-08, 8.65532538723e-08, 8.73542180955e-08, 8.8156298059e-08,
  8.89597752521e-08, 8.97649321908e-08, 9.05720531451e-08, 9.138142487e-08,
  9.21933373471e-08, 9.30080845407e-08, 9.38259651738e-08, 9.46472835298e-08,
  9.54723502847e-08, 9.63014833769e-08, 9.71350089201e-08, 9.79732621669e-08,
  9.88165885297e-08, 9.96653446693e-08, 1.00519899658e-07, 1.0138063623e-07,
  1.02247952126e-07, 1.03122261554e-07, 1.04003996769e-07, 1.04893609795e-07,
  1.05791574313e-07, 1.06698387725e-07, 1.07614573423e-07, 1.08540683296e-07,
  1.09477300508e-07, 1.1042504257e-07, 1.11384564771e-07, 1.12356564007e-07,
  1.13341783071e-07, 1.14341015475e-07, 1.15355110887e-07, 1.16384981291e-07,
  1.17431607977e-07, 1.18496049514e-07, 1.19579450872e-07, 1.20683053909e-07,
  1.21808209468e-07, 1.2295639141e-07, 1.24129212952e-07, 1.25328445797e-07,
  1.26556042658e-07, 1.27814163916e-07, 1.29105209375e-07, 1.30431856341e-07,
  1.31797105598e-07, 1.3320433736e-07, 1.34657379914e-07, 1.36160594606e-07,
  1.37718982103e-07, 1.39338316679e-07, 1.41025317971e-07, 1.42787873535e-07,
  1.44635331499e-07, 1.4657889173e-07, 1.48632138436e-07, 1.50811780719e-07,
  1.53138707402e-07, 1.55639532047e-07, 1.58348931426e-07, 1.61313325908e-07,
  1.64596952856e-07, 1.68292495203e-07, 1.72541128694e-07, 1.77574279496e-07,
  1.83813550477e-07, 1.92166040885e-07, 2.05295471952e-07, 2.22600839893e-07
};

// -------------------------------------------------------------------------
// GenZigMsg00: generate random normal variates by implementation of the
//  ziggurat algorithm from
//  Marsaglia G, Tsang WW. (2000) Journal of Statistical Software 5(8), 1-7.
//  (Usually the Zigurrat algorithm is the fastest one)
// 
// GNU SL modified the algorithm:
// 
// 1) use 128 steps instead of 256 to decrease the amount of static
// data necessary.  
// 
// 2) use an acceptance sampling from an exponential wedge
// exp(-R*(x-R/2)) for the tail of the base strip to simplify the
// implementation.  The area of exponential wedge is used in
// calculating 'v' and the coefficients in ziggurat table, so the
// coefficients differ slightly from those in the Marsaglia and Tsang
// paper.
// 
// For r.n.g. issues see
// Leong et al. (2005) Journal of Statistical Software 5(7), 1-4.
// For issues of using r.n. in the rejection method see
// Doornik et al. (2005) Tech report, University of Oxford.
// 
int RVNorm::GenZigMsg00( double* rv )
{
    if( rv == NULL )
        return PSL_ERR_ADDRESS;

    const int   maxn = 100;
    unsigned long   i, j;
    double  x, y;
    int n, sign;

    const unsigned long range = GetRng().GetMax() - GetRng().GetMin();
    const unsigned long offset = GetRng().GetMin();

    for( n = 0; n < maxn; n++ ) {
        if( 0xFFFFFFFF <= range ) {
            unsigned long   k = GetRng().Get() - offset;
            i = k & 0xFF;
            j = ( k >> 8 ) & 0xFFFFFF;
        }
        else if( 0x00FFFFFF <= range ) {
            unsigned long   k1 = GetRng().Get() - offset;
            unsigned long   k2 = GetRng().Get() - offset;
            i = (k1 & 0xFF);
            j = (k2 & 0x00FFFFFF);
        }
        else {
            i = GetRng().GetI( 256 );//choose the step
            j = GetRng().GetI( 16777216 );//sample from 2^24
        }

        sign = ( i & 0x80 )? 1: -1;
        i &= 0x7f;

        x = j * grvn_wtab[i];

        if( j < grvn_ktab[i])
            break;

        if( i < 127 ) {
            double  y0, y1, U1;
            y0 = grvn_ytab[i];
            y1 = grvn_ytab[i+1];
            U1 = GetRng().GetDouble();
            y = y1 + ( y0 - y1 )* U1;
        }
        else {
            double  U1, U2;
            U1 = 1.0 - GetRng().GetDouble0();
            U2 = GetRng().GetDouble();
            x = grvn_ZM_R - log( U1 ) * grvn_1oZM_R;
            y = exp( -grvn_ZM_R *( x - 0.5 * grvn_ZM_R )) * U2;
        }

        if( y < exp( -0.5 * x * x ))
            break;
    }
    if( sign < 0 )
        *rv = -x;
    else
        *rv = x;
    if( maxn <= n )
        return PSL_MAXITERATS;
    return PSL_OK;
}
