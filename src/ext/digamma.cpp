/* Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author: G. Jungman */
// *** the code adopted from GSL ***

#include <math.h>
#include "psl.h"
#include "pslcodes.h"
#include "cheb.h"
#include "exp.h"
#include "gamma.h"
#include "zeta.h"
#include "digamma.h"


// Private-Section----------------------------------------------------------
//

/* Chebyshev fits from SLATEC code for psi(x)

 Series for PSI        on the interval  0.         to  1.00000D+00
                                       with weighted error   2.03E-17
                                        log weighted error  16.69
                              significant figures required  16.39
                                   decimal places required  17.37

 Series for APSI       on the interval  0.         to  2.50000D-01
                                       with weighted error   5.54E-17
                                        log weighted error  16.26
                              significant figures required  14.42
                                   decimal places required  16.86

*/

static double psics_data[23] = {
  -.038057080835217922,
   .491415393029387130, 
  -.056815747821244730,
   .008357821225914313,
  -.001333232857994342,
   .000220313287069308,
  -.000037040238178456,
   .000006283793654854,
  -.000001071263908506,
   .000000183128394654,
  -.000000031353509361,
   .000000005372808776,
  -.000000000921168141,
   .000000000157981265,
  -.000000000027098646,
   .000000000004648722,
  -.000000000000797527,
   .000000000000136827,
  -.000000000000023475,
   .000000000000004027,
  -.000000000000000691,
   .000000000000000118,
  -.000000000000000020
};

static Tcheb_series psi_cs = {
    psics_data,
    22,
    -1, 1,
    17
};

static double apsics_data[16] = {    
  -.0204749044678185,
  -.0101801271534859,
   .0000559718725387,
  -.0000012917176570,
   .0000000572858606,
  -.0000000038213539,
   .0000000003397434,
  -.0000000000374838,
   .0000000000048990,
  -.0000000000007344,
   .0000000000001233,
  -.0000000000000228,
   .0000000000000045,
  -.0000000000000009,
   .0000000000000002,
  -.0000000000000000 
};    

static Tcheb_series apsi_cs = {
    apsics_data,
    15,
    -1, 1,
    9
};

#define PSI_TABLE_NMAX 100
static double psi_table[PSI_TABLE_NMAX+1] = {
  0.0,  // Infinity                 // psi(0)
 -SLC_EULER,                        // psi(1)
  0.42278433509846713939348790992,  // ...
  0.92278433509846713939348790992,
  1.25611766843180047272682124325,
  1.50611766843180047272682124325,
  1.70611766843180047272682124325,
  1.87278433509846713939348790992,
  2.01564147795560999653634505277,
  2.14064147795560999653634505277,
  2.25175258906672110764745616389,
  2.35175258906672110764745616389,
  2.44266167997581201673836525479,
  2.52599501330914535007169858813,
  2.60291809023222227314862166505,
  2.67434666166079370172005023648,
  2.74101332832746036838671690315,
  2.80351332832746036838671690315,
  2.86233685773922507426906984432,
  2.91789241329478062982462539988,
  2.97052399224214905087725697883,
  3.02052399224214905087725697883,
  3.06814303986119666992487602645,
  3.11359758531574212447033057190,
  3.15707584618530734186163491973,
  3.1987425128519740085283015864,
  3.2387425128519740085283015864,
  3.2772040513135124700667631249,
  3.3142410883505495071038001619,
  3.3499553740648352213895144476,
  3.3844381326855248765619282407,
  3.4177714660188582098952615740,
  3.4500295305349872421533260902,
  3.4812795305349872421533260902,
  3.5115825608380175451836291205,
  3.5409943255438998981248055911,
  3.5695657541153284695533770196,
  3.5973435318931062473311547974,
  3.6243705589201332743581818244,
  3.6506863483938174848844976139,
  3.6763273740348431259101386396,
  3.7013273740348431259101386396,
  3.7257176179372821503003825420,
  3.7495271417468059598241920658,
  3.7727829557002943319172153216,
  3.7955102284275670591899425943,
  3.8177324506497892814121648166,
  3.8394715810845718901078169905,
  3.8607481768292527411716467777,
  3.8815815101625860745049801110,
  3.9019896734278921969539597029,
  3.9219896734278921969539597029,
  3.9415975165651470989147440166,
  3.9608282857959163296839747858,
  3.9796962103242182164764276160,
  3.9982147288427367349949461345,
  4.0163965470245549168131279527,
  4.0342536898816977739559850956,
  4.0517975495308205809735289552,
  4.0690389288411654085597358518,
  4.0859880813835382899156680552,
  4.1026547480502049565823347218,
  4.1190481906731557762544658694,
  4.1351772229312202923834981274,
  4.1510502388042361653993711433,
  4.1666752388042361653993711433,
  4.1820598541888515500147557587,
  4.1972113693403667015299072739,
  4.2121367424746950597388624977,
  4.2268426248276362362094507330,
  4.2413353784508246420065521823,
  4.2556210927365389277208378966,
  4.2697055997787924488475984600,
  4.2835944886676813377364873489,
  4.2972931188046676391063503626,
  4.3108066323181811526198638761,
  4.3241399656515144859531972094,
  4.3372978603883565912163551041,
  4.3502848733753695782293421171,
  4.3631053861958823987421626300,
  4.3757636140439836645649474401,
  4.3882636140439836645649474401,
  4.4006092930563293435772931191,
  4.4128044150075488557724150703,
  4.4248526077786331931218126607,
  4.4367573696833950978837174226,
  4.4485220755657480390601880108,
  4.4601499825424922251066996387,
  4.4716442354160554434975042364,
  4.4830078717796918071338678728,
  4.4942438268358715824147667492,
  4.5053549379469826935258778603,
  4.5163439489359936825368668713,
  4.5272135141533849868846929582,
  4.5379662023254279976373811303,
  4.5486045001977684231692960239,
  4.5591308159872421073798223397,
  4.5695474826539087740464890064,
  4.5798567610044242379640147796,
  4.5900608426370772991885045755,
  4.6001618527380874001986055856
};

#define PSI_1_TABLE_NMAX 100
static double psi_1_table[PSI_1_TABLE_NMAX+1] = {
  0.0,  // Infinity                 // psi(1,0)
  SLC_PI * SLC_PI / 6.0,            // psi(1,1)
  0.644934066848226436472415,       // ... 
  0.394934066848226436472415,
  0.2838229557371153253613041,
  0.2213229557371153253613041,
  0.1813229557371153253613041,
  0.1535451779593375475835263,
  0.1331370146940314251345467,
  0.1175120146940314251345467,
  0.1051663356816857461222010,
  0.0951663356816857461222010,
  0.0869018728717683907503002,
  0.0799574284273239463058557,
  0.0740402686640103368384001,
  0.0689382278476838062261552,
  0.0644937834032393617817108,
  0.0605875334032393617817108,
  0.0571273257907826143768665,
  0.0540409060376961946237801,
  0.0512708229352031198315363,
  0.0487708229352031198315363,
  0.0465032492390579951149830,
  0.0444371335365786562720078,
  0.0425467743683366902984728,
  0.0408106632572255791873617,
  0.0392106632572255791873617,
  0.0377313733163971768204978,
  0.0363596312039143235969038,
  0.0350841209998326909438426,
  0.0338950603577399442137594,
  0.0327839492466288331026483,
  0.0317433665203020901265817,
  0.03076680402030209012658168,
  0.02984853037475571730748159,
  0.02898347847164153045627052,
  0.02816715194102928555831133,
  0.02739554700275768062003973,
  0.02666508681283803124093089,
  0.02597256603721476254286995,
  0.02531510384129102815759710,
  0.02469010384129102815759710,
  0.02409521984367056414807896,
  0.02352832641963428296894063,
  0.02298749353699501850166102,
  0.02247096461137518379091722,
  0.02197713745088135663042339,
  0.02150454765882086513703965,
  0.02105185413233829383780923,
  0.02061782635456051606003145,
  0.02020133322669712580597065,
  0.01980133322669712580597065,
  0.01941686571420193164987683,
  0.01904704322899483105816086,
  0.01869104465298913508094477,
  0.01834810912486842177504628,
  0.01801753061247172756017024,
  0.01769865306145131939690494,
  0.01739086605006319997554452,
  0.01709360088954001329302371,
  0.01680632711763538818529605,
  0.01652854933985761040751827,
  0.01625980437882562975715546,
  0.01599965869724394401313881,
  0.01574770606433893015574400,
  0.01550356543933893015574400,
  0.01526687904880638577704578,
  0.01503731063741979257227076,
  0.01481454387422086185273411,
  0.01459828089844231513993134,
  0.01438824099085987447620523,
  0.01418415935820681325171544,
  0.01398578601958352422176106,
  0.01379288478501562298719316,
  0.01360523231738567365335942,
  0.01342261726990576130858221,
  0.01324483949212798353080444,
  0.01307170929822216635628920,
  0.01290304679189732236910755,
  0.01273868124291638877278934,
  0.01257845051066194236996928,
  0.01242220051066194236996928,
  0.01226978472038606978956995,
  0.01212106372098095378719041,
  0.01197590477193174490346273,
  0.01183418141592267460867815,
  0.01169577311142440471248438,
  0.01156056489076458859566448,
  0.01142844704164317229232189,
  0.01129931481023821361463594,
  0.01117306812421372175754719,
  0.01104961133409026496742374,
  0.01092885297157366069257770,
  0.01081070552355853781923177,
  0.01069508522063334415522437,
  0.01058191183901270133041676,
  0.01047110851491297833872701,
  0.01036260157046853389428257,
  0.01025632035036012704977199,  // ... 
  0.01015219706839427948625679,  // psi(1,99) 
  0.01005016666333357139524567   // psi(1,100)
};

//
//
//

// -------------------------------------------------------------------------
// * digamma for x both positive and negative; we do both
// * cases here because of the way we use even/odd parts
// * of the function
//
static int psi_x( const double x, double* result, double* err )
{
    const double y = fabs( x );

    if( x == 0.0 || x == -1.0 || x == -2.0 ) {
        return PSL_ERR_DOMAIN;
    }
    else if( y >= 2.0 ) {
        const double t = 8.0/( y * y ) - 1.0;
        double  pres, perr;
        cheb_eval_e( &apsi_cs, t, &pres, &perr );
        if( x < 0.0 ) {
            const double s = sin( SLC_PI * x );
            const double c = cos( SLC_PI * x );
            if( fabs( s ) < 2.0 * SLC_SQRT_DBL_MIN ) {
                return PSL_ERR_DOMAIN;
            }
            else {
                if( result ) {
                    *result = log( y ) - 0.5 / x + pres - SLC_PI * c / s;
                    if( err ) {
                        *err  = SLC_PI * fabs( x ) * SLC_DBL_EPSILON /( s * s );
                        *err += perr;
                        *err += SLC_DBL_EPSILON * fabs( *result );
                    }
                }
                return PSL_OK;
            }
        }
        else {
            if( result ) {
                *result = log( y ) - 0.5 / x + pres;
                if( err ) {
                    *err  = perr;
                    *err += SLC_DBL_EPSILON * fabs( *result );
                }
            }
            return PSL_OK;
        }
    }
    else { // -2 < x < 2
        double  pres, perr;

        if( x < -1.0 ) { // x = -2 + v
            const double v  = x + 2.0;
            const double t1 = 1.0 / x;
            const double t2 = 1.0 /( x + 1.0 );
            const double t3 = 1.0 / v;
            cheb_eval_e( &psi_cs, 2.0 * v - 1.0, &pres, &perr );

            if( result ) {
                *result = -( t1 + t2 + t3 ) + pres;
                if( err ) {
                    *err  = SLC_DBL_EPSILON * ( fabs( t1 ) + fabs( x /( t2 * t2 )) + fabs( x /( t3 * t3 )));
                    *err += perr;
                    *err += SLC_DBL_EPSILON * fabs( *result );
                }
            }
            return PSL_OK;
        }
        else if( x < 0.0 ) { // x = -1 + v
            const double v  = x + 1.0;
            const double t1 = 1.0 / x;
            const double t2 = 1.0 / v;
            cheb_eval_e( &psi_cs, 2.0 * v - 1.0, &pres, &perr );

            if( result ) {
                *result = -( t1 + t2 ) + pres;
                if( err ) {
                    *err  = SLC_DBL_EPSILON * ( fabs( t1 ) + fabs( x /( t2 * t2 )));
                    *err += perr;
                    *err += SLC_DBL_EPSILON * fabs( *result );
                }
            }
            return PSL_OK;
        }
        else if( x < 1.0 ) { // x = v
            const double t1 = 1.0 / x;
            cheb_eval_e( &psi_cs, 2.0 * x - 1.0, &pres, &perr );

            if( result ) {
                *result = -t1 + pres;
                if( err ) {
                    *err  = SLC_DBL_EPSILON * t1;
                    *err += perr;
                    *err += SLC_DBL_EPSILON * fabs( *result );
                }
            }
            return PSL_OK;
        }
        else { // x = 1 + v
            const double v = x - 1.0;
            return cheb_eval_e( &psi_cs, 2.0 * v - 1.0, result, err );
        }
    }
}


// -------------------------------------------------------------------------
// generic polygamma, (d/dx)^n [psi(x)];
// assumes n >= 0 and x > 0
//
static int psi_n_xg0( const int n, const double x, double* result, double* err )
{
    if( n == 0 ) {
        return psl_psi_e( x, result, err );
    }
    else {
        // Abramowitz + Stegun 6.4.10 
        double  ln_nf, ln_nferr;
        double  hzeta, hzetaerr;
        int     hzstatus = psl_hzeta_e( n + 1.0, x, &hzeta, &hzetaerr );
        int     nfstatus = psl_lnfact_e(( unsigned int ) n, &ln_nf, &ln_nferr );
        int     estatus  = psl_exp_mult_err_e( ln_nf, ln_nferr, hzeta, hzetaerr, result, err );
        if( result )
            if( SLC_EVEN( n ))
                *result = -*result;
        if( estatus != PSL_SUCCESS )
            return estatus;
        if( nfstatus != PSL_SUCCESS )
            return nfstatus;
        if( hzstatus != PSL_SUCCESS )
            return hzstatus;
        return PSL_SUCCESS;
    }
}

// =========================================================================
// Psi of natural numbers
//
int psl_psi_int_e( const int n, double* result, double* err )
{
    if( n <= 0 ) {
        return PSL_ERR_DOMAIN;
    }
    else if( n <= PSI_TABLE_NMAX ) {
        if( result ) {
            *result = psi_table[n];
            if( err )
                *err = SLC_DBL_EPSILON * fabs( *result );
        }
        return PSL_OK;
    }
    else {
        // Abramowitz+Stegun 6.3.18
        const double c2 = -1.0 / 12.0;
        const double c3 =  1.0 / 120.0;
        const double c4 = -1.0 / 252.0;
        const double c5 =  1.0 / 240.0;
        const double ni2 = ( 1.0 / n ) * ( 1.0 / n );
        const double ser = ni2 * ( c2 + ni2 * ( c3 + ni2 * ( c4 + ni2 * c5 )));
        if( result ) {
            *result = log( n ) - 0.5 / n + ser;
            if( err ) {
                *err  = SLC_DBL_EPSILON * ( fabs( log( n )) + fabs( 0.5 / n ) + fabs( ser ));
                *err += SLC_DBL_EPSILON * fabs( *result );
            }
        }
        return PSL_OK;
    }
}

// -------------------------------------------------------------------------
// Psi of real numbers
//
int psl_psi_e( const double x, double* result, double* err )
{
    return psi_x( x, result, err );
}


// -------------------------------------------------------------------------
// Trigamma function, psi'(n), n>0
//
int psl_psi_1_int_e( const int n, double* result, double* err )
{
    if( n <= 0 ) {
        return PSL_ERR_DOMAIN;
    }
    else if( n <= PSI_1_TABLE_NMAX ) {
        if( result ) {
            *result = psi_1_table[n];
            if( err )
                *err = SLC_DBL_EPSILON * *result;
        }
        return PSL_SUCCESS;
    }
    else {
        // Abramowitz+Stegun 6.4.12; double-precision for n > 100
        const double c0 = -1.0 / 30.0;
        const double c1 =  1.0 / 42.0;
        const double c2 = -1.0 / 30.0;
        const double ni2 = ( 1.0 / n ) * ( 1.0 / n );
        const double ser =  ni2 * ni2 *( c0 + ni2 *( c1 + c2 * ni2 ));
        if( result ) {
            *result = ( 1.0 + 0.5 / n + 1.0 /( 6.0 * n * n ) + ser ) / n;
            if( err )
                *err = SLC_DBL_EPSILON * *result;
        }
        return PSL_SUCCESS;
    }
}


// -------------------------------------------------------------------------
// Trigamma function, psi'(n) for general x
//
int psl_psi_1_e( const double x, double* result, double* err )
{
    if( x == 0.0 || x == -1.0 || x == -2.0 ) {
        return PSL_ERR_DOMAIN;
    }
    else if( 0.0 < x ) {
        return psi_n_xg0( 1, x, result, err );
    }
    else if( -5.0 < x ) {
        // Abramowitz + Stegun 6.4.6
        int     M = -floor( x );
        double  fx = x + M;
        double  sum = 0.0;
        int     m;

        if( fx == 0.0 )
            return PSL_ERR_DOMAIN;

        for( m = 0; m < M; m++ )
            sum += 1.0 /(( x + m )*( x + m ));

        int     psistatus = psi_n_xg0( 1, fx, result, err );
        if( result ) {
            *result += sum;
            if( err )
                *err += M * SLC_DBL_EPSILON * sum;
        }
        return psistatus;
    }
    else {
        // Abramowitz + Stegun 6.4.7
        const double sin_px = sin( SLC_PI * x );
        const double d = SLC_PI * SLC_PI /( sin_px * sin_px );
        double  r, rerr;
        int     psistatus = psi_n_xg0( 1, 1.0-x, &r, &rerr );
        if( result ) {
            *result = d - r;
            if( err )
                *err = rerr + 2.0 * SLC_DBL_EPSILON * d;
        }
        return psistatus;
    }
}


// -------------------------------------------------------------------------
// Polygamma function, psi^(n)(x) for n>=0, x>0
//
int psl_psi_n_e( const int n, const double x, double* result, double* err )
{
    if( n == 0 ) {
        return psl_psi_e( x, result, err );
    }
    else if( n == 1 ) {
        return psl_psi_1_e( x, result, err );
    }
    else if( n < 0 || x <= 0.0 ) {
        return PSL_ERR_DOMAIN;
    }
    else {
        double  ln_nf, ln_nferr;
        double  hzeta, hzetaerr;
        int     hzstatus = psl_hzeta_e( n+1.0, x, &hzeta, &hzetaerr );
        int     nfstatus = psl_lnfact_e(( unsigned int )n, &ln_nf, &ln_nferr );
        int     estatus  = psl_exp_mult_err_e( ln_nf, ln_nferr, hzeta, hzetaerr, result, err );
        if( result )
            if( SLC_EVEN( n ))
                *result = -*result;
        if( estatus != PSL_SUCCESS )
            return estatus;
        if( nfstatus != PSL_SUCCESS )
            return nfstatus;
        if( hzstatus != PSL_SUCCESS )
            return hzstatus;
        return PSL_SUCCESS;
    }
}


// -------------------------------------------------------------------------
// Functions w/ Natural Prototypes
//

// double psl_psi_int( const int n )
// {
//     double  result, err;
//     if( psl_psi_int_e( n, &result, &err ) != PSL_OK )
//         ;
//     return result;
// }
// 
// double psl_psi( const double x )
// {
//     double  result, err;
//     if( psl_psi_e( x, &result, &err ) != PSL_OK )
//         ;
//     return result;
// }

