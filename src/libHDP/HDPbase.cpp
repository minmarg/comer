/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "data.h"
#include "ext/psl.h"
#include "ext/gamma.h"
#include "logitnormal.h"
#include "BinarySearchStructure.h"
#include "HDPbase.h"


// =========================================================================
//global object for application only
HDPbase  HDPBASE;
HDPbase  HDPctBASE;
//
void SetHDPBASEhlp( HDPbase& obj, const char* filename )
{
    mystring msg = "Reading ";
    msg += filename;
    message( msg.c_str());
    obj.ReadParameters( filename );
    obj.SetPriorProbs();
}
void SetHDPBASE( const char* filename )
{
    SetHDPBASEhlp( HDPBASE, filename );
}
void SetHDPctBASE( const char* filename )
{
    SetHDPBASEhlp( HDPctBASE, filename );
}

// -------------------------------------------------------------------------

double HDPbase::s_defdpmtau_ = 1.0;//determines creation rate of new table
double HDPbase::s_defdpmgamma_ = 1.0;//determines creation rate of new dish

//default dimensionality of samples
int HDPbase::s_defsampledim_ = NUMAA;
//default dish size
int HDPbase::s_defdishsize_ = 50;
//default table size
int HDPbase::s_deftablesize_ = 3;
//default basin size
int HDPbase::s_defbasinsize_ = 1000;

double HDPbase::s_defkappa_pp_a_ = 1.0;//prior param. a for kappa
double HDPbase::s_defkappa_pp_b_ = 1.0;//prior param. b for kappa
//default value of the NIW kappa parameter (no. measurements)
double HDPbase::s_defkappa_ = 1.0;
//default value of the NIW nu parameter (deg. of freedom)
double HDPbase::s_defnu_ = 1.0;

// -------------------------------------------------------------------------
// constructor: initialization
//
HDPbase::HDPbase()
:   uninfprior_( true ),
    clst4each( false ),
    basin_( NULL ),
    menu_( NULL ),
    chain_( NULL ),
    ctxtsize_( 0 ),
    rdnorsts_( 0 ),
    totalnos_( 0 ),
    restype_( TRT_novalues ),
    mhupdate_( 0 ),
    gibbsit_( 0 ),
    restarted_( false ),
    mixweight_( 1.0 ),
    nosupclsts_( -1 ),
    adjweight_( 0.1 ),
    scores_( NULL ),
    tau_( 1.0 ),
    gamma_( 1.0 ),
    iterread_( 0 ),
    mlpd_( -1.e9 ),
    lastlpd_( -1.e9 ),
    S0scalefac_( 1.0 ),
    adjdegf_( 0.0 )
{
}

// -------------------------------------------------------------------------
// destructor:
//
HDPbase::~HDPbase()
{
    DestroyBasin();
    DestroyMenu();
    DestroyChain();
    ctxtsize_ = 0;
    totalnos_ = 0;
    adjdegf_ = 0.0;
}

// =========================================================================
// Table2Mtx: copy data from table to matrix
//
void HDPbase::Table2Mtx( const Table* table, Pslmatrix* mtx ) const
{
    if( GetMenu() == NULL )
        throw myruntime_error( "Table2Mtx: Null menu." );

    if( table == NULL || mtx == NULL )
        throw myruntime_error( "Table2Mtx: Memory access error." );

    const int   dim = GetMenu()->GetDim();
    const int   nov = table->GetActualSize();//number of vectors at the table
    const Pslvector* vec = NULL;
    int n, i, j;
    mtx->Reserve( dim, nov );

    for( n = 0, j = 0; n < table->GetSize(); n++ ){
        if( table->GetVectorNIndAt( n ) < 0 )
            continue;
        vec = table->GetVectorNAt( n );
        if( vec == NULL )
            continue;
        for( i = 0; i < vec->GetSize(); i++ )
            mtx->SetValueAt( i, j, vec->GetValueAt( i ));
        j++;
    }
}





// -------------------------------------------------------------------------
// SetPriorProbs: calculate and set prior probabilities for each dish
//
void HDPbase::SetPriorProbs()
{
    mystring        preamb = "HDPbase: SetPriorProbs: ";
    Menu*           menu = GetMenu();
    const double    gamma = GetDPMGamma();

    if( menu == NULL )
        throw myruntime_error( preamb + "Null HDP menu.");
    if( menu->GetSize() != menu->GetActualSize() || menu->GetSize() < 1 )
        throw myruntime_error( preamb + "Invalid HDP menu size.");
    if( GetDPMGamma() <= 0.0 )
        throw myruntime_error( preamb + "Invalid HDP parameter gamma value.");

    const Dish* dish = NULL;
    double  sum;
    int nodshs = menu->GetSize();
    int k, nk;

    //normalize prior probabilities
    sum = gamma;
    for( k = 0; k < nodshs; k++ ) {
        dish = menu->GetDishAt(k);
        if( dish == NULL )
            throw myruntime_error( preamb + "Null HDP dish.");
        nk = dish->GetDishSize();
        if( nk < 1 )
            throw myruntime_error( preamb + "Invalid HDP dish size.");
        sum += ( double )nk;
    }
    for( k = 0; k < nodshs; k++ ) {
        dish = menu->GetDishAt(k);
        nk = dish->GetDishSize();
        menu->SetPriorProbAt( k, (double)nk / sum );
    }
    menu->SetPriorProbNewDish( gamma / sum );
}

// -------------------------------------------------------------------------
// PPrCompare: compare posterior probabilities
static inline
int PPrCompare( const void* key1, const void* key2, void* pars )
{
    void**  params = (void**)pars;
    if( params == NULL )
        throw myruntime_error("PPrCompare: Null parameters.");
    const Menu* menu = (const Menu*)params[0];
    const double* lprior = (const double*)params[1];
    if( menu == NULL || lprior == NULL )
        throw myruntime_error("PPrCompare: Null parameters.");
    Dish*   dish = NULL;
    int     nodshs = menu->GetSize();
    int     k1 = (int)(ssize_t)key1;
    int     k2 = (int)(ssize_t)key2;
    double  pr1 = *lprior;
    double  pr2 = *lprior;
    double  diff;
    if( 0 < k1 ) {
        if( nodshs <= k1 || ( dish = menu->GetDishAt(k1)) == NULL )
            throw myruntime_error("PPrCompare: Null dish.");
        pr1 = dish->GetTmpValue();
    }
    if( 0 < k2 ) {
        if( nodshs <= k2 || ( dish = menu->GetDishAt(k2)) == NULL )
            throw myruntime_error("PPrCompare: Null dish.");
        pr2 = dish->GetTmpValue();
    }
    //key2-key1 -- to sort in descending order
    diff = pr2 - pr1;
    return ( diff < 0.0 )? -1: (( 0.0 < diff )? 1: 0 );
}

// =========================================================================
// CalcPPProbs: calculate posterior predictive probabilities of central 
//  vector of logit-normal variable;
//  lnvar -- input logit-normal variable
//  ppr -- background posterior probability of vector `lnvar'
//  ppprobs -- output posterior probabilities
//  cindcs -- indices of clusters 
//  NOTE: input value `lnvar' will be changed on exit
//
void HDPbase::CalcPPProbs( Pslvector& lnvar, double* bppr, Pslvector* ppprobs, Ivector* cindcs,
                           double lpfact, bool usepriors, bool tonormal ) const
{
    mystring errstr;
    mystring preamb = "HDPbase: CalcPPProbs: ";
    if( GetDPMGamma() <= 0.0 )
        throw myruntime_error( preamb + "Invalid HDP parameter gamma value.");

    const Menu*     menu = GetMenu();
    const double    tau0 = GetDPMTau();
    const double    gamma = GetDPMGamma();
    int             nospcls = GetNoSupClusters();
//     const double    cdLPFACT = 0.02;//log-space factor in calculating probability threshold; TRAINING value
    const double    cdLPFACT = lpfact;//0.02;//log-space factor in calculating probability threshold; WORKING value
    const double    cdMINPRB = 0.1;//minimum probability threshold
    const bool      cbPRIORS = usepriors;//if false, posteriors are approximated to likelihoods

    if( menu == NULL )
        throw myruntime_error( preamb + "Null HDP menu.");
    if( menu->GetSize() != menu->GetActualSize() || menu->GetSize() < 1 )
        throw myruntime_error( preamb + "Invalid HDP menu size.");

    int     nogrps = GetReadNoGroups();
    int     nodshs = menu->GetSize();
    int     dim = menu->GetDim();
    int     adim = dim + 1;
    int     ctxtsz = menu->GetCtx();

    Dish*       dish = NULL;
    Pslvector   tvec( dim*ctxtsz );//normal transform
    double* pv, *ptv;
    double  thlprob;//threshold log probability
    double  lprob, lprior, maxlp;
    double  pnorm;
    bool    pset;
    size_t  c;
    int     n, d, k, nk, mk, mtot;
    int     err;

    void*   params[] = { (void*)menu, (void*)&lprior };
    BinarySearchStructure dshnxs( PPrCompare, nodshs+1, true/*keep*/, (void*)params );

    if( dim < 1 )
        throw myruntime_error( preamb + "Invalid HDP dimensions.");
    if( ctxtsz < 1 )
        throw myruntime_error( preamb + "Wrong HDP context size.");
    if( nogrps < 1 )
        throw myruntime_error( preamb + "Invalid no. groups.");
    if( nospcls == -1 )
        nospcls = nodshs + 1;//1 for prior
    if( nospcls < 2 || nodshs+1 < nospcls )
        throw myruntime_error( preamb + "Wrong number of HDP support clusters.");

    if( bppr == NULL || ppprobs == NULL || cindcs == NULL )
        throw myruntime_error( preamb + "Null output data and/or vectors.");

    if( tonormal ) {
        if( lnvar.GetSize() != adim * ctxtsz )
            throw myruntime_error( preamb + "Invalid size of input vector.");
        //set normal vector
        pv = lnvar.GetVector();
        ptv = tvec.GetVector();
        for( n = 0; n < ctxtsz; n++, pv += adim, ptv += dim ) {
            //TRANSFORM to multivariate normal
            ::LogitNormal2Normal( pv, adim, 1.e-1, false );
            //reduce dimenions of resulting vector
            for( d = 0; d < dim; d++ )
                ptv[d] = pv[d];
        }
    }
    else {
        if( lnvar.GetSize() != dim * ctxtsz )
            throw myruntime_error( preamb + "Invalid size of input vector.");
        tvec = lnvar;
    }

    ppprobs->Allocate( nospcls );
    cindcs->Allocate( nospcls );
    ppprobs->Clear();
    cindcs->Clear();
    *bppr = 0.0;

    mtot = 0; //total number of tables (local clusters)
    maxlp = -LOC_DBL_MAX;
    //calculate and save log prior probability to lprior
    PriorProbVec( &tvec, &lprior );
    if( !isfinite( lprior ))
        throw myruntime_error( preamb + "Invalid HDP prior prob. value.");
    if( maxlp < lprior )
        maxlp = lprior;
    //calculate max exponent in log posterior predictive probabilities over all dishes
    for( k = 0; k < nodshs; k++ ) {
        dish = menu->GetDishAt(k);
        if( dish == NULL )
            throw myruntime_error( preamb + "Null HDP dish.");
        nk = dish->GetDishSize();
        mk = dish->GetReadNoTables();
        if( nk < 1 )
            throw myruntime_error( preamb + "Invalid HDP dish size.");
        if( mk < 1 || nk < mk )
            throw myruntime_error( preamb + "Invalid no. tables for HDP dish.");
        mtot += mk;
        ProbVecOfDish( &tvec, k, &lprob );
        //density values can be >1
        if( !isfinite( lprob ))
            throw myruntime_error( preamb + "Invalid HDP posterior predictive prob. value.");
        if( maxlp < lprob )
            maxlp = lprob;
        dish->SetTmpValue( lprob );
    }
    //prior probability and set threshold probability
    pnorm = 0.0;
    lprior -= maxlp;
    if( lprior < SLC_LOG_DBL_MIN || !isfinite( lprior ))
        lprior = 0.0;
    else {
        lprior = exp( lprior );
        if( cbPRIORS )///
            ///lprior *= gamma;//apprx.
            lprior *= ((double)nogrps)*tau0*gamma;
        if( lprior )
            pnorm += lprior;
    }
    //exponentiate log probabilities
    for( k = 0; k < nodshs; k++ ) {
        pset = true;
        dish = menu->GetDishAt(k);
        if( dish == NULL )
            throw myruntime_error( preamb + "Null HDP dish.");
        nk = dish->GetDishSize();
        mk = dish->GetReadNoTables();
        if( nk < 1 )
            throw myruntime_error( preamb + "Invalid HDP dish size.");
        if( mk < 1 || nk < mk || mtot < mk )
            throw myruntime_error( preamb + "Invalid no. tables for HDP dish.");
        lprob = dish->GetTmpValue();
        lprob -= maxlp;
        if( lprob < SLC_LOG_DBL_MIN || !isfinite( lprob ))
            lprob = 0.0;
        else {
            lprob = exp( lprob );
            if( cbPRIORS )///
                ///lprob *= ( double )nk;//apprx.
                lprob *= (double)nk *((double)mtot+gamma ) + (double)(nogrps)*(double)mk*tau0;
        }
        if( pset )
            dish->SetTmpValue( lprob );
        if( lprob )
            pnorm += lprob;
    }

    if( pnorm <= 0.0 )
        return;

    //normalize probabilities
    if( pnorm <= 0.0 )
        throw myruntime_error( preamb + "Failed to normalize HDP posterior probabilities.");
    if( lprior )
        lprior /= pnorm;
    dshnxs.Push((const void*)-1 );
    thlprob = pow( lprior, cdLPFACT );
//     thlprob = cdMINPRB;

    for( k = 0; k < nodshs; k++ ) {
        dish = menu->GetDishAt(k);
        if( dish == NULL )
            throw myruntime_error( preamb + "Null HDP dish.");
        if( dish->GetTmpValue()) {
            lprob = dish->GetTmpValue() / pnorm;
//             if( lprob < thlprob ) {
//                 dish->SetTmpValue( 0.0 );
//                 continue;
//             }
            dish->SetTmpValue( lprob );
            dshnxs.Push((const void*)k );
        }
    }

    if( dshnxs.GetSize() < 1 )
        throw myruntime_error( preamb + "Too few clusters selected.");

    //write top probabilities and calculate background posterior probability
    *bppr = 0.0;
    for( c = 0; c < dshnxs.GetSize() && c < (size_t)(nospcls-1); c++ ) {
        k = (int)(ssize_t)dshnxs.GetValueAt(c);
        if( k < 0 )
            continue;//prior will be placed at the end
        if( nodshs < k || ( dish = menu->GetDishAt(k)) == NULL )
            throw myruntime_error( preamb + "Invalid HDP dish.");
        if( c == 0 && dish->GetTmpValue() < thlprob )
            break;//the greatest probability is less than the threshold
        cindcs->Push( k );
        ppprobs->Push( dish->GetTmpValue());
        *bppr += dish->GetTmpValue() * ( cbPRIORS? 1.0: menu->GetPriorProbAt(k));//already mixed
    }
    cindcs->Push( -1 );
    ppprobs->Push( lprior );
    *bppr += lprior * ( cbPRIORS? 1.0: menu->GetPriorProbNewDish());//already mixed
}

// =========================================================================
// MixLNVar: mix central vector of logit-normal variable using clusters 
//  encapsulated in this class object;
//  lnvar -- input logit-normal variable
//  lnmixed -- output mixed logit-normal variable
//  NOTE: input value `lnvar' will be changed on exit
  //
void HDPbase::MixCLNVar( Pslvector& lnvar, Pslvector* lnmixed ) const
{
    mystring errstr;
    mystring preamb = "HDPbase: MixLNVar: ";
    if( GetDPMGamma() <= 0.0 )
        throw myruntime_error( preamb + "Invalid HDP parameter gamma value.");

    const Menu*     menu = GetMenu();
    const double    gamma = GetDPMGamma();
    const double    mixwgt = GetMixWeight();
    const double    lpfact = 0.1;//log-space factor in calculating probability threshold
    const bool      sngdsh = true;//use single dish of highest probability in mixing

    if( menu == NULL )
        throw myruntime_error( preamb + "Null HDP menu.");
    if( menu->GetSize() != menu->GetActualSize() || menu->GetSize() < 1 )
        throw myruntime_error( preamb + "Invalid HDP menu size.");

    if( mixwgt <= 0.0 )
        return;
    if( 1.0 < mixwgt )
        throw myruntime_error( preamb + "Invalid HDP weight of mixing.");

    int     nodshs = menu->GetSize();
    int     dim = menu->GetDim();
    int     adim = dim + 1;
    int     ctxtsz = menu->GetCtx();
    int     parity = ( ctxtsz & 1 ) ^ 1;
    int     hlf = ctxtsz >> 1;
    int     mid = hlf - parity;

    const Pslvector* mu = NULL;
    Dish*       dish = NULL;
    Pslvector   tvec( dim*ctxtsz );//normal transform
    Pslvector   subv;
    Pslvector   mxdv;
    double* pv, *ptv;
    double  thlprob;//threshold log probability
    double  lprob, lprior, maxlp;
    double  pnorm, sum, lnv;
    int     n, d, k, kmaxp, nk;
    int     sgn, err;

    if( dim < 1 )
        throw myruntime_error( preamb + "Invalid HDP dimensions.");
    if( ctxtsz < 1 )
        throw myruntime_error( preamb + "Wrong HDP context size.");

    if( lnvar.GetSize() != adim * ctxtsz )
        throw myruntime_error( preamb + "Invalid size of input vector.");
    if( lnmixed == NULL || lnmixed->GetSize() != adim )
        throw myruntime_error( preamb + "Invalid size of output vector.");

    //copy input central vector to output vector
    subv = lnvar.SubVector( mid*adim, adim );
    lnmixed->Copy( subv );

    //set normal vector
    pv = lnvar.GetVector();
    ptv = tvec.GetVector();
    for( n = 0; n < ctxtsz; n++, pv += adim, ptv += dim ) {
        //TRANSFORM to multivariate normal
        ::LogitNormal2Normal( pv, adim, 1.e-1, false );
        //reduce dimenions of resulting vector
        for( d = 0; d < dim; d++ )
            ptv[d] = pv[d];
    }

    maxlp = -LOC_DBL_MAX;
    //calculate and save log prior probability to lprior
    PriorProbVec( &tvec, &lprior );
    if( !isfinite( lprior ))
        throw myruntime_error( preamb + "Invalid HDP prior prob. value.");
    if( maxlp < lprior )
        maxlp = lprior;
    //get threshold probability
    thlprob = lprior * lpfact;
    //calculate log posterior predictive probabilities for each dish
    for( k = 0, kmaxp = -1; k < nodshs; k++ ) {
        dish = menu->GetDishAt(k);
        if( dish == NULL )
            throw myruntime_error( preamb + "Null HDP dish.");
        nk = dish->GetDishSize();
        if( nk < 1 )
            throw myruntime_error( preamb + "Invalid HDP dish size.");
        ProbVecOfDish( &tvec, k, &lprob );
        //density values can be >1
        if( !isfinite( lprob ))
            throw myruntime_error( preamb + "Invalid HDP posterior predictive prob. value.");
        if( maxlp < lprob && thlprob <= lprob ) {
            maxlp = lprob;
            kmaxp = k;
        }
        dish->SetTmpValue( lprob );
    }
    //exponentiate log probabilities
    for( k = 0; k < nodshs; k++ ) {
        dish = menu->GetDishAt(k);
        if( dish == NULL )
            throw myruntime_error( preamb + "Null HDP dish.");
        lprob = dish->GetTmpValue();
        if( lprob < thlprob ) {
            dish->SetTmpValue(0.0);
            continue;
        }
        if( sngdsh && kmaxp != k ) {
            dish->SetTmpValue(0.0);
            continue;
        }
        lprob -= maxlp;
        if( lprob < SLC_LOG_DBL_MIN || !isfinite( lprob ))
            lprob = 0.0;
        else
            lprob = exp( lprob );
        dish->SetTmpValue( lprob );
    }
    //check threshold of prior probability
    if( lprior < thlprob || ( sngdsh && 0 < kmaxp ))
        lprior = 0.0;
    else {
        lprior -= maxlp;
        if( lprior < SLC_LOG_DBL_MIN || !isfinite( lprior ))
            lprior = 0.0;
        else
            lprior = exp( lprior );
    }
    //normalize probabilities
    pnorm = 0.0;
    for( k = 0; k < nodshs; k++ ) {
        dish = menu->GetDishAt(k);
        if( dish == NULL )
            throw myruntime_error( preamb + "Null HDP dish.");
        nk = dish->GetDishSize();
        if( nk < 1 )
            throw myruntime_error( preamb + "Invalid HDP dish size.");
        if( dish->GetTmpValue()) {
            dish->SetTmpValue(( double )nk * dish->GetTmpValue());
            pnorm += dish->GetTmpValue();
        }
    }
    if( lprior ) {
        lprior *= gamma;
        pnorm += lprior;
    }

    if( pnorm <= 0.0 )
        return;

    for( k = 0; k < nodshs; k++ ) {
        dish = menu->GetDishAt(k);
        if( dish == NULL )
            throw myruntime_error( preamb + "Null HDP dish.");
        if( dish->GetTmpValue())
            dish->SetTmpValue( dish->GetTmpValue() / pnorm );
    }
    if( lprior )
        lprior /= pnorm;

    //mix central vector with mean vectors of dishes
    mxdv = lnmixed->SubVector( 0, dim );
    mxdv.Zero();
    for( k = 0; k < nodshs; k++ ) {
        dish = menu->GetDishAt(k);
        mu = menu->GetMuVectorAt(k);
        if( mu == NULL || dish == NULL )
            throw myruntime_error( preamb + "Null HDP dish mean vector.");
        subv = mu->SubVector( mid*dim, dim );
        if(( err = mxdv.Superposition( dish->GetTmpValue(), subv )) != 0 )
            throw myruntime_error( preamb + TranslatePSLError( err ));
    }
    subv = tvec.SubVector( mid*dim, dim );
    if(( err = mxdv.Superposition( lprior, subv )) != 0 )
        throw myruntime_error( preamb + TranslatePSLError( err ));

    //apply weight of mixing
    if( mixwgt < 1.0 ) {
        if(( err = mxdv.MultiplyBy( mixwgt )) != 0 )
            throw myruntime_error( preamb + TranslatePSLError( err ));
        if(( err = mxdv.Superposition( 1.0 - mixwgt, subv )) != 0 )
            throw myruntime_error( preamb + TranslatePSLError( err ));
    }

    //back-transform to logit-normal space
    if(( err = mxdv.Exp()) != 0 )
        throw myruntime_error( preamb + TranslatePSLError( err ));
    pnorm = 1.0 + mxdv.Sum();
    if( !isfinite( pnorm ))
        throw myruntime_error( preamb + "Failed to back-transform HDP-mixed vector.");
    sum = 0.0;
    for( d = 0; d < dim; d++ ) {
        lnv = mxdv.GetValueAt(d) / pnorm;
        lnmixed->SetValueAt( d, lnv );
        sum += lnv;
    }
    lnmixed->SetValueAt( dim, lnv = 1.0 - sum );
}





// =========================================================================
// CalcPriorParams: calculate prior parameters of NIW distribution;
//  return flag of addition of noise
//
bool HDPbase::CalcPriorParams( bool addnoise, double stdfact )
{
    SPDmatrix*  S0 = NULL;
    Pslvector   stds;
    double  val;
    int novs;
    int n, err;

    if( GetBasin() == NULL )
        throw myruntime_error("CalcPriorParams: Null basin.");
    if( GetMenu() == NULL )
        throw myruntime_error("CalcPriorParams: Null menu.");

    if(( err = CalcPriorParams( stds )) == 0 )
        return false;

    if( !addnoise )
        throw myruntime_error( TranslatePSLError( err ));

    novs = GetBasin()->GetSize();
    S0 = GetMenu()->GetS0();
    if( S0 == NULL )
        throw myruntime_error("CalcPriorParams: Null scatter matrix.");
    if( novs < 1 )
        throw myruntime_error("CalcPriorParams: No vectors in basin.");
    if( stdfact < 0 )
        throw myruntime_error("CalcPriorParams: Negative std. deviation factor.");
    //set standard deviations for each dimension
    stds.Reserve( S0->GetNoRows());
    for( n = 0; n < stds.GetSize(); n++ ) {
        val = S0->GetValueAt( n, n ) /( double )novs;
        if( val < 0 )
            throw myruntime_error("CalcPriorParams: Negative standard deviation.");
        stds.SetValueAt( n, stdfact * sqrt( val ));
    }

    if(( err = CalcPriorParams( stds )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));
    return true;
}

// =========================================================================
// CalcPriorParams: calculate prior parameters of NIW distribution
//
int HDPbase::CalcPriorParams( const Pslvector& stds )
{
    if( GetBasin() == NULL )
        return PSL_ERR_ADDRESS;
    if( GetMenu() == NULL )
        return PSL_ERR_ADDRESS;
    if( GetS0ScaleFac() <= 0.0 )
        throw myruntime_error("HDPbase: SetUninfPriorParams: Invalid scale factor.");

    const Pslvector*    frq = NULL;
    Pslvector*          mu = NULL;
    SPDmatrix*          S0 = NULL;
    SPDmatrix*          S0inv = NULL;
    mystring    preamb = "CalcPriorParams: ";
    mystring    errstr;
    double      ldet = 0.0;
    double      inos = 0.0;
    int         err, dim = 0;
    int         n, nos;

    try {
        for( n = 0, nos = 0; n < GetBasin()->GetSize(); n++ ) {
            frq = GetBasin()->GetValueAt( n );
            if( frq == NULL )
                continue;
            if( mu == NULL ) {
                mu = new Pslvector( dim = frq->GetSize());
                if( mu == NULL )
                    throw myruntime_error( "Not enough memory." );
                if( stds.GetSize() && frq->GetSize() != stds.GetSize())
                    throw myruntime_error( "Inconsistent size of vector of std. devs." );
            }
            if( frq->GetSize() != dim )
                throw myruntime_error( "Inconsistent vector lengths." );
            if(( err = mu->Superposition( 1.0, *frq )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
            nos++;
        }
        if( mu == NULL || nos == 0 || dim == 0 )
            return 0;

        inos = 1.0 /( double )nos;

        //mean vector
        if(( err = mu->MultiplyBy( inos )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        Pslvector   dev( dim );
        Pslvector   adev( dim );
        Pslmatrix   Sdev( dim, dim );

        S0 = new SPDmatrix( dim );
        if( S0 == NULL )
            throw myruntime_error( "Not enough memory." );

        S0inv = new SPDmatrix( dim );
        if( S0inv == NULL )
            throw myruntime_error( "Not enough memory." );

        for( n = 0; n < GetBasin()->GetSize(); n++ ) {
            frq = GetBasin()->GetValueAt( n );
            if( frq == NULL )
                continue;

            //deviation of sample from the mean vector
            dev = *frq;
            if( stds.GetSize())
                dev.AddGNoise( stds );
            if(( err = dev.Superposition( -1.0, *mu )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
            //sum of deviations
            if(( err = adev.Superposition( 1.0, dev )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            if(( err = Sdev.Mul( dev, dev )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            if(( err = S0->Add( Sdev )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
        }

        //take square of sum of deviations
        if(( err = Sdev.Mul( adev, adev )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        //divide by number of observations
        if(( err = Sdev.Scale( inos )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        //correct the round off error of covariance matrix S0
        if(( err = S0->Sub( Sdev )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        //S0 is now covariance matrix multiplied by the number of observations;
        //multiplication below is valid only if S0 is cov. matrix;
        ////scale covariance matrix S0 to obtain valid NIW dist. parameter of scale matrix
        ////if(( err = S0->Scale( GetDefKappaParam())) != 0 )
        ////    throw myruntime_error( TranslatePSLError( err ));

            //NOTE:divide scale matrix by number of observations
            if(( err = S0->Scale( inos )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

        if(( err = S0->Scale( GetS0ScaleFac())) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        //SAVE scatter matrix here
        GetMenu()->SetS0ScaleFac( GetS0ScaleFac());
        GetMenu()->SetS0( S0 );
        //solve for inverse matrix and calculate determinant
        *S0inv = *S0;
        if(( err = S0inv->CholeskyDecompose()) == 0 ) {
            if(( err = S0inv->CDedLogDet( &ldet )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
            if(( err = S0inv->CDedInvert()) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
        }
        else {
            return err;//return for adding noise
            SPDmatrix   S0invhlp( dim );//for LU-decomposition
            S0invhlp = *S0;
            if(( err = S0invhlp.LUDecompose()) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
            if(( err = S0invhlp.LUDedLogDet( &ldet )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
            if(( err = S0invhlp.LUDedInvert( *S0inv )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
        }


        //SAVE
        /*GetMenu()->SetS0( S0 );*/ S0 = NULL;//saved above
        GetMenu()->SetInvS0( S0inv ); S0inv = NULL;
        GetMenu()->SetLDetS0( ldet );
        GetMenu()->SetMu0( mu ); mu = NULL;
        GetMenu()->SetKappa0_pp_a( GetDefKappa_pp_a());
        GetMenu()->SetKappa0_pp_b( GetDefKappa_pp_b());
        if( !GetMenu()->GetKappa0())
            GetMenu()->SetKappa0( inos / GetS0ScaleFac()/*1*//*dim*//*GetDefKappaParam()*/);
        if( !GetMenu()->GetNu0())
            GetMenu()->SetNu0( GetBasin()->GetSize()/*1*//*dim*//*GetDefNuParam()*/);
        GetMenu()->SetDim( dim );
        GetMenu()->SetCtx( GetCtxtSize());

        //probability factors
        CalcPriorProbFact();

    } catch( myexception const& ex ) {
        errstr = preamb + ex.what();
    }

    if( mu ) { delete mu; mu = NULL; }
    if( S0 ) { delete S0; S0 = NULL; }
    if( S0inv ) { delete S0inv; S0inv = NULL; }

    if( !errstr.empty())
        throw myruntime_error( errstr );
    return 0;
}

// -------------------------------------------------------------------------
// SetUninfPriorParamsS0I: set uninformative prior parameters of NIW 
//     distribution; scale matrix S0 is proportional to I
//
int HDPbase::SetUninfPriorParamsS0I( int dim, int ctx, const Pslvector& mvec, 
                                     double kp0preset, double nu0preset )
{
    if( GetMenu() == NULL )
        throw myruntime_error("HDPbase: SetUninfPriorParamsS0I: Null menu.");
    if( GetS0ScaleFac() <= 0.0 )
        throw myruntime_error("HDPbase: SetUninfPriorParamsS0I: Invalid scale factor.");
    if( dim < 1 || ctx < 1 )
        throw myruntime_error("HDPbase: SetUninfPriorParamsS0I: Invalid dimensionality.");
    if( dim != mvec.GetSize())
        throw myruntime_error("HDPbase: SetUninfPriorParamsS0I: Inconsistent dimensionality.");

    Pslvector*          mu0 = NULL;
    SPDmatrix*          S0 = NULL;
    SPDmatrix*          S0inv = NULL;
    mystring    preamb = "SetUninfPriorParamsS0I: ";
    mystring    errstr;
    double      scale = GetS0ScaleFac();
    double      ldet = 0.0;
    int         err;

    try {
        mu0 = new Pslvector( mvec );
        S0 = new SPDmatrix( dim );
        S0inv = new SPDmatrix( dim );
        if( mu0 == NULL || S0 == NULL || S0inv == NULL )
            throw myruntime_error( "Not enough memory." );

        S0->SetIdentity();
        S0inv->SetIdentity();

        if(( err = S0->Scale( scale )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = S0inv->Scale( 1.0 / scale )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        ldet = dim * log( scale );

        GetMenu()->SetS0ScaleFac( GetS0ScaleFac());
        GetMenu()->SetS0( S0 ); S0 = NULL;
        GetMenu()->SetInvS0( S0inv ); S0inv = NULL;
        GetMenu()->SetLDetS0( ldet );
        GetMenu()->SetMu0( mu0 ); mu0 = NULL;
        GetMenu()->SetKappa0_pp_a( GetDefKappa_pp_a());
        GetMenu()->SetKappa0_pp_b( GetDefKappa_pp_b());
        if( 0.0 < kp0preset )
                GetMenu()->SetKappa0( kp0preset );
        else    GetMenu()->SetKappa0( 1.0 / GetS0ScaleFac());
        if( 0.0 < nu0preset )
                GetMenu()->SetNu0( nu0preset );
        else    GetMenu()->SetNu0( GetS0ScaleFac());
        GetMenu()->SetDim( dim );
        GetMenu()->SetCtx( ctx );

        //probability factors
        CalcPriorProbFact();

    } catch( myexception const& ex ) {
        errstr = preamb + ex.what();
    }

    if( mu0 ) { delete mu0; mu0 = NULL; }
    if( S0 ) { delete S0; S0 = NULL; }
    if( S0inv ) { delete S0inv; S0inv = NULL; }

    if( !errstr.empty())
        throw myruntime_error( errstr );
    return 0;
}

// =========================================================================
// CalcPriorProbFact: recalculate probability factors
//
void HDPbase::CalcPriorProbFact()
{
    if( GetMenu() == NULL )
        throw myruntime_error("HDPbase: CalcPriorProbFact: Null menu.");
    if( GetMenu()->GetKappa0() <= 0.0 || GetMenu()->GetNu0() <= 0.0/*( double )dim*/ )
        throw myruntime_error( "Invalid number of degrees of freedom." );

    double  pfact = 0.0;
    double  term1, term2;
    double  operr;
    int     err;
    const double    do2 = ( double)GetMenu()->GetDim() * 0.5;
    const double    nu0p1o2 = ( GetMenu()->GetNu0() + 1.0 )* 0.5;
    const double    kp0 = GetMenu()->GetKappa0();
    const double    kp0p1 = kp0 + 1.0;

    const double    nu_t_o2 = nu0p1o2;//student's t nu over 2

//     const double    dofdom_o2 = nu_t_o2 - do2;
    const double    dofdom_o2 = NU_t_o2;
    const double    dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2

    if(( err = psl_lngamma_e( dofdom_pdim_o2, &term1, &operr )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));
    if(( err = psl_lngamma_e( dofdom_o2, &term2, &operr )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    pfact = term1 - term2;
    pfact -= do2 * SLC_LNPI;
    pfact += do2 * ( log( kp0 ) - log( kp0p1 ));

    GetMenu()->SetLogPriorProbFact( pfact );
}

// =========================================================================
// RecalcMenuParams: recalculate parameters for each dish of the menu
//
void HDPbase::RecalcMenuParams()
{
    if( GetMenu() == NULL )
        return;

    Dish*   dish = NULL;
    int     k;

    for( k = 0; k < GetMenu()->GetSize(); k++ )
    {
        dish = GetMenu()->GetDishAt( k );
        if( dish == NULL )
            continue;
        RecalcDishParams( k );
    }
}

// -------------------------------------------------------------------------
// RecalcDishParams: recalculate parameters for the given dish of the menu
//
void HDPbase::RecalcDishParams( int k )
{
    if( GetBasin() == NULL )
        throw myruntime_error( "RecalcDishParams: Null basin." );
    if( GetMenu() == NULL )
        throw myruntime_error( "RecalcDishParams: Null menu." );

    const int           dim = GetMenu()->GetDim();
    const double        kp0 = GetMenu()->GetKappa0();
    const Pslvector*    mu0 = GetMenu()->GetMu0();
    const SPDmatrix*    L0 = GetMenu()->GetS0();

    if( 0.0 < kp0 && ( mu0 == NULL || L0 == NULL ))
        throw myruntime_error( "RecalcDishParams: Null prior parameters." );
    if( dim <= 0 )
        throw myruntime_error( "RecalcDishParams: Invalid dimensionality." );

    const Pslvector*    val = NULL;
    Pslvector*          mu = NULL;
    SPDmatrix*          L = NULL;
    SPDmatrix*          Linv = NULL;
    Dish*       dish = NULL;
    mystring    preamb = "RecalcDishParams: ";
    mystring    errstr;
    char        bufstr[KBYTE];
    double      ldet = 0.0;
    double      inos = 0.0;
    int         err;
    int         n, nos;

    if( k < 0 || GetMenu()->GetSize() <= k )
        throw myruntime_error( "RecalcDishParams: Invalid dish index." );

    dish = GetMenu()->GetDishAt( k );
    if( dish == NULL )
        throw myruntime_error( "RecalcDishParams: Null dish at index." );

    try {
        mu = new Pslvector( dim );
        if( mu == NULL )
            throw myruntime_error( "Not enough memory." );
        L = new SPDmatrix( dim );
        if( L == NULL )
            throw myruntime_error( "Not enough memory." );
#ifdef PROBMTXINVSM
        Linv = new SPDmatrix( dim );
        if( Linv == NULL )
            throw myruntime_error( "Not enough memory." );
#endif
        for( n = 0, nos = 0; n < dish->GetSize(); n++ ) {
            if( dish->GetVectorNIndAt( n ) < 0 )
                continue;
            val = dish->GetVectorNAt( n );
            if( val == NULL )
                continue;
            if( val->GetSize() != dim )
                throw myruntime_error( "Inconsistent vector dimensionality." );
            if(( err = mu->Superposition( 1.0, *val )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
            nos++;
        }
        if( dish->GetActualSize() != nos ) {
            sprintf( bufstr, "Dish %d size inconsistency.", k );
            throw myruntime_error( bufstr );
        }
        if( nos != 0 &&
            nos != dish->GetDishSize()) {
            //when proposal size is defined
            if(( err = mu->MultiplyBy(( double )dish->GetDishSize()/( double )nos )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
            nos = dish->GetDishSize();
        }
        if( nos == 0 ) {
            sprintf( bufstr, "Dish %d has no members.", k );
            throw myruntime_error( bufstr );
        }
        inos = 1.0 /( double )nos;

        //mean vector
        if(( err = mu->MultiplyBy( inos )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        Pslvector   eta( *mu );//unweighted mean vector
        double      tot = kp0 + ( double )nos;
        double      tpr = kp0 * ( double )nos;

        if( 0.0 < kp0 ) {
            Pslvector   m0w( *mu0 );
            //weigh prior mean vector
            if(( err = m0w.MultiplyBy( kp0 / tot )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
            //weigh dish mean vector
            if(( err = mu->MultiplyBy(( double )nos / tot )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
            //weighted mean vector
            if(( err = mu->Superposition( 1.0, m0w )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
        }

        //proceed with SCALE matrix...
        Pslvector   dev( dim );
        Pslmatrix   Sdev( dim, dim );
        SPDmatrix   S( dim );

        //1.calculate S
        for( n = 0; n < dish->GetSize(); n++ ) {
            if( dish->GetVectorNIndAt( n ) < 0 )
                continue;
            val = dish->GetVectorNAt( n );
            if( val == NULL )
                continue;

            //deviation of sample from the mean vector
            dev = *val;
            if(( err = dev.Superposition( -1.0, eta )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            if(( err = Sdev.Mul( dev, dev )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            if(( err = S.Add( Sdev )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
        }

        *L = S;
        if( 0.0 < kp0 ) {
            //2.term of deviation of means
            dev = eta;
            if(( err = dev.Superposition( -1.0, *mu0 )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            if(( err = Sdev.Mul( dev, dev )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            if(( err = Sdev.Scale( tpr / tot )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            //3.final parameter scale matrix
            if(( err = L->Add( *L0 )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
            if(( err = L->Add( Sdev )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
        }

        //4.solve for inverse matrix and calculate determinant
#ifndef PROBMTXINVSM//if NOT defined
        Linv = &S;
#endif
        *Linv = *L;
        if(( err = Linv->CholeskyDecompose()) == 0 ) {
            if(( err = Linv->CDedLogDet( &ldet )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
#ifdef PROBMTXINVSM
            if(( err = Linv->CDedInvert()) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
#endif
        }
        else {
            SPDmatrix   Linvhlp( dim );//for LU-decomposition
            Linvhlp = *L;
            if(( err = Linvhlp.LUDecompose()) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
            if(( err = Linvhlp.LUDedLogDet( &ldet )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
#ifdef PROBMTXINVSM
            if(( err = Linvhlp.LUDedInvert( *Linv )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
#endif
        }

        //SAVE
        GetMenu()->SetSMatrixAt( k, L ); L = NULL;
#ifdef PROBMTXINVSM
        GetMenu()->SetInvSMatrixAt( k, Linv ); Linv = NULL;
#endif
        GetMenu()->SetLDetSMAt( k, ldet );
        GetMenu()->SetMuVectorAt( k, mu ); mu = NULL;

    } catch( myexception const& ex ) {
        errstr = preamb + ex.what();
    }

    if( mu ) { delete mu; mu = NULL; }
    if( L ) { delete L; L = NULL; }
#ifdef PROBMTXINVSM
    if( Linv ) { delete Linv; Linv = NULL; }
#endif
    if( !errstr.empty())
        throw myruntime_error( errstr );
}





// =========================================================================
// CalcProbOfData: calculate log probability of all data given current 
//  configuration of clusters
//
void HDPbase::CalcProbOfData( double* lprob )
{
    if( GetMenu() == NULL )
        throw myruntime_error( "CalcProbOfData: Null menu." );
    if( lprob == NULL )
        throw myruntime_error( "CalcProbOfData: Memory access error." );

    Dish*   dish = NULL;
    double  lp;
    int     k;

    *lprob = 0.0;
    for( k = 0; k < GetMenu()->GetSize(); k++ ) {
        dish = GetMenu()->GetDishAt( k );
        if( dish == NULL )
            continue;
        CalcProbOfDish( k, &lp );
        GetMenu()->SetProbAt( k, lp );
        *lprob += lp;
    }
    SetLastLProbData( *lprob );
}

// =========================================================================
// CalcProbOfDish: calculate log probability of dish k
//
void HDPbase::CalcProbOfDish( int k, double* lprob )
{
    mystring    preamb = "CalcProbOfDish: ";
    mystring    errstr;
    int         err;

    if( GetMenu() == NULL )
        throw myruntime_error( preamb + "Null menu." );

    *lprob = -1.e6;

    if( k < 0 || GetMenu()->GetSize() <= k )
        throw myruntime_error( preamb + "Invalid dish index." );

    const Dish*         dishk = GetMenu()->GetDishAt( k );
    const Pslvector*    mu = GetMenu()->GetMuVectorAt( k );
    const SPDmatrix*    L = GetMenu()->GetSMatrixAt( k );
    const double        ldet = GetMenu()->GetLDetSMAt( k );
    const double        L0ldet = GetMenu()->GetLDetS0();

    if( dishk == NULL || mu == NULL || L == NULL )
        throw myruntime_error( preamb + "Memory access error." );

    const int           dim = GetMenu()->GetDim();
    const double        kp0 = GetMenu()->GetKappa0();
    const double        nu0 = GetMenu()->GetNu0();
    const int           nk = dishk->GetActualSize();
    const double        kp = kp0 + ( double )nk;
    const double        nu = nu0 + ( double )nk;

    if( dim <= 0 )
        throw myruntime_error( preamb + "Invalid dimensionality." );
    if( nk <= 0 )
        throw myruntime_error( preamb + "Dish has no members." );
    if( kp0 <= 0.0 || nu0 <= 0.0/*( double )dim*/ )
        throw myruntime_error( preamb + "Invalid prior parameters kappa and nu." );

    double  lnsum = 0.0;
    double  term1, term2;
    double  operr;
    double  j;
    const double    do2 = ( double )dim * 0.5;
    const double    nuo2 = nu * 0.5;
    const double    nko2 = ( double )nk * 0.5;

    const double    nu_t_o2 = nuo2;//student's t nu over 2

//     const double    dofdom_o2 = nu_t_o2 - do2;
    const double    dofdom_o2 = NU_t_o2;
    const double    dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2
    const double    dofdom_pdim1_o2 = dofdom_pdim_o2 + 0.5;//(deg.of freedom+1 + dim.) over 2
    const double    dofdom_pdim1_mv_o2 = dofdom_pdim1_o2 - nko2;//(deg.of freedom+1 + dim.- L) / 2
    const double    dofdom_pdim_mv_o2 = dofdom_pdim_o2 - nko2;//(deg.of freedom + dim.- L) over 2

    //calculate probability next
    if( nk == 1 ) {
        if(( err = psl_lngamma_e( dofdom_pdim_o2, &term1, &operr )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));
        if(( err = psl_lngamma_e( dofdom_o2, &term2, &operr )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));
        lnsum += term1 - term2;
    }
    else
        for( j = 0.5; j <= do2; j += 0.5 )
        {
            if(( err = psl_lngamma_e( dofdom_pdim1_o2 - j, &term1, &operr )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
            if(( err = psl_lngamma_e( dofdom_pdim1_mv_o2 - j, &term2, &operr )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            lnsum += term1 - term2;
        }
    lnsum -= ( double )nk * do2 * SLC_LNPI;
    lnsum += do2 * ( log( kp0 ) - log( kp ));
    lnsum += dofdom_pdim_mv_o2 * L0ldet;
    lnsum -= dofdom_pdim_o2 * ldet;

    *lprob = lnsum;
}





// =========================================================================
// TestMenu: test menu structure
//
bool HDPbase::TestMenu() const
{
    if( GetBasin() == NULL || GetMenu() == NULL || GetChain() == NULL )
        return false;

    const Restaurant*   rest = NULL;
    const Table*        tbl = NULL;
    const Dish*         dsh = NULL;
    const Pslvector*    frq = NULL;
    mystring    preamb = "TestMenu: ";
    mystring    merror;
    const int   nods = GetMenu()->GetSize();
    int*    dishs = NULL;
    int*    vecs = NULL;
    int     notbls = 0, novecs = 0;
    int     nts = 0, nvs = 0;
    int     ddx;
    int     r, t, k, v;

    try {
        dishs = ( int* )malloc( nods * sizeof( int ));
        vecs = ( int* )malloc( nods * sizeof( int ));
        if( dishs == NULL || vecs == NULL )
            throw myruntime_error("Not enough memory.");

        memset( dishs, 0, nods * sizeof( int ));
        memset( vecs, 0, nods * sizeof( int ));

        for( r = 0; r < GetChain()->GetSize(); r++ ) {
            rest = GetChain()->GetRestaurantAt( r );
            if( rest == NULL )
                throw myruntime_error("Null restaurant.");
            for( t = 0; t < rest->GetSize(); t++ ) {
                tbl = rest->GetTableAt( t );
                if( tbl == NULL )
                    continue;
                k = tbl->GetDishIndex();
                dsh = tbl->GetDish();
                if( k < 0 || nods <= k || dsh == NULL )
                    throw myruntime_error("Invalid table's dish index.");
                dishs[k]++;
                vecs[k] += tbl->GetActualSize();
                notbls++;
                novecs += tbl->GetActualSize();
                for( v = 0; v < tbl->GetSize(); v++ ) {
                    if( tbl->GetVectorNIndAt( v ) < 0 )
                        continue;
                    frq = tbl->GetVectorNAt( v );
                    ddx = tbl->GetVectorNDishIndAt( v );
                    if( ddx < 0 || dsh->GetSize() <= ddx )
                        throw myruntime_error("Invalid dish index saved in table vector.");
                    if( dsh->GetVectorNAt( ddx ) != frq )
                        throw myruntime_error("Table vector not found in dish.");
                }
            }
        }
        for( k = 0; k < nods; k++ ) {
            dsh = GetMenu()->GetDishAt( k );
            if( dsh == NULL )
                continue;
            if( dsh->GetNoTables() < 1 || dishs[k] != dsh->GetNoTables())
                throw myruntime_error("Dish has invalid number of tables.");
            if( dsh->GetActualSize() < 1 || vecs[k] != dsh->GetActualSize())
                throw myruntime_error("Dish has invalid number of tables.");
            nts += dsh->GetNoTables();
            nvs += dsh->GetActualSize();
        }
        if( notbls != nts )
            throw myruntime_error("Inconsistent number of tables.");
        if( novecs != nvs )
            throw myruntime_error("Inconsistent number of vectors.");
    } catch( myexception const& ex ) {
        merror = preamb + ex.what();
    }

    if( dishs ) {
        free( dishs );
        dishs = NULL;
    }
    if( vecs ) {
        free( vecs );
        vecs = NULL;
    }
    if( !merror.empty()) {
        error( merror.c_str());
        return false;
    }
    return true;
}





// =========================================================================
// LogitNormalErrCorr: correct error of logit-normal variable to sum 
//     exactly to 1;
//  result is written back to `f'
//
void HDPbase::LogitNormalErrCorr( Pslvector* f, double tol )
{
    int             dim;
    double          vv, vsum, corr;
    int             r, non;
    char            strbuf[KBYTE];
    const double    accuracy = 1.0e-6;

    if( f == NULL )
        return;

    dim = f->GetSize();
    if( dim <= 1 )
        throw myruntime_error("LogitNormalErrCorr: Invalid dimensionality.");

    vsum = 0.0;
    non = 0;
    for( r = 0; r < dim; r++ ) {
        vv = f->GetValueAt( r );
        if( vv < 0.0 )
            throw myruntime_error("LogitNormalErrCorr: Invalid logit-normal value.");
        if( vv != 0.0 )
            non++;
        vsum += vv;
    }

    corr = vsum - 1.0;
    if( corr == 0.0 )
        return;

    if( tol < corr ) {
        sprintf( strbuf, "LogitNormalErrCorr: Error tolerance exceeded: %g.", corr );
        throw myruntime_error( strbuf );
    }

    if( fabs( corr ) <= accuracy ) {
        for( r = 0; r < dim; r++ ) {
            vv = f->GetValueAt( r );
            //correct small error with the first value met
            if(( 0.0 < corr && corr <= vv ) ||
               ( corr < 0.0 && (( non < dim && vv == 0.0 ) || 
                                ( dim <= non && vv - corr <= 1.0 ))))
            {
                vv -= corr;
                f->SetValueAt( r, vv );
                break;
            }
        }
    }
    else if( accuracy < fabs( corr ) && non ) {
        corr /= ( double )non;
        for( r = 0; r < dim; r++ ) {
            vv = f->GetValueAt( r );
            if( vv == 0.0 )
                continue;
            vv -= corr;
            non--;
            //correct overflow/undeflow error
            if( vv < 0.0 ) {
                if( non )
                    corr -= vv /( double )non;
                vv = 0.0;
            }
            else if( 1.0 < vv ) {
                if( non )
                    corr += ( 1.0 - vv )/( double )non;
                vv = 1.0;
            }
            f->SetValueAt( r, vv );
        }
    }
}

// -------------------------------------------------------------------------
// LogitNormal2Normal: tranform logit-normally distributed variable to 
//      normally distributed variable;
//  tranformed variable is written back to `f' and 
//  NOTE: dimensionality is reduced by 1
//
void HDPbase::LogitNormal2Normal( Pslvector* f, double tol )
{
    const bool      bCH0 = false;//check 0s
    const bool      bCHACC = false;//check accuracy
    int             dim;
    double          vv, vlast, vsum, corr;
    int             r, noz, non;
    double          fake = 1.0e-4;//fake value
    const double    accuracy = 1.0e-6;

    if( f == NULL )
        return;

    dim = f->GetSize();
    if( dim <= 1 )
        throw myruntime_error("LogitNormal2Normal: Invalid dimensionality.");

    vsum = f->Sum();
    if( bCHACC ) {
        if( vsum < 1.0 - accuracy || 1.0 + accuracy < vsum )
            LogitNormalErrCorr( f, tol );
    } else {
        if( vsum < 1.0 - tol || 1.0 + tol < vsum )
            throw myruntime_error("LogitNormal2Normal: Not a logit-normal variable.");
    }

    vsum = 0.0;
    non = noz = 0;
    for( r = 0; r < dim; r++ ) {
        vv = f->GetValueAt( r );
        if( vv < 0.0 )
            throw myruntime_error("LogitNormal2Normal: Invalid logit-normal value.");
        if( vv == 0.0 ) {
            noz++;
            continue;
        }
        non++;
        vsum += vv;
    }

    if( bCHACC )
        if( vsum < 1.0 - accuracy || 1.0 + accuracy < vsum )
            throw myruntime_error("LogitNormal2Normal: Not a logit-normal variable.");
    if(( bCH0? noz < 0: 0 < noz )|| non <= 0 || noz + non != dim )
        throw myruntime_error("LogitNormal2Normal: Invalid logit-normal vector.");

    if( !bCH0 ) {
        //make sure there are no null frequencies
        corr = fake * ( double )noz /( double )non;
        if( 0.0 < corr ) {
            vsum = 0.0;
            for( r = 0; r < dim; r++ ) {
                vv = f->GetValueAt( r );
                if( vv == 0.0 ) {
                    vsum += fake;
                    //null frequency, add fake value
                    f->SetValueAt( r, fake );
                    continue;
                }
                vv -= corr;
                non--;
                if( vv <= 0.0 ) {
                    vsum += fake;
                    //less than correction value; recalc. `corr'
                    f->SetValueAt( r, fake );
                    corr += ( fake - vv ) /( double )non;
                    continue;
                }
                vsum += vv;
                f->SetValueAt( r, vv );
            }
        }
    }
    vlast = f->GetValueAt( dim - 1 );
    if( vlast <= 0.0 )
        throw myruntime_error("LogitNormal2Normal: Not a logit-normal variable: last el.");

    if( bCHACC )
        if( vsum < 1.0 - accuracy || 1.0 + accuracy < vsum || vlast <= 0.0 )
            throw myruntime_error("LogitNormal2Normal: Not a logit-normal variable: correct. 2.");

    //apply transformation to normal distribution
    for( r = 0; r < dim - 1; r++ ) {
        vv = f->GetValueAt( r );
        if( vv <= 0.0 )
            throw myruntime_error("LogitNormal2Normal: Invalid logit-normal value: correct. 2.");
        vv = log( vv / vlast );
        f->SetValueAt( r, vv );
    }
    f->SetValueAt( dim - 1, 0.0 );
}
