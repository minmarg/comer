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
#include "HDPbase.h"



// =========================================================================
// PriorProbVec: compute prior probability of t-distributed vector
//
void HDPbase::PriorProbVec( const Pslvector* vec, double* lprob ) const
{
    PriorProbVecHlpObs( vec, lprob );
//     PriorProbVecHlp( vec, lprob );
}

// -------------------------------------------------------------------------
// PriorProbVecHlpObs: compute prior probability of t-distributed vector;
//  a version requiring inverse of scale matrix
//
void HDPbase::PriorProbVecHlpObs( const Pslvector* vec, double* lprob ) const
{
    mystring    preamb = "PriorProbVecHlpObs: ";
    mystring    errstr;
    char        bufstr[KBYTE];
    int         err;

//     if( GetBasin() == NULL )
//         throw myruntime_error( preamb + "Nil Basin." );
    if( GetMenu() == NULL )
        throw myruntime_error( preamb + "Nil Menu." );
    if( lprob == NULL )
        throw myruntime_error( preamb + "Memory access error." );
    if( vec == NULL )
        throw myruntime_error( preamb + "Invalid vector address." );

    *lprob = -1.e6;

    const Pslvector*    mu0 = GetMenu()->GetMu0();
    const SPDmatrix*    L0 = GetMenu()->GetS0();
    const SPDmatrix*    L0inv = GetMenu()->GetInvS0();
    const double        L0ldet = GetMenu()->GetLDetS0();
    const double        lpfact = GetMenu()->GetLogPriorProbFact();

    if( mu0 == NULL || L0 == NULL || L0inv == NULL )
        throw myruntime_error( preamb + "Memory access error." );

    const int           dim = GetMenu()->GetDim();
    const double        kp0 = GetMenu()->GetKappa0();
    const double        nu0 = GetMenu()->GetNu0();

    if( dim <= 0 )
        throw myruntime_error( preamb + "Invalid dimensionality." );
    if( kp0 <= 0.0 || nu0 <= 0.0/*( double )dim*/ )
        throw myruntime_error( preamb + "Number of degrees of freedom is too small." );

    double  lnsum = 0.0;
    const double    do2 = ( double)dim * 0.5;
    const double    no2 = nu0 * 0.5;
    const double    nu0p1o2 = no2 + 0.5;//(nu+1)/2
    const double    kp01 = kp0 + 1.0;

    const double    nu_t_o2 = nu0p1o2;//student's t nu over 2

//     const double    dofdom_o2 = nu_t_o2 - do2;
    const double    dofdom_o2 = NU_t_o2;
    const double    dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2

    Pslvector   dev, scl;
    Pslmatrix   trv, pdt;

    dev = *vec;
    if(( err = dev.Superposition( -1.0, *mu0 )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = dev.Transpose( trv )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = pdt.Mul( trv, *L0inv )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = scl.Mul( pdt, dev )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if( scl.GetValueAt( 0 ) <= -kp01 / kp0 )
        throw myruntime_error( preamb + "Invalid determinant." );

    lnsum = lpfact;
    lnsum -= 0.5 * L0ldet;
    lnsum -= dofdom_pdim_o2 * log( 1.0 + kp0 / kp01  * scl.GetValueAt( 0 ));

    *lprob = lnsum;
//     *prob = exp( lnsum );
}

// -------------------------------------------------------------------------
// PriorProbVecHlp: compute prior probability of t-distributed vector;
//
void HDPbase::PriorProbVecHlp( const Pslvector* vec, double* lprob ) const
{
    mystring    preamb = "PriorProbVecHlp: ";
    mystring    errstr;
    char        bufstr[KBYTE];
    int         err;

//     if( GetBasin() == NULL )
//         throw myruntime_error( preamb + "Null Basin." );
    if( GetMenu() == NULL )
        throw myruntime_error( preamb + "Nil Menu." );
    if( lprob == NULL )
        throw myruntime_error( preamb + "Memory access error." );
    if( vec == NULL )
        throw myruntime_error( preamb + "Invalid vector address." );

    *lprob = -1.e6;

    const Pslvector*    mu0 = GetMenu()->GetMu0();
    const SPDmatrix*    L0 = GetMenu()->GetS0();
    const double        L0ldet = GetMenu()->GetLDetS0();
    const double        lpfact = GetMenu()->GetLogPriorProbFact();

    if( mu0 == NULL || L0 == NULL )
        throw myruntime_error( preamb + "Memory access error." );

    const int           dim = GetMenu()->GetDim();
    const double        kp0 = GetMenu()->GetKappa0();
    const double        nu0 = GetMenu()->GetNu0();

    if( dim <= 0 )
        throw myruntime_error( preamb + "Invalid dimensionality." );
    if( kp0 <= 0.0 || nu0 <= 0.0/*( double )dim*/ )
        throw myruntime_error( preamb + "Number of degrees of freedom is too small." );

    double  lnsum = 0.0;
    const double    do2 = ( double)dim * 0.5;
    const double    no2 = nu0 * 0.5;
    const double    nu0p1o2 = no2 + 0.5;//(nu+1)/2
    const double    kp01 = kp0 + 1.0;

    const double    nu_t_o2 = nu0p1o2;//student's t nu over 2

//     const double    dofdom_o2 = nu_t_o2 - do2;
    const double    dofdom_o2 = NU_t_o2;
    const double    dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2
    const double    dofdom_pdim_m1_o2 = dofdom_pdim_o2 - 0.5;//(deg.of freedom-1 + dim.) over 2

    Pslvector   dev;
    SPDmatrix   MCM( dim );
    double      ldetMCM = 0.0;

    dev = *vec;
    if(( err = dev.Superposition( -1.0, *mu0 )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = MCM.Mul( dev, dev )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = MCM.Scale( kp0 / kp01 )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = MCM.Add( *L0 )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    //factorize by Cholesky decomp. and calculate determinant
    if(( err = MCM.CholeskyDecompose()) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = MCM.CDedLogDet( &ldetMCM )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));


    lnsum = lpfact;
    lnsum += dofdom_pdim_m1_o2 * L0ldet;
    lnsum -= dofdom_pdim_o2 * ldetMCM;

    *lprob = lnsum;
//     *prob = exp( lnsum );
}





// =========================================================================
// ProbVecOfDish: compute probability of t-distributed vector to get dish k
//
void HDPbase::ProbVecOfDish( const Pslvector* vec, int k, double* lprob ) const
{
#ifdef PROBMTXINVSM
    ProbVecOfDishHlpObs( vec, k, lprob );
#else
    ProbVecOfDishHlp( vec, k, lprob );
#endif
}

#ifdef PROBMTXINVSM

// -------------------------------------------------------------------------
// ProbVecOfDishHlpObs: compute probability of t-distributed vector to get 
//     dish k; requires inverse of scale matrix
//
void HDPbase::ProbVecOfDishHlpObs( const Pslvector* vec, int k, double* lprob ) const
{
    mystring    preamb = "ProbVecOfDishHlpObs: ";
    mystring    errstr;
    char        bufstr[KBYTE];
    int         err;

//     if( GetBasin() == NULL )
//         throw myruntime_error( preamb + "Nil Basin." );
    if( GetMenu() == NULL )
        throw myruntime_error( preamb + "Nil Menu." );
    if( lprob == NULL )
        throw myruntime_error( preamb + "Memory access error." );
    if( vec == NULL )
        throw myruntime_error( preamb + "Invalid vector address." );

    *lprob = -1.e6;

    if( k < 0 || GetMenu()->GetSize() <= k )
        throw myruntime_error( preamb + "Invalid cluster index." );

    const Dish*         dishk = GetMenu()->GetDishAt( k );
    const Pslvector*    mu = GetMenu()->GetMuVectorAt( k );
    const SPDmatrix*    L = GetMenu()->GetSMatrixAt( k );
    const SPDmatrix*    Linv = GetMenu()->GetInvSMatrixAt( k );
    const double        ldet = GetMenu()->GetLDetSMAt( k );

    if( dishk == NULL || 
        mu == NULL || L == NULL || Linv == NULL )
        throw myruntime_error( preamb + "Memory access error." );

    const int           dim = GetMenu()->GetDim();
    const double        kp0 = GetMenu()->GetKappa0();
    const double        nu0 = GetMenu()->GetNu0();
    const int           nk = dishk->GetDishSize();
    const double        kp = kp0 + nk;
    const double        nu = nu0 + nk;

    if( dim <= 0 )
        throw myruntime_error( preamb + "Invalid dimensionality." );
    if( nk <= 0 )
        throw myruntime_error( preamb + "Cluster has no members." );
    if( kp <= 0.0 || nu <= 0.0/*( double )dim*/ )
        throw myruntime_error( preamb + "Cluster has too few members." );

    double  lnsum = 0.0;
    double  term1, term2;
    double  operr;
    const double    do2 = ( double)dim * 0.5;
    const double    no2 = nu * 0.5;
    const double    nup1o2 = no2 + 0.5;//(nu+1)/2
    const double    kp1 = kp + 1.0;

    const double    nu_t_o2 = nup1o2;//student's t nu over 2

//     const double    dofdom_o2 = nu_t_o2 - do2;
    const double    dofdom_o2 = NU_t_o2;
    const double    dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2

    if(( err = psl_lngamma_e( dofdom_pdim_o2, &term1, &operr )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));
    if(( err = psl_lngamma_e( dofdom_o2, &term2, &operr )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    Pslvector   dev, scl;
    Pslmatrix   trv, pdt;

    dev = *vec;
    if(( err = dev.Superposition( -1.0, *mu )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = dev.Transpose( trv )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = pdt.Mul( trv, *Linv )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = scl.Mul( pdt, dev )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if( scl.GetValueAt( 0 ) <= -kp1 / kp )
        throw myruntime_error( preamb + "Invalid determinant." );

    lnsum = term1 - term2;
    lnsum -= do2 * SLC_LNPI;
    lnsum += do2 * ( log( kp ) - log( kp1 ));
    lnsum -= 0.5 * ldet;
    lnsum -= dofdom_pdim_o2 * log( 1.0 + kp / kp1  * scl.GetValueAt( 0 ));

    *lprob = lnsum;
//     *prob = exp( lnsum );
}

#endif

// -------------------------------------------------------------------------
// ProbVecOfDishHlp: compute probability of t-distributed vector to get 
//     dish k;
//
void HDPbase::ProbVecOfDishHlp( const Pslvector* vec, int k, double* lprob ) const
{
    mystring    preamb = "ProbVecOfDishHlp: ";
    mystring    errstr;
    char        bufstr[KBYTE];
    int         err;

//     if( GetBasin() == NULL )
//         throw myruntime_error( preamb + "Nil Basin." );
    if( GetMenu() == NULL )
        throw myruntime_error( preamb + "Nil Menu." );
    if( lprob == NULL )
        throw myruntime_error( preamb + "Memory access error." );
    if( vec == NULL )
        throw myruntime_error( preamb + "Invalid vector address." );

    *lprob = -1.e6;

    if( k < 0 || GetMenu()->GetSize() <= k )
        throw myruntime_error( preamb + "Invalid cluster index." );

    const Dish*         dishk = GetMenu()->GetDishAt( k );
    const Pslvector*    mu = GetMenu()->GetMuVectorAt( k );
    const SPDmatrix*    L = GetMenu()->GetSMatrixAt( k );
    const double        ldet = GetMenu()->GetLDetSMAt( k );

    if( dishk == NULL || mu == NULL || L == NULL )
        throw myruntime_error( preamb + "Memory access error." );

    const int           dim = GetMenu()->GetDim();
    const double        kp0 = GetMenu()->GetKappa0();
    const double        nu0 = GetMenu()->GetNu0();
    const int           nk = dishk->GetDishSize();
    const double        kp = kp0 + nk;
    const double        nu = nu0 + nk;

    if( dim <= 0 )
        throw myruntime_error( preamb + "Invalid dimensionality." );
//     if( nk <= 0 )
//         throw myruntime_error( preamb + "Cluster has no members." );
    if( kp <= 0.0 || nu <= 0.0/*( double )dim*/ )
        throw myruntime_error( preamb + "Cluster has too few members." );

    double  lnsum = 0.0;
    double  term1, term2;
    double  operr;
    const double    do2 = ( double)dim * 0.5;
    const double    no2 = nu * 0.5;
    const double    nup1o2 = no2 + 0.5;//(nu+1)/2
    const double    kp1 = kp + 1.0;

    const double    nu_t_o2 = nup1o2;//student's t nu over 2

//     const double    dofdom_o2 = nu_t_o2 - do2;
    const double    dofdom_o2 = NU_t_o2;
    const double    dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2
    const double    dofdom_pdim_m1_o2 = dofdom_pdim_o2 - 0.5;//(deg.of freedom-1 + dim.) over 2

    if(( err = psl_lngamma_e( dofdom_pdim_o2, &term1, &operr )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));
    if(( err = psl_lngamma_e( dofdom_o2, &term2, &operr )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    Pslvector   dev;
    SPDmatrix   MCM( dim );
    double      ldetMCM = 0.0;

    dev = *vec;
    if(( err = dev.Superposition( -1.0, *mu )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = MCM.Mul( dev, dev )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = MCM.Scale( kp / kp1 )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = MCM.Add( *L )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    //factorize by Cholesky decomp. and calculate determinant
    if(( err = MCM.CholeskyDecompose()) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = MCM.CDedLogDet( &ldetMCM )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));


    lnsum = term1 - term2;
    lnsum -= do2 * SLC_LNPI;
    lnsum += do2 * ( log( kp ) - log( kp1 ));
    lnsum += dofdom_pdim_m1_o2 * ldet;
    lnsum -= dofdom_pdim_o2 * ldetMCM;

    *lprob = lnsum;
//     *prob = exp( lnsum );
}





// =========================================================================
// ProbVecOfDishExc: compute probability of t-distributed vector to get 
//  dish k excluded that vector
//
void HDPbase::ProbVecOfDishExc( const Pslvector* vec, int k, double* lprob ) const
{
#ifdef PROBMTXINVSM
    ProbVecOfDishExcHlpObs( vec, k, lprob );
#else
    ProbVecOfDishExcHlp( vec, k, lprob );
#endif
}

#ifdef PROBMTXINVSM

// -------------------------------------------------------------------------
// ProbVecOfDishExcHlpObs: compute probability of t-distributed member 
//  vector to get dish k excluded the given member vector;
//  requires inverse of scale matrix
// NOTE: the member vector n is supposed to belong to dish k already
//
void HDPbase::ProbVecOfDishExcHlpObs( const Pslvector* vec, int k, double* lprob ) const
{
    mystring    preamb = "ProbVecOfDishExcHlpObs: ";
    mystring    errstr;
    char        bufstr[KBYTE];
    int         err;

//     if( GetBasin() == NULL )
//         throw myruntime_error( preamb + "Nil Basin." );
    if( GetMenu() == NULL )
        throw myruntime_error( preamb + "Nil Menu." );
    if( lprob == NULL )
        throw myruntime_error( preamb + "Memory access error." );
    if( vec == NULL )
        throw myruntime_error( preamb + "Invalid vector address." );

    *lprob = -1.e6;

    if( k < 0 || GetMenu()->GetSize() <= k )
        throw myruntime_error( preamb + "Invalid cluster index." );

    const Dish*         dishk = GetMenu()->GetDishAt( k );
    const Pslvector*    mu = GetMenu()->GetMuVectorAt( k );
    const SPDmatrix*    L = GetMenu()->GetSMatrixAt( k );
    const SPDmatrix*    Linv = GetMenu()->GetInvSMatrixAt( k );
    const double        ldet = GetMenu()->GetLDetSMAt( k );

    if( dishk == NULL || mu == NULL || L == NULL || Linv == NULL )
        throw myruntime_error( preamb + "Memory access error." );

    const int           dim = GetMenu()->GetDim();
    const double        kp0 = GetMenu()->GetKappa0();
    const double        nu0 = GetMenu()->GetNu0();
    const int           nk = dishk->GetDishSize();
    const double        kp = kp0 + nk;
    const double        nu = nu0 + nk;

    if( dim <= 0 )
        throw myruntime_error( preamb + "Invalid dimensionality." );
    if( nk <= 0 )
        throw myruntime_error( preamb + "Cluster has no members." );
    if( kp <= 0.0 || nu <= 0.0/*( double )dim*/ )
        throw myruntime_error( preamb + "Cluster has too few members." );

    if( nk == 1 ) {
        //dish has only one member; probability equals to prior
        PriorProbVec( vec, lprob );
        return;
    }

    double  lnsum = 0.0;
    double  term1, term2;
    double  operr;
    const double    do2 = ( double)dim * 0.5;
    const double    no2 = nu * 0.5;
    const double    num1o2 = no2 - 0.5;//(nu-1)/2
    const double    km1 = kp - 1.0;

    const double    nu_t_o2 = no2;//student's t nu over 2

//     const double    dofdom_o2 = nu_t_o2 - do2;
    const double    dofdom_o2 = NU_t_o2;
    const double    dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2
    const double    dofdom_pdim_m1_o2 = dofdom_pdim_o2 - 0.5;//(deg.of freedom-1 + dim.) over 2

    if(( err = psl_lngamma_e( dofdom_pdim_o2, &term1, &operr )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));
    if(( err = psl_lngamma_e( dofdom_o2, &term2, &operr )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    Pslvector   muexc;//mean vector of dish k excluded vec
    Pslvector   dev, scl;
    Pslmatrix   trv, pdt;//scale matrix of dish k excluded vec

    //mu_{n_k-1}
    muexc = *mu;
    if(( err = muexc.MultiplyBy( kp / km1 )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = muexc.Superposition( -1.0 / km1, *vec )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    //vec-mu_{n_k-1}
    dev = *vec;
    if(( err = dev.Superposition( -1.0, muexc )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    //|L_{n_k-1}|
    if(( err = dev.Transpose( trv )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = pdt.Mul( trv, *Linv )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = scl.Mul( pdt, dev )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if( scl.GetValueAt( 0 ) <= kp / km1 )
        throw myruntime_error( preamb + "Invalid determinant." );

    lnsum = term1 - term2;
    lnsum -= do2 * SLC_LNPI;
    lnsum += do2 * ( log( km1 ) - log( kp ));
    lnsum -= dofdom_pdim_o2 * log( 1.0 + km1 / kp * scl.GetValueAt( 0 ));
    lnsum -= 0.5 * ldet;

    *lprob = lnsum;
//     *prob = exp( lnsum );
}

#endif

// -------------------------------------------------------------------------
// ProbVecOfDishExcHlp: compute probability of t-distributed member 
//  vector to get dish k excluded the given member vector;
// NOTE: the member vector n is supposed to belong to dish k already
//
void HDPbase::ProbVecOfDishExcHlp( const Pslvector* vec, int k, double* lprob ) const
{
    mystring    preamb = "ProbVecOfDishExcHlp: ";
    mystring    errstr;
    char        bufstr[KBYTE];
    int         err;

//     if( GetBasin() == NULL )
//         throw myruntime_error( preamb + "Nil Basin." );
    if( GetMenu() == NULL )
        throw myruntime_error( preamb + "Nil Menu." );
    if( lprob == NULL )
        throw myruntime_error( preamb + "Memory access error." );
    if( vec == NULL )
        throw myruntime_error( preamb + "Invalid vector address." );

    *lprob = -1.e6;

    if( k < 0 || GetMenu()->GetSize() <= k )
        throw myruntime_error( preamb + "Invalid cluster index." );

    const Dish*         dishk = GetMenu()->GetDishAt( k );
    const Pslvector*    mu = GetMenu()->GetMuVectorAt( k );
    const SPDmatrix*    L = GetMenu()->GetSMatrixAt( k );
    const double        ldet = GetMenu()->GetLDetSMAt( k );

    if( dishk == NULL || mu == NULL || L == NULL )
        throw myruntime_error( preamb + "Memory access error." );

    const int           dim = GetMenu()->GetDim();
    const double        kp0 = GetMenu()->GetKappa0();
    const double        nu0 = GetMenu()->GetNu0();
    const int           nk = dishk->GetDishSize();
    const double        kp = kp0 + nk;
    const double        nu = nu0 + nk;

    if( dim <= 0 )
        throw myruntime_error( preamb + "Invalid dimensionality." );
    if( nk <= 0 )
        throw myruntime_error( preamb + "Cluster has no members." );
    if( kp <= 0.0 || nu <= 0.0/*( double )dim*/ )
        throw myruntime_error( preamb + "Cluster has too few members." );

    if( nk == 1 ) {
        //dish has only one member; probability equals to prior
        PriorProbVec( vec, lprob );
        return;
    }

    double  lnsum = 0.0;
    double  term1, term2;
    double  operr;
    const double    do2 = ( double)dim * 0.5;
    const double    no2 = nu * 0.5;
    const double    num1o2 = no2 - 0.5;//(nu-1)/2
    const double    km1 = kp - 1.0;

    const double    nu_t_o2 = no2;//student's t nu over 2

//     const double    dofdom_o2 = nu_t_o2 - do2;
    const double    dofdom_o2 = NU_t_o2;
    const double    dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2
    const double    dofdom_pdim_m1_o2 = dofdom_pdim_o2 - 0.5;//(deg.of freedom-1 + dim.) over 2

    if(( err = psl_lngamma_e( dofdom_pdim_o2, &term1, &operr )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));
    if(( err = psl_lngamma_e( dofdom_o2, &term2, &operr )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    Pslvector   muexc;//mean vector of dish k excluded vec
    Pslvector   dev;
    SPDmatrix   Lexc( dim );//scale matrix of dish k excluded vec
    double      Lexcldet = 0.0;

    //mu_{n_k-1}
    muexc = *mu;
    if(( err = muexc.MultiplyBy( kp / km1 )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = muexc.Superposition( -1.0 / km1, *vec )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    //vec-mu_{n_k-1}
    dev = *vec;
    if(( err = dev.Superposition( -1.0, muexc )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    //L_{n_k-1}
    if(( err = Lexc.Mul( dev, dev )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = Lexc.Scale( -km1 / kp )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = Lexc.Add( *L )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    //factor by Cholesky decomp. and calculate determinant
    if(( err = Lexc.CholeskyDecompose()) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = Lexc.CDedLogDet( &Lexcldet )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));


    lnsum = term1 - term2;
    lnsum -= do2 * SLC_LNPI;
    lnsum += do2 * ( log( km1 ) - log( kp ));
    lnsum += dofdom_pdim_m1_o2 * Lexcldet;
    lnsum -= dofdom_pdim_o2 * ldet;

    *lprob = lnsum;
//     *prob = exp( lnsum );
}
