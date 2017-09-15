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


// PRIOR PROBABILITY SECTION
// =========================================================================
// PriorProbMtx: computation of prior probability of t-distributed 
//     matrix (set of vectors: table)
//
void HDPbase::PriorProbMtx( const Table* tbl, double* lprob ) const
{
    if( GetMenu() == NULL )
        throw myruntime_error( "PriorProbMtx: Null Menu." );
    if( tbl == NULL )
        throw myruntime_error( "PriorProbMtx: Invalid matrix address." );

    const int       dd = GetMenu()->GetDim();
    const int       vv = tbl->GetActualSize();//number of vectors at the table
    const double    dim = ( double )dd;
    const double    nov = ( double )vv;
    const double    nov2 = SQUARE( nov );
    const double    dim1pvv = dim * ( 1.0 + 2.5 * nov );

    if( nov2 < dim1pvv ) {
        Pslmatrix   mtx;
        Table2Mtx( tbl, &mtx );
        PriorProbMtx( &mtx, lprob );//comp.: //A*L^2+A^2*L + A^3/3+A
    }
    else
        PriorProbMtxHlp( tbl, lprob );//comp.: A^2*L+A^2+(approx.:2.5*A*L) + A^3/3+A
}

// -------------------------------------------------------------------------
// PriorProbMtx: compute prior probability of t-distributed matrix
//
void HDPbase::PriorProbMtx( const Pslmatrix* mtx, double* lprob ) const
{
    if( GetMenu() == NULL )
        throw myruntime_error( "PriorProbMtx: Null Menu." );
    if( mtx == NULL )
        throw myruntime_error( "PriorProbMtx: Invalid matrix address." );

    const int       dd = GetMenu()->GetDim();
    const int       vv = mtx->GetNoCols();//number of vectors in matrix
    const double    dim = ( double )dd;
    const double    nov = ( double )vv;
    const double    nov313 = ( double )CUBE( nov ) * 1.33 + nov;
    const double    dim303 = ( double )CUBE( dim ) * 0.33 + dim;

    //NOTE: calculations of LxL matrices fixed
    //NOTE: LU decomposition should be used since a matrix may be not >0
    //NOTE: calculation complexity by this approach becomes unfavourable
    if( 0&&nov313 < dim303 )
        PriorProbMtxHlpObs( mtx, lprob );//comp.: A^2*L+A*L^2+L^3 + L^3/3+L
    else
        PriorProbMtxHlp( mtx, lprob );//comp.: //A*L^2+A^2*L + A^3/3+A
}


// -------------------------------------------------------------------------
// PriorProbMtxHlpObs: compute prior probability of t-distributed matrix
//  4-matrix mult. version; requires inverse of scale matrix
//  matrix mul. complexity: A^2*L+A*L^2+L^3 (L<=A); A*L^2+A^2*L+A*L^2 (L>A)
//  complexity of determinant calc.: L^3/3+L
//
void HDPbase::PriorProbMtxHlpObs( const Pslmatrix* mtx, double* lprob ) const
{
    mystring    preamb = "PriorProbMtxHlpObs: ";
    mystring    errstr;
    char        bufstr[KBYTE];
    int         err;

//     if( GetBasin() == NULL )
//         throw myruntime_error( preamb + "Null Basin." );
    if( GetMenu() == NULL )
        throw myruntime_error( preamb + "Null Menu." );
    if( lprob == NULL )
        throw myruntime_error( preamb + "Memory access error." );
    if( mtx == NULL )
        throw myruntime_error( preamb + "Invalid matrix address." );

    *lprob = -1.e6;

    const Pslvector*    mu0 = GetMenu()->GetMu0();
    const SPDmatrix*    L0 = GetMenu()->GetS0();
    const SPDmatrix*    L0inv = GetMenu()->GetInvS0();
    const double        L0ldet = GetMenu()->GetLDetS0();

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
    double  term1, term2;
    double  operr;
    double  j;
    const int       nov = mtx->GetNoCols();//number of vectors in matrix
    const double    novv = ( double )nov;
    const double    do2 = ( double)dim * 0.5;
    const double    nu0o2 = nu0 * 0.5;
    const double    vvo2 = ( double)nov * 0.5;
    const double    nu0p1o2 = nu0o2 + 0.5;//(nu0+1)/2
    const double    nu0pvo2 = nu0o2 + vvo2;//(nu0+L)/2
    const double    nu0pv1o2 = nu0pvo2 + 0.5;//(nu0+L+1)/2
    const double    kp0pv = kp0 + novv;

    const double    nu_t_o2 = nu0pvo2;//student t's nu over 2

//     const double    dofdom_o2 = nu_t_o2 - do2;
    const double    dofdom_o2 = NU_t_o2;
    const double    dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2
    const double    dofdom_pdim1_o2 = dofdom_pdim_o2 + 0.5;//(deg.of freedom+1 + dim.) over 2
    const double    dofdom_pdim1_mv_o2 = dofdom_pdim1_o2 - vvo2;//(deg.of freedom+1 + dim.- L) / 2

    Pslvector   unit( nov );
    Pslmatrix   C( nov, nov );//centering matrix
    Pslmatrix   m01TmH;
    Pslmatrix   SSM;//scaled scatter matrix
    Pslmatrix   MCM;//combined matrices
    double      ldetMCM = 0.0;

    if(( err = C.SetCentering( kp0pv )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    unit.SetUnity();

    if(( err = m01TmH.Mul( *mu0, unit )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = m01TmH.Sub( *mtx )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = m01TmH.Transpose( SSM )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if( dim < nov ) {//NOTE: favourable chain multiplication complexity
        if(( err = MCM.Mul( C, SSM )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = SSM.Mul( MCM, *L0inv )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = MCM.Mul( SSM, m01TmH )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));
    }
    else {//nov <= dim
        if(( err = MCM.Mul( *L0inv, m01TmH )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = m01TmH.Mul( SSM, MCM )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = MCM.Mul( C, m01TmH )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));
    }
    if(( err = MCM.AddIdentity()) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if( MCM.GetNoRows() <= 2 ) {
        if( MCM.GetNoRows() == 1 )
            ldetMCM = MCM.GetValueAt( 0,0 );
        else if( MCM.GetNoRows() == 2 )
            ldetMCM = MCM.GetValueAt( 0,0 ) * MCM.GetValueAt( 1,1 ) - 
                      MCM.GetValueAt( 1,0 ) * MCM.GetValueAt( 0,1 );
        if( ldetMCM <= 0.0 )
            throw myruntime_error( preamb + "Invalid determinant." );
        ldetMCM = log( ldetMCM );
    }
    else {
        //factorize MCM by LU-decomp. and calculate determinant of MCM
        NSmatrix    MCMLU(nov);//MCM LU-decomposed
        MCMLU = MCM;
        if(( err = MCMLU.LUDecompose()) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = MCMLU.LUDedLogDet( &ldetMCM )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));
    }

    //calculate probability next
    if( nov == 1 ) {
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
    lnsum -= novv * do2 * SLC_LNPI;
    lnsum += do2 * ( log( kp0 ) - log( kp0pv ));
    lnsum -= vvo2 * L0ldet;
    lnsum -= dofdom_pdim_o2 * ldetMCM;

    *lprob = lnsum;
//     *prob = exp( lnsum );
}

// -------------------------------------------------------------------------
// PriorProbMtxHlp: compute prior probability of t-distributed matrix;
//  3-matrix mult. version
//  matrix mul. complexity: A*L^2+A^2*L
//  complexity of determinant calc.: A^3/3+A
//
void HDPbase::PriorProbMtxHlp( const Pslmatrix* mtx, double* lprob ) const
{
    mystring    preamb = "PriorProbMtxHlp: ";
    mystring    errstr;
    char        bufstr[KBYTE];
    int         err;

//     if( GetBasin() == NULL )
//         throw myruntime_error( preamb + "Null Basin." );
    if( GetMenu() == NULL )
        throw myruntime_error( preamb + "Null Menu." );
    if( lprob == NULL )
        throw myruntime_error( preamb + "Memory access error." );
    if( mtx == NULL )
        throw myruntime_error( preamb + "Invalid matrix address." );

    *lprob = -1.e6;

    const Pslvector*    mu0 = GetMenu()->GetMu0();
    const SPDmatrix*    L0 = GetMenu()->GetS0();
    const double        L0ldet = GetMenu()->GetLDetS0();

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
    double  term1, term2;
    double  operr;
    double  j;
    const int       nov = mtx->GetNoCols();//number of vectors in matrix
    const double    novv = ( double )nov;
    const double    do2 = ( double)dim * 0.5;
    const double    nu0o2 = nu0 * 0.5;
    const double    vvo2 = ( double)nov * 0.5;
    const double    nu0p1o2 = nu0o2 + 0.5;//(nu0+1)/2
    const double    nu0pvo2 = nu0o2 + vvo2;//(nu0+L)/2
    const double    nu0pv1o2 = nu0pvo2 + 0.5;//(nu0+L+1)/2
    const double    kp0pv = kp0 + novv;

    const double    nu_t_o2 = nu0pvo2;//student's t nu over 2

//     const double    dofdom_o2 = nu_t_o2 - do2;
    const double    dofdom_o2 = NU_t_o2;
    const double    dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2
    const double    dofdom_pdim1_o2 = dofdom_pdim_o2 + 0.5;//(deg.of freedom+1 + dim.) over 2
    const double    dofdom_pdim1_mv_o2 = dofdom_pdim1_o2 - vvo2;//(deg.of freedom+1 + dim.- L) / 2
    const double    dofdom_pdim_mv_o2 = dofdom_pdim_o2 - vvo2;//(deg.of freedom + dim.- L) over 2

    Pslvector   unit( nov );
    Pslmatrix   C( nov, nov );//centering matrix
    Pslmatrix   m01TmHT;
    Pslmatrix   SSM;//scaled scatter matrix
    Pslmatrix   MCM;//combined matrices
    double      ldetMCM = 0.0;

    if(( err = C.SetCentering( kp0pv )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    unit.SetUnity();

    if(( err = MCM.Mul( *mu0, unit )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = MCM.Sub( *mtx )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = MCM.Transpose( m01TmHT )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = SSM.Mul( MCM, C )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = MCM.Mul( SSM, m01TmHT )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = MCM.Add( *L0 )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    //factorize MCM by LU-decomp. and calculate determinant of MCM
    if(( err = ( reinterpret_cast<SPDmatrix&>( MCM )).CholeskyDecompose()) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = ( reinterpret_cast<SPDmatrix&>( MCM )).CDedLogDet( &ldetMCM )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));


    //calculate probability next
    if( nov == 1 ) {
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
    lnsum -= novv * do2 * SLC_LNPI;
    lnsum += do2 * ( log( kp0 ) - log( kp0pv ));
    lnsum += dofdom_pdim_mv_o2 * L0ldet;
    lnsum -= dofdom_pdim_o2 * ldetMCM;

    *lprob = lnsum;
//     *prob = exp( lnsum );
}


// -------------------------------------------------------------------------
// PriorProbMtxHlp: alternative computation of prior probability of 
//     t-distributed matrix (set of vectors: table)
//
void HDPbase::PriorProbMtxHlp( const Table* tbl, double* lprob ) const
{
    mystring    preamb = "PriorProbMtxHlp: ";
    mystring    errstr;
    char        bufstr[KBYTE];
    int         err, n;

//     if( GetBasin() == NULL )
//         throw myruntime_error( preamb + "Null Basin." );
    if( GetMenu() == NULL )
        throw myruntime_error( preamb + "Null Menu." );
    if( lprob == NULL )
        throw myruntime_error( preamb + "Memory access error." );
    if( tbl == NULL )
        throw myruntime_error( preamb + "Invalid matrix address." );

    *lprob = -1.e6;

    const Pslvector*    mu0 = GetMenu()->GetMu0();
    const SPDmatrix*    L0 = GetMenu()->GetS0();
    const double        L0ldet = GetMenu()->GetLDetS0();

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
    double  term1, term2;
    double  operr;
    double  j;
    const int       nov = tbl->GetActualSize();//number of vectors in matrix
    const double    novv = ( double )nov;
    const double    do2 = ( double)dim * 0.5;
    const double    nu0o2 = nu0 * 0.5;
    const double    vvo2 = ( double)nov * 0.5;
    const double    nu0p1o2 = nu0o2 + 0.5;//(nu0+1)/2
    const double    nu0pvo2 = nu0o2 + vvo2;//(nu0+L)/2
    const double    nu0pv1o2 = nu0pvo2 + 0.5;//(nu0+L+1)/2
    const double    kp0pv = kp0 + novv;

    const double    nu_t_o2 = nu0pvo2;//student's t nu over 2

//     const double    dofdom_o2 = nu_t_o2 - do2;
    const double    dofdom_o2 = NU_t_o2;
    const double    dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2
    const double    dofdom_pdim1_o2 = dofdom_pdim_o2 + 0.5;//(deg.of freedom+1 + dim.) over 2
    const double    dofdom_pdim1_mv_o2 = dofdom_pdim1_o2 - vvo2;//(deg.of freedom+1 + dim.- L) / 2
    const double    dofdom_pdim_mv_o2 = dofdom_pdim_o2 - vvo2;//(deg.of freedom + dim.- L) over 2

    const Pslvector*    vec;
    Pslvector   mnew( dim );//mean of table
    Pslvector   dev;
    SPDmatrix   LSL( dim );//lambda_{L}
    SPDmatrix   Sdev( dim );
    double      ldetLSL = 0.0;

    for( n = 0; n < tbl->GetSize(); n++ )
    {
        if( tbl->GetVectorNIndAt( n ) < 0 )
            continue;
        vec = tbl->GetVectorNAt( n );
        if( vec == NULL )
            continue;
        if(( err = mnew.Superposition( 1.0, *vec )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        dev = *vec;
        if(( err = dev.Superposition( -1.0, *mu0 )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = Sdev.Mul( dev, dev )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = LSL.Add( Sdev )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));
    }

    dev = mnew;
    if(( err = dev.Superposition( -novv, *mu0 )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = Sdev.Mul( dev, dev )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = Sdev.Scale( -1.0 / kp0pv )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = LSL.Add( Sdev )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = LSL.Add( *L0 )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    //factorize LSL by Cholesky decomp. (LU-decomp.) and calculate determinant of LSL
    if(( err = LSL.CholeskyDecompose()) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = LSL.CDedLogDet( &ldetLSL )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));


    //calculate probability next
    if( nov == 1 ) {
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
    lnsum -= novv * do2 * SLC_LNPI;
    lnsum += do2 * ( log( kp0 ) - log( kp0pv ));
    lnsum += dofdom_pdim_mv_o2 * L0ldet;
    lnsum -= dofdom_pdim_o2 * ldetLSL;

    *lprob = lnsum;
//     *prob = exp( lnsum );
}





// POSTERIOR PREDICITVE SECTION
// =========================================================================
// ProbMtxOfDish: computation of probability of t-distributed 
//     matrix (set of vectors: table) to get dish k
//
void HDPbase::ProbMtxOfDish( const Table* tbl, int k, double* lprob ) const
{
    if( GetMenu() == NULL )
        throw myruntime_error( "ProbMtxOfDish: Null Menu." );
    if( tbl == NULL )
        throw myruntime_error( "ProbMtxOfDish: Invalid matrix address." );

    const int       dd = GetMenu()->GetDim();
    const int       vv = tbl->GetActualSize();//number of vectors at the table
    const double    dim = ( double )dd;
    const double    nov = ( double )vv;
    const double    nov2 = SQUARE( nov );
    const double    dim1pvv = dim * ( 1.0 + 2.5 * nov );

    if( nov2 < dim1pvv ) {
        Pslmatrix   mtx;
        Table2Mtx( tbl, &mtx );
        ProbMtxOfDish( &mtx, k, lprob );//comp.: //A*L^2+A^2*L + A^3/3+A
    }
    else
        ProbMtxOfDishHlp( tbl, k, lprob );//comp.: A^2*L+A^2+(approx.:2.5*A*L) + A^3/3+A
}

// -------------------------------------------------------------------------
// ProbMtxOfDish: compute probability of t-distributed matrix to get 
//  dish k
//
void HDPbase::ProbMtxOfDish( const Pslmatrix* mtx, int k, double* lprob ) const
{
#ifdef PROBMTXINVSM
    if( GetMenu() == NULL )
        throw myruntime_error( "ProbMtxOfDish: Null Menu." );
    if( mtx == NULL )
        throw myruntime_error( "ProbMtxOfDish: Invalid matrix address." );

    const int       dd = GetMenu()->GetDim();
    const int       vv = mtx->GetNoCols();//number of vectors in matrix
    const double    dim = ( double )dd;
    const double    nov = ( double )vv;
    const double    nov313 = ( double )CUBE( nov ) * 1.33 + nov;
    const double    dim303 = ( double )CUBE( dim ) * 0.33 + dim;

    //NOTE: LU decomposition should be used since a matrix may be not >0
    //NOTE: calculation complexity by this approach becomes unfavourable
    if( 0&&nov313 < dim303 )
        ProbMtxOfDishHlpObs( mtx, k, lprob );//comp.: A^2*L+A*L^2+L^3 + L^3/3+L
    else
#endif
        ProbMtxOfDishHlp( mtx, k, lprob );//comp.: //A*L^2+A^2*L + A^3/3+A
}


#ifdef PROBMTXINVSM

// -------------------------------------------------------------------------
// ProbMtxOfDishHlpObs: compute probability of t-distributed matrix to get 
//  dish k; 4-matrix mult. version; requires inverse of scale matrix
//  matrix mul. complexity: A^2*L+A*L^2+L^3 (L<=A); A*L^2+A^2*L+A*L^2 (L>A)
//  complexity of determinant calc.: L^3/3+L
//
void HDPbase::ProbMtxOfDishHlpObs( const Pslmatrix* mtx, int k, double* lprob ) const
{
    mystring    preamb = "ProbMtxOfDishHlpObs: ";
    mystring    errstr;
    char        bufstr[KBYTE];
    int         err;

    if( GetBasin() == NULL )
        throw myruntime_error( preamb + "Null Basin." );
    if( GetMenu() == NULL )
        throw myruntime_error( preamb + "Null Menu." );
    if( lprob == NULL )
        throw myruntime_error( preamb + "Memory access error." );
    if( mtx == NULL )
        throw myruntime_error( preamb + "Invalid matrix address." );

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
    double  j;
    const int       nov = mtx->GetNoCols();//number of vectors in matrix
    const double    novv = ( double )nov;
    const double    do2 = ( double)dim * 0.5;
    const double    nuo2 = nu * 0.5;
    const double    vvo2 = ( double)nov * 0.5;
    const double    nup1o2 = nuo2 + 0.5;//(nu+1)/2
    const double    nupvo2 = nuo2 + vvo2;//(nu+L)/2
    const double    nupv1o2 = nupvo2 + 0.5;//(nu+L+1)/2
    const double    kppv = kp + novv;

    const double    nu_t_o2 = nupvo2;//student's t nu over 2

//     const double    dofdom_o2 = nu_t_o2 - do2;
    const double    dofdom_o2 = NU_t_o2;
    const double    dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2
    const double    dofdom_pdim1_o2 = dofdom_pdim_o2 + 0.5;//(deg.of freedom+1 + dim.) over 2
    const double    dofdom_pdim1_mv_o2 = dofdom_pdim1_o2 - vvo2;//(deg.of freedom+1 + dim.- L) / 2

    Pslvector   unit( nov );
    Pslmatrix   C( nov, nov );//centering matrix
    Pslmatrix   mu1TmH;
    Pslmatrix   SSM;//scaled scatter matrix
    Pslmatrix   MCM;//combined matrices
    double      ldetMCM = 0.0;

    if(( err = C.SetCentering( kppv )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    unit.SetUnity();
//     unit.SetAllToValue( -1.0 );

    if(( err = mu1TmH.Mul( *mu, unit )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = mu1TmH.Sub( *mtx )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = mu1TmH.Transpose( SSM )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if( dim < nov ) {//NOTE: favourable chain multiplication complexity
        if(( err = MCM.Mul( C, SSM )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = SSM.Mul( MCM, *Linv )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = MCM.Mul( SSM, mu1TmH )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));
    }
    else {//nov <= dim
        if(( err = MCM.Mul( *Linv, mu1TmH )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = mu1TmH.Mul( SSM, MCM )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = MCM.Mul( C, mu1TmH )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));
    }
    if(( err = MCM.AddIdentity()) != 0 )
        throw myruntime_error( TranslatePSLError( err ));


    if( MCM.GetNoRows() <= 2 ) {
        if( MCM.GetNoRows() == 1 )
            ldetMCM = MCM.GetValueAt( 0,0 );
        else if( MCM.GetNoRows() == 2 )
            ldetMCM = MCM.GetValueAt( 0,0 ) * MCM.GetValueAt( 1,1 ) - 
                      MCM.GetValueAt( 1,0 ) * MCM.GetValueAt( 0,1 );
        if( ldetMCM <= 0.0 )
            throw myruntime_error( preamb + "Invalid determinant." );
        ldetMCM = log( ldetMCM );
    }
    else {
        //factorize MCM by LU-decomp. and calculate determinant of MCM
        NSmatrix    MCMLU(nov);//MCM LU-decomposed
        MCMLU = MCM;
        if(( err = MCMLU.LUDecompose()) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = MCMLU.LUDedLogDet( &ldetMCM )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));
    }

    //calculate probability next
    if( nov == 1 ) {
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
    //NOTE:INVALID!
//     for( j = 0.5; j <= vvo2; j += 0.5 )
//     {
//         if(( err = psl_lngamma_e( nupv1o2 - j, &term1, &operr )) != 0 )
//             throw myruntime_error( TranslatePSLError( err ));
//         if(( err = psl_lngamma_e( nupv1o2 - do2 - j, &term2, &operr )) != 0 )
//             throw myruntime_error( TranslatePSLError( err ));
// 
//         lnsum += term1 - term2;
//     }
    lnsum -= novv * do2 * SLC_LNPI;
    lnsum += do2 * ( log( kp ) - log( kppv ));
    lnsum -= vvo2 * ldet;
    lnsum -= dofdom_pdim_o2 * ldetMCM;

    *lprob = lnsum;
//     *prob = exp( lnsum );
}

#endif//PROBMTXINVSM

// -------------------------------------------------------------------------
// ProbMtxOfDishHlp: compute probability of t-distributed matrix to get 
//  dish k; 3-matrix mult. version
//  matrix mul. complexity: A*L^2+A^2*L
//  complexity of determinant calc.: A^3/3+A
//
void HDPbase::ProbMtxOfDishHlp( const Pslmatrix* mtx, int k, double* lprob ) const
{
    mystring    preamb = "ProbMtxOfDishHlp: ";
    mystring    errstr;
    char        bufstr[KBYTE];
    int         err;

    if( GetBasin() == NULL )
        throw myruntime_error( preamb + "Null Basin." );
    if( GetMenu() == NULL )
        throw myruntime_error( preamb + "Null Menu." );
    if( lprob == NULL )
        throw myruntime_error( preamb + "Memory access error." );
    if( mtx == NULL )
        throw myruntime_error( preamb + "Invalid matrix address." );

    *lprob = -1.e6;

    if( k < 0 || GetMenu()->GetSize() <= k )
        throw myruntime_error( preamb + "Invalid cluster index." );

    const Dish*         dishk = GetMenu()->GetDishAt( k );
    const Pslvector*    mu = GetMenu()->GetMuVectorAt( k );
    const SPDmatrix*    L = GetMenu()->GetSMatrixAt( k );
    const double        ldet = GetMenu()->GetLDetSMAt( k );

    if( dishk == NULL || 
        mu == NULL || L == NULL )
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
    double  j;
    const int       nov = mtx->GetNoCols();//number of vectors in matrix
    const double    novv = ( double )nov;
    const double    do2 = ( double)dim * 0.5;
    const double    nuo2 = nu * 0.5;
    const double    vvo2 = ( double )nov * 0.5;
    const double    nup1o2 = nuo2 + 0.5;//(nu+1)/2
    const double    nupvo2 = nuo2 + vvo2;//(nu+L)/2
    const double    nupv1o2 = nupvo2 + 0.5;//(nu+L+1)/2
    const double    kppv = kp + novv;

    const double    nu_t_o2 = nupvo2;//student's t nu over 2

//     const double    dofdom_o2 = nu_t_o2 - do2;
    const double    dofdom_o2 = NU_t_o2;
    const double    dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2
    const double    dofdom_pdim1_o2 = dofdom_pdim_o2 + 0.5;//(deg.of freedom+1 + dim.) over 2
    const double    dofdom_pdim1_mv_o2 = dofdom_pdim1_o2 - vvo2;//(deg.of freedom+1 + dim.- L) / 2
    const double    dofdom_pdim_mv_o2 = dofdom_pdim_o2 - vvo2;//(deg.of freedom + dim.- L) over 2

    Pslvector   unit( nov );
    Pslmatrix   C( nov, nov );//centering matrix
    Pslmatrix   mu1TmHT;
    Pslmatrix   SSM;//scaled scatter matrix
    Pslmatrix   MCM;//combined matrices
    double      ldetMCM = 0.0;

    if(( err = C.SetCentering( kppv )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    unit.SetUnity();

    if(( err = MCM.Mul( *mu, unit )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = MCM.Sub( *mtx )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = MCM.Transpose( mu1TmHT )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = SSM.Mul( MCM, C )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = MCM.Mul( SSM, mu1TmHT )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = MCM.Add( *L )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    //factorize MCM by LU-decomp. and calculate determinant of MCM
    if(( err = ( reinterpret_cast<SPDmatrix&>( MCM )).CholeskyDecompose()) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = ( reinterpret_cast<SPDmatrix&>( MCM )).CDedLogDet( &ldetMCM )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));


    //calculate probability next
    if( nov == 1 ) {
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
    lnsum -= novv * do2 * SLC_LNPI;
    lnsum += do2 * ( log( kp ) - log( kppv ));
    lnsum += dofdom_pdim_mv_o2 * ldet;
    lnsum -= dofdom_pdim_o2 * ldetMCM;

    *lprob = lnsum;
//     *prob = exp( lnsum );
}


// -------------------------------------------------------------------------
// ProbMtxOfDishHlp: alternative computation of probability of t-distributed 
//     matrix (set of vectors: table) to get dish k
//
void HDPbase::ProbMtxOfDishHlp( const Table* tbl, int k, double* lprob ) const
{
    mystring    preamb = "ProbMtxOfDishHlp: ";
    mystring    errstr;
    char        bufstr[KBYTE];
    int         err, n;

    if( GetBasin() == NULL )
        throw myruntime_error( preamb + "Null Basin." );
    if( GetMenu() == NULL )
        throw myruntime_error( preamb + "Null Menu." );
    if( lprob == NULL )
        throw myruntime_error( preamb + "Memory access error." );
    if( tbl == NULL )
        throw myruntime_error( preamb + "Invalid matrix address." );

    *lprob = -1.e6;

    if( k < 0 || GetMenu()->GetSize() <= k )
        throw myruntime_error( preamb + "Invalid cluster index." );

    const Dish*         dishk = GetMenu()->GetDishAt( k );
    const Pslvector*    mu = GetMenu()->GetMuVectorAt( k );
    const SPDmatrix*    L = GetMenu()->GetSMatrixAt( k );
    const double        ldet = GetMenu()->GetLDetSMAt( k );

    if( dishk == NULL || 
        mu == NULL || L == NULL )
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
    double  j;
    const int       nov = tbl->GetActualSize();//number of vectors in matrix
    const double    novv = ( double )nov;
    const double    do2 = ( double)dim * 0.5;
    const double    nuo2 = nu * 0.5;
    const double    vvo2 = ( double )nov * 0.5;
    const double    nup1o2 = nuo2 + 0.5;//(nu+1)/2
    const double    nupvo2 = nuo2 + vvo2;//(nu+L)/2
    const double    nupv1o2 = nupvo2 + 0.5;//(nu+L+1)/2
    const double    kppv = kp + novv;

    const double    nu_t_o2 = nupvo2;//student's t nu over 2

//     const double    dofdom_o2 = nu_t_o2 - do2;
    const double    dofdom_o2 = NU_t_o2;
    const double    dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2
    const double    dofdom_pdim1_o2 = dofdom_pdim_o2 + 0.5;//(deg.of freedom+1 + dim.) over 2
    const double    dofdom_pdim1_mv_o2 = dofdom_pdim1_o2 - vvo2;//(deg.of freedom+1 + dim.- L) / 2
    const double    dofdom_pdim_mv_o2 = dofdom_pdim_o2 - vvo2;//(deg.of freedom + dim.- L) over 2

    const Pslvector*    vec;
    Pslvector   mnew( dim );//mean of table
    Pslvector   dev;
    SPDmatrix   LSL( dim );//lambda_{n_k+L}
    SPDmatrix   Sdev( dim );
    double      ldetLSL = 0.0;

    for( n = 0; n < tbl->GetSize(); n++ )
    {
        if( tbl->GetVectorNIndAt( n ) < 0 )
            continue;
        vec = tbl->GetVectorNAt( n );
        if( vec == NULL )
            continue;
        if(( err = mnew.Superposition( 1.0, *vec )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        dev = *vec;
        if(( err = dev.Superposition( -1.0, *mu )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = Sdev.Mul( dev, dev )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = LSL.Add( Sdev )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));
    }

    dev = mnew;
    if(( err = dev.Superposition( -novv, *mu )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = Sdev.Mul( dev, dev )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = Sdev.Scale( -1.0 / kppv )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = LSL.Add( Sdev )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = LSL.Add( *L )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    //factorize LSL by Cholesky decomp. (LU-decomp.) and calculate determinant of LSL
    if(( err = LSL.CholeskyDecompose()) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = LSL.CDedLogDet( &ldetLSL )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));


    //calculate probability next
    if( nov == 1 ) {
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
    lnsum -= novv * do2 * SLC_LNPI;
    lnsum += do2 * ( log( kp ) - log( kppv ));
    lnsum += dofdom_pdim_mv_o2 * ldet;
    lnsum -= dofdom_pdim_o2 * ldetLSL;

    *lprob = lnsum;
//     *prob = exp( lnsum );
}





// EXCLUSIVE POSTERIOR PREDICITVE SECTION
// =========================================================================
// ProbMtxOfDishExc: computation of exclusive probability of t-distributed 
//     matrix (set of vectors: table) to get dish k
//
void HDPbase::ProbMtxOfDishExc( const Table* tbl, int k, double* lprob ) const
{
    if( GetMenu() == NULL )
        throw myruntime_error( "ProbMtxOfDishExc: Null Menu." );
    if( tbl == NULL )
        throw myruntime_error( "ProbMtxOfDishExc: Invalid matrix address." );

    const int       dd = GetMenu()->GetDim();
    const int       vv = tbl->GetActualSize();//number of vectors at the table
    const double    dim = ( double )dd;
    const double    nov = ( double )vv;
    const double    nov2 = SQUARE( nov );
    //auxiliary sum over columns compensate for 2.5 factor in compl. expr.
    const double    dim1pvv = dim * ( 1.0 + 1.0 * nov );

    if( nov2 < dim1pvv ) {
        Pslmatrix   mtx;
        Table2Mtx( tbl, &mtx );
        ProbMtxOfDishExc( &mtx, k, lprob );//comp.: //A*L^2+A^2*L + A^3/3+A
    }
    else
        ProbMtxOfDishExcHlp( tbl, k, lprob );//comp.: A^2*L+A^2+A*L + A^3/3+A
}

// -------------------------------------------------------------------------
// ProbMtxOfDishExc: compute exclusive probability of t-distributed 
//  matrix to get dish k
//
void HDPbase::ProbMtxOfDishExc( const Pslmatrix* mtx, int k, double* lprob ) const
{
#ifdef PROBMTXINVSM
    if( GetMenu() == NULL )
        throw myruntime_error( "ProbMtxOfDishExc: Null Menu." );
    if( mtx == NULL )
        throw myruntime_error( "ProbMtxOfDishExc: Invalid matrix address." );

    const int       dd = GetMenu()->GetDim();
    const int       vv = mtx->GetNoCols();//number of vectors in matrix
    const double    dim = ( double )dd;
    const double    nov = ( double )vv;
    const double    nov313 = ( double )CUBE( nov ) * 1.33 + nov;
    const double    dim303 = ( double )CUBE( dim ) * 0.33 + dim;

    //NOTE: LU decomposition should be used since a matrix may be not >0
    //NOTE: calculation complexity by this approach becomes unfavourable
    if( 0&&nov313 < dim303 )
        ProbMtxOfDishExcHlpObs( mtx, k, lprob );//comp.: A^2*L+A*L^2+L^3 + L^3/3+L
    else
#endif
        ProbMtxOfDishExcHlp( mtx, k, lprob );//comp.: //A*L^2+A^2*L + A^3/3+A
}


#ifdef PROBMTXINVSM

// -------------------------------------------------------------------------
// ProbMtxOfDishExcHlpObs: compute exclusive probability of t-distributed 
//  matrix to get dish k
// NOTE: the member matrix mtx is supposed to belong to dish k already
//  4-matrix mult. version; requires inverse of scale matrix
//  matrix mul. complexity: A^2*L+A*L^2+L^3 (L<=A); A*L^2+A^2*L+A*L^2 (L>A)
//  complexity of determinant calc.: L^3/3+L
//
void HDPbase::ProbMtxOfDishExcHlpObs( const Pslmatrix* mtx, int k, double* lprob ) const
{
    mystring    preamb = "ProbMtxOfDishExcHlpObs: ";
    mystring    errstr;
    char        bufstr[KBYTE];
    int         err;

    if( GetBasin() == NULL )
        throw myruntime_error( preamb + "Null Basin." );
    if( GetMenu() == NULL )
        throw myruntime_error( preamb + "Null Menu." );
    if( lprob == NULL )
        throw myruntime_error( preamb + "Memory access error." );
    if( mtx == NULL )
        throw myruntime_error( preamb + "Invalid matrix address." );

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
    const int           nov = mtx->GetNoCols();//number of vectors in matrix

    if( dim <= 0 )
        throw myruntime_error( preamb + "Invalid dimensionality." );
    if( nk <= 0 )
        throw myruntime_error( preamb + "Cluster has no members." );
    if( kp <= 0.0 || nu <= 0.0/*( double )dim*/ )
        throw myruntime_error( preamb + "Cluster has too few members." );
    if( nk < nov )
        throw myruntime_error( preamb + "Invalid matrix size." );
    if( nk == nov ) {
        //matrix size coincides with the dish size; predictive is equal to prior
        PriorProbMtx( mtx, lprob );
        return;
    }

    double  lnsum = 0.0;
    double  term1, term2;
    double  operr;
    double  j;
    const double    novv = ( double )nov;
    const double    do2 = ( double)dim * 0.5;
    const double    nuo2 = nu * 0.5;
    const double    vvo2 = ( double)nov * 0.5;
    const double    nup1o2 = nuo2 + 0.5;//(nu+1)/2
    const double    nupvo2 = nuo2 - vvo2;//(nu-L)/2
    const double    nupv1o2 = nupvo2 + 0.5;//(nu-L+1)/2
    const double    kppv = kp - novv;

    const double    nu_t_o2 = nuo2;//student's t nu over 2

//     const double    dofdom_o2 = nu_t_o2 - do2;
    const double    dofdom_o2 = NU_t_o2;
    const double    dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2
    const double    dofdom_pdim1_o2 = dofdom_pdim_o2 + 0.5;//(deg.of freedom+1 + dim.) over 2
    const double    dofdom_pdim1_mv_o2 = dofdom_pdim1_o2 - vvo2;//(deg.of freedom+1 + dim.- L) / 2
    const double    dofdom_pdim_mv_o2 = dofdom_pdim_o2 - vvo2;//(deg.of freedom + dim.- L) over 2

    Pslvector   muexc;
    Pslvector   unit( nov );
    Pslmatrix   C( nov, nov );//centering matrix
    Pslmatrix   mu1TmH;
    Pslmatrix   SSM;//scaled scatter matrix
    Pslmatrix   MCM;//combined matrices
    double      ldetMCM = 0.0;

    if(( err = C.SetCentering( kp )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    //minus because Lambda_{n_k-L} is expressed via Lambda_{n_k}
    if(( err = C.Scale( -1.0 )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    unit.SetUnity();

    //mu_{n_k-L}
    if(( err = mtx->SumColumns( muexc )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = muexc.Superposition( -kp, *mu )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = muexc.MultiplyBy( -1.0 / kppv )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    //mu_{n_k-L} 1^T - H
    if(( err = mu1TmH.Mul( muexc, unit )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = mu1TmH.Sub( *mtx )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = mu1TmH.Transpose( SSM )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    //- {mu_{n_k-L} 1^T - H}^T Lambda^{-1}_{n_k} {mu_{n_k-L} 1^T - H} C
    if( dim < nov ) {//NOTE: favourable chain multiplication complexity
        if(( err = MCM.Mul( C, SSM )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = SSM.Mul( MCM, *Linv )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = MCM.Mul( SSM, mu1TmH )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));
    }
    else {//nov <= dim
        if(( err = MCM.Mul( *Linv, mu1TmH )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = mu1TmH.Mul( SSM, MCM )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = MCM.Mul( C, mu1TmH )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));
    }
    if(( err = MCM.AddIdentity()) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if( MCM.GetNoRows() <= 2 ) {
        if( MCM.GetNoRows() == 1 )
            ldetMCM = MCM.GetValueAt( 0,0 );
        else if( MCM.GetNoRows() == 2 )
            ldetMCM = MCM.GetValueAt( 0,0 ) * MCM.GetValueAt( 1,1 ) - 
                      MCM.GetValueAt( 1,0 ) * MCM.GetValueAt( 0,1 );
        if( ldetMCM <= 0.0 )
            throw myruntime_error( preamb + "Invalid determinant." );
        ldetMCM = log( ldetMCM );
    }
    else {
        //factorize MCM by LU-decomp. and calculate determinant of MCM
        NSmatrix    MCMLU(nov);//MCM LU-decomposed
        MCMLU = MCM;
        if(( err = MCMLU.LUDecompose()) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = MCMLU.LUDedLogDet( &ldetMCM )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));
    }

    //calculate probability next
    if( nov == 1 ) {
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
    lnsum -= novv * do2 * SLC_LNPI;
    lnsum += do2 * ( log( kppv ) - log( kp ));
    lnsum += dofdom_pdim_mv_o2 * ldetMCM;
    lnsum -= vvo2 * ldet;

    *lprob = lnsum;
//     *prob = exp( lnsum );
}

#endif//PROBMTXINVSM

// -------------------------------------------------------------------------
// ProbMtxOfDishExcHlp: compute exclusive probability of t-distributed 
//  matrix to get dish k
// NOTE: the member matrix mtx is supposed to belong to dish k already
//  3-matrix mult. version
//  matrix mul. complexity: A*L^2+A^2*L
//  complexity of determinant calc.: A^3/3+A
//
void HDPbase::ProbMtxOfDishExcHlp( const Pslmatrix* mtx, int k, double* lprob ) const
{
    mystring    preamb = "ProbMtxOfDishExcHlp: ";
    mystring    errstr;
    char        bufstr[KBYTE];
    int         err;

    if( GetBasin() == NULL )
        throw myruntime_error( preamb + "Null Basin." );
    if( GetMenu() == NULL )
        throw myruntime_error( preamb + "Null Menu." );
    if( lprob == NULL )
        throw myruntime_error( preamb + "Memory access error." );
    if( mtx == NULL )
        throw myruntime_error( preamb + "Invalid matrix address." );

    *lprob = -1.e6;

    if( k < 0 || GetMenu()->GetSize() <= k )
        throw myruntime_error( preamb + "Invalid cluster index." );

    const Dish*         dishk = GetMenu()->GetDishAt( k );
    const Pslvector*    mu = GetMenu()->GetMuVectorAt( k );
    const SPDmatrix*    L = GetMenu()->GetSMatrixAt( k );
    const double        ldet = GetMenu()->GetLDetSMAt( k );

    if( dishk == NULL || 
        mu == NULL || L == NULL )
        throw myruntime_error( preamb + "Memory access error." );

    const int           dim = GetMenu()->GetDim();
    const double        kp0 = GetMenu()->GetKappa0();
    const double        nu0 = GetMenu()->GetNu0();
    const int           nk = dishk->GetDishSize();
    const double        kp = kp0 + nk;
    const double        nu = nu0 + nk;
    const int           nov = mtx->GetNoCols();//number of vectors in matrix

    if( dim <= 0 )
        throw myruntime_error( preamb + "Invalid dimensionality." );
    if( nk <= 0 )
        throw myruntime_error( preamb + "Cluster has no members." );
    if( kp <= 0.0 || nu <= 0.0/*( double )dim*/ )
        throw myruntime_error( preamb + "Cluster has too few members." );
    if( nk < nov )
        throw myruntime_error( preamb + "Invalid matrix size." );
    if( nk == nov ) {
        //matrix size coincides with the dish size; predictive is equal to prior
        PriorProbMtx( mtx, lprob );
        return;
    }

    double  lnsum = 0.0;
    double  term1, term2;
    double  operr;
    double  j;
    const double    novv = ( double )nov;
    const double    do2 = ( double)dim * 0.5;
    const double    nuo2 = nu * 0.5;
    const double    vvo2 = ( double )nov * 0.5;
    const double    nup1o2 = nuo2 + 0.5;//(nu+1)/2
    const double    nupvo2 = nuo2 - vvo2;//(nu-L)/2
    const double    nupv1o2 = nupvo2 + 0.5;//(nu-L+1)/2
    const double    kppv = kp - novv;

    const double    nu_t_o2 = nuo2;//student's t nu over 2

//     const double    dofdom_o2 = nu_t_o2 - do2;
    const double    dofdom_o2 = NU_t_o2;
    const double    dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2
    const double    dofdom_pdim1_o2 = dofdom_pdim_o2 + 0.5;//(deg.of freedom+1 + dim.) over 2
    const double    dofdom_pdim1_mv_o2 = dofdom_pdim1_o2 - vvo2;//(deg.of freedom+1 + dim.- L) / 2
    const double    dofdom_pdim_mv_o2 = dofdom_pdim_o2 - vvo2;//(deg.of freedom + dim.- L) over 2

    Pslvector   muexc;//mtx-excluded mean vector of dish k
    Pslvector   unit( nov );
    Pslmatrix   C( nov, nov );//centering matrix
    Pslmatrix   mu1TmHT;
    Pslmatrix   SSM;//scaled scatter matrix
    Pslmatrix   MCM;//combined matrices
    double      ldetMCM = 0.0;

    if(( err = C.SetCentering( kp )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    unit.SetUnity();

    //mu_{n_k-L}
    if(( err = mtx->SumColumns( muexc )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = muexc.Superposition( -kp, *mu )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = muexc.MultiplyBy( -1.0 / kppv )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    //mu_{n_k-L} 1^T - H
    if(( err = MCM.Mul( muexc, unit )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = MCM.Sub( *mtx )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = MCM.Transpose( mu1TmHT )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    //{mu_{n_k-L} 1^T - H} C {mu_{n_k-L} 1^T - H}^T
    if(( err = SSM.Mul( MCM, C )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = MCM.Mul( SSM, mu1TmHT )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = MCM.SubFrom( *L )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    //factorize MCM by LU-decomp. and calculate determinant of MCM
    if(( err = ( reinterpret_cast<SPDmatrix&>( MCM )).CholeskyDecompose()) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = ( reinterpret_cast<SPDmatrix&>( MCM )).CDedLogDet( &ldetMCM )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));


    //calculate probability next
    if( nov == 1 ) {
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
    lnsum -= novv * do2 * SLC_LNPI;
    lnsum += do2 * ( log( kppv ) - log( kp ));
    lnsum += dofdom_pdim_mv_o2 * ldetMCM;
    lnsum -= dofdom_pdim_o2 * ldet;

    *lprob = lnsum;
//     *prob = exp( lnsum );
}

// -------------------------------------------------------------------------
// ProbMtxOfDishExcHlp: alternative computation of exclusive probability of 
//     t-distributed matrix (set of vectors: table) to get dish k
// NOTE: the table `tbl' is supposed to belong to dish k already
//
void HDPbase::ProbMtxOfDishExcHlp( const Table* tbl, int k, double* lprob ) const
{
    mystring    preamb = "ProbMtxOfDishExcHlp: ";
    mystring    errstr;
    char        bufstr[KBYTE];
    int         err, n;

    if( GetBasin() == NULL )
        throw myruntime_error( preamb + "Null Basin." );
    if( GetMenu() == NULL )
        throw myruntime_error( preamb + "Null Menu." );
    if( lprob == NULL )
        throw myruntime_error( preamb + "Memory access error." );
    if( tbl == NULL )
        throw myruntime_error( preamb + "Invalid matrix address." );

    *lprob = -1.e6;

    if( k < 0 || GetMenu()->GetSize() <= k )
        throw myruntime_error( preamb + "Invalid cluster index." );

    const Dish*         dishk = GetMenu()->GetDishAt( k );
    const Pslvector*    mu = GetMenu()->GetMuVectorAt( k );
    const SPDmatrix*    L = GetMenu()->GetSMatrixAt( k );
    const double        ldet = GetMenu()->GetLDetSMAt( k );

    if( dishk == NULL || 
        mu == NULL || L == NULL )
        throw myruntime_error( preamb + "Memory access error." );

    const int           dim = GetMenu()->GetDim();
    const double        kp0 = GetMenu()->GetKappa0();
    const double        nu0 = GetMenu()->GetNu0();
    const int           nk = dishk->GetDishSize();
    const double        kp = kp0 + nk;
    const double        nu = nu0 + nk;
    const int           nov = tbl->GetActualSize();//number of vectors in matrix

    if( dim <= 0 )
        throw myruntime_error( preamb + "Invalid dimensionality." );
    if( nk <= 0 )
        throw myruntime_error( preamb + "Cluster has no members." );
    if( kp <= 0.0 || nu <= 0.0/*( double )dim*/ )
        throw myruntime_error( preamb + "Cluster has too few members." );
    if( nk < nov )
        throw myruntime_error( preamb + "Invalid matrix size." );
    if( nk == nov ) {
        //table size coincides with the dish size; predictive is equal to prior
        PriorProbMtx( tbl, lprob );
        return;
    }

    double  lnsum = 0.0;
    double  term1, term2;
    double  operr;
    double  j;
    const double    novv = ( double )nov;
    const double    do2 = ( double)dim * 0.5;
    const double    nuo2 = nu * 0.5;
    const double    vvo2 = ( double )nov * 0.5;
    const double    nup1o2 = nuo2 + 0.5;//(nu+1)/2
    const double    nupvo2 = nuo2 - vvo2;//(nu-L)/2
    const double    nupv1o2 = nupvo2 + 0.5;//(nu-L+1)/2
    const double    kppv = kp - novv;

    const double    nu_t_o2 = nuo2;//student's t nu over 2

//     const double    dofdom_o2 = nu_t_o2 - do2;
    const double    dofdom_o2 = NU_t_o2;
    const double    dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2
    const double    dofdom_pdim1_o2 = dofdom_pdim_o2 + 0.5;//(deg.of freedom+1 + dim.) over 2
    const double    dofdom_pdim1_mv_o2 = dofdom_pdim1_o2 - vvo2;//(deg.of freedom+1 + dim.- L) / 2
    const double    dofdom_pdim_mv_o2 = dofdom_pdim_o2 - vvo2;//(deg.of freedom + dim.- L) over 2

    const Pslvector*    vec;
    Pslvector   mnew( dim );//mean of table
    Pslvector   muexc;//tbl-excluded mean vector of dish k
    Pslvector   dev;
    SPDmatrix   LSL( dim );//lambda_{n_k-L}
    SPDmatrix   Sdev( dim );
    double      ldetLSL = 0.0;

    //mu_{n_k-L}
    for( n = 0; n < tbl->GetSize(); n++ )
    {
        if( tbl->GetVectorNIndAt( n ) < 0 )
            continue;
        vec = tbl->GetVectorNAt( n );
        if( vec == NULL )
            continue;
        if(( err = mnew.Superposition( 1.0, *vec )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));
    }
    muexc = mnew;
    if(( err = muexc.Superposition( -kp, *mu )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = muexc.MultiplyBy( -1.0 / kppv )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    //-SUM (eta - mu_{n_k-L}) (eta - mu_{n_k-L})^T
    for( n = 0; n < tbl->GetSize(); n++ )
    {
        if( tbl->GetVectorNIndAt( n ) < 0 )
            continue;
        vec = tbl->GetVectorNAt( n );
        if( vec == NULL )
            continue;
        dev = *vec;
        if(( err = dev.Superposition( -1.0, muexc )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = Sdev.Mul( dev, dev )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = LSL.Sub( Sdev )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));
    }

    //(bar{eta} - mu_{n_k-L}) (bar{eta} - mu_{n_k-L})^T
    dev = mnew;
    if(( err = dev.Superposition( -novv, muexc )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = Sdev.Mul( dev, dev )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = Sdev.Scale( 1.0 / kp )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    ////lambda_{n_k-L}
    if(( err = LSL.Add( Sdev )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = LSL.Add( *L )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    //factorize LSL by Cholesky decomp. (LU-decomp.) and calculate determinant of LSL
    if(( err = LSL.CholeskyDecompose()) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    if(( err = LSL.CDedLogDet( &ldetLSL )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));


    //calculate probability next
    if( nov == 1 ) {
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
    lnsum -= novv * do2 * SLC_LNPI;
    lnsum += do2 * ( log( kppv ) - log( kp ));
    lnsum += dofdom_pdim_mv_o2 * ldetLSL;
    lnsum -= dofdom_pdim_o2 * ldet;

    *lprob = lnsum;
//     *prob = exp( lnsum );
}
