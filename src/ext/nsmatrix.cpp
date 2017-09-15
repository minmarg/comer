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
#include "psl.h"
#include "pslcodes.h"
#include "nsmatrix.h"

// -------------------------------------------------------------------------
// constructor: initialization
//
NSmatrix::NSmatrix( int nrc )
:   Pslmatrix( nrc, nrc ),
    rowperm_( NULL ),
    szoperm_( 0 ),
    moneton_( 1 )
{
}

// -------------------------------------------------------------------------
// copy constructor
//
NSmatrix::NSmatrix( const NSmatrix& right )
:   Pslmatrix(),
    rowperm_( NULL ),
    szoperm_( 0 ),
    moneton_( 1 )
{
    *this = right;
}

// -------------------------------------------------------------------------
// constructor: default
//
NSmatrix::NSmatrix()
:   Pslmatrix(),
    rowperm_( NULL ),
    szoperm_( 0 ),
    moneton_( 1 )
{
}

// -------------------------------------------------------------------------
// destructor:
// -------------------------------------------------------------------------

NSmatrix::~NSmatrix()
{
    DestroyPermuts();
}

// -------------------------------------------------------------------------
// operator=: unmastered assignment
//
NSmatrix& NSmatrix::operator=( const NSmatrix& right )
{
    moneton_ = 1;
    DestroyPermuts();
    Pslmatrix::operator=( right );

    CopyPermuts( right );
    moneton_ = right.moneton_;
    return *this;
}

// -------------------------------------------------------------------------
// operator=: unmastered assignment
//
NSmatrix& NSmatrix::operator=( const Pslmatrix& right )
{
    if( right.GetNoRows() != right.GetNoCols()) {
        priverror( "Failed to assign: Not square matrix." );
        return *this;
    }
    moneton_ = 1;
    DestroyPermuts();
    Pslmatrix::operator=( right );
    return *this;
}

// -------------------------------------------------------------------------
// priverror: trivial error handler
// -------------------------------------------------------------------------

void NSmatrix::priverror( const char* errstr ) const
{
    fprintf( stderr, "ERROR: NSmatrix: %s\n", errstr );
    abort();
}

// -------------------------------------------------------------------------
// Singular: check the LU-decomposed matrix for singularity
//
bool NSmatrix::LUDedSingular() const
{
    int norows = GetNoRows();
    int n;
    for( n = 0; n < norows; n++ ) {
        if( GetValueAt( n, n ) == 0.0 )
            return true;
    }
    return false;
}

// =========================================================================
// LUDecompose: factorize the matrix into LU Decomposition: P A = LU. 
//  L and U are lower and upper triangular matrices, respectively, on output
//  overwritten with in this matrix, A. On output, row permutations are 
//  written in rowperm_ and sign (-1)^n in moneton_.
//
int NSmatrix::LUDecompose()
{
    const int   norows = GetNoRows();
    const int   nocols = GetNoCols();

    if( norows != nocols ) {
#ifdef NSmatrixTESTPRINT
        fprintf( stderr, "PSL_ERR_DIM: LUDecompose.\n");
#endif
        return PSL_ERR_DIM;
    }

    double  A_jj, A_ij, max;
    double  A_ik, A_jk;
    int     i, j, k;
    int     piv;//pivot row
    int     err;

    moneton_ = 1;
    InitPermuts( nocols );

    for( j = 0; j < nocols - 1; j++ ) {
        //find maximum in the jth column
        max = fabs( GetValueAt( j, j ));
        piv = j;

        for( i = j + 1; i < norows; i++ ) {
            A_ij = fabs( GetValueAt( i, j ));

            if( max < A_ij ) {
                max = A_ij;
                piv = i;
            }
        }

        if( piv != j ) {
            if(( err = SwapRows( j, piv )) != PSL_SUCCESS ) {
#ifdef NSmatrixTESTPRINT
                fprintf( stderr, "LUDecompose: Swap of rows failed.\n");
#endif
                return err;
            }
            if(( err = SwapPermPoss( j, piv )) != PSL_SUCCESS ) {
#ifdef NSmatrixTESTPRINT
                fprintf( stderr, "LUDecompose: Swap in permutation vector failed.\n");
#endif
                return err;
            }
            moneton_ = -moneton_;
        }

        A_jj = GetValueAt( j, j );

        if( A_jj == 0.0 ) {
#ifdef NSmatrixTESTPRINT
            fprintf( stderr, "PSL_ERR_ILLEGAL: LUDecompose: Division by 0.\n");
#endif
            return PSL_ERR_ILLEGAL;
        }

        for( i = j + 1; i < norows; i++ ) {
            A_ij = GetValueAt( i, j ) / A_jj;
            SetValueAt( i, j, A_ij );

            for( k = j + 1; k < nocols; k++ ) {
                A_ik = GetValueAt( i, k );
                A_jk = GetValueAt( j, k );
                SetValueAt( i, k, A_ik - A_ij * A_jk );
            }
        }
    }
    return PSL_SUCCESS;
}

// =========================================================================
// LUDedSolve: solve system of equations: Ax = b, where PA = LU 
//  decomposed into two matrices L and U written in the lower and upper
//  triangular parts of A, respectively. On input xb is b vector, on output
//  it holds the solution.
//  First, Ly = P^T b is solved for y; second, U x = y is solved for x.
//
int NSmatrix::LUDedSolve( Pslvector& xb ) const
{
    int status;
    if( GetNoRows() != GetNoCols()) {
#ifdef NSmatrixTESTPRINT
        fprintf( stderr, "PSL_ERR_DIM: LUDedSolve.\n");
#endif
        return PSL_ERR_DIM;
    }
    if( GetNoCols() != xb.GetSize()) {
#ifdef NSmatrixTESTPRINT
        fprintf( stderr, "PSL_ERR_DIM: LUDedSolve.\n");
#endif
        return PSL_ERR_DIM;
    }
    if( LUDedSingular()) {
#ifdef NSmatrixTESTPRINT
        fprintf( stderr, "PSL_ERR_ILLEGAL: LUDedSolve: Singular decomposition.\n");
#endif
        return PSL_ERR_ILLEGAL;
    }

    if(( status = PermuteVector( xb )) != PSL_SUCCESS ) {
#ifdef NSmatrixTESTPRINT
        fprintf( stderr, "LUDedSolve: Failed to permute vector.\n");
#endif
        return status;
    }

    if(( status = SolveLower( xb, true )) != PSL_SUCCESS ) {
#ifdef NSmatrixTESTPRINT
        fprintf( stderr, "LUDedSolve: %s\n", TranslatePSLError( status ));
#endif
        return status;
    }
    if(( status = SolveUpper( xb )) != PSL_SUCCESS ) {
#ifdef NSmatrixTESTPRINT
        fprintf( stderr, "LUDedSolve: %s\n", TranslatePSLError( status ));
#endif
        return status;
    }
    return PSL_SUCCESS;
}

// =========================================================================
// LUDedInvert: invert the matrix factorized by LU decomposition.
//  On input this matrix is assumed to have already decomposed: PA = L U.
//  Inversion is accomplished solving the system Ax = I, where I is the 
//  identity matrix.
//
int NSmatrix::LUDedInvert( Pslmatrix& inv ) const
{
    const int   norows = GetNoRows();
    const int   nocols = GetNoCols();

    if( norows != nocols ||
        inv.GetNoRows() != inv.GetNoCols() ||
        norows != inv.GetNoRows()) {
#ifdef NSmatrixTESTPRINT
        fprintf( stderr, "PSL_ERR_DIM: LUDedInvert.\n");
#endif
        return PSL_ERR_DIM;
    }
    if( LUDedSingular()) {
#ifdef NSmatrixTESTPRINT
        fprintf( stderr, "PSL_ERR_ILLEGAL: LUDedInvert: Singular decomposition.\n");
#endif
        return PSL_ERR_ILLEGAL;
    }

    int j, status;

    inv.SetIdentity();

    for( j = 0; j < nocols; j++ ) {
        Pslvector   cv = inv.ColVector( j );
        if(( status = LUDedSolve( cv )) != PSL_SUCCESS ) {
#ifdef NSmatrixTESTPRINT
            fprintf( stderr, "LUDedInvert: %s\n", TranslatePSLError( status ));
#endif
            return status;
        }
    }
    return PSL_SUCCESS;
}

// =========================================================================
// LUDedDet: calculate determinant of the matrix factorized by  LU
//  decomposition. Since the matrix is: PA = LU, and the diagonal 
//  elements of L are all 1, the determinant is equal to the product of the 
//  diagonal elements of U.
//
int NSmatrix::LUDedDet( double* det ) const
{
    if( !det )
        return 1;

    const int   norows = GetNoRows();
    const int   nocols = GetNoCols();
    double      res = moneton_;
    int         i;

    if( norows != nocols ) {
#ifdef NSmatrixTESTPRINT
        fprintf( stderr, "NSmatrix: LUDedDet: PSL_ERR_DIM.\n");
#endif
        return PSL_ERR_DIM;
    }
    for( i = 0; i < norows; i++ )
        res *= GetValueAt( i, i );

    *det = res;
    return PSL_SUCCESS;
}

// =========================================================================
// LUDedLogDet: calculate log determinant of the matrix factorized by  LU
//  decomposition. Since the matrix is: PA = LU, and the diagonal 
//  elements of L are all 1, the determinant is equal to the product of the 
//  diagonal elements of U.
//
int NSmatrix::LUDedLogDet( double* ldet ) const
{
    if( !ldet )
        return 1;

    const int   norows = GetNoRows();
    const int   nocols = GetNoCols();
    double      lres = 0.0;moneton_;
    double      val;
    int         i, nn = 0;

    if( moneton_ < 0 )
        nn = 1;

    if( norows != nocols ) {
#ifdef NSmatrixTESTPRINT
        fprintf( stderr, "NSmatrix: LUDedLogDet: PSL_ERR_DIM.\n");
#endif
        return PSL_ERR_DIM;
    }
    for( i = 0; i < norows; i++ ) {
        if( GetValueAt(i,i) < 0.0 )
            nn++;
        val = fabs( GetValueAt(i,i));
        if( val == 0.0 ) {
#ifdef SPDmatrixTESTPRINT
            fprintf( stderr, "NSmatrix: CDedLogDet: PSL_ERR_ILLEGAL.\n");
#endif
            return PSL_ERR_ILLEGAL;
        }
        lres += log( val );
    }
    *ldet = lres;
    if( nn & 1 ) {
#ifdef SPDmatrixTESTPRINT
        fprintf( stderr, "NSmatrix: CDedLogDet: PSL_ERR_ILLEGAL: Log of negative determinant.\n");
#endif
        return PSL_ERR_ILLEGAL;
    }
    return PSL_SUCCESS;
}
