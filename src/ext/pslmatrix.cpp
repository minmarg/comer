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
#include "pslvector.h"
#include "pslmatrix.h"

// -------------------------------------------------------------------------
// constructor: initialization
//
Pslmatrix::Pslmatrix( int nr, int nc )
:   values_( NULL ),
    norows_( 0 ),
    nocols_( 0 ),
    capacity_( 0 ),
    rowdim_( 0 ),
    master_( true )
{
    Reserve( nr, nc );
}

// -------------------------------------------------------------------------
// constructor: unmastered copy
//
Pslmatrix::Pslmatrix( const Pslmatrix& right )
:   values_( NULL ),
    norows_( 0 ),
    nocols_( 0 ),
    capacity_( 0 ),
    rowdim_( 0 ),
    master_( true )
{
    *this = right;
}

// -------------------------------------------------------------------------
// constructor: default
//
Pslmatrix::Pslmatrix()
:   values_( NULL ),
    norows_( 0 ),
    nocols_( 0 ),
    capacity_( 0 ),
    rowdim_( 0 ),
    master_( true )
{
}

// -------------------------------------------------------------------------
// destructor:
// -------------------------------------------------------------------------

Pslmatrix::~Pslmatrix()
{
    DestroyValues();
}

// -------------------------------------------------------------------------
// operator=: unmastered assignment
//
Pslmatrix& Pslmatrix::operator=( const Pslmatrix& right )
{
    DestroyValues();
    rowdim_ = 1;
    master_ = true;

    if( right.GetMaster()) {
        if( right.GetNoRows() < 1 || right.GetNoCols() < 1 )
            return *this;
        rowdim_ = right.GetPhysRowDim();
        Reserve( right.GetNoRows(), right.GetNoCols());
        Copy( right );
    }
    else {
        master_ = false;
        values_ = right.values_;
        norows_ = right.norows_;
        nocols_ = right.nocols_;
        rowdim_ = right.rowdim_;
    }
    return *this;
}

// -------------------------------------------------------------------------
// priverror: trivial error handler
// -------------------------------------------------------------------------

void Pslmatrix::priverror( const char* errstr ) const
{
    fprintf( stderr, "ERROR: Pslmatrix: %s.\n", errstr );
    abort();
}

// -------------------------------------------------------------------------
// Realloc: memory allocation
// -------------------------------------------------------------------------

void Pslmatrix::Realloc( int newcap )
{
    if( !GetMaster())
        priverror( "Only master is allowed to control memory." );

    double* tmp_values = NULL;

    if( newcap <= capacity_ )
        return;

    if( capacity_ == 0 ) {
        tmp_values = ( double* )malloc( sizeof( double ) * newcap );
    } else {
        tmp_values = ( double* )realloc( values_, sizeof( double ) * newcap );
    }

    if( !tmp_values )
        priverror( "Not enough memory" );

    values_ = tmp_values;

    // fill uninitialized memory with zeros
    memset( values_ + capacity_, 0, sizeof( double ) * ( newcap - capacity_ ));
    capacity_ = newcap;
}

// -------------------------------------------------------------------------
// Copy: copy elements of matrix to this matrix
// -------------------------------------------------------------------------

void Pslmatrix::Copy( const Pslmatrix& mtx )
{
    if( mtx.GetNoRows() < 1 || mtx.GetNoCols() < 1 )
        return;

    int noelms = GetNoRows() * GetNoCols();
    int mtxnoelems = mtx.GetNoRows() * mtx.GetNoCols();

    if( GetMaster() && GetCapacity() < mtxnoelems )
        return;

    if( GetNoRows() < mtx.GetNoRows() || GetNoCols() < mtx.GetNoCols()) {
        if( GetMaster()) {
            SetPhysRowDim( mtx.GetPhysRowDim());
            SetDims( mtx.GetNoRows(), mtx.GetNoCols());
        }
        else
            priverror( "Unable to copy: Not enough slave's space." );
    }

    int norows = GetNoRows();
    int nocols = GetNoCols();

    if( mtx.GetNoRows() < norows )
        norows = mtx.GetNoRows();
    if( mtx.GetNoCols() < nocols )
        nocols = mtx.GetNoCols();

    for( int n = 0; n < norows; n++ )
        for( int m = 0; m < nocols; m++ )
            SetValueAt( n, m, mtx.GetValueAt( n, m ));
}

// -------------------------------------------------------------------------
// SetIdentity: set this matrix equal to the identity matrix
//
void Pslmatrix::SetIdentity()
{
    int norows = GetNoRows();
    int nocols = GetNoCols();
    int i, j;
    for( i = 0; i < norows; i++ ) {
        for( j = 0; j < i && j < nocols; j++ )
            SetValueAt( i, j, 0.0 );
        if( j < nocols )
            SetValueAt( i, i, 1.0 );
        for( ++j; j < nocols; j++ )
            SetValueAt( i, j, 0.0 );
    }
}

// -------------------------------------------------------------------------
// SetCentering: set this matrix equal to the centering matrix
//
int Pslmatrix::SetCentering( double nrm )
{
    if( nrm <= 0.0 )
        return PSL_ERR_DOMAIN;
    int norows = GetNoRows();
    int nocols = GetNoCols();
    int i, j;
    double  val = -1.0 / nrm;
    double  dval = 1.0 + val;
    for( i = 0; i < norows; i++ ) {
        for( j = 0; j < i && j < nocols; j++ )
            SetValueAt( i, j, val );
        if( j < nocols )
            SetValueAt( i, i, dval );
        for( ++j; j < nocols; j++ )
            SetValueAt( i, j, val );
    }
    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// Zero: assign all elements to zero
// -------------------------------------------------------------------------

void Pslmatrix::Zero()
{
    for( int n = 0; n < GetNoRows(); n++ )
        for( int m = 0; m < GetNoCols(); m++ )
            SetValueAt( n, m, 0.0 );
}

// -------------------------------------------------------------------------
// Clear: clears all elements
// -------------------------------------------------------------------------

void Pslmatrix::Clear()
{
    if( !GetMaster())
        priverror( "Only master is allowed to clear data." );

    if( GetValues())
        Zero();
    SetDims( 0, 0 );
}

// -------------------------------------------------------------------------
// Print: print matrix to file
// -------------------------------------------------------------------------

void Pslmatrix::Print( FILE* fp, const char* format ) const
{
    if( fp == NULL )
        return;

    int n, m;
    if( format == NULL )
        format = " %g";

    for( n = 0; n < GetNoRows(); n++ ) {
        for( m = 0; m < GetNoCols(); m++ )
            fprintf( fp, format, GetValueAt( n, m ));

        fprintf( fp, "\n");
    }
}

// =========================================================================
// MATRIX OPERATORS
// =========================================================================
// Transpose: apply transpose operator in-place to this (square) matrix
//
int Pslmatrix::Transpose()
{
    int norows = GetNoRows();
    int nocols = GetNoCols();
    int     n, m;
    double  tmp;
    if( norows < 1 || norows != nocols )
        return PSL_ERR_DIM;

    for( n = 0; n < norows; n++ )
        for( m = n + 1; m < nocols; m++ ) {
            tmp = GetValueAt( n, m );
            SetValueAt( n, m, GetValueAt( m, n ));
            SetValueAt( m, n, tmp );
        }
    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// Transpose: apply transpose operator to the matrix
//
int Pslmatrix::Transpose( Pslmatrix& tr ) const
{
    int norows = GetNoRows();
    int nocols = GetNoCols();
    int n, m;
    if( norows < 1 || nocols < 1 )
        return PSL_ERR_DIM;

    tr.Reserve( nocols, norows );
    for( n = 0; n < nocols; n++ )
        for( m = 0; m < norows; m++ )
            tr.SetValueAt( n, m, GetValueAt( m, n ));

    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// Add: add matrix to this one
//
int Pslmatrix::Add( const Pslmatrix& mt )
{
    int norows = GetNoRows();
    int nocols = GetNoCols();
    int     n, m;
    double  val;
    if( norows != mt.GetNoRows() || nocols != mt.GetNoCols())
        return PSL_ERR_DIM;

    for( n = 0; n < norows; n++ )
        for( m = 0; m < nocols; m++ ) {
            val = mt.GetValueAt( n, m );
            if( val )
                AddValueAt( n, m, val );
        }
    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// Add: add constant double to each of the elements
//
int Pslmatrix::Add( double value )
{
    int norows = GetNoRows();
    int nocols = GetNoCols();
    int n, m;

    if( value == 0.0 )
        return PSL_SUCCESS;

    for( n = 0; n < norows; n++ )
        for( m = 0; m < nocols; m++ )
            AddValueAt( n, m, value );

    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// AddToDiag: add value to each element of the diagonal of the matrix
//
int Pslmatrix::AddToDiag( double value )
{
    int norows = GetNoRows();
    int nocols = GetNoCols();
    int i;
    for( i = 0; i < norows; i++ ) {
        if( nocols <= i )
            break;
        AddValueAt( i, i, value );
    }
    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// Sub: subtract matrix from this one
//
int Pslmatrix::Sub( const Pslmatrix& mt )
{
    int norows = GetNoRows();
    int nocols = GetNoCols();
    int n, m;
    if( norows != mt.GetNoRows() || nocols != mt.GetNoCols())
        return PSL_ERR_DIM;

    for( n = 0; n < norows; n++ )
        for( m = 0; m < nocols; m++ )
            AddValueAt( n, m, -mt.GetValueAt( n, m ));

    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// SubFrom: subtract this matrix from the one given by `mt'
//
int Pslmatrix::SubFrom( const Pslmatrix& mt )
{
    int norows = GetNoRows();
    int nocols = GetNoCols();
    int n, m;
    double  mtval, val;
    if( norows != mt.GetNoRows() || nocols != mt.GetNoCols())
        return PSL_ERR_DIM;

    for( n = 0; n < norows; n++ )
        for( m = 0; m < nocols; m++ ) {
            mtval = mt.GetValueAt( n, m );
            val = GetValueAt( n, m );
            SetValueAt( n, m, mtval - val );
        }

    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// SumColumns: sum all columns and write result to vector
//
int Pslmatrix::SumColumns( Pslvector& vv ) const
{
    int norows = GetNoRows();
    int nocols = GetNoCols();
    int n, m;

    vv.Reserve( norows );
    vv.Zero();

    for( m = 0; m < nocols; m++ ) {
        const Pslvector col = ColVector( m ) ;
        vv.Superposition( 1.0, col );
    }
    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// Scale: scale matrix by the given value
//
int Pslmatrix::Scale( double value )
{
    int norows = GetNoRows();
    int nocols = GetNoCols();
    int n, m;

    if( value == 1.0 )
        return PSL_SUCCESS;
    else if( value == -1.0 ) {
        for( n = 0; n < norows; n++ )
            for( m = 0; m < nocols; m++ )
                ChangeSignAt( n, m );
    }
    else if( value == 0.0 ) {
        for( n = 0; n < norows; n++ )
            for( m = 0; m < nocols; m++ )
                SetValueAt( n, m, 0.0 );
    }
    else {
        for( n = 0; n < norows; n++ )
            for( m = 0; m < nocols; m++ )
                MulValueAt( n, m, value );
    }
    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// Mul: mutliply two matrices and save the result
//
int Pslmatrix::Mul( const Pslmatrix& m1, const Pslmatrix& m2 )
{
    int m1rows = m1.GetNoRows();
    int m1cols = m1.GetNoCols();
    int m2rows = m2.GetNoRows();
    int m2cols = m2.GetNoCols();
    int     n, m, k;
    double  sum;
    if( m1cols != m2rows )
        return PSL_ERR_DIM;

    Reserve( m1rows, m2cols );
    for( n = 0; n < m1rows; n++ )
        for( m = 0; m < m2cols; m++ ) {
            sum = 0.0;
            for( k = 0; k < m1cols; k++ )
                sum += m1.GetValueAt( n, k ) * m2.GetValueAt( k, m );
            SetValueAt( n, m, sum );
        }
    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// Mul: multiply vector by transposed vector and save the result
//
int Pslmatrix::Mul( const Pslvector& v1, const Pslvector& v2 )
{
    int v1sz = v1.GetSize();
    int v2sz = v2.GetSize();
    int n, m;
    if( v1sz < 1 || v2sz < 1 )
        return PSL_SUCCESS;

    Reserve( v1sz, v2sz );
    for( n = 0; n < v1sz; n++ )
        for( m = 0; m < v2sz; m++ )
            SetValueAt( n, m, v1.GetValueAt( n ) * v2.GetValueAt( m ));
    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// KronProduct: Kronecker product of two matrices; result is saved in this 
//  matrix
//
int Pslmatrix::KronProduct( const Pslmatrix& mm1, const Pslmatrix& mm2 )
{
    int m1rows = mm1.GetNoRows();
    int m1cols = mm1.GetNoCols();
    int m2rows = mm2.GetNoRows();
    int m2cols = mm2.GetNoCols();
    int myrows = m1rows * m2rows;
    int mycols = m1cols * m2cols;
    int     n1, m1, n2, m2, n, m;
    double  val1, val2;

    Reserve( myrows, mycols );
    for( n1 = 0; n1 < m1rows; n1++ )
        for( m1 = 0; m1 < m1cols; m1++ ) {
            n = n1 * m2rows;
            m = m1 * m2cols;
            val1 = mm1.GetValueAt( n1, m1 );
            for( n2 = 0; n2 < m2rows; n2++ )
                for( m2 = 0; m2 < m2cols; m2++ ) {
                    val2 = mm2.GetValueAt( n2, m2 );
                    SetValueAt( n + n2, m + m2, val1 * val2 );
                }
        }
    return PSL_SUCCESS;
}

// =========================================================================
// SOLVE methods
// =========================================================================
// SolveLower: solve lower triangular part of the matrix assuming the upper
//  part being 0: Ax = b. On input xb is assumed to be equal to b, on output
//  xb contains the solution
//
int Pslmatrix::SolveLower( Pslvector& xb, bool unitdiag ) const
{
    if( GetNoRows() != GetNoCols()) {
        return PSL_ERR_DIM;
    }
    if( GetNoCols() != xb.GetSize()) {
        return PSL_ERR_DIM;
    }

    int     noelems = xb.GetSize();
    double  A_ii, A_ij, tmp;
    int     i, j;

    for( i = 0; i < noelems; i++ ) {
        tmp = xb.GetValueAt( i );
        for( j = 0; j < i; j++ ) {
            A_ij = GetValueAt( i, j );
            tmp -= A_ij * xb.GetValueAt( j );
        }
        if( !unitdiag ) {
            A_ii = GetValueAt( i, i );
            if( A_ii == 0.0 )
                return PSL_ERR_ILLEGAL;
            tmp /= A_ii;
        }
        xb.SetValueAt( i, tmp );
    }
    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// SolveUpper: solve upper triangular part of the matrix assuming the lower
//  part being 0: Ax = b. On input xb is assumed to be equal to b, on output
//  xb contains the solution
//
int Pslmatrix::SolveUpper( Pslvector& xb, bool unitdiag ) const
{
    if( GetNoRows() != GetNoCols()) {
        return PSL_ERR_DIM;
    }
    if( GetNoCols() != xb.GetSize()) {
        return PSL_ERR_DIM;
    }

    int     noelems = xb.GetSize();
    double  A_ii, A_ij, tmp;
    int     i, j;

    for( i = noelems - 1; 0 <= i; i-- ) {
        tmp = xb.GetValueAt( i );
        for( j = i + 1; j < noelems; j++ ) {
            A_ij = GetValueAt( i, j );
            tmp -= A_ij * xb.GetValueAt( j );
        }
        if( !unitdiag ) {
            A_ii = GetValueAt( i, i );
            if( A_ii == 0.0 )
                return PSL_ERR_ILLEGAL;
            tmp /= A_ii;
        }
        xb.SetValueAt( i, tmp );
    }
    return PSL_SUCCESS;
}

// =========================================================================
// MultiplyLower: multiply lower triangular part of the matrix by vector b
//  assuming the upper part of the matrix being 0: Ab = x. On input bx is 
//  assumed to be equal to b, on output bx contains the result
//
int Pslmatrix::MultiplyLower( Pslvector& bx ) const
{
    if( GetNoRows() != GetNoCols()) {
        return PSL_ERR_DIM;
    }
    if( GetNoCols() != bx.GetSize()) {
        return PSL_ERR_DIM;
    }

    int     noelems = bx.GetSize();
    double  A_ij, tmp;
    int     i, j;

    for( i = noelems - 1; 0 <= i; i-- ) {
        tmp = 0.0;
        for( j = 0; j <= i; j++ ) {
            A_ij = GetValueAt( i, j );
            tmp += A_ij * bx.GetValueAt( j );
        }
        bx.SetValueAt( i, tmp );
    }
    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// MultiplyUpper: multiply upper triangular part of the matrix by vector b
//  assuming the lower part of the matrix being 0: Ab = x. On input bx is 
//  assumed to be equal to b, on output bx contains the result
//
int Pslmatrix::MultiplyUpper( Pslvector& bx ) const
{
    if( GetNoRows() != GetNoCols()) {
        return PSL_ERR_DIM;
    }
    if( GetNoCols() != bx.GetSize()) {
        return PSL_ERR_DIM;
    }

    int     noelems = bx.GetSize();
    double  A_ij, tmp;
    int     i, j;

    for( i = 0; i < noelems; i++ ) {
        tmp = 0.0;
        for( j = i; j < noelems; j++ ) {
            A_ij = GetValueAt( i, j );
            tmp += A_ij * bx.GetValueAt( j );
        }
        bx.SetValueAt( i, tmp );
    }
    return PSL_SUCCESS;
}

// =========================================================================
// SWAP operations
// =========================================================================
// SwapRows: swap rows
//
int Pslmatrix::SwapRows( int i1, int i2 )
{
    double  tmp;
    int norows = GetNoRows();
    int nocols = GetNoCols();
    int n;

    if( norows <= i1 || norows <= i2 )
        return PSL_ERR_ADDRESS;
    if( i1 == i2 )
        return PSL_SUCCESS;

    for( n = 0; n < nocols; n++ ) {
        tmp = GetValueAt( i1, n );
        SetValueAt( i1, n, GetValueAt( i2, n ));
        SetValueAt( i2, n, tmp );
    }
    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// SwapCols: swap columns
//
int Pslmatrix::SwapCols( int j1, int j2 )
{
    double  tmp;
    int norows = GetNoRows();
    int nocols = GetNoCols();
    int n;

    if( nocols <= j1 || nocols <= j2 )
        return PSL_ERR_ADDRESS;
    if( j1 == j2 )
        return PSL_SUCCESS;

    for( n = 0; n < norows; n++ ) {
        tmp = GetValueAt( n, j1 );
        SetValueAt( n, j1, GetValueAt( n, j2 ));
        SetValueAt( n, j2, tmp );
    }
    return PSL_SUCCESS;
}

// =========================================================================
// SUB- vector, matrix operations
// =========================================================================
// Submatrix: make submatrix of this matrix
//
const Pslmatrix Pslmatrix::SubMatrix( int i, int j, int nr, int nc ) const
{
    Pslmatrix   sub;

    if( i < 0 || GetNoRows() <= i )
        priverror( "Submatrix row index is out of range." );
    if( j < 0 || GetNoCols() <= j )
        priverror( "Submatrix column index is out of range." );
    if( nr < 1 || GetNoRows() < i + nr )
        priverror( "Submatrix number of rows overflows." );
    if( nc < 1 || GetNoCols() < j + nc )
        priverror( "Submatrix number of columns overflows." );

    sub.values_ = values_ + i*rowdim_ + j;
    sub.norows_ = nr;
    sub.nocols_ = nc;
    sub.rowdim_ = rowdim_;
    sub.master_ = false;
    return sub;
}

// =========================================================================
// Stack: stack rows of matrix into single vector (private matrix 
//  representation)
//
const Pslvector Pslmatrix::Stack() const
{
    Pslvector   stack;

    stack.SetVector( values_, norows_ * nocols_ );
    stack.SetStride( 1 );
    return stack;
}

// -------------------------------------------------------------------------
// RowVector: make new vector from the matrix row
//
const Pslvector Pslmatrix::RowVector( int i ) const
{
    Pslvector   row;

    if( i < 0 || GetNoRows() <= i )
        priverror( "RowVector row index is out of range." );

    row.SetVector( values_ + i*rowdim_, nocols_ );
    row.SetStride( 1 );
    return row;
}

// -------------------------------------------------------------------------
// SubRowVector: make new vector from the matrix subrow elements
//
const Pslvector Pslmatrix::SubRowVector( int i, int offset, int n ) const
{
    Pslvector   subrow;

    if( i < 0 || GetNoRows() <= i )
        priverror( "SubRowVector row index is out of range." );
    if( n < 1 || offset < 0 )
        priverror( "SubRowVector number of elements is invalid." );
    if( GetNoCols() < offset + n )
        priverror( "SubRowVector number of elements overflows." );

    subrow.SetVector( values_ + i*rowdim_ + offset, n );
    subrow.SetStride( 1 );
    return subrow;
}

// =========================================================================
// ColVector: make new vector from the matrix column
//
const Pslvector Pslmatrix::ColVector( int j ) const
{
    Pslvector   col;

    if( j < 0 || GetNoCols() <= j )
        priverror( "ColVector column index is out of range." );

    col.SetVector( values_ + j, norows_ );
    col.SetStride( rowdim_ );
    return col;
}

// -------------------------------------------------------------------------
// SubColVector: make new vector from the matrix subcolumn elements
//
const Pslvector Pslmatrix::SubColVector( int j, int offset, int n ) const
{
    Pslvector   subcol;

    if( j < 0 || GetNoCols() <= j )
        priverror( "SubColVector column index is out of range." );
    if( n < 1 || offset < 0 )
        priverror( "SubColVector number of elements is invalid." );
    if( GetNoRows() < offset + n )
        priverror( "SubColVector number of elements overflows." );

    subcol.SetVector( values_ + offset*rowdim_ + j, n );
    subcol.SetStride( rowdim_ );
    return subcol;
}

// =========================================================================
// DiagVector: make new vector from the matrix diagonal elements
//
const Pslvector Pslmatrix::DiagVector() const
{
    Pslvector   diag;

    diag.SetVector( values_, SLC_MIN( norows_, nocols_ ));
    diag.SetStride( rowdim_ + 1 );
    return diag;
}

// -------------------------------------------------------------------------
// SubDiagVector: make new vector from the matrix subdiagonal elements
//
const Pslvector Pslmatrix::SubDiagVector( int k ) const
{
    Pslvector   subdiag;

    if( k < 0 || GetNoRows() <= k )
        priverror( "SubDiagVector index is out of range." );

    subdiag.SetVector( values_ + k*rowdim_, SLC_MIN( norows_ - k, nocols_ ));
    subdiag.SetStride( rowdim_ + 1 );
    return subdiag;
}

// -------------------------------------------------------------------------
// SuperDiagVector: make new vector from the upper matrix subdiagonal 
//     elements
//
const Pslvector Pslmatrix::SuperDiagVector( int k ) const
{
    Pslvector   superdiag;

    if( k < 0 || GetNoCols() <= k )
        priverror( "SuperDiagVector index is out of range." );

    superdiag.SetVector( values_ + k, SLC_MIN( norows_, nocols_ - k ));
    superdiag.SetStride( rowdim_ + 1 );
    return superdiag;
}
