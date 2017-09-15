/***************************************************************************
 *   Copyright (C) 2010 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __nsmatrix__
#define __nsmatrix__

#include <stdlib.h>
#include <string.h>
#include "pslcodes.h"
#include "pslmatrix.h"

#ifndef NSmatrixTESTPRINT
#define NSmatrixTESTPRINT
#endif

// -------------------------------------------------------------------------
// Invertibale square matrix
//
class NSmatrix: public Pslmatrix
{
public:
    NSmatrix( int nrc );
    NSmatrix( const NSmatrix& );
    virtual ~NSmatrix();

    virtual NSmatrix&   operator=( const NSmatrix& );
    virtual NSmatrix&   operator=( const Pslmatrix& );
    virtual void        priverror( const char* ) const;

    int     LUDecompose();
    int     LUDedSolve( Pslvector& xb ) const;
    int     LUDedInvert( Pslmatrix& inv ) const;
    int     LUDedDet( double* ) const;
    int     LUDedLogDet( double* ) const;

protected:
    explicit    NSmatrix();

    void    InitPermuts( int size );
    void    CopyPermuts( const NSmatrix& );
    void    DestroyPermuts();
    int     GetPermutsAt( int pos ) const;
    int     SwapPermPoss( int i1, int i2 );
    int     PermuteVector( Pslvector& v ) const;
    int     GetSizeOfRowPerms() const {return szoperm_; }

private:
    bool    LUDedSingular() const;

private:
    int*    rowperm_;//row permutation vector containing indices of permuted positions
    int     szoperm_;//length of rowperm_
    int     moneton_;//(-1)^n; positive if number of interchanges is even; negative orherwise
};


// -------------------------------------------------------------------------
// InitPermuts: initialize vector of row permutations
//
inline
void NSmatrix::InitPermuts( int size )
{
    if( !GetMaster()) {
        priverror( "Slave attempt of initialization of vector of permutations." );
        return;
    }
    if( size < 1 )
        return;
    if( size != szoperm_ ) {
        DestroyPermuts();
        rowperm_ = ( int* )malloc( size * sizeof( int ));
        if( rowperm_ == NULL ) {
            priverror( "Not enough memory." );
            return;
        }
        szoperm_ = size;
    }
    for( int n = 0; n < szoperm_; n++ )
        rowperm_[n] = n;
}

// -------------------------------------------------------------------------
// CopyPermuts: copy vector of row permutations
//
inline
void NSmatrix::CopyPermuts( const NSmatrix& right )
{
    DestroyPermuts();
    if( GetMaster()) {
        rowperm_ = NULL;
        szoperm_ = right.szoperm_;
        if( szoperm_ < 1 || right.rowperm_ == NULL )
            return;
        rowperm_ = ( int* )malloc( szoperm_ * sizeof( int ));
        if( rowperm_ == NULL ) {
            priverror( "Not enough memory." );
            return;
        }
        memcpy( rowperm_, right.rowperm_, szoperm_ * sizeof( int ));
        return;
    }
    else {
        rowperm_ = right.rowperm_;
        szoperm_ = right.szoperm_;
    }
}

// -------------------------------------------------------------------------
// GetPermutsAt: get value of row permutations
//
inline
int NSmatrix::GetPermutsAt( int pos ) const
{
    if( rowperm_ == NULL || GetSizeOfRowPerms() <= pos ) {
        priverror( "Memory access error." );
        return pos;
    }
    return rowperm_[pos];
}

// -------------------------------------------------------------------------
// SwapPermPoss: swap position values in the permutation vector
//
inline
int NSmatrix::SwapPermPoss( int i1, int i2 )
{
    if( !GetMaster()) {
        priverror( "Slave attempt to swap permutation vector values." );
        return PSL_ERR_ILLEGAL;
    }
    int size = GetSizeOfRowPerms();
    int tmp;
    if( size <= i1 || size <= i2 || rowperm_ == NULL ) {
        priverror( "Memory access error." );
        return PSL_ERR_ADDRESS;
    }
    if( i1 == i2 )
        return PSL_SUCCESS;
    tmp = rowperm_[i1];
    rowperm_[i1] = rowperm_[i2];
    rowperm_[i2] = tmp;
    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// PermuteVector: permute vector using information of row permutations
//
inline
int NSmatrix::PermuteVector( Pslvector& v ) const
{
    int size = GetSizeOfRowPerms();
    int n, k, pk;
    double  tmp;

    if( rowperm_ == NULL ) {
        priverror( "Memory access error." );
        return PSL_ERR_ADDRESS;
    }
    if( v.GetSize() != GetSizeOfRowPerms()) {
        priverror( "Unable to permute vector; inconsistent dimensions." );
        return PSL_ERR_DIM;
    }
    for( n =  0; n < size; n++ ) {
        for( k = GetPermutsAt( n ); n < k; k = GetPermutsAt( k ));
        if( k < n )
            continue;
        //k == n
        pk = GetPermutsAt( k );
        if( pk == n )
            continue;

        tmp = v.GetValueAt( n );
        while( pk != n ) {
            v.SetValueAt( k, v.GetValueAt( pk ));
            k = pk;
            pk = GetPermutsAt( k );
        }
        v.SetValueAt( k, tmp );
    }
    return PSL_SUCCESS;
}

// -------------------------------------------------------------------------
// DestroyPermuts: destroy vector of row permutations
//
inline
void NSmatrix::DestroyPermuts()
{
    if( !GetMaster())
        return;
    if( rowperm_ ) {
        free( rowperm_ );
        rowperm_ = NULL;
    }
    szoperm_ = 0;
}

#endif//__nsmatrix__
