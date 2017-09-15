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
#include "rvdisc.h"

// =========================================================================
// RVDisc
//
RVDisc::RVDisc( Rng& rng, TRVDsc met )
:   RVar( rng ),
    met_( met ),
    probs_( NULL ),
    psize_( 0 ),
    changed_( true )
{
}

RVDisc::~RVDisc()
{
}

// -------------------------------------------------------------------------
// GenIAlias: generate discrete random variable given probability vector;
//  using alias method proposed by 
//      Walker AJ. (1977) ACM Trans Math Software 3(3), 253-6.
//  Alias tables are constructed using O(n) algorithm by
//      Kronmal and Peterson. (1979) The American Statistician 33(4), 214-8.
//  The Kronmal&Peterson method is also present in 
//      L.Devroye. (1986) Non Uniform Random Variate Generation, p.109;
//  the general walker algorithm may also be found in there as well as in
//      Knuth TAOCP v2 p.120 and p.585 (solution to exercise 7)
//
int RVDisc::GenIAlias( int* ind )
{
    const double* probs = GetProbs();
    int plen = GetPSize();
    bool chn = GetChanged();
    int err;

    if( ind == NULL )
        return PSL_ERR_ADDRESS;

    if( probs == NULL || plen < 1 )
        return PSL_ERR_INVALID;

    if( chn && ( err = MakeAliasTables()) != PSL_OK )
        return err;

    int     x = ( int )GetRng().GetI(( unsigned long )plen );
    double  v = GetRng().GetDouble();

    if( x < 0 || plen < x )
        return PSL_ERR_ILLEGAL;

    if( plen <= x )
        x--;

    if( q_.GetSize() <= x || j_.GetSize() <= x )
        return PSL_ERR_ILLEGAL;

    if( v <= q_.GetValueAt( x ))
        *ind = x;
    else
        *ind = j_.GetValueAt( x );

    return PSL_OK;
}

// -------------------------------------------------------------------------
// MakeAliasTables: make alias table and associated vector of probabilities
// NOTE: given probabilities must be normalized!
//
int RVDisc::MakeAliasTables()
{
    Ivector grt, sml;//greater, smaller sets (stacks) of indices
    double  qv;
    int plen = GetPSize();
    int hplen = plen >> 1 + 1;
    const double* probs = GetProbs();
    static const double eps = 1.e-12 * ( double )plen;
    static const double qeps = eps * ( double )plen;
    int n, k, l, err;

    q_.Clear();
    j_.Clear();

    if( probs == NULL || plen < 1 )
        return PSL_ERR_INVALID;

    grt.Allocate( hplen );
    sml.Allocate( hplen );

    q_.Allocate( plen );
    //preset number of elements;
    //for those unused values, q_ will be 1
    j_.Reserve( plen );

    for( n = 0; n < plen; n++ ) {
        qv = ( double )plen * probs[n];
        q_.Push( qv );
        if( qv < 1.0 )
            sml.Push( n );
        else
            grt.Push( n );
    }

    //check if normalized
    qv = q_.Sum();
    if( qv < plen - qeps || plen + qeps < qv )
        return PSL_ERR_INVALID;

    //there should always be some elements in `grt'
    if( grt.GetSize() < 1 )
        return PSL_ERR_ILLEGAL;
  
    for( ; sml.GetSize(); ) {
        k = grt.GetLast();
        l = sml.GetLast();//-- q_[l] finalized
// fprintf(stdout,"  %d %d; k=%d l=%d q[k]=%.14g q[l]=%.14g",grt.GetSize(),sml.GetSize(),k,l,q_.GetValueAt(k),q_.GetValueAt(l));
        j_.SetValueAt( l, k );//-- j_[l] finalized
        q_.AddValueAt( k, q_.GetValueAt(l) - 1.0 );
// fprintf(stdout," q[k]=%g\n",q_.GetValueAt(k));
        qv = q_.GetValueAt(k);
        //remove last element from smaller set;
        //NOTE: at this point! not to mix the order
        sml.DecDim();
        if( !sml.GetSize()) {
            //manage fp errors
            if( qv < 1.0 - eps || 1.0 + eps < qv )
                return PSL_ERR_ILLEGAL;
            q_.SetValueAt( k, 1.0 );
            break;
        }
        if( qv < 1.0 ) {
            grt.DecDim();//remove last element from greater set
            sml.Push( k );//add index k to smaller set
        }
    }

    ResetChanged();
    return PSL_OK;
}

// =========================================================================
// GenILinear: generate discrete random variable by inverse search; 
//  the algorithm may be useful when probability vector changes often and 
//  alias method looses its effectiveness;
//  NOTE: this is the simplest implementation of linear search; more 
//  elaborate implementation may be constructed using binary trees or other
//  effective data structures (e.g. Knuth,Devroye,etc.); however when 
//  probability vector does not change frequently, the alias method 
//  should be preferable
//
int RVDisc::GenIInvLin( int* ind )
{
    const double* probs = GetProbs();
    int plen = GetPSize();
    int n, err;

    if( ind == NULL )
        return PSL_ERR_ADDRESS;

    if( probs == NULL || plen < 1 )
        return PSL_ERR_INVALID;

    double  v = GetRng().GetDouble1();
    double  F = 0.0, pf;

    for( n = 0; n < plen; n++ ) {
        pf = F;
        F += probs[n];
        if( pf <= v && v < F ) {
            *ind = n;
            return PSL_OK;
        }
    }
    *ind = plen;
    return PSL_OK;
}
