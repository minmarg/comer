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
#include <time.h>
#include "ext/psl.h"
#include "ext/rng.h"
#include "ext/rv/rvexp.h"
#include "slicesmpl.h"

//FILE GLOBALS
//random number generators
MTRng   _slicesmplRNGss_;
MTRng   _slicesmplRNGtst_;

//default size of vectors
const unsigned int SliceSmpl::s_defszpnts_ = 1000; //default vector size of sampled points
const unsigned int SliceSmpl::s_defszpntscur_ = 100; //default vector size of currently sampled points

// -------------------------------------------------------------------------
// constructor
//
SliceSmpl::SliceSmpl()
:   feval_( NULL ),
    feparams_( NULL ),

    w_( 1.0 ),
    m_( 0 ),
    mset_( false ),

    xlbset_( false ),
    xrbset_( false ),
    xlb_( 0.0 ),
    xrb_( 0.0 ),
    xinit_( 0.0 ),
    lfxinit_( 0.0 ),

    thin_( 3 ),
    burnin_( 100 ),
    noss_( 1 ),
    accpnts_( NULL ),
    accpntscur_( NULL ),
    nofevals_( 0 ),
    nofevalscur_( 0 ),

    run_( false )
{
    time_t tm;
    clock_t cl = clock();
    time( &tm );
    _slicesmplRNGss_.Set((unsigned long)(size_t)this +(unsigned long)(size_t)(&_slicesmplRNGss_) +
            (unsigned long)tm +(unsigned long)cl);
    _slicesmplRNGtst_.Set((unsigned long)(size_t)this +(unsigned long)(size_t)(&_slicesmplRNGtst_) +
            (unsigned long)tm +(unsigned long)cl);
    InitAccpnts();
    InitCurAccpnts();
}

// destructor
//
SliceSmpl::~SliceSmpl()
{
    DestroyCurAccpnts();
    DestroyAccpnts();
}

// =========================================================================
// Run: start sampling from f(x)
//
int SliceSmpl::Run()
{
    const int   lcn_mlimit = 1e6;
    double  x1, lfx1;
    double  value;
    int n, t, err;

    if( GetAllSamples() == NULL || GetLastSamples() == NULL )
        return PSL_ERR_ADDRESS;

    if( GetNoReqSamples() <= 0 )
        return 0;

    if( GetBurnin() < 0 )
        return SliceSmpl_ERR_BININVALID;

    if( GetThin() <= 0 )
        return SliceSmpl_ERR_THININVALID;

    if( GetW() <= gscd_SLICE_XEPS )
        return SliceSmpl_ERR_WTOOSMALL;

    if( IsMSet()) {
        if( GetM() <= 0 )
            return SliceSmpl_ERR_MINVALID;
        if( lcn_mlimit < GetM())
            return SliceSmpl_ERR_MTOOLARGE;
    }

    if( IsLeftXSet() && IsRightXSet() && GetRightX() <= GetLeftX())
        return SliceSmpl_ERR_BNDSINVALID;

    if( IsLeftXSet() && GetX0() < GetLeftX())
        return SliceSmpl_ERR_INITXDOMAIN;

    if( IsRightXSet() && GetRightX() < GetX0())
        return SliceSmpl_ERR_INITXDOMAIN;

    if( GetFEval() == NULL )
        return SliceSmpl_ERR_FEVALADDRESS;


    GetLastSamples()->Clear();
    //initialise current number of function evaluations
    SetNoCurFevals( 0 );

    //get log f(x) at initial point
    if(( err = ( *GetFEval())( GetX0(), &value, GetParams())) != PSL_OK )
        return err;
    SetLFvalAtX0( value );
    IncNoFevals();
    IncNoCurFevals();

    //omit number of iterations in burn-in
    for( n = 0; n < GetBurnin(); n++ ) {
        //sample a new point
        if(( err = Sample( GetX0(), GetLFvalAtX0(), &x1, &lfx1 )) != PSL_OK )
            return err;
        SetX0( x1 );
        SetLFvalAtX0( lfx1 );
    }
    //execute slice sampling
    for( n = 0; n < GetNoReqSamples(); n++ ) {
        for( t = 0; t < GetThin(); t++ ) {
            //sample a new point
            if(( err = Sample( GetX0(), GetLFvalAtX0(), &x1, &lfx1 )) != PSL_OK )
                return err;
        }
        //push accepted point
        GetAllSamples()->Push( x1 );
        GetLastSamples()->Push( x1 );
        SetX0( x1 );
        SetLFvalAtX0( lfx1 );
    }

    return 0;
}

// -------------------------------------------------------------------------
// Sample: sample from slice
//
int SliceSmpl::Sample( double x0, double lfx0, double* x1, double* lfx1 )
{
    const int   lcn_maxit = 1000;
    bool    mset = IsMSet() && 1 < GetM();
    double  logy, val;
    double  u1, u2, L, R;
    int     J = 1, K = 1;
    int n, err;

    if( x1 == NULL || lfx1 == NULL )
        return PSL_ERR_ADDRESS;

    RVExp   re( _slicesmplRNGss_ );//exponential random variable
    re.SetScale( 1.0 );
    if(( err = re.Gen( &val )) != 0 )
        return err;

    //determine the slice, density y, level, in log terms.
    logy = lfx0 - val;

    //find the initial interval to sample from;
    //guarantee x0 is in [L,R], even with roundoff
    u1 = GetW() * _slicesmplRNGtst_.GetDouble();
    L = x0 - u1;
    R = x0 + ( GetW() - u1 );

    if( mset ) {
        //set limits on number of steps
        J = ( int )floor( GetM() * _slicesmplRNGtst_.GetDouble());
        K = ( GetM() - 1 ) - J;
    }

    //expand interval until its ends are outside the slice,
    //or until the limit on steps is reached
    for( n = 0; 0 < J &&( IsLeftXSet()? GetLeftX() < L: 1 ); n++, L -= GetW()) {
        if( lcn_maxit <= n )
            return SliceSmpl_ERR_STEPOUTITER;
        //log f(L)
        if(( err = ( *GetFEval())( L, &val, GetParams())) != PSL_OK )
            return err;
        IncNoFevals();
        IncNoCurFevals();
        if( val <= logy )
            break;
        if( mset ) J--;
    }
    for( n = 0; 0 < K &&( IsRightXSet()? R < GetRightX(): 1 ); n++, R += GetW()) {
        if( lcn_maxit <= n )
            return SliceSmpl_ERR_STEPOUTITER;
        //log f(R)
        if(( err = ( *GetFEval())( R, &val, GetParams())) != PSL_OK )
            return err;
        IncNoFevals();
        IncNoCurFevals();
        if( val <= logy )
            break;
        if( mset ) K--;
    }

    //shrink interval to left and right bounds
    if( IsLeftXSet() && L < GetLeftX())
        L = GetLeftX();
    if( IsRightXSet() && GetRightX() < R )
        R = GetRightX();

    //sample from the interval, shrinking it on each rejection
    for( n = 0;; n++ ) { 
        if( lcn_maxit <= n )
            return SliceSmpl_ERR_SHRINKINITER;

        *x1 = L + ( R - L ) * _slicesmplRNGtst_.GetDouble();
        //log f(L)
        if(( err = ( *GetFEval())( *x1, lfx1, GetParams())) != PSL_OK )
            return err;
        IncNoFevals();
        IncNoCurFevals();

        if( logy <= *lfx1 )
            break;
        if( x0 < *x1 )
            R = *x1;
        else
            L = *x1;
    }

    return 0;
}
