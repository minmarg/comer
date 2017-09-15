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
#include "cars.h"

//FILE GLOBALS
//random number generators
MTRng   _carsRNGss_;
MTRng   _carsRNGuh_;

//default size of vectors
const int CARS::s_defsz_ = 100;

// -------------------------------------------------------------------------
// constructor
//
CARS::CARS()
:   noss_( 0 ),
    heval_( NULL ),
    heparams_( NULL ),
    accpnts_( NULL ),
    xv_( NULL ),
    hxv_( NULL ),
    hdxv_( NULL ),
    srtxv_( NULL ),
    xlb_( 0.0 ),
    xub_( 0.0 ),
    uhlb_( 0.0 ),
    uhub_( 0.0 ),
    zv_( NULL ),
    uhzv_( NULL ),
    maxuh_( 0.0 ),
    euhcdf_( NULL ),
    euharea_( 0.0 ),
    lowb_( false ),
    uppb_( false ),
    run_( false )
{
    time_t  tm;
    clock_t cl = clock();
    time( &tm );
    _carsRNGss_.Set((unsigned long)(size_t)this +(unsigned long)(size_t)(&_carsRNGss_) +
            (unsigned long)tm +(unsigned long)cl);
    _carsRNGuh_.Set((unsigned long)(size_t)this +(unsigned long)(size_t)(&_carsRNGuh_) +
            (unsigned long)tm +(unsigned long)cl);
    InitVectors();
}

// destructor
//
CARS::~CARS()
{
    DestroyVectors();
}

// =========================================================================
// Run: start sampling distribution f(x)
//
int CARS::Run()
{
    int s;

    SetRun( true );
    Initialize();
    for( s = 0; s < GetNoIterats(); s++ )
        Sample();

    return 0;
}

// -------------------------------------------------------------------------
// Initialize: initialization step in sampling by the ARS algorithm
//
void CARS::Initialize()
{
    int     ndx;
    double  x0, uhval, hx0, hdx0;
    double  intexpuh; //integral of exponentiated upper hull
    int n;

    CheckVectors();

    ndx = ( int )( size_t )GetSrtX()->GetValueAt( 0 );
    if( ndx < 0 || GetX()->GetSize() <= ndx )
        throw myruntime_error( "CARS: Initialize: Invalid index.");

    x0 = GetX()->GetValueAt( ndx );
    hx0 = GetHX()->GetValueAt( ndx );
    hdx0 = GetHDX()->GetValueAt( ndx );
    if( GetLowerSet()) {
        uhval = ( GetLowerX() - x0 ) * hdx0 + hx0;
        SetUHatLower( uhval );
    }
    if( GetUpperSet()) {
        uhval = ( GetUpperX() - x0 ) * hdx0 + hx0;
        SetUHatUpper( uhval );
    }

    if( GetUpperSet()) {
        if( GetLowerSet()) {
            SetMaxUHull( SLC_MAX( GetUHatUpper(), GetUHatLower()));
            if( fabs( hdx0 ) < gpsc_EPSILON )
                //derivative at lower bound is (in limit) flat
                // u_k(x) = h(x_j); INT = e^{h(x_j)} (xub-xlb)
                intexpuh = prexp(( GetUHatUpper() + GetUHatLower()) * 0.5 - GetMaxUHull()) * 
                          ( GetUpperX() - GetLowerX());
            else
                //INT_{xlb}^{xub} e^{u_k(x)} dx 
                //      = e^{h(xub)} (1 - e^{(xlb-xub)h'(xub)}) / h'(xub)
                //      = e^{h(xub)} (1 - e^{uh(xlb)-uh(xub)}) / h'(xub)
                intexpuh = prexp( GetUHatUpper() - GetMaxUHull()) * 
                          ( 1 - prexp( GetUHatLower() - GetUHatUpper())) / hdx0;
        }
        else {
            SetMaxUHull( GetUHatUpper());
            if( fabs( hdx0 ) < gpsc_EPSILON )
                throw myruntime_error(
                    "CARS: Unable to calculate integral: Try other starting point.");
            //INT_{-inf}^{xub} e^{u_k(x)} dx = e^{h(xub)} / h'(xub)
            intexpuh = 1.0 / hdx0; //e^{h(xub)-h(xub)}
        }
    }
    else {
        if( GetLowerSet()) {
            SetMaxUHull( GetUHatLower());
            if( fabs( hdx0 ) < gpsc_EPSILON )
                throw myruntime_error(
                    "CARS: Unable to calculate integral: Try other starting point.");
            //INT_{xlb}^{inf} e^{u_k(x)} dx = -e^{h(xlb)} / h'(xlb)
            intexpuh = -1.0 / hdx0;
        }
        else {
            intexpuh = 0.0;
            if( GetX()->GetSize() < 2 )
                throw myruntime_error(
                    "CARS: Unbounded sampling requires >=2 initial arguments.");
        }
    }

    SetEUHArea( intexpuh );
    GetEUHCDF()->Push( 1.0 );

    //update upper hull
    for( n = 0; n+1 < GetSrtX()->GetSize(); n++ ) {
        UpdateUH( n );
        if( 0 < n )
            UpdateUH( n-1, false/*update rather than insert*/ );
    }
    if( n < GetSrtX()->GetSize())
        Update( n );

    ndx = ( int )( size_t )GetSrtX()->GetValueAt( 0 );
    if( ndx < 0 || GetHDX()->GetSize() <= ndx )
        throw myruntime_error( "CARS: Initialize: Invalid index.");
    if( !GetLowerSet() && GetHDX()->GetValueAt( ndx ) < gpsc_EPSILON )
        throw myruntime_error("CARS: derivative value of the beginning point is invalid.");

    ndx = ( int )( size_t )GetSrtX()->GetValueAt( GetSrtX()->GetSize()-1 );
    if( ndx < 0 || GetHDX()->GetSize() <= ndx )
        throw myruntime_error( "CARS: Initialize: Invalid index.");
    if( !GetUpperSet() && -gpsc_EPSILON < GetHDX()->GetValueAt( ndx ))
        throw myruntime_error("CARS: derivative value of the end point is invalid.");
}

// -------------------------------------------------------------------------
// Sample: step of sampling by the ARS algorithm
//
void CARS::Sample()
{
    CheckVectors();

    double  urv1, logurv1;
    double  x0, xn;
    double  x, hx, hdx;
    double  x1, hx1, hdx1;
    double  smpl, hsmpl;
    double  uhsmpl, lhsmpl;
    bool    sampled;
    int     in, nd, loc;

    for( sampled = false; !sampled; ) {
        SampleUH( &smpl, &in );
#ifdef CARSTESTPRINT
        fprintf(stderr,"Sample: smpl=%g in=%d\n",smpl,in);
#endif
        urv1 = _carsRNGss_.GetDouble0();
        if( urv1 <= 0.0 || 1.0 < urv1 )
            throw myruntime_error("CARS: Sample: Invalid RNG value.");
        logurv1 = log( urv1 );

        if( in < 0 || GetSrtX()->GetSize() <= in )
            throw myruntime_error("CARS: Sample: Invalid index.");
        nd = ( int )( size_t )GetSrtX()->GetValueAt( in );
        if( nd < 0 || GetX()->GetSize() <= nd )
            throw myruntime_error("CARS: Sample: Memory access error.");

        x = GetX()->GetValueAt( nd );
        hx = GetHX()->GetValueAt( nd );
        hdx = GetHDX()->GetValueAt( nd );

        if( in == 0 )
            x0 = x;
        else {
            nd = ( int )( size_t )GetSrtX()->GetValueAt( 0 );
            if( nd < 0 || GetX()->GetSize() <= nd )
                throw myruntime_error("CARS: Sample: Memory access error.");
            x0 = GetX()->GetValueAt( nd );
        }
        if( in == GetSrtX()->GetSize()-1 )
            xn = x;
        else {
            nd = ( int )( size_t )GetSrtX()->GetValueAt( GetSrtX()->GetSize()-1 );
            if( nd < 0 || GetX()->GetSize() <= nd )
                throw myruntime_error("CARS: Sample: Memory access error.");
            xn = GetX()->GetValueAt( nd );
        }        

        //upper hull value at smpl
        uhsmpl = hx - GetMaxUHull() + ( smpl - x ) * hdx;

        if( x0 < smpl && smpl < xn ) {
            if( x < smpl ) {
                //smpl is above x(i); needing x(i+1) for lower hull value
                if( GetSrtX()->GetSize() <= in+1 )
                    throw myruntime_error("CARS: Sample: Index error.");
                nd = ( int )( size_t )GetSrtX()->GetValueAt( in+1 );
                if( nd < 0 || GetX()->GetSize() <= nd )
                    throw myruntime_error("CARS: Sample: Memory access error.");
                x1 = GetX()->GetValueAt( nd );
                hx1 = GetHX()->GetValueAt( nd );
            }
            else {
                x1 = x; hx1 = hx;
                //smpl is below x(i); needing x(i-1) for lower hull value
                if( in < 1 )
                    throw myruntime_error("CARS: Sample: Index error.");
                nd = ( int )( size_t )GetSrtX()->GetValueAt( in-1 );
                if( nd < 0 || GetX()->GetSize() <= nd )
                    throw myruntime_error("CARS: Sample: Memory access error.");
                x = GetX()->GetValueAt( nd );
                hx = GetHX()->GetValueAt( nd );
            }

            if( x1 <= x )
                throw myruntime_error("CARS: Sample: Unsorted arguments.");
            //lower hull value at smpl
            lhsmpl = hx1 - GetMaxUHull() + ( smpl - x1 ) * ( hx1 - hx ) / ( x1 - x );

            if( logurv1 < lhsmpl - uhsmpl ) {
                //urv < e^{lh(smpl)-uh(smpl)}
                sampled = true;
                continue;
            }
        }
        //include smpl into vector of x's
        NewX( smpl, false/*do not clear*/, &loc, &hsmpl );

        if( logurv1 < hsmpl - GetMaxUHull() - uhsmpl )
            //urv < e^{h(smpl)-uh(smpl)}
            sampled = true;
        Update( loc );
    }
    //push accepted point
    GetResPoints()->Push( smpl );
}

// -------------------------------------------------------------------------
// SampleUH: sampling from by the ARS algorithm
//
void CARS::SampleUH( double* smpl, int* in )
{
    CheckVectors();
    if( smpl == NULL || in == NULL )
        throw myruntime_error("CARS: SampleUH: Invalid parameters.");

    double  urv = _carsRNGuh_.GetDouble0();
    double  intexpuh = GetEUHArea(); //integral of exponentiated upper hull
    double  loginteuh = 0.0; //log of integral of exp. upper hull
    double  x, hx, hdx, z, uhz;
    double  x1, hx1, hdx1, z1, uhz1;
    double  cdf, cdf1, luha;
    int sgnhdx1;
    int i, nd;

    if( urv <= 0.0 || 1.0 < urv )
        throw myruntime_error("CARS: SampleUH: Invalid RNG value.");
    if( GetEUHCDF()->GetSize() != GetX()->GetSize() || GetEUHCDF()->GetSize() < 1 )
        throw myruntime_error("CARS: SampleUH: Invalid cdf size.");
    if( GetZ()->GetSize() != GetX()->GetSize())
        throw myruntime_error("CARS: SampleUH: Invalid size of upper hull.");

    //find cdf piece to sample from
    for( i = 0; i < GetEUHCDF()->GetSize(); i++ )
    {
        nd = ( int )( size_t )GetSrtX()->GetValueAt( i );
        if( nd < 0 || GetX()->GetSize() <= nd )
            throw myruntime_error("CARS: SampleUH: Memory access error.");
        x1 = GetX()->GetValueAt( nd );
        hx1 = GetHX()->GetValueAt( nd );
        hdx1 = GetHDX()->GetValueAt( nd );
        sgnhdx1 = signbit( hdx1 )? -1: 1;

        if( i < GetEUHCDF()->GetSize()) {
            z1 = GetZ()->GetValueAt( nd );
            uhz1 = GetUHZ()->GetValueAt( nd );
        }

        cdf1 = GetEUHCDF()->GetValueAt( i );
        if( urv <= cdf1 )
            break;

        x = x1; hx = hx1; hdx = hdx1;
        z = z1; uhz = uhz1;
        cdf = cdf1;
    }
#ifdef CARSTESTPRINT
    fprintf(stderr,"SampleUH: urv=%g cdf=%g i=%d; EUHCDF:...\n",urv,cdf,i);
    GetEUHCDF()->Print(stderr," %.6g");
#endif

    if( GetEUHCDF()->GetSize() <= i )
        throw myruntime_error("CARS: SampleUH: Invalid cdf.");
    if( uhz == gpsc_SMALLPW )
        throw myruntime_error("CARS: SampleUH: Invalid upper hull value");

    if( 0.0 < intexpuh )
        loginteuh = log( intexpuh );

    if(( *in = i ) <= 0 ) {
        //sample below the lowest point
        if( GetLowerSet()) {
            luha = GetUHatLower() - GetMaxUHull() - loginteuh;
            if( fabs( hdx1 ) < gpsc_EPSILON )
                //from (see below):
                //smpl = xlb + log(1 + S*(urv-0) h'(x_0) / e^{uh(xlb)}) / h'(x_0)
                //smpl = xlb + S urv / e^{uh(xlb)}
                *smpl = GetLowerX() + prexp( -luha ) * urv;
            else {
                //from (see below):
                //smpl = xlb + log(1 + S*(urv-0) h'(x_0) / e^{uh(xlb)}) / h'(x_0)
                luha = log( urv ) + log( fabs( hdx1 )) - luha;
#ifdef CARSTESTPRINT
                fprintf(stderr,"SampleUH: uh(xlb)=%g maxuh=%g logint=%g luha=%g hdx1=%g sgn=%d\n",
                    GetUHatLower(),GetMaxUHull(),loginteuh,luha,hdx1,sgnhdx1);
#endif
                if( luha < gpsc_MAXEPOW ) {
                    luha = sgnhdx1 * prexp( luha );
                    if( luha <= -1.0 )
                        throw myruntime_error("CARS: SampleUH: Failed to sample from upper hull.");
                    *smpl = GetLowerX() + log( 1.0 + luha ) / hdx1;
                }
                else
                    *smpl = GetLowerX() + luha / hdx1;
            }
        }
        else {
            //From equation:
            //INT_{-inf}^{smpl} e^{h(x_0)+(x-x_0)h'(x_0)} dx = S*(urv-0)
            //e^{h(x_0)} e^{-x_0 h'(x_0)} 1/h'(x_0) (e^{smpl h'(x_0)}-e^{-inf h'(x_0)})
            //      = S urv
            //e^{h(x_0)} e^{(smpl-x_0)h'(x_0)} / h'(x_0) = S urv
            //smpl = x_0 + log(S urv h'(x_0) / e^{h(x_0)}) / h'(x_0)
            if( hdx1 <= 0.0 )
                throw myruntime_error("CARS: SampleUH: Non-positive derivative at x0.");
            *smpl = ( x1 * hdx1 + loginteuh + log( urv * hdx1 ) - hx1 + GetMaxUHull()) / hdx1;
        }
    }
    else {
        //sample between i-1 and i
        luha = uhz - GetMaxUHull() - loginteuh;
        if( fabs( hdx1 ) < gpsc_EPSILON )
            //from (see below): log(1+x) = x as x->0
            //smpl = z_{i-1} + log(1 + S*(urv-cdf_{i-1}) h'(x_i) / e^{uh(z_{i-1})}) / h'(x_i)
            //smpl = z_{i-1} + S*(urv-cdf_{i-1}) / e^{uh(z_{i-1})}
            *smpl = z + prexp( -luha ) * ( urv - cdf );
        else {
            //From equation:
            //INT_{z_{i-1}}^{smpl} e^{uh(z_{i-1})+(x-z_{i_1})h'(x_i)} dx = S*(urv-cdf_{i-1})
            //e^{uh(z_{i-1})} e^{-z_{i-1}h'(x_i)} 1/h'(x_i) (e^{smpl h'(x_i)}-e^{z_{i-1}h'(x_i)})
            //      = S*(urv-cdf_{i-1})
            //e^{uh(z_{i-1})} (e^{(smpl-z_{i-1})h'(x_i)} - 1) / h'(x_i) = S*(urv-cdf_{i-1})
            //smpl = z_{i-1} + log(1 + S*(urv-cdf_{i-1}) h'(x_i) / e^{uh(z_{i-1})}) / h'(x_i)
            if( cdf < urv  ) {
                luha = log( urv - cdf ) + log( fabs( hdx1 )) - luha;
                if( luha < gpsc_MAXEPOW ) {
                    luha = sgnhdx1 * prexp( luha );
                    if( luha <= -1.0 )
                        throw myruntime_error("CARS: SampleUH: Failed to sample from upper hull.");
                    *smpl = z + log( 1.0 + luha ) / hdx1;
#ifdef CARSTESTPRINT
                    fprintf(stderr,"SampleUH: SMPL=%g z=%g hdx1=%g sgn=%d\n",*smpl,z,hdx1,sgnhdx1);
#endif
                }
                else
                    *smpl = z + luha / hdx1;
            }
            else
                *smpl = z + luha / hdx1;
        }
    }
}

// -------------------------------------------------------------------------
// Update: update step in sampling by the ARS algorithm
//
void CARS::Update( int n )
{
    CheckVectors();

    double  uhval, mval;
    int i;

    //update upper hull
    UpdateUH( n );
    if( 0 < n )
        UpdateUH( n-1, false/*update rather than insert*/ );

    //update max value of upper hull
    for( i = 0; i < GetUHZ()->GetSize(); i++ ) {
        uhval = GetUHZ()->GetValueAt( i );
        if( i == 0 || GetMaxUHull() < uhval )
            SetMaxUHull( uhval );
    }
    if( GetLowerSet() && GetMaxUHull() < GetUHatLower())
        SetMaxUHull( GetUHatLower());
    if( GetUpperSet() && GetMaxUHull() < GetUHatUpper())
        SetMaxUHull( GetUHatUpper());

    //update integral of exponentiated upper hull
    UpdateEUHCDF();
}

// -------------------------------------------------------------------------
// UpdateUH: update upper hull at point `n'; `insert' indicates insertion or
//  replacement of the new calculated value
//
void CARS::UpdateUH( int n, bool insert )
{
    CheckVectors();

    int     nd = -1, nd1 = -1;
    double  x, x1;
    double  hx, hx1;
    double  hdx, hdx1;
    double  z, uhz, uhval;

    if( n < 0 || GetSrtX()->GetSize() <= n )
        throw myruntime_error("CARS: UpdateUH: Memory access error.");
    if( GetZ()->GetSize() < GetX()->GetSize()) {
        GetZ()->Allocate( GetX()->GetCapacity());
        GetZ()->Reserve( GetX()->GetSize());
        GetUHZ()->Allocate( GetX()->GetCapacity());
        GetUHZ()->Reserve( GetX()->GetSize());
    }

    nd = ( int )( size_t )GetSrtX()->GetValueAt( n );
    if( nd < 0 || GetX()->GetSize() <= nd )
        throw myruntime_error("CARS: UpdateUH: Memory access error.");
    x = GetX()->GetValueAt( nd );
    hx = GetHX()->GetValueAt( nd );
    hdx = GetHDX()->GetValueAt( nd );
    if( n + 1 < GetSrtX()->GetSize()) {
        nd1 = ( int )( size_t )GetSrtX()->GetValueAt( n+1 );
        if( nd1 < 0 || GetX()->GetSize() <= nd1 )
            throw myruntime_error("CARS: UpdateUH: Memory access error.");
        x1 = GetX()->GetValueAt( nd1 );
        hx1 = GetHX()->GetValueAt( nd1 );
        hdx1 = GetHDX()->GetValueAt( nd1 );
        if( x1 < x )
            throw myruntime_error("CARS: UpdateUH: Unsorted arguments.");
    }

    if( nd1 < 0 ) {
        GetZ()->SetValueAt( nd, 0.0 );
        GetUHZ()->SetValueAt( nd, gpsc_SMALLPW );//as mark of the last value
        //this is the last point
        if( GetUpperSet()) {
            uhval = ( GetUpperX() - x ) * hdx + hx;
            SetUHatUpper( uhval );
        }
        return;
    }

    if( hdx < hdx1 )
        throw myruntime_error("CARS: UpdateUH: Non log-concave function.");
    Intersection( nd, nd1, &z, &uhz );

    //update values
    GetZ()->SetValueAt( nd, z );
    GetUHZ()->SetValueAt( nd, uhz );

    if( n == 0 ) {
        if( GetLowerSet()) {
            uhval = ( GetLowerX() - x ) * hdx + hx;
            SetUHatLower( uhval );
        }
    }
}

// -------------------------------------------------------------------------
// Intersection: calculate intersection point betwween 2 tangents
//
void CARS::Intersection( int n1, int n2, double* pz, double* puhz )
{
    CheckVectors();

    double  x1, hx1, hdx1;
    double  x2, hx2, hdx2;
    double  hx12, hx21;
    double  diffhdx;
    double  z = 0.0, uhz = 0.0;

    if( n1 < 0 || GetX()->GetSize() <= n1 || 
        n2 < 0 || GetX()->GetSize() <= n2 )
        throw myruntime_error("CARS: Intersection: Memory access error.");

    x1 = GetX()->GetValueAt( n1 );
    hx1 = GetHX()->GetValueAt( n1 );
    hdx1 = GetHDX()->GetValueAt( n1 );

    x2 = GetX()->GetValueAt( n2 );
    hx2 = GetHX()->GetValueAt( n2 );
    hdx2 = GetHDX()->GetValueAt( n2 );

    if( x2 < x1 )
        throw myruntime_error("CARS: Intersection: Unsorted arguments.");

    //test for non-concavity                                          
    hx12 = hx1 + ( x2 - x1 ) * hdx1;
    hx21 = hx2 + ( x1 - x2 ) * hdx2;
    if( hx21 < hx1 || hx12 < hx2 )
        throw myruntime_error(
          "CARS: Intersection: Fragmentally non log-concave function.");

    diffhdx = hdx2 - hdx1;

    if( fabs( diffhdx ) <= gpsc_EPSILON ) {
        //intersection is at the midpoint, if nearly parallel
        z = ( x1 + x2 ) * 0.5;
        uhz = ( hx1 + hx2 ) * 0.5;
    }
    //find crosspoint from the equality of two tangent lines.
    //if... for greater numerical precision
    else if( fabs( hdx1 ) < fabs( hdx2 )) {
        z = x2 + ( hx1 - hx2 + ( x2 - x1 ) * hdx1 ) / diffhdx;
        uhz = hx1 + ( z - x1 ) * hdx1;
    }
    else {
        z = x1 + ( hx1 - hx2 + ( x2 - x1 ) * hdx2 ) / diffhdx;
        uhz = hx2 + ( z - x2 ) * hdx2;
    }

    if( z < x1 || x2 < z )
        throw myruntime_error("CARS: Intersection: Intersection at invalid point.");

    if( pz ) *pz = z;
    if( puhz ) *puhz = uhz;
}

// -------------------------------------------------------------------------
// UpdateEUHCDF: calculate CDF of exponentiated upper hull
//
void CARS::UpdateEUHCDF()
{
    CheckVectors();

    if( GetX()->GetSize() != GetZ()->GetSize())
        throw myruntime_error( "CARS: UpdateEUHCDF: Inconsistent vector sizes.");
    if( GetX()->GetSize() < 2 )
        //no upper hull
        return;

    const double accuracy = 1.e-6;

    int     ndx;
    double  x, hx, hdx, z, uhz;
    double  x1, hx1, hdx1, z1, uhz1;
    double  diffuhz, cdf;
    double  intexpuh; //integral of exponentiated upper hull
    char    strbuf[KBYTE];
    int i;

    GetEUHCDF()->Clear();

#ifdef CARSTESTPRINT
    fprintf(stderr,"UpdateEUHCDF: SrtX:...\n");
    for(i=0;i<GetSrtX()->GetSize();i++)
        fprintf(stderr," %d",(int)GetSrtX()->GetValueAt(i));
    fprintf(stderr,"\nUpdateEUHCDF: X:...\n");
    GetX()->Print(stderr," %.6f");
    fprintf(stderr,"UpdateEUHCDF: Z:...\n");
    GetZ()->Print(stderr," %.6f");
#endif

    ndx = ( int )( size_t )GetSrtX()->GetValueAt( i = 0 );
    if( ndx < 0 || GetX()->GetSize() <= ndx )
        throw myruntime_error( "CARS: UpdateEUHCDF: Invalid index.");

    x1 = x = GetX()->GetValueAt( ndx );
    hx1 = hx = GetHX()->GetValueAt( ndx );
    hdx1 = hdx = GetHDX()->GetValueAt( ndx );
    z = GetZ()->GetValueAt( ndx );
    uhz = GetUHZ()->GetValueAt( ndx );

    if( uhz == gpsc_SMALLPW )
        throw myruntime_error( "CARS: UpdateEUHCDF: Invalid upper hull value.");

    //lower bound processing
    if( GetLowerSet()) {
        if( fabs( hdx ) < gpsc_EPSILON )
            //derivative at the first x is (in limit) flat
            // u_k(x) = uh(xlb); INT = e^{uh(xlb)} (z-xlb)
            intexpuh = prexp( GetUHatLower() - GetMaxUHull()) * ( z - GetLowerX());
        else {
            diffuhz = GetUHatLower() - uhz;
            if( gpsc_MAXEPOW < diffuhz )
                //uh(xlb) >> uh(z)
                intexpuh = -prexp( GetUHatLower() - GetMaxUHull()) / hdx;
            else
                //INT_{xlb}^{z} e^{u_k(x)} dx = e^{uh(z)} (1 - e^{uh(xlb)-uh(z)}) / uh'(z)
                intexpuh = prexp( uhz - GetMaxUHull()) * ( 1.0 - prexp( diffuhz )) / hdx;
        }
    }
    else {
        if( fabs( hdx ) < gpsc_EPSILON )
            intexpuh = 0.0;
        else
            //INT_{-inf}^{z} e^{u_k(x)} dx = e^{uh(z)} / uh'(z)
            intexpuh = prexp( uhz - GetMaxUHull()) / hdx;
    }

    GetEUHCDF()->Push( intexpuh );

    //processing of interior of upper hull
    for( i = 1; i < GetSrtX()->GetSize(); i++ )
    {
        ndx = ( int )( size_t )GetSrtX()->GetValueAt( i );
        if( ndx < 0 || GetX()->GetSize() <= ndx )
            throw myruntime_error( "CARS: UpdateEUHCDF: Invalid index.");

        x1 = GetX()->GetValueAt( ndx );
        hx1 = GetHX()->GetValueAt( ndx );
        hdx1 = GetHDX()->GetValueAt( ndx );

        if( GetSrtX()->GetSize() <= i+1 )
            break;

        z1 = GetZ()->GetValueAt( ndx );
        uhz1 = GetUHZ()->GetValueAt( ndx );
        if( uhz1 == gpsc_SMALLPW )
            throw myruntime_error( "CARS: UpdateEUHCDF: Invalid upper hull value.");

        if( fabs( hdx1 ) < gpsc_EPSILON )
            //derivative at the first x is (in limit) flat
            // u_k(z) = uh(z); INT = e^{uh(z)} (z1-z)
            intexpuh += prexp(( uhz1 + uhz ) * 0.5 - GetMaxUHull()) * ( z1 - z );
        else {
            diffuhz = uhz - uhz1;
            if( gpsc_MAXEPOW < diffuhz )
                //uh(z) >> uh(z1)
                intexpuh -= prexp( uhz - GetMaxUHull()) / hdx1;
            else
                //INT_{z}^{z1} e^{u_k(x)} dx = INT e^{uh(z1)+(x-z1)uh'(z1)}
                //  = e^{uh(z1)} e^{-z1 uh'(z1)} 1/uh'(z1) (e^{z1 uh'(z1)}-e^{z uh'(z1)})
                //  = e^{uh(z1)} (1 - e^{(z-z1)uh'(z1)}) / uh'(z1)
                //  = e^{uh(z1)} (1 - e^{uh(z)-uh(z1)}) / uh'(z1)
                intexpuh += prexp( uhz1 - GetMaxUHull()) * ( 1.0 - prexp( diffuhz )) / hdx1;
        }
        GetEUHCDF()->Push( intexpuh );
        x = x1; hx = hx1; hdx = hdx1;
        z = z1; uhz = uhz1;
    }

    //upper bound processing
    if( GetUpperSet()) {
        if( fabs( hdx1 ) < gpsc_EPSILON )
            //derivative at the first x is (in limit) flat
            // u_k(x) = uh(xub); INT = e^{uh(xub)} (xub-z)
            intexpuh += prexp(( GetUHatUpper() + uhz ) * 0.5 - GetMaxUHull()) * 
                      ( GetUpperX() - z );
        else {
            diffuhz = uhz - GetUHatUpper();
            if( gpsc_MAXEPOW < diffuhz )
                //uh(z) >> uh(xub)
                intexpuh -= prexp( uhz - GetMaxUHull()) / hdx1;
            else
                //INT_{z}^{xub} e^{u_k(x)} dx = e^{uh(xub)} (1 - e^{uh(z)-uh(xub)}) / uh'(xub)
                intexpuh += prexp( GetUHatUpper() - GetMaxUHull()) * 
                      ( 1.0 - prexp( diffuhz )) / hdx1;
        }
    }
    else {
        if( gpsc_EPSILON <= fabs( hdx1 ))
            //INT_{z}^{inf} e^{u_k(x)} dx = INT e^{uh(z)+(x-z)uh'(x1)}
            //  = e^{uh(z)} e^{-z uh'(x1)} 1/uh'(x1) (e^{-inf uh'(x1)}-e^{z uh'(x1)})
            //  = -e^{uh(z)} / uh'(x1)
            intexpuh -= prexp( uhz - GetMaxUHull()) / hdx1;
    }

    GetEUHCDF()->Push( intexpuh );
    SetEUHArea( intexpuh );

    //normalize cdf
    for( i = 0, cdf = 0.0; i < GetEUHCDF()->GetSize(); i++ ) {
        cdf = GetEUHCDF()->GetValueAt( i );
        if( cdf < 0 )
            throw myruntime_error( "CARS: UpdateEUHCDF: Negative CDF value.");
        if( 0 < cdf ) {
            cdf /= intexpuh;
            GetEUHCDF()->SetValueAt( i, cdf );
        }
    }

    if( cdf < 1.0 - accuracy || 1.0 + accuracy < cdf ) {
        sprintf( strbuf, "CARS: UpdateEUHCDF: Not conserved CDF: %.6g", cdf );
        throw myruntime_error( strbuf );
    }
    if( GetEUHCDF()->GetSize() && cdf < 1.0 )
        GetEUHCDF()->SetValueAt( GetEUHCDF()->GetSize()-1, 1.0 );
}
