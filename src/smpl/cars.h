/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __cars__
#define __cars__

#include <math.h>
#include "myexcept.h"
#include "BinarySearchStructure.h"
#include "ext/pslvector.h"

// #ifndef CARSTESTPRINT
// #define CARSTESTPRINT
// #endif

static const double gpsc_SMALLPW = -65536.0; //extremely small expoenent

static const double gpsc_MINEPOW = -64.0; //minimum exponent
static const double gpsc_MAXEPOW = 64.0; //maximum exponent
static const double gpsc_EPSILON = exp( gpsc_MINEPOW ); //the lowest positive value used

//IMPORTANT: the function should evaluate log f(x) and d/dx log f(x) at x
typedef int ( *THEvaluator )( double x, double* hx, double* hdx, void* );

// -------------------------------------------------------------------------
// The implementation of Adaptive Rejection Sampling, adapted from the 
// original Gilks code; Gilks & Wild. (1992) Applied Statistics 41, 337-348.
//
class CARS
{
public:
    CARS();
    ~CARS();

    int         GetNoIterats() const { return noss_; }
    void        SetNoIterats( int value ) { noss_ = value; }

    //IMPORTANT: the function should evaluate log f(x) and d/dx log f(x) at x
    THEvaluator GetHEval() const { return heval_; }
    void        SetHEval( THEvaluator faddr )  { heval_ = faddr; }

    void*       GetParams() const { return heparams_; }
    void        SetParams( void* value ) { heparams_ = value; }

    void        PushInitialX( double xi ) { NewX( xi, true ); }

    double      GetLowerX() const { return xlb_; }
    void        SetLowerX( double value ) { xlb_ = value; SetLowerSet( true ); }
    double      GetUpperX() const { return xub_; }
    void        SetUpperX( double value ) { xub_ = value; SetUpperSet( true ); }

    double      GetUHatLower() const { return uhlb_; }
    double      GetUHatUpper() const { return uhub_; }
    double      GetMaxUHull() const { return maxuh_; }

    bool        GetLowerSet() const { return lowb_; }
    void        SetLowerSet( bool value ) { lowb_ = value; }
    bool        GetUpperSet() const { return uppb_; }
    void        SetUpperSet( bool value ) { uppb_ = value; }

    int         Run();

    const Pslvector*    GetResPoints() const { return accpnts_; }
    const Pslvector*    GetX() const { return xv_; }
    const Pslvector*    GetHX() const { return hxv_; }
    const Pslvector*    GetHDX() const { return hdxv_; }
    const BinarySearchStructure* GetSrtX() const { return srtxv_; }
    const Pslvector*    GetZ() const { return zv_; }
    const Pslvector*    GetUHZ() const { return uhzv_; }
    const Pslvector*    GetEUHCDF() const { return euhcdf_; }

protected:
    Pslvector*      GetResPoints() { return accpnts_; }
    Pslvector*      GetX() { return xv_; }
    Pslvector*      GetHX() { return hxv_; }
    Pslvector*      GetHDX() { return hdxv_; }
    BinarySearchStructure*  GetSrtX() { return srtxv_; }
    Pslvector*      GetZ() { return zv_; }
    Pslvector*      GetUHZ() { return uhzv_; }
    Pslvector*      GetEUHCDF() { return euhcdf_; }

    void            NewX( double xi, bool check, int* = NULL, double* = NULL, double* = NULL );

    double          GetEUHArea() const { return euharea_; }
    void            SetEUHArea( double value ) { euharea_ = value; }

    void            SetUHatLower( double value ) { uhlb_ = value; }
    void            SetUHatUpper( double value ) { uhub_ = value; }
    void            SetMaxUHull( double value ) { maxuh_ = value; }

    void            Initialize();
    void            Sample();
    void            SampleUH( double* smpl, int* in );
    void            Update( int n );
    void            UpdateUH( int n, bool insert = true );
    void            UpdateEUHCDF();
    void            Intersection( int n1, int n2, double* pz, double* puhz );

    void            InitVectors();
    void            ClearVectors();
    void            DestroyVectors();

    bool            GetRun() const { return run_; }
    void            SetRun( bool value ) { run_ = value; }

    void            CheckVectors() const;

    static int      GetDefVecSize() { return s_defsz_; }
    static double   prexp( double arg );

protected:
    int             noss_; //number of points to sample
    THEvaluator     heval_; //evaluator of hx and hdx
    void*           heparams_; //parameters to be passed to evaluator
    Pslvector*      accpnts_; //accepted sample points
    Pslvector*      xv_; //argument vector (x-values)
    Pslvector*      hxv_; //vector of values of ln f(x) (paper notation)
    Pslvector*      hdxv_; //vector of derivative values: d/dx ln f(x)
    BinarySearchStructure* srtxv_; //indices of sorted x_
    double          xlb_; //lower bound of x if given
    double          xub_; //upper bound of x if given
    double          uhlb_; //upper hull value at lower bound of x
    double          uhub_; //upper hull value at upper bound of x
    Pslvector*      zv_; //vector of upper hull intersection points
    Pslvector*      uhzv_; //values of upper hull at intersection points
    double          maxuh_; //maximum value of upper hull
    Pslvector*      euhcdf_; //cdf of the exponentiated upper hull
    double          euharea_; //area of the exponentiated upper hull
    bool            lowb_; //whether lower bound is set
    bool            uppb_; //whether upper bound is set
    bool            run_; //recently `Run' has been called
    static const int    s_defsz_; //default size of vectors
};

// -------------------------------------------------------------------------
// BSSCompare: compare elements of argument vector `x_'
//
inline
int BSSCompare( const void* key1, const void* key2, void* pars )
{
    if( pars == NULL )
        throw myruntime_error("CARS: BSSCompare: Memory access error.");

    const CARS* pthis = ( const CARS* )pars;
    int ndx1 = ( int )( size_t )key1;
    int ndx2 = ( int )( size_t )key2;

    if( pthis->GetX() == NULL )
        throw myruntime_error("CARS: BSSCompare: Memory access error.");

    if( ndx1 < 0 || pthis->GetX()->GetSize() <= ndx1 ||
        ndx2 < 0 || pthis->GetX()->GetSize() <= ndx2 )
        throw myruntime_error("CARS: BSSCompare: Memory access error.");

    double  diff = pthis->GetX()->GetValueAt( ndx1 ) - pthis->GetX()->GetValueAt( ndx2 );
    return ( diff < 0.0 )? -1: (( 0.0 < diff )? 1: 0 );
}

// -------------------------------------------------------------------------
// prexp: private exp function
//
inline
double CARS::prexp( double arg )
{
    if( arg < gpsc_MINEPOW )
        return 0.0;
    return exp( arg );
}

// -------------------------------------------------------------------------
// CheckVectors: test vectors for errors
//
inline
void CARS::CheckVectors() const
{
    if( !GetResPoints() ||
        !GetX() || !GetHX() || !GetHDX() ||
        !GetSrtX() || !GetZ() || !GetUHZ() || !GetEUHCDF() ||
        !GetHEval())
        throw myruntime_error("CARS: CheckVectors: Null vectors.");
    if( GetX()->GetSize() != GetHX()->GetSize() || 
        GetHDX()->GetSize() != GetSrtX()->GetSize() ||
        GetHDX()->GetSize() != GetX()->GetSize() ||
        GetX()->GetSize() < GetZ()->GetSize() || 
        GetZ()->GetSize() != GetUHZ()->GetSize())
        throw myruntime_error("CARS: CheckVectors: Inconsistent vectors.");
    if( GetX()->GetSize() < 1 )
        throw myruntime_error("CARS: No initial arguments (>=1).");
}

// -------------------------------------------------------------------------
// InitVectors: initialize member vectors
//
inline
void CARS::InitVectors()
{
    DestroyVectors();
    accpnts_ = new Pslvector();
    xv_ = new Pslvector();
    hxv_ = new Pslvector();
    hdxv_ = new Pslvector();
    srtxv_ = new BinarySearchStructure( BSSCompare, GetDefVecSize(), true/*dupl.*/, this );
    zv_ = new Pslvector();
    uhzv_ = new Pslvector();
    euhcdf_ = new Pslvector();

    if( accpnts_ == NULL ||
        xv_ == NULL || hxv_ == NULL || hdxv_ == NULL ||
        srtxv_ == NULL || zv_ == NULL || uhzv_ == NULL || euhcdf_ == NULL )
        throw myruntime_error("CARS: Not enough memory.");

    accpnts_->Allocate( GetDefVecSize());
    xv_->Allocate( GetDefVecSize());
    hxv_->Allocate( GetDefVecSize());
    hdxv_->Allocate( GetDefVecSize());
    zv_->Allocate( GetDefVecSize());
    uhzv_->Allocate( GetDefVecSize());
    euhcdf_->Allocate( GetDefVecSize());
}

// -------------------------------------------------------------------------
// ClearVectors: clear member vectors
//
inline
void CARS::ClearVectors()
{
    if( accpnts_ == NULL ||
        xv_ == NULL || hxv_ == NULL || hdxv_ == NULL ||
        srtxv_ == NULL || zv_ == NULL || uhzv_ == NULL || euhcdf_ == NULL )
        throw myruntime_error("CARS: Memory access error.");
    accpnts_->Clear();
    xv_->Clear();
    hxv_->Clear();
    hdxv_->Clear();
    srtxv_->Clear();
    zv_->Clear();
    uhzv_->Clear();
    euhcdf_->Clear();
}

// -------------------------------------------------------------------------
// DestroyVectors: destroy member vectors
//
inline
void CARS::DestroyVectors()
{
    if( accpnts_ ) delete accpnts_; accpnts_ = NULL;
    if( xv_ ) delete xv_; xv_ = NULL;
    if( hxv_ ) delete hxv_; hxv_ = NULL;
    if( hdxv_ ) delete hdxv_; hdxv_ = NULL;
    if( srtxv_ ) delete srtxv_; srtxv_ = NULL;
    if( zv_ ) delete zv_; zv_ = NULL;
    if( uhzv_ ) delete uhzv_; uhzv_ = NULL;
    if( euhcdf_ ) delete euhcdf_; euhcdf_ = NULL;
}

// -------------------------------------------------------------------------
// NewX: push new value of x; clear all values in vectors if
//  recently `Run' has been called
//
inline
void CARS::NewX( double xi, bool check, int* ploc, double* phxi, double* phdxi )
{
    if( !GetX() || !GetHX() || !GetHDX() ||
        !GetSrtX() || !GetHEval())
        throw myruntime_error("CARS: Memory access error.");
    if( check && GetRun()) {
        ClearVectors();
        SetRun( false );
    }
    double  hx, hdx;
    int len = GetX()->GetSize();
    int loc;
    int err;

    GetX()->Push( xi );
    GetSrtX()->Push(( const void* )len, &loc );
    if(( err = ( *GetHEval())( xi, &hx, &hdx, GetParams())) != 0 )
        throw myruntime_error("CARS: Evaluation of function values failed.");
    GetHX()->Push( hx );
    GetHDX()->Push( hdx );

    if( ploc ) *ploc = loc;
    if( phxi ) *phxi = hx;
    if( phdxi ) *phdxi = hdx;
}

#endif//__cars__
