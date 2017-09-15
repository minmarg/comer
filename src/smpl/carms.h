/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __carms__
#define __carms__

#include <math.h>
#include "ext/pslcodes.h"
#include "ext/pslvector.h"
#include "BinarySearchStructure.h"
#include "myexcept.h"

// #ifndef CARMSTESTPRINT
// #define CARMSTESTPRINT
// #endif

static const double gscd_ARMS_XEPS = 1.e-5;     // critical relative x-value difference
static const double gscd_ARMS_YEPS = 0.1;       // critical y-value difference
static const double gscd_ARMS_EYEPS = 0.001;    // critical relative exp(y) difference
static const double gscd_ARMS_YCEIL = 50.0;     // maximum exponent

typedef struct SARMSpoint_ {    // point in the x,y plane
  double    x, y;               // x and y coordinates
  double    ey;                 // exp(y-ymax+YCEIL)
  double    cdf;                // integral up to x of rejection envelope
  int       f;                  // is y an evaluated point of log-density (or intersection)
  struct SARMSpoint_* pl, *pr;  // envelope points to left and right of x
} SARMSpoint;

// -------------------------------------------------------------------------

//NOTE: function should evaluate log f(x) at x
typedef int ( *TARMSEvaluator )( double x, double* fx, void* );

// -------------------------------------------------------------------------
// The implementation of Adaptive Rejection Metropolis Sampling;
// adapted from the original Gilks code; 
// Gilks, Best & Tan. (1992) Applied Statistics 44, 455-72.
//
class CARMS
{
public:
    CARMS();
    ~CARMS();

    int         GetNoReqSamples() const { return noss_; }
    void        SetNoReqSamples( int value ) { noss_ = value; }

    //NOTE: function should evaluate log f(x) at x
    TARMSEvaluator  GetFEval() const { return feval_; }
    void            SetFEval( TARMSEvaluator faddr )  { feval_ = faddr; }

    void*       GetParams() const { return feparams_; }
    void        SetParams( void* value ) { feparams_ = value; }

    double      GetLeftX() const { return xlb_; }
    void        SetLeftX( double value ) { xlb_ = value; }
    double      GetRightX() const { return xrb_; }
    void        SetRightX( double value ) { xrb_ = value; }

    int         PushInitialX( double xi );

    double      GetAdjCnvx() const { return adjconv_; }
    void        SetAdjCnvx( double value ) { adjconv_ = value; }

    int         GetMaxNoEnvPoints() const { return maxnoenvpnts_; }
    void        SetMaxNoEnvPoints( int value ) { maxnoenvpnts_ = value; }

    bool        GetDoMHupdate() const { return meton_; }
    void        SetDoMHupdate( bool value ) { meton_ = value; }

    double      GetPrevMCvalue() const { return metxprev_; }
    void        SetPrevMCvalue( double value ) { metxprev_ = value; }

    double      GetFvalAtPrevMCvalue() const { return metyprev_; }
    void        SetFvalAtPrevMCvalue( double value ) { metyprev_ = value; }

    int         Run();

    const Pslvector*    GetInitialXs() const { return xinit_; }

    const Pslvector*    GetPerCents() const { return percents_; }
    void                SetPerCents( const Pslvector* value ) { percents_ = value; }

    const Pslvector*    GetCentiles() const { return centiles_; }

    const Pslvector*    GetLastSamples() const { return accpntscur_; }
    const Pslvector*    GetAllSamples() const { return accpnts_; }

    int             GetNoFevals() const { return nofevals_; }
    int             GetNoCurFevals() const { return nofevalscur_; }

protected:
    int     Initialize();
    int     Sample( SARMSpoint* p );
    int     Invert( double prob, SARMSpoint* p );
    int     Test( SARMSpoint* p, int* );
    int     Update( SARMSpoint* p );
    int     CalcCDF();
    int     Intersection( SARMSpoint* q );
    int     Area( SARMSpoint* q, double* area );

    int     InitEnvelope();
    void    DestroyEnvelope() { if( envelope_ ) free( envelope_ ); envelope_ = NULL; }

    int     InitXinit();
    void    DestroyXinit(){ if( xinit_ ) delete xinit_; xinit_ = NULL; 
                            if( xinitpos_ ) delete xinitpos_; xinitpos_ = NULL; 
    }
    int     InitAccpnts();
    void    DestroyAccpnts() { if( accpnts_ ) delete accpnts_; accpnts_ = NULL; }

    int     InitCurAccpnts();
    void    DestroyCurAccpnts() { if( accpntscur_ ) delete accpntscur_; accpntscur_ = NULL; }

    int     InitCentiles( int size );
    void    DestroyCentiles() { if( centiles_ ) delete centiles_; centiles_ = NULL; }

    static int  expshift( double y, double y0, double* result );
    static int  logshift( double y, double y0, double* result );

    static int  XvalCompare( const void* key1, const void* key2, void* pars );

    static unsigned int GetDefInitPntVecSize() { return s_definitpnts_; }
    static unsigned int GetDefPntVecSize() { return s_defszpnts_; }
    static unsigned int GetDefCurPntVecSize() { return s_defszpntscur_; }

    const SARMSpoint*   GetEnvelope() const { return envelope_; }
    SARMSpoint*         GetEnvelope() { return envelope_; }

    int             GetNoEnvPoints() const { return noenvpnts_; }
    void            SetNoEnvPoints( int value ) { noenvpnts_ = value; }
    void            IncNoEnvPoints() { noenvpnts_++; }

    double          GetEnvMax() const { return envmax_; }
    void            SetEnvMax( double value ) { envmax_ = value; }

    int             GenInitialXs();
    int             GetInitialXat( int n, double* value );

    Pslvector*      GetInitialXs() { return xinit_; }
    const BinarySearchStructure* GetInitialXposs() const { return xinitpos_; }
    BinarySearchStructure* GetInitialXposs() { return xinitpos_; }
    Pslvector*      GetLastSamples() { return accpntscur_; }
    Pslvector*      GetAllSamples() { return accpnts_; }
    Pslvector*      GetCentiles() { return centiles_; }

    void            SetNoFevals( int value ) { nofevals_ = value; }
    void            IncNoFevals() { nofevals_++; }

    void            SetNoCurFevals( int value ) { nofevalscur_ = value; }
    void            IncNoCurFevals() { nofevalscur_++; }

    bool            GetRun() const { return run_; }
    void            SetRun( bool value ) { run_ = value; }

protected:
    bool    meton_;     // whether metropolis is to be used
    double  metxprev_;  // previous Markov chain iterate
    double  metyprev_;  // current log density at metxprev_

    TARMSEvaluator  feval_;     //evaluator of fx
    void*           feparams_;  //parameters for evaluator

    SARMSpoint*     envelope_; // rejection envelope
    int             noenvpnts_; // number of points in current envelope
    int             maxnoenvpnts_; //maximum number of envelope points
    double          envmax_; // the maximum y-value in the current envelope
    double          adjconv_; //adjustment for convexity

    double          xlb_; //left bound of x
    double          xrb_; //right bound of x
    Pslvector*      xinit_; //initial x points
    BinarySearchStructure* xinitpos_; //positions of initial x points

    int             noss_; //number of points to sample
    Pslvector*      accpnts_; //accepted sample points
    Pslvector*      accpntscur_; //accepted sample points by last run
    int             nofevals_; //number of function evaluations performed totally
    int             nofevalscur_; //number of function evaluations performed by last run

    const Pslvector* percents_; //percentages for envelope centiles
    Pslvector*      centiles_; //requested centiles

    bool            run_; //recently `Run' has been called

    static const unsigned int   s_definitpnts_; //default vector size of initial x points
    static const unsigned int   s_defszpnts_; //default vector size of sampled points
    static const unsigned int   s_defszpntscur_; //default vector size of currently sampled points
};

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// InitEnvelope: initialize vector of envelope points
//
inline
int CARMS::InitEnvelope()
{
    DestroyEnvelope();
    envelope_ = ( SARMSpoint* )malloc( GetMaxNoEnvPoints() * sizeof( SARMSpoint ));
    if( envelope_ == NULL )
        return PSL_ERR_MEMORY;
    memset( envelope_, 0, GetMaxNoEnvPoints() * sizeof( SARMSpoint ));
    return 0;
}

// -------------------------------------------------------------------------
// InitXinit: initialize vector of initial x points
//
inline
int CARMS::InitXinit()
{
    DestroyXinit();
    xinit_ = new Pslvector();
    xinitpos_ = new BinarySearchStructure( XvalCompare, GetDefInitPntVecSize(), false/*no dups*/, this );
    if( xinit_ == NULL || xinitpos_ == NULL )
        return PSL_ERR_MEMORY;
    try {
        xinit_->Allocate( GetDefInitPntVecSize());
    } catch( myexception const& ex ) {
        return PSL_ERR_MEMORY;
    }
    return 0;
}

// XvalCompare: static comparator of x points (for sorted array of values)
inline
int CARMS::XvalCompare( const void* key1, const void* key2, void* pars )
{
    if( pars == NULL )
        throw myruntime_error("CARMS: XvalCompare: Memory access error.");

    const CARMS* pthis = ( const CARMS* )pars;
    int ndx1 = ( int )( size_t )key1;
    int ndx2 = ( int )( size_t )key2;

    if( pthis->GetInitialXs() == NULL )
        throw myruntime_error("CARMS: XvalCompare: Memory access error.");

    if( ndx1 < 0 || pthis->GetInitialXs()->GetSize() <= ndx1 ||
        ndx2 < 0 || pthis->GetInitialXs()->GetSize() <= ndx2 )
        throw myruntime_error("CARMS: XvalCompare: Memory access error.");

    double  diff = pthis->GetInitialXs()->GetValueAt( ndx1 ) - pthis->GetInitialXs()->GetValueAt( ndx2 );
    return ( diff < 0.0 )? -1: (( 0.0 < diff )? 1: 0 );
}

// -------------------------------------------------------------------------
// InitAccpnts: initialize vector of sampled points
//
inline
int CARMS::InitAccpnts()
{
    DestroyAccpnts();
    accpnts_ = new Pslvector();
    if( accpnts_ == NULL )
        return PSL_ERR_MEMORY;
    try {
        accpnts_->Allocate( GetDefPntVecSize());
    } catch( myexception const& ex ) {
        return PSL_ERR_MEMORY;
    }
    return 0;
}

// -------------------------------------------------------------------------
// InitCurAccpnts: initialize vector of currently sampled points
//
inline
int CARMS::InitCurAccpnts()
{
    DestroyCurAccpnts();
    accpntscur_ = new Pslvector();
    if( accpntscur_ == NULL )
        return PSL_ERR_MEMORY;
    try {
        accpntscur_->Allocate( GetDefCurPntVecSize());
    } catch( myexception const& ex ) {
        return PSL_ERR_MEMORY;
    }
    return 0;
}

// -------------------------------------------------------------------------
// InitCentiles: initialize centiles vector
//
inline
int CARMS::InitCentiles( int size )
{
    if( size <= 0 )
        return PSL_ERR_ADDRESS;
    DestroyCentiles();
    centiles_ = new Pslvector();
    if( centiles_ == NULL )
        return PSL_ERR_MEMORY;
    try {
        centiles_->Allocate( size );
    } catch( myexception const& ex ) {
        return PSL_ERR_MEMORY;
    }
    return 0;
}

// -------------------------------------------------------------------------
// GetInitialXat: get initial point value
//
inline
int CARMS::GetInitialXat( int n, double* value )
{
    if( value == NULL )
        return PSL_ERR_ADDRESS;
    if( xinit_ == NULL || xinitpos_ == NULL ||
        xinit_->GetSize() != xinitpos_->GetSize() ||
        n < 0 || xinit_->GetSize() <= n )
        return PSL_ERR_ADDRESS;
    int loc = ( int )( size_t )xinitpos_->GetValueAt( n );
    if( loc < 0 || xinit_->GetSize() <= loc )
        return PSL_ERR_ADDRESS;
    *value = xinit_->GetValueAt( loc );
    return 0;
}

// -------------------------------------------------------------------------
// expshift: exponentiate shifted y without underflow
//
inline
int CARMS::expshift( double y, double y0, double* result )
{
    if( result == NULL )
        return 0;
    *result = 0.0;
    if( -gscd_ARMS_YCEIL < y - y0 )
        *result = exp( y - y0 );
//     if( -gscd_ARMS_YCEIL - gscd_ARMS_YCEIL < y - y0 )
//         *result = exp( y - y0 + gscd_ARMS_YCEIL );
    return 0;
}

// -------------------------------------------------------------------------
// logshift: inverse of expshift
//
inline
int CARMS::logshift( double y, double y0, double* result )
{
    if( result == NULL )
        return PSL_ERR_ADDRESS;
    if( y < 0.0 )
        return PSL_ERR_DOMAIN;
    *result = -gscd_ARMS_YCEIL;
    if( y )
        *result = log( y ) + y0;
//     *result = log( y ) + y0 - gscd_ARMS_YCEIL;
    return 0;
}

// -------------------------------------------------------------------------
// Error codes
//
#define CARMS_ERR_FEVALADDRESS  10001
#define CARMS_ERR_PERCDOMAIN    10003
#define CARMS_ERR_MCVALDOMAIN   10005
#define CARMS_ERR_BNDSINVALID   10007
#define CARMS_ERR_INITXDOMAIN   10009
#define CARMS_ERR_INITXINSUFF   10011
#define CARMS_ERR_INITXTOOMANY  10013
#define CARMS_ERR_INITXSET      10015
#define CARMS_ERR_NEGCNVX       10017
#define CARMS_ERR_CNVXVLTED     10019

// -------------------------------------------------------------------------
//
inline const char* TranslateCARMSError( int code )
{
    switch( code ) {
        case CARMS_ERR_FEVALADDRESS:    return "No function evaluator";
        case CARMS_ERR_PERCDOMAIN:      return "Percentage requesting centile is out of range";
        case CARMS_ERR_MCVALDOMAIN:     return "Previous markov chain iterate out of range";
        case CARMS_ERR_BNDSINVALID:     return "Invalid function domain (bounds)";
        case CARMS_ERR_INITXDOMAIN:     return "Initial point is out of domain";
        case CARMS_ERR_INITXINSUFF:     return "Insufficient number of initial points";
        case CARMS_ERR_INITXTOOMANY:    return "Number of initial points exceeds limit for envelope";
        case CARMS_ERR_INITXSET:        return "Initial points set already";
        case CARMS_ERR_NEGCNVX:         return "Negative convexity parameter";
        case CARMS_ERR_CNVXVLTED:       return "Not a convex density in interval supplied";
    }
    return TranslatePSLError( code );
}

#endif//__carms__
