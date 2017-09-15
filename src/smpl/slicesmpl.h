/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __slicesmpl__
#define __slicesmpl__

#include <math.h>
#include "ext/pslcodes.h"
#include "ext/pslvector.h"
#include "myexcept.h"

// #ifndef SliceSmplTESTPRINT
// #define SliceSmplTESTPRINT
// #endif

static const double gscd_SLICE_XEPS = 1.e-35;   // critical relative x-value difference

static const double gscd_SLICE_YEPS = 0.1;       // critical y-value difference
static const double gscd_SLICE_EYEPS = 0.001;    // critical relative exp(y) difference
static const double gscd_SLICE_YCEIL = 50.0;     // maximum exponent

// -------------------------------------------------------------------------

//NOTE: function should evaluate log f(x) at x
typedef int ( *TSSLogEvaluator )( double x, double* lfx, void* );

// -------------------------------------------------------------------------
// The implementation of Slice Sampling;
// Neal. (2003) Annals of Statistics 31(3), 705-67.
//
class SliceSmpl
{
public:
    SliceSmpl();
    ~SliceSmpl();

    int         GetNoReqSamples() const { return noss_; }
    void        SetNoReqSamples( int value ) { noss_ = value; }

    //NOTE: function should evaluate log f(x) at x
    TSSLogEvaluator GetFEval() const { return feval_; }
    void            SetFEval( TSSLogEvaluator faddr )  { feval_ = faddr; }

    void*       GetParams() const { return feparams_; }
    void        SetParams( void* value ) { feparams_ = value; }

    double      GetW() const { return w_; }
    void        SetW( double value ) { w_ = value; }

    int         GetM() const { return m_; }
    void        SetM( int value ) { m_ = value; mset_ = true; }
    bool        IsMSet() const { return mset_; }
    void        ResetM() { mset_ = false; }

    double      GetLeftX() const { return xlb_; }
    void        SetLeftX( double value ) { xlb_ = value; xlbset_ = true; }
    double      GetRightX() const { return xrb_; }
    void        SetRightX( double value ) { xrb_ = value; xrbset_ = true; }

    bool        IsLeftXSet() const { return xlbset_; }
    void        ResetLeftX() { xlbset_ = false; }
    bool        IsRightXSet() const { return xrbset_; }
    void        ResetRightX() { xrbset_ = false; }

    double      GetX0() const { return xinit_; }
    void        SetX0( double value ) { xinit_ = value; }

    int         GetThin() const { return thin_; }
    void        SetThin( int value ) { thin_ = value; }

    int         GetBurnin() const { return burnin_; }
    void        SetBurnin( int value ) { burnin_ = value; }

    int         Run();

    const Pslvector*    GetLastSamples() const { return accpntscur_; }
    const Pslvector*    GetAllSamples() const { return accpnts_; }

    int             GetNoFevals() const { return nofevals_; }
    int             GetNoCurFevals() const { return nofevalscur_; }

protected:
    int     Sample( double x0, double lfx0, double* x1, double* lfx1 );

    int     InitAccpnts();
    void    DestroyAccpnts() { if( accpnts_ ) delete accpnts_; accpnts_ = NULL; }

    int     InitCurAccpnts();
    void    DestroyCurAccpnts() { if( accpntscur_ ) delete accpntscur_; accpntscur_ = NULL; }

    static int  expshift( double y, double y0, double* result );
    static int  logshift( double y, double y0, double* result );

    static unsigned int GetDefPntVecSize() { return s_defszpnts_; }
    static unsigned int GetDefCurPntVecSize() { return s_defszpntscur_; }

    double          GetLFvalAtX0() const { return lfxinit_; }
    void            SetLFvalAtX0( double value ) { lfxinit_ = value; }

    Pslvector*      GetLastSamples() { return accpntscur_; }
    Pslvector*      GetAllSamples() { return accpnts_; }

    void            SetNoFevals( int value ) { nofevals_ = value; }
    void            IncNoFevals() { nofevals_++; }

    void            SetNoCurFevals( int value ) { nofevalscur_ = value; }
    void            IncNoCurFevals() { nofevalscur_++; }

    bool            GetRun() const { return run_; }
    void            SetRun( bool value ) { run_ = value; }

protected:
    TSSLogEvaluator feval_;     //evaluator of log f(x)
    void*           feparams_;  //parameters for evaluator

    double          w_; //width of steps for creating interval in 'stepping out' procedure
    int             m_; //limit on number of steps in 'stepping out' procedure
    bool            mset_; //limit on number of steps set

    bool            xlbset_; //left bound of x set
    bool            xrbset_; //right bound of x set
    double          xlb_; //left bound of x
    double          xrb_; //right bound of x
    double          xinit_; //initial x point
    double          lfxinit_; //log f(x) at initial x point

    int             thin_; //generate each serial `thin' point (others are omiited in generated sequence)
    int             burnin_; //number of burn-in iterations
    int             noss_; //number of points to sample
    Pslvector*      accpnts_; //accepted sample points
    Pslvector*      accpntscur_; //accepted sample points by last run
    int             nofevals_; //number of function evaluations performed totally
    int             nofevalscur_; //number of function evaluations performed by last run

    bool            run_; //recently `Run' has been called

    static const unsigned int   s_defszpnts_; //default vector size of sampled points
    static const unsigned int   s_defszpntscur_; //default vector size of currently sampled points
};

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// InitAccpnts: initialize vector of sampled points
//
inline
int SliceSmpl::InitAccpnts()
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
int SliceSmpl::InitCurAccpnts()
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

// // -------------------------------------------------------------------------
// // expshift: exponentiate shifted y without underflow
// //
// inline
// int SliceSmpl::expshift( double y, double y0, double* result )
// {
//     if( result == NULL )
//         return 0;
//     *result = 0.0;
//     if( -gscd_ARMS_YCEIL < y - y0 )
//         *result = exp( y - y0 );
// //     if( -gscd_ARMS_YCEIL - gscd_ARMS_YCEIL < y - y0 )
// //         *result = exp( y - y0 + gscd_ARMS_YCEIL );
//     return 0;
// }
// 
// // -------------------------------------------------------------------------
// // logshift: inverse of expshift
// //
// inline
// int SliceSmpl::logshift( double y, double y0, double* result )
// {
//     if( result == NULL )
//         return PSL_ERR_ADDRESS;
//     if( y < 0.0 )
//         return PSL_ERR_DOMAIN;
//     *result = -gscd_ARMS_YCEIL;
//     if( y )
//         *result = log( y ) + y0;
// //     *result = log( y ) + y0 - gscd_ARMS_YCEIL;
//     return 0;
// }

// -------------------------------------------------------------------------
// Error codes
//
#define SliceSmpl_ERR_BININVALID    10101
#define SliceSmpl_ERR_THININVALID   10103
#define SliceSmpl_ERR_WTOOSMALL     10105
#define SliceSmpl_ERR_MINVALID      10107
#define SliceSmpl_ERR_MTOOLARGE     10109
#define SliceSmpl_ERR_FEVALADDRESS  10111
#define SliceSmpl_ERR_BNDSINVALID   10113
#define SliceSmpl_ERR_INITXDOMAIN   10115
#define SliceSmpl_ERR_STEPOUTITER   10117
#define SliceSmpl_ERR_SHRINKINITER  10119

// -------------------------------------------------------------------------
//
inline const char* TranslateSliceSmplError( int code )
{
    switch( code ) {
        case SliceSmpl_ERR_BININVALID:      return "Invalid burn-in number of iterations";
        case SliceSmpl_ERR_THININVALID:     return "Invalid thinness parameter";
        case SliceSmpl_ERR_WTOOSMALL:       return "Too small Width parameter w";
        case SliceSmpl_ERR_MINVALID:        return "Invalid step limit m";
        case SliceSmpl_ERR_MTOOLARGE:       return "Too large step limit m";
        case SliceSmpl_ERR_FEVALADDRESS:    return "No function evaluator";
        case SliceSmpl_ERR_BNDSINVALID:     return "Invalid function domain (bounds)";
        case SliceSmpl_ERR_INITXDOMAIN:     return "Initial point is out of domain";
        case SliceSmpl_ERR_STEPOUTITER:     return "Too many iterations in stepping-out procedure";
        case SliceSmpl_ERR_SHRINKINITER:    return "Too many iterations in shrinking-in procedure";
    }
    return TranslatePSLError( code );
}

#endif//__slicesmpl__
