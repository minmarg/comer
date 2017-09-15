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
#include "carms.h"

//FILE GLOBALS
//random number generators
MTRng   _carmsRNGss_;
MTRng   _carmsRNGtst_;

//default size of vectors
const unsigned int CARMS::s_definitpnts_ = 10; //default vector size of initial x points
const unsigned int CARMS::s_defszpnts_ = 1000; //default vector size of sampled points
const unsigned int CARMS::s_defszpntscur_ = 100; //default vector size of currently sampled points

// -------------------------------------------------------------------------
// constructor
//
CARMS::CARMS()
:   meton_( true ),
    metxprev_( 0.0 ),
    metyprev_( 0.0 ),

    feval_( NULL ),
    feparams_( NULL ),

    envelope_( NULL ),
    noenvpnts_( 0 ),
    maxnoenvpnts_( 100 ),
    envmax_( 0.0 ),
    adjconv_( 1.0 ),

    xlb_( 0.0 ),
    xrb_( 0.0 ),
    xinit_( NULL ),
    xinitpos_( NULL ),

    noss_( 1 ),
    accpnts_( NULL ),
    accpntscur_( NULL ),
    nofevals_( 0 ),
    nofevalscur_( 0 ),

    percents_( NULL ),
    centiles_( NULL ),

    run_( false )
{
    time_t tm;
    clock_t cl = clock();
    time( &tm );
    _carmsRNGss_.Set((unsigned long)(size_t)this +(unsigned long)(size_t)(&_carmsRNGss_) +
            (unsigned long)tm +(unsigned long)cl);
    _carmsRNGtst_.Set((unsigned long)(size_t)this +(unsigned long)(size_t)(&_carmsRNGtst_) +
            (unsigned long)tm +(unsigned long)cl);
    InitEnvelope();
    InitXinit();
    InitAccpnts();
    InitCurAccpnts();
}

// destructor
//
CARMS::~CARMS()
{
    DestroyCurAccpnts();
    DestroyAccpnts();
    DestroyXinit();
    DestroyEnvelope();
    DestroyCentiles();
}

// =========================================================================
// Run: start sampling from f(x)
//
int CARMS::Run()
{
    SARMSpoint  pwork; //working point, not yet incorporated in envelope
    double  value;
    int ns = 0; //number of values currently sampled
    int n, err;

    const Pslvector* prcents = GetPerCents();

    if( GetNoReqSamples() <= 0 )
        return 0;

    if( GetAllSamples() == NULL || GetLastSamples() == NULL )
        return PSL_ERR_ADDRESS;

    if( GetFEval() == NULL )
        return CARMS_ERR_FEVALADDRESS;

    GetLastSamples()->Clear();
    //initialise current number of function evaluations
    SetNoCurFevals( 0 );

    if( prcents ) {
        //check requested envelope centiles
        for( n = 0; n < prcents->GetSize(); n++ ) {
            value = prcents->GetValueAt( n );
            if( value < 0.0 || 100.0 < value )
                //percentage requesting centile is out of range
                return CARMS_ERR_PERCDOMAIN;
        }
        if(( err = InitCentiles( prcents->GetSize())) != PSL_OK )
            return err;
    }

    //set up initial envelope
    if(( err = Initialize()) != PSL_OK )
        return err;

    //previous MC point
    if( GetDoMHupdate()) {
        if( GetPrevMCvalue() < GetLeftX() || GetRightX() < GetPrevMCvalue()) {
            //previous markov chain iterate out of range
            return CARMS_ERR_MCVALDOMAIN;
        }
        if(( err = ( *GetFEval())( GetPrevMCvalue(), &value, GetParams())) != PSL_OK )
            return err;
        SetFvalAtPrevMCvalue( value );
    }

    //do adaptive rejection
    for( ns = 0; ns < GetNoReqSamples(); ) {
        //sample a new point
        if(( err = Sample( &pwork )) != PSL_OK )
            return err;

        //perform rejection (and perhaps metropolis) tests
        if(( err = Test( &pwork, &n )) != PSL_OK )
            return err;
        if( n == 1 ) {
            //point accepted
            GetAllSamples()->Push( pwork.x );
            GetLastSamples()->Push( pwork.x );
            ns++;
        }
    }

    if( prcents ) {
        //calculate requested envelope centiles
        for( n = 0; n < prcents->GetSize(); n++ ) {
            Invert( prcents->GetValueAt(n)/100.0, &pwork );
            if( GetCentiles())
                GetCentiles()->Push( pwork.x );
        }
    }

    return 0;
}

// -------------------------------------------------------------------------
// PushInitialX: save initial point; at least three of them are required
//
int CARMS::PushInitialX( double xi )
{
    if( GetInitialXs() == NULL || GetInitialXposs() == NULL )
        return PSL_ERR_ADDRESS;

    if( GetRightX() <= GetLeftX())
        return CARMS_ERR_BNDSINVALID;

    if( xi <= GetLeftX() || GetRightX() <= xi )
        return CARMS_ERR_INITXDOMAIN;

    int len = GetInitialXs()->GetSize();

    if(( GetMaxNoEnvPoints() >> 1 ) <= len )
        return CARMS_ERR_INITXTOOMANY;

    if( GetRun())
        return CARMS_ERR_INITXSET;

    try {
        GetInitialXs()->Push( xi );
        if( !GetInitialXposs()->Push(( const void* )len ))
            GetInitialXs()->DecDim();
    } catch( myexception const& ex ) {
        return PSL_ERR_ADDRESS;
    }

    return 0;
}

// -------------------------------------------------------------------------
// GenInitialXs: generate initial points; at least three of them required
//
int CARMS::GenInitialXs()
{
    if( GetInitialXs() == NULL || GetInitialXposs() == NULL )
        return PSL_ERR_ADDRESS;

    if( GetRightX() <= GetLeftX())
        return CARMS_ERR_BNDSINVALID;

    double  intvlen = GetRightX() - GetLeftX();
    double  value, delta = 1.0;
    int maxinpnts = GetMaxNoEnvPoints() >> 1;
    int npnts = 5;
    int n, len;

    if( maxinpnts <= npnts || ( int )intvlen <= npnts )
        npnts = 5;

    if( maxinpnts <= npnts || ( int )intvlen <= npnts )
        npnts = 3;

    delta = intvlen /( double )( npnts + 1 );

    for( n = 1; n <= npnts; n++ ) {
        value = GetLeftX() + ( double )n * delta;
        len = GetInitialXs()->GetSize();
        try {
            GetInitialXs()->Push( value );
            if( !GetInitialXposs()->Push(( const void* )len ))
                GetInitialXs()->DecDim();
        } catch( myexception const& ex ) {
            return PSL_ERR_ADDRESS;
        }
    }

    return 0;
}

// -------------------------------------------------------------------------
// Initialize: set up initial envelope
//
int CARMS::Initialize()
{
    SARMSpoint* q = NULL;
    double  xval, value;
    int n, k, nn;
    int err;

    if( GetFEval() == NULL )
        return CARMS_ERR_FEVALADDRESS;

    if( GetInitialXs() == NULL || GetInitialXposs() == NULL )
        return PSL_ERR_ADDRESS;

    if( GetRightX() <= GetLeftX())
        return CARMS_ERR_BNDSINVALID;

    if( GetInitialXs()->GetSize() < 3 )
        if(( err = GenInitialXs()) != PSL_OK )
            return err;

    if( GetInitialXs()->GetSize() < 3 )
        //too few initial points
        return CARMS_ERR_INITXINSUFF;

    nn = TIMES2( GetInitialXs()->GetSize()) + 1;
    if( GetMaxNoEnvPoints() < nn )
        //too many initial points
        return CARMS_ERR_INITXTOOMANY;

    if( GetAdjCnvx() < 0.0 )
        //negative convexity parameter
        return CARMS_ERR_NEGCNVX;

    if( GetRun())
        return 0;
    
    InitEnvelope();

    if(( q = GetEnvelope()) == NULL )
        return PSL_ERR_ADDRESS;

    //set up envelope point
    //first position is for left pnt bound; start with next log density point;
    //last position if for right pnt bound
    q->x = GetLeftX();
    q->pl = NULL;
    q->pr = q + 1;
    q++;
    for( n = 1, k = 0; n < nn-1; n++, q++ ) {
        if( q == NULL )
            return PSL_ERR_ADDRESS;
        if( n & 1 ) {
            //point on log density (otherwise intersection point)
            if(( err = GetInitialXat( k++, &xval )) != PSL_OK )
                return err;
            if(( err = ( *GetFEval())( xval, &value, GetParams())) != PSL_OK )
                return err;
            IncNoFevals();
            IncNoCurFevals();
            q->x = xval;
            q->y = value;
            q->f = 1;
        }
        q->pl = q - 1;
        q->pr = q + 1;
    }
    if( q == NULL )
        return PSL_ERR_ADDRESS;
    q->x = GetRightX();
    q->pl = q - 1;
    q->pr = NULL;

    //calculate intersection points
    q = GetEnvelope();
    for( n = 0; n < nn; n += 2, q += 2 ) {
        if(( err = Intersection( q )) != PSL_OK )
          //envelope violation without metropolis
          return err;
    }

    //exponentiate and integrate envelope
    if(( err = CalcCDF()) != PSL_OK )
        return err;

    //save number of pnts currently in envelope
    SetNoEnvPoints( nn );

    return 0;
}


// -------------------------------------------------------------------------
// Intersection: find intersection of two chords;
//     q is point to store point of intersection
//
int CARMS::Intersection( SARMSpoint* q )
{
    double  gl, gr, grl, dl, dr;
    int il, ir, irl;
    int err;

    if( q == NULL )
        return PSL_ERR_ADDRESS;
    if( q->f )
      //this is not an intersection point
      return PSL_ERR_INVALID;

    //calculate coordinates of point of intersection
    if( q->pl != NULL ) {
        if( q->pl->pl == NULL )
            return PSL_ERR_ADDRESS;
        if( q->pl->pl->pl != NULL ) {
            //chord gradient can be calculated at left end of interval
            gl = ( q->pl->y - q->pl->pl->pl->y )/( q->pl->x - q->pl->pl->pl->x );
            il = 1;
        }
    } else
        //no chord gradient on left
        il = 0;

    if( q->pr != NULL ) {
        //chord gradient can be calculated at right end of interval
        if( q->pr->pr == NULL )
            return PSL_ERR_ADDRESS;
        if( q->pr->pr->pr != NULL ) {
            gr = ( q->pr->y - q->pr->pr->pr->y )/( q->pr->x - q->pr->pr->pr->x );
            ir = 1;
        }
    } else
        //no chord gradient on right
        ir = 0;

    if( q->pl != NULL && q->pr != NULL ){
        //chord gradient can be calculated across interval
        grl = ( q->pr->y - q->pl->y )/( q->pr->x - q->pl->x );
        irl = 1;
    } else
        irl = 0;

    if( irl && il && gl < grl ) {
        //convexity on left exceeds current threshold
        if( !GetDoMHupdate())
            //envelope violation without metropolis
            return CARMS_ERR_CNVXVLTED;
        //adjust left gradient
        gl = gl + ( 1.0 + GetAdjCnvx()) * ( grl - gl );
    }

    if( irl && ir && grl < gr ) {
        //convexity on right exceeds current threshold
        if( !GetDoMHupdate()){
            //envelope violation without metropolis
            return CARMS_ERR_CNVXVLTED;
        }
        //adjust right gradient
        gr = gr + ( 1.0 + GetAdjCnvx()) * ( grl - gr );
    }

    if( il && irl ) {
        //diff. of derivatives at a single point (delta y on right side)
        dr = ( gl - grl ) * ( q->pr->x - q->pl->x );
        if( dr < gscd_ARMS_YEPS )
            //adjust dr to avoid numerical problems
            dr = gscd_ARMS_YEPS;
    }

    if( ir && irl ) {
        //diff. of derivatives at a single point (delta y on left side)
        dl = ( grl - gr ) * ( q->pr->x - q->pl->x );
        if( dl < gscd_ARMS_YEPS )
            //adjust dl to avoid numerical problems
            dl = gscd_ARMS_YEPS;
    }

    if( il && ir && irl ) {
        //gradients on both sides
        //NOTE:MM:from syst. of eqtns: dr/(xr-xl)=(yp-yr)/(xp-xl) and dl/(xr-xl)=(yp-yr)/(xr-xp)
        q->x = ( dl * q->pr->x + dr * q->pl->x )/( dl + dr );
        ///q->y = ( dl * q->pr->y + dr * q->pl->y + dl * dr )/( dl + dr );//MM:possibly error in org.code
        q->y = ( dl * q->pr->y + dr * q->pr->y + dl * dr )/( dl + dr );//NOTE:MM
    } else if( il && irl ) {
        //gradient only on left side
        q->x = q->pr->x;
        q->y = q->pr->y + dr;
    } else if( ir && irl ) {
        //gradient only on right side
        q->x = q->pl->x;
        q->y = q->pl->y + dl;
    } else if( il ) {
        //right hand bound
        q->y = q->pl->y + gl * ( q->x - q->pl->x );
    } else if( ir ) {
        //left hand bound
        q->y = q->pr->y - gr * ( q->pr->x - q->x );
    } else
        //gradient on neither side - should be impossible
        return PSL_ERR_INVALID;

    if(( q->pl && q->x < q->pl->x )||( q->pr && q->pr->x < q->x ))
        //intersection point outside interval (through imprecision)
        return PSL_ERR_INVALID;

    return 0;
}

// -------------------------------------------------------------------------
// Sample: sample from piecewise exponential envelope
//
int CARMS::Sample( SARMSpoint* p )
{
    double prob;
    int err;

    //sample a uniform
    prob = _carmsRNGss_.GetDouble();
    //get x-value correponding to a cumulative probability prob
    if(( err = Invert( prob, p )) != PSL_OK )
        return err;

    return 0;
}

// -------------------------------------------------------------------------
// Invert: obtain a point corresponding to a qiven cumulative probability;
//     p is working POINT holding the sampled value
//
int CARMS::Invert( double prob, SARMSpoint* p )
{
    SARMSpoint* q = GetEnvelope();
    double  u, xl, xr;
    double  yl, yr;
    double  eyl, eyr;
    double  frac, z;
    double  value = 0.0;
    int err;

    if( q == NULL || p == NULL )
        return PSL_ERR_ADDRESS;
    //find rightmost point in envelope
    for( ; q->pr; q = q->pr );

    //find exponential piece containing point implied by prob
    u = prob * q->cdf;
    for( ; q->pl && u < q->pl->cdf; q = q->pl );

    if( q->pl == NULL )
        return PSL_ERR_ADDRESS;

    //set left and right pnts of p, etc.
    p->pl = q->pl;
    p->pr = q;
    p->f = 0;
    p->cdf = u;

    //calculate fraction of integral within this piece
    frac = ( u - q->pl->cdf )/( q->cdf - q->pl->cdf );

    //get the required x-value
    if( q->pl->x == q->x ) {
        //interval is of zero length
        p->x = q->x;
        p->y = q->y;
        p->ey = q->ey;
    } else {
        xl = q->pl->x;
        xr = q->x;
        yl = q->pl->y;
        yr = q->y;
        eyl = q->pl->ey;
        eyr = q->ey;
        if( fabs( yr - yl ) < gscd_ARMS_YEPS ) {
            //linear approximation was used in integration in function cumulate
            if( gscd_ARMS_EYEPS * fabs( eyr + eyl ) < fabs( eyr - eyl )) {
                p->x = xl + ( xr - xl )/( eyr - eyl ) *
                            ( -eyl + sqrt(( 1. - frac ) * eyl*eyl + frac*eyr*eyr ));
            } else
                p->x = xl + ( xr - xl ) * frac;
            p->ey = ( p->x - xl )/( xr - xl ) * ( eyr - eyl ) + eyl;
            if(( err = logshift( p->ey, GetEnvMax(), &value )) != PSL_OK )
                return err;
            p->y = value;
        } else {
            //piece was integrated exactly in function cumulate
            if(( err = logshift(( 1. - frac ) * eyl + frac * eyr, GetEnvMax(), &value )) != PSL_OK )
                return err;
            p->x = xl + ( xr - xl )/( yr - yl ) * ( -yl + value );
            p->y = ( p->x - xl )/( xr - xl ) * (yr - yl) + yl;
            if(( err = expshift( p->y, GetEnvMax(), &value )) != PSL_OK )
                return err;
            p->ey = value;
        }
    }

    //guard against imprecision yielding point outside interval
    if( p->x < xl || xr < p->x )
        return PSL_ERR_ROUNDOFF;

    return 0;
}

// -------------------------------------------------------------------------
// Test: perform rejection, squeezing, and metropolis tests;
//     p is point to be tested
//
int CARMS::Test( SARMSpoint* p, int* acc )
{
    SARMSpoint* ql, *qr;
    double  u, y;
    double  ysqueez, ynew, yold;
    double  znew, zold, w;
    double  value = 0.0;
    int err;

    if( GetEnvelope() == NULL || GetFEval() == NULL || 
        p == NULL || acc == NULL )
        return PSL_ERR_ADDRESS;

    *acc = 0;
    //for rejection test
    u = _carmsRNGtst_.GetDouble0() * p->ey;
    if(( err = logshift( u, GetEnvMax(), &value )) != PSL_OK )
        return err;
    y = value;

    if( p->pl == NULL || p->pr == NULL )
        return PSL_ERR_ADDRESS;

    if( !GetDoMHupdate() && p->pl->pl && p->pr->pr ) {
        //perform squeezing test
        if( p->pl->f )
            ql = p->pl;
        else
            ql = p->pl->pl;

        if( p->pr->f )
            qr = p->pr;
        else
            qr = p->pr->pr;

        //avg area height
        ysqueez = ( qr->y * ( p->x - ql->x ) + ql->y * ( qr->x - p->x )) /
                  ( qr->x - ql->x );
        if( y <= ysqueez ) {
            //accept point at squeezing step
            *acc = 1;
            return 0;
        }
    }

    //evaluate log density at point to be tested
    if(( err = ( *GetFEval())( p->x, &value, GetParams())) != PSL_OK )
        return err;
    IncNoFevals();
    IncNoCurFevals();
    ynew = value;

    //perform rejection test
    if( !GetDoMHupdate() || ( GetDoMHupdate() && ynew <= y )) {
        //update envelope
        p->y = ynew;
        if(( err = expshift( p->y, GetEnvMax(), &value )) != PSL_OK )
            return err;
        p->ey = value;
        p->f = 1;
        if(( err = Update( p )) != PSL_OK )
            //envelope violation without metropolis
            return err;

        //perform rejection test
        if( ynew <= y )
            //reject point at rejection step
            return 0;
        else {
            //accept point at rejection step
            *acc = 1;
            return 0;
        }
    }

    //continue with metropolis step
    yold = GetFvalAtPrevMCvalue();

    //find envelope piece containing metxprev_
    ql = GetEnvelope();
    for( ; ql->pl; ql = ql->pl );

    if( ql->pr == NULL )
        return PSL_ERR_ADDRESS;

    for( ; ql->pr->x < GetPrevMCvalue(); ql = ql->pr );

    if( ql->pr == NULL )
        return PSL_ERR_ADDRESS;

    qr = ql->pr;
    //calculate height of envelope at metxprev_
    w = ( GetPrevMCvalue() - ql->x )/( qr->x - ql->x );
    zold = ql->y + w * ( qr->y - ql->y );//envelope value at metxprev_
    znew = p->y;
    if( yold < zold ) zold = yold;
    if( ynew < znew ) znew = ynew;
    //M-H test weight in log-space
    w = ynew - znew - yold + zold;
    if( 0.0 < w )
        w = 0.0;

    if( -gscd_ARMS_YCEIL < w )
        w = exp( w );
    else
        w = 0.0;

    u = _carmsRNGtst_.GetDouble();
    if( w < u ) {
        //M-H test violated, replace current point with previous MC iterate value
        p->x = GetPrevMCvalue();
        p->y = GetFvalAtPrevMCvalue();
        if(( err = expshift( p->y, GetEnvMax(), &value )) != PSL_OK )
            return err;
        p->ey = value;
        p->f = 1;
        p->pl = ql;
        p->pr = qr;
    } else {
        //trial point accepted by metropolis, update previous MC iterate value
        SetPrevMCvalue( p->x );
        SetFvalAtPrevMCvalue( ynew );
    }
    *acc = 1;
    return 0;
}

// -------------------------------------------------------------------------
// Update: update envelope to incorporate new point on log density;
//     p is point to be incorporated
//
int CARMS::Update( SARMSpoint* p )
{
    SARMSpoint* m;
    SARMSpoint* ql, *qr, *q;
    double  value;
    int err;

    if( GetEnvelope() == NULL || GetFEval() == NULL ||
        p == NULL || p->pl == NULL || p->pr == NULL )
        return PSL_ERR_ADDRESS;

    if( !p->f || GetMaxNoEnvPoints() - 2 < GetNoEnvPoints())
        //y-value has not been evaluated or no room for further points;
        //ignore this point
        return 0;

    //copy working point p to a new point q
    q = GetEnvelope() + GetNoEnvPoints();
    if( q == NULL )
        return PSL_ERR_ADDRESS;
    IncNoEnvPoints();
    q->x = p->x;
    q->y = p->y;
    q->f = 1;

    //allocate an unused point for a new intersection
    m = GetEnvelope() + GetNoEnvPoints();
    if( m == NULL )
        return PSL_ERR_ADDRESS;
    IncNoEnvPoints();
    m->f = 0;
    if( p->pl->f && !p->pr->f ) {
        //left end of piece is on log density;
        //set up new intersection in interval between p->pl and p
        m->pl = p->pl;
        m->pr = q;
        q->pl = m;
        q->pr = p->pr;
        m->pl->pr = m;
        q->pr->pl = q;
    } else if( !p->pl->f && p->pr->f ) {
        //right end of interval is on log density;
        //set up new intersection in interval between p and p->pr
        m->pr = p->pr;
        m->pl = q;
        q->pr = m;
        q->pl = p->pl;
        m->pr->pl = m;
        q->pl->pr = q;
    } else
        //should be impossible: adjacent pnts are both density or intersection pnts
        return PSL_ERR_INVALID;

    //adjust position of q within interval if too close to an endpoint
    if( q->pl->pl != NULL )
        ql = q->pl->pl;
    else
        ql = q->pl;

    if( q->pr->pr != NULL )
        qr = q->pr->pr;
    else
        qr = q->pr;

    if( q->x < ( 1. - gscd_ARMS_XEPS ) * ql->x + gscd_ARMS_XEPS * qr->x ) {
        //q too close to left end of interval
        q->x = ( 1. - gscd_ARMS_XEPS ) * ql->x + gscd_ARMS_XEPS * qr->x;
        if(( err = ( *GetFEval())( q->x, &value, GetParams())) != PSL_OK )
            return err;
        IncNoFevals();
        IncNoCurFevals();
        q->y = value;
    } else if ( gscd_ARMS_XEPS * ql->x + ( 1. - gscd_ARMS_XEPS ) * qr->x < q->x ) {
        //q too close to right end of interval
        q->x = gscd_ARMS_XEPS * ql->x + ( 1. - gscd_ARMS_XEPS ) * qr->x;
        if(( err = ( *GetFEval())( q->x, &value, GetParams())) != PSL_OK )
            return err;
        IncNoFevals();
        IncNoCurFevals();
        q->y = value;
    }

    //revise intersection points
    if(( err = Intersection( q->pl )) != PSL_OK )
        return err;

    if(( err = Intersection( q->pr )) != PSL_OK )
        return err;

    if( q->pl->pl != NULL )
        if(( err = Intersection( q->pl->pl->pl )) != PSL_OK )
            return err;

    if( q->pr->pr != NULL )
        if(( err = Intersection( q->pr->pr->pr )) != PSL_OK )
            return err;

    //exponentiate and integrate new envelope
    if(( err = CalcCDF()) != PSL_OK )
        return err;

    return 0;
}

// -------------------------------------------------------------------------
// CalcCDF: exponentiate and integrate envelope
//
int CARMS::CalcCDF()
{
    SARMSpoint* q, *qlmost;
    double  value = 0.0;
    int err;

    if( GetEnvelope() == NULL )
        return PSL_ERR_ADDRESS;

    //find left end of envelope
    for( qlmost = GetEnvelope(); qlmost->pl; qlmost = qlmost->pl );

    //find maximum y-value: search envelope
    value = qlmost->y;
    for( q = qlmost->pr; q; q = q->pr )
        if( value < q->y )
            value = q->y;
    SetEnvMax( value );

    //exponentiate envelope
    for( q = qlmost; q; q = q->pr ) {
        if(( err = expshift( q->y, GetEnvMax(), &value )) != PSL_OK )
            return err;
        q->ey = value;
    }

    //integrate exponentiated envelope
    qlmost->cdf = 0.;
    for( q = qlmost->pr; q; q = q->pr ) {
        if(( err = Area( q, &value )) != PSL_OK )
            return err;
        q->cdf = q->pl->cdf + value;
// fprintf(stderr,"CDF: x=%f cdf=%f\n",q->x,q->cdf);
    }
    return 0;
}

// -------------------------------------------------------------------------
// Area: integrate piece of exponentiated envelope to left of point q
//
int CARMS::Area( SARMSpoint* q, double* value )
{
    if( q == NULL || value == NULL )
        return PSL_ERR_ADDRESS;

    *value = 0.0;

    if( q->pl == NULL )
      //leftmost point in envelope
      return 0;

    if( q->pl->x == q->x )
      //interval is zero length
      return 0;

    if( fabs( q->y - q->pl->y ) < gscd_ARMS_YEPS )
        //integrate straight line piece
        *value = 0.5 *( q->ey + q->pl->ey )*( q->x - q->pl->x );
    else
        //integrate exponential piece
        *value = ( q->ey - q->pl->ey )/( q->y - q->pl->y )*( q->x - q->pl->x );
// fprintf(stderr,"Area: xl=%f x=%f yl=%f y=%f  eyl=%f ey=%f  value=%f\n",
//         q->pl->x,q->x,q->pl->y,q->y,q->pl->ey,q->ey,*value);

    return 0;
}

