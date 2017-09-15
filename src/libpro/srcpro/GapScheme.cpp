/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <stdlib.h>
#include <string.h>

#include "mystring.h"
#include "myexcept.h"

#include "GapScheme.h"
#include "Serializer.h"

//local constants
const double    lc_tpm = 1.0;
const double    lc_tpi = 0.4;
//static members
//transition power indices
double GapScheme::s_trnpowind_[P_NSTATES] = 
    {lc_tpm,lc_tpm,lc_tpm, lc_tpi,lc_tpi,lc_tpi, lc_tpi,lc_tpi,lc_tpi};


// -------------------------------------------------------------------------
// default constructor: gap opening costs are initialized to zero, no memory
//     allocation is performed
// -------------------------------------------------------------------------

GapScheme::GapScheme( double gap_open, double gap_extend )
:
    fixed( false ),
    posvec( NULL ),
    posext( NULL ),

    orgtrobs_( NULL ),
    trnprobs_( NULL ),
    logtrnps_( NULL ),
    nsltrnps_( NULL ),
    vector( NULL ),
    extvec( NULL ),
    weights( NULL ),
    deleteInt( NULL ),
    length( 0 ),
    openCost( gap_open ),
    extendCost( gap_extend ),
    scaledopenCost( gap_open ),
    scaledextendCost( gap_extend ),
    aacids( NULL ),
    allocated( 0 ),
    scalefactor( 1 ),
    scoremult_( 1.0 ),
    configuration( NULL ),

    positionalaccorr( NULL ),

    gapprobfactevalue( 1.0 ),
    gapprobfactweight( 0.0 ),
    gapprobfactshift( 0.0 ),

    acorrnumerator1st( -1.0 ),
    acorrnumerator2nd( -1.0 ),
    acorrlogscale( -1.0 ),
    acorrdenomscale( -1.0 ),

    thickness( 0 ),
    contextevalue( -1.0 ),
    contextadjusted( false ),
    useposaccorrect( false )
{
    ResetContextEvalue();
    for( size_t n = 0; n < DCnt; n++ )
        deletes[n] = NULL;
}

// -------------------------------------------------------------------------
// destructor: deallocates memory if it was allocated before
// -------------------------------------------------------------------------

GapScheme::~GapScheme()
{
    destroy();
}

// -------------------------------------------------------------------------
// IsCompatible: verifies whether the vector containing gap opening costs is 
//     compositionally identical to another one
// -------------------------------------------------------------------------

bool GapScheme::IsComposIdentical( const GapScheme& one ) const
{
    if( GetColumns() != one.GetColumns())
        return false;

    for( int n = 0; n < GetColumns(); n++ )
        if( aacids[n] != one.aacids[n] )
            return false;

    return true;
}

// -------------------------------------------------------------------------
// destroy: deallocate memory and reset values
// -------------------------------------------------------------------------

void GapScheme::destroy()
{
    if( orgtrobs_ ) { free( orgtrobs_ ); orgtrobs_ = NULL; }
    if( trnprobs_ ) { free( trnprobs_ ); trnprobs_ = NULL; }
    if( logtrnps_ ) { free( logtrnps_ ); logtrnps_ = NULL; }
    if( nsltrnps_ ) { free( nsltrnps_ ); nsltrnps_ = NULL; }

    if( posext ) { free( posext ); posext = NULL; }
    if( posvec ) { free( posvec ); posvec = NULL; }

    if( weights) { free( weights); weights= NULL; }
    if( extvec ) { free( extvec ); extvec = NULL; }
    if( vector ) { free( vector ); vector = NULL; }
    if( aacids ) { free( aacids ); aacids = NULL; }

    if( positionalaccorr ) { free( positionalaccorr ); positionalaccorr = NULL; }

    if( deleteInt ) { free( deleteInt ); deleteInt = NULL; }

    for( size_t n = 0; n < DCnt; n++ )
        if( deletes[n] ) {
            free( deletes[n] ); deletes[n] = NULL;
        }

    length = 0;
    allocated = 0;
}

// -------------------------------------------------------------------------
// Clear: erase memory and reset values
// -------------------------------------------------------------------------

void GapScheme::Clear()
{
    if( allocated <= 0 )
        return;

    if( orgtrobs_ ) memset( orgtrobs_, 0, sizeof( double ) * P_NSTATES * ( allocated + 1 ));
    if( trnprobs_ ) memset( trnprobs_, 0, sizeof( double ) * P_NSTATES * ( allocated + 1 ));
    if( logtrnps_ ) memset( logtrnps_, 0, sizeof( double ) * P_NSTATES * ( allocated + 1 ));
    if( nsltrnps_ ) memset( nsltrnps_, 0, sizeof( double ) * P_NSTATES * ( allocated + 1 ));

    if( posext )    memset( posext, 0, sizeof( double ) * allocated );
    if( posvec )    memset( posvec, 0, sizeof( double ) * allocated );

    if( weights)    memset( weights, 0, sizeof( double ) * allocated );
    if( extvec )    memset( extvec, 0, sizeof( double ) * allocated );
    if( vector )    memset( vector, 0, sizeof( double ) * allocated );
    if( aacids )    memset( aacids, 0, sizeof( char ) * allocated );

    if( positionalaccorr )  memset( positionalaccorr, 0, sizeof( double ) * allocated );

    if( deleteInt ) memset( deleteInt, 0, sizeof( int ) * allocated );

    for( size_t n = 0; n < DCnt; n++ )
        if( deletes[n] )
            memset( deletes[n], 0, sizeof( double ) * allocated );

    length = 0;
}

// -------------------------------------------------------------------------
// reallocate: memory reallocation for the vector of gap costs and vector
//     'aacids'
// -------------------------------------------------------------------------

void GapScheme::reallocate( int howmuch )
{
    double(*tmp_orgtrobs )[P_NSTATES];
    double(*tmp_trnprobs )[P_NSTATES];
    double(*tmp_logtrnps )[P_NSTATES];
    double(*tmp_nsltrnps )[P_NSTATES];
    int*    tmp_deleteInt;
    double* tmp_deletes[DCnt];
    double* tmp_weights;
    double* tmp_posvec;
    double* tmp_posext;
    double* tmp_vector;
    double* tmp_extvec;
    char*   tmp_aacids;
    double* tmp_positionalaccorr;

    if( howmuch <= allocated )
        return;

    if( allocated <= 0 ) {
        tmp_orgtrobs = ( double(*)[P_NSTATES] )malloc( sizeof( double ) * P_NSTATES * ( howmuch + 1 ));
        tmp_trnprobs = ( double(*)[P_NSTATES] )malloc( sizeof( double ) * P_NSTATES * ( howmuch + 1 ));
        tmp_logtrnps = ( double(*)[P_NSTATES] )malloc( sizeof( double ) * P_NSTATES * ( howmuch + 1 ));
        tmp_nsltrnps = ( double(*)[P_NSTATES] )malloc( sizeof( double ) * P_NSTATES * ( howmuch + 1 ));
        tmp_deleteInt = ( int* )malloc( sizeof( int ) * howmuch );

        tmp_deletes[DBeg]   = ( double* )malloc( sizeof( double ) * howmuch );
        tmp_deletes[DEnd]   = ( double* )malloc( sizeof( double ) * howmuch );
        tmp_deletes[DSlope] = ( double* )malloc( sizeof( double ) * howmuch );

        tmp_weights = ( double* )malloc( sizeof( double ) * howmuch );

        tmp_posext  = ( double* )malloc( sizeof( double ) * howmuch );
        tmp_posvec  = ( double* )malloc( sizeof( double ) * howmuch );

        tmp_extvec  = ( double* )malloc( sizeof( double ) * howmuch );
        tmp_vector  = ( double* )malloc( sizeof( double ) * howmuch );
        tmp_aacids  = ( char*   )malloc( sizeof( char   ) * howmuch );

        tmp_positionalaccorr = ( double* )malloc( sizeof( double ) * howmuch );

    } else {
        tmp_orgtrobs = ( double(*)[P_NSTATES] )realloc( orgtrobs_, sizeof( double ) * P_NSTATES * ( howmuch + 1 ));
        tmp_trnprobs = ( double(*)[P_NSTATES] )realloc( trnprobs_, sizeof( double ) * P_NSTATES * ( howmuch + 1 ));
        tmp_logtrnps = ( double(*)[P_NSTATES] )realloc( logtrnps_, sizeof( double ) * P_NSTATES * ( howmuch + 1 ));
        tmp_nsltrnps = ( double(*)[P_NSTATES] )realloc( nsltrnps_, sizeof( double ) * P_NSTATES * ( howmuch + 1 ));
        tmp_deleteInt = ( int* )realloc( deleteInt, sizeof( int ) * howmuch );

        tmp_deletes[DBeg]   = ( double* )realloc( deletes[DBeg],    sizeof( double ) * howmuch );
        tmp_deletes[DEnd]   = ( double* )realloc( deletes[DEnd],    sizeof( double ) * howmuch );
        tmp_deletes[DSlope] = ( double* )realloc( deletes[DSlope],  sizeof( double ) * howmuch );

        tmp_weights = ( double* )realloc( weights,  sizeof( double  ) * howmuch );

        tmp_posext  = ( double* )realloc( posext,   sizeof( double  ) * howmuch );
        tmp_posvec  = ( double* )realloc( posvec,   sizeof( double  ) * howmuch );

        tmp_extvec  = ( double* )realloc( extvec,   sizeof( double  ) * howmuch );
        tmp_vector  = ( double* )realloc( vector,   sizeof( double  ) * howmuch );
        tmp_aacids  = ( char*   )realloc( aacids,   sizeof( char    ) * howmuch );

        tmp_positionalaccorr = ( double* )realloc( positionalaccorr, sizeof( double ) * howmuch );
    }

    if( !tmp_orgtrobs ||
        !tmp_trnprobs ||
        !tmp_logtrnps ||
        !tmp_nsltrnps ||
        !tmp_deleteInt ||
        !tmp_deletes[DBeg] ||
        !tmp_deletes[DEnd] ||
        !tmp_deletes[DSlope] ||
        !tmp_weights ||
        !tmp_posext || !tmp_posvec ||
        !tmp_extvec || !tmp_vector || !tmp_aacids ||
        !tmp_positionalaccorr
    )
        throw myruntime_error( mystring( "GapScheme: Not enough memory." ));

    orgtrobs_ = tmp_orgtrobs;
    trnprobs_ = tmp_trnprobs;
    logtrnps_ = tmp_logtrnps;
    nsltrnps_ = tmp_nsltrnps;
    deleteInt = tmp_deleteInt;
    deletes[DBeg]   = tmp_deletes[DBeg];
    deletes[DEnd]   = tmp_deletes[DEnd];
    deletes[DSlope] = tmp_deletes[DSlope];
    weights = tmp_weights;
    posext  = tmp_posext;
    posvec  = tmp_posvec;
    extvec  = tmp_extvec;
    vector  = tmp_vector;
    aacids  = tmp_aacids;
    positionalaccorr = tmp_positionalaccorr;

    if( allocated <= 0 ) {
        memset( tmp_orgtrobs, 0, sizeof( double ) * P_NSTATES );
        memset( tmp_trnprobs, 0, sizeof( double ) * P_NSTATES );
        memset( tmp_logtrnps, 0, sizeof( double ) * P_NSTATES );
        memset( tmp_nsltrnps, 0, sizeof( double ) * P_NSTATES );
    }
    tmp_orgtrobs++;
    tmp_trnprobs++;
    tmp_logtrnps++;
    tmp_nsltrnps++;

    // fill uninitialized memory with zeros
    memset( tmp_orgtrobs + allocated,  0, sizeof( double ) * P_NSTATES *( howmuch - allocated ));
    memset( tmp_trnprobs + allocated,  0, sizeof( double ) * P_NSTATES *( howmuch - allocated ));
    memset( tmp_logtrnps + allocated,  0, sizeof( double ) * P_NSTATES *( howmuch - allocated ));
    memset( tmp_nsltrnps + allocated,  0, sizeof( double ) * P_NSTATES *( howmuch - allocated ));
    memset( deleteInt + allocated,  0, sizeof( int ) * ( howmuch - allocated ));

    memset( deletes[DBeg]   + allocated,  0, sizeof( double ) * ( howmuch - allocated ));
    memset( deletes[DEnd]   + allocated,  0, sizeof( double ) * ( howmuch - allocated ));
    memset( deletes[DSlope] + allocated,  0, sizeof( double ) * ( howmuch - allocated ));

    memset( weights + allocated,  0, sizeof( double ) * ( howmuch - allocated ));

    memset( posext  + allocated,  0, sizeof( double ) * ( howmuch - allocated ));
    memset( posvec  + allocated,  0, sizeof( double ) * ( howmuch - allocated ));

    memset( extvec  + allocated,  0, sizeof( double ) * ( howmuch - allocated ));
    memset( vector  + allocated,  0, sizeof( double ) * ( howmuch - allocated ));
    memset( aacids  + allocated,  0, sizeof( char   ) * ( howmuch - allocated ));

    memset( positionalaccorr + allocated,  0, sizeof( double ) * ( howmuch - allocated ));

    allocated = howmuch;
}

// -------------------------------------------------------------------------
// Initialize: initializes data dependent on the insertion- and deletion-
//     state information
// -------------------------------------------------------------------------

void GapScheme::Initialize()
{
    InitializePosGapCosts();
    ComputeDeleteSlopes();
    ReevaluateTransProbs();
}

// -------------------------------------------------------------------------
// Initialize: initializes the gap opening and extension vectors
// -------------------------------------------------------------------------

void GapScheme::InitializePosGapCosts()
{
    InitializePosGapOpenCosts();
    InitializePosGapExtendCosts();
}

// InitializePosGapOpenCosts: positional gap open costs
//
void GapScheme::InitializePosGapOpenCosts()
{
    for( int n = 0; n < GetColumns(); n++ )
        SetPosOpenAt( GetScaledOpenCost(), n );
}

// InitializePosGapExtendCosts: positional gap extend costs
//
void GapScheme::InitializePosGapExtendCosts()
{
    for( int n = 0; n < GetColumns(); n++ )
        SetPosExtendAt( GetScaledExtendCost(), n );
}

// -------------------------------------------------------------------------
// ReevaluateTransProbs: reevaluate transition probabilities
//
void GapScheme::ReevaluateTransProbs()
{
    int n, s;
    double  val, logval, logfct;
    int     scale = GetScaleFactor();
    //multiplier of scores in the procedure of obtaining required lambda (cbs)
    double  mult = GetScoresMultiplier();

    if( scale <= 0 )
        scale = 1;
    //if( mult <= 0.0 )///cbs does not affect log values when commented
        mult = 1.0;

    for( n = -1; n < GetColumns(); n++ ) {
        for( s = 0; s < P_NSTATES; s++ ) {
            val = GetOrgTrProbsAt( s, n );
            logval = LOG_PROB_MIN;
            if( 0.0 < val ) {
                logval = log( val ) * mult;
                logfct = GetTransPowIndex( s );
                if( P_MM + gTPTRANS_NTPS <= s && logfct != 1.0 )
                    logval *= logfct;
                val = exp( logval );
            }
            else
                val = 0.0;
            SetTransProbsAt( val, s, n );
            SetLogTransProbsAt( logval * scale, s, n );
            SetNSLogTransProbsAt( logval, s, n );
        }
    }
}

// -------------------------------------------------------------------------
// InitializePosACcorrection: initialize positional values of correction for
//     autocorrelation function given array of information contents
// -------------------------------------------------------------------------

void GapScheme::InitializePosACcorrection( const double* posinfcontents, double correct )
{
    double  numerator = GetAutocorrectionNumerator1st();
    double  defaultac = DEFAULT_AC_CORRECTION * 10.0;
    double  deviation = 0.25;
    double  value, entropy, sysrelent;

//     if( posinfcontents == NULL )
//         throw myruntime_error( mystring( "GapScheme: InitializePosACcorrection: Null address of entropies." ));

    defaultac += deviation;

    if( GetConfiguration()) {
        sysrelent = GetConfiguration()[ProcomUngapped].GetH();
        if( 0.0 < sysrelent )
            defaultac = numerator / sqrt( sysrelent );
    }

    for( int n = 0; n < GetColumns(); n++ ) {
        entropy = -1.0;
        if( posinfcontents )
            entropy = posinfcontents[n];
        if( 0.0 < entropy ) {
            if( 0.0 < correct )
                value = correct / sqrt( entropy );
            else
                value = numerator / sqrt( entropy );
        } else
            value = defaultac;
        SetPosACcorrectionAt( value * ( double )GetScaleFactor(), n );
    }
}

// -------------------------------------------------------------------------
// ResetACcorrection: resets correction for autocorrelation function to its
//     default value
// -------------------------------------------------------------------------

void GapScheme::ResetACcorrection()
{
    double  correct = DEFAULT_AC_CORRECTION;
    double  sysrelent = 1.0;
    double  deviation = 0.25;

    if( GetAutocorrectionNumerator1st() < 0.0 )
        throw myruntime_error( mystring( "GapScheme: Invalid numerator for autocorrection." ));

    if( GetConfiguration()) {
        sysrelent = GetConfiguration()[ProcomUngapped].GetH();
        correct = sysrelent + deviation;
        if( 0.0 < sysrelent )
            correct = GetAutocorrectionNumerator1st() / sqrt( sysrelent );
    }

    SetACcorrection( correct );
}

// -------------------------------------------------------------------------
// ResetPosACcorrection: resets positional values of correction for
//     autocorrelation function
// -------------------------------------------------------------------------

void GapScheme::ResetPosACcorrection()
{
    double  numerator = GetAutocorrectionNumerator1st();
    double  defaultac = DEFAULT_AC_CORRECTION * 10.0;
    double  deviation = 0.25;
    double  sysrelent;

    defaultac += deviation;

    if( GetConfiguration()) {
        sysrelent = GetConfiguration()[ProcomUngapped].GetH();
        if( 0.0 < sysrelent )
            defaultac = numerator / sqrt( sysrelent );
    }

    defaultac *= ( double )GetScaleFactor();

    for( int n = 0; n < GetColumns(); n++ ) {
        SetPosACcorrectionAt( defaultac, n );
    }
}

// -------------------------------------------------------------------------
// SetACcorrection: sets correction for autocorrelation function
// -------------------------------------------------------------------------

void GapScheme::SetACcorrection( double value )
{
    autcorrection = value;
    SetScaledACcorrection( autcorrection * ( double )GetScaleFactor());
}

// -------------------------------------------------------------------------
// SetACcorrectionByEval: sets correction for autocorrelation function
//     according to the given evalue
// -------------------------------------------------------------------------

void GapScheme::SetACcorrectionByEval( double evalue, double relent, const double* posinfcontents )
{
    double  correct = 0.0;
    double  sysrelent = 1.0;
    double  evalbound = 1e-5;   //used to be computed
    double  logscale = GetAutocorrectionDenomScale();
    double  evalscale = GetAutocorrectionLogScale();    //not a typo!
    double  upper = GetAutocorrectionNumerator2nd();

    if( GetAutocorrectionNumerator2nd() < 0.0 ||
        GetAutocorrectionLogScale() < 0.0 ||
        GetAutocorrectionDenomScale() < 0.0 )
        throw myruntime_error( mystring( "GapScheme: Invalid parameter values for autocorrection on 2nd pass." ));

    if( GetConfiguration())
        sysrelent = GetConfiguration()[ProcomUngapped].GetH();

    if( 0.0 < sysrelent )
        upper /= sqrt( sysrelent );
//     if( 0.0 < relent )
//         upper /= sqrt( relent );


    evalbound = exp( -evalscale - 1.0 / ( logscale * upper ));
    correct = upper;

    if( evalbound < evalue )
        //avoid phase transition of function
        correct = upper;
    else
        if( 0.0 < evalue ) {
            correct = -1.0 / (( log( evalue ) + evalscale ) * logscale );
            if( upper < correct )
                correct = upper;
            else
                if( posinfcontents )
                    InitializePosACcorrection( posinfcontents, correct );
        }

    SetACcorrection( correct );
}

// -------------------------------------------------------------------------
// ResetExtendCost: resets extend cost to its initial value
// -------------------------------------------------------------------------

void GapScheme::ResetExtendCost()
{
    SetScaledExtendCost( GetExtendCost() * GetScaleFactor());
    InitializePosGapExtendCosts();
}

// -------------------------------------------------------------------------
// SetExtendCostByEval: determines extend cost according to e-value given
// TODO: requires more investigation
// -------------------------------------------------------------------------

void GapScheme::SetExtendCostByEval( double evalue, double relent )
{
    if( evalue < 0.0 )
        return;
/// see todo; constants should be reestimated
    return;

    double  numerator = -30.0;
    double  maxcost = -TIMES2( GetExtendCost());
    double  mincost = 0.0;//-GetExtendCost();
    double  cost = mincost;

    if( 0.0 < evalue ) {
        cost += numerator / log( evalue );
        if( maxcost < cost || cost < 0.0 )
            cost = maxcost;
    }

    SetScaledExtendCost( -cost * GetScaleFactor());
    InitializePosGapExtendCosts();
}

// -------------------------------------------------------------------------
// AdjustContextByEval: adjusts context parameters according to e-value
//     given
// -------------------------------------------------------------------------

void GapScheme::AdjustContextByEval( double evalue, double relent, const double* posinfcontents )
{
    SetACcorrectionByEval( evalue, relent, posinfcontents );
    SetExtendCostByEval( evalue, relent );
    SetContextOn();
    SetContextEvalue( evalue );
}

// -------------------------------------------------------------------------
// ResetContext: resets context parameters
// -------------------------------------------------------------------------

void GapScheme::ResetContext()
{
    if( GetContextOn())
        ResetPosACcorrection();
    ResetACcorrection();
    ResetExtendCost();
    ResetThickness();
    SetContextOff();
}

// -------------------------------------------------------------------------
// Prepare: prepares gap costs for using for alignments
// -------------------------------------------------------------------------

void GapScheme::Prepare( int scfact )
{
    SetScaleFactor( scfact );
    SetScaledOpenCost( GetOpenCost() * GetScaleFactor());
    if( !GetContextOn()) {
        SetScaledExtendCost( GetExtendCost() * GetScaleFactor());
        SetACcorrection( GetACcorrection());
        InitializePosGapExtendCosts();
        ResetContextEvalue();
    }
    InitializePosGapOpenCosts();
//     InitializePosGapCosts();
    ComputeCosts();
}

// Prepare: overloaded
//
void GapScheme::Prepare( int scfact, const double* posinfcontents )
{
    Prepare( scfact );
    if( !GetContextOn())
        InitializePosACcorrection( posinfcontents );
    ReevaluateTransProbs();
}

// Prepare: overloaded
//
void GapScheme::Prepare( double opencost, double extncost, bool fixedcosts )
{
    SetFixed( fixedcosts );
    SetOpenCost( opencost );
    SetExtendCost( extncost );
    Prepare(( int ) 1 );
    ReevaluateTransProbs();
}

// -------------------------------------------------------------------------
// ComputeDeleteSlopes: computes delete-state slopes for all positions
// -------------------------------------------------------------------------

void GapScheme::ComputeDeleteSlopes()
{
    for( int n = 0; n < GetColumns(); n++ )
        ComputeDeletesSlopeAt( n );
}

// -------------------------------------------------------------------------
// Push: push values related to the gap processing scheme
// -------------------------------------------------------------------------

void GapScheme::Push( const double( *tps )[P_NSTATES], double w, double delbeg, double delend, int delta, char aa )
{
    if( allocated < GetColumns() + 1 ) {
        reallocate( allocated + allocated + 1 );
    }
    PushAt( tps, w, delbeg, delend, delta, aa, GetColumns());
}

// -------------------------------------------------------------------------
// PushAt: pushes values at the given position
// -------------------------------------------------------------------------

void GapScheme::PushAt( const double( *tps )[P_NSTATES], double w, double delbeg, double delend, int delta, 
                        char aa, int position )
{
    if( allocated < position + 1 )
        throw myruntime_error( "GapScheme: Memory access error." );

    weights[position] = w;
    aacids [position] = aa;

    if( GetColumns() < position + 1 )
        SetColumns( position + 1 );

    SetOrgTrProbsAt( tps, position );
    SetDeleteStateAt( delbeg, delend, delta, position );
}

// -------------------------------------------------------------------------
// Autocorrelation: computes normalized sum of autocorrelations given
//     vector of values; i.e. computes all autocorrelation function values
//     by varying autocorrelation window size.
// -------------------------------------------------------------------------

double GapScheme::Autocorrelation( double* maxs, int wsize, const double* info, const double* posco )
{
    if( maxs == NULL || info == NULL || posco == NULL )
        throw myruntime_error( mystring( "GapScheme: Memory access error." ));

#ifdef __DEBUG__
    if( wsize <= 0 )
        throw myruntime_error( mystring( "GapScheme: Invalid window paramaters." ));
#endif

    double  ro = 0.0;                       //sum of autocorrelations
    double  ac = GetScaledACcorrection();   //correction for autocorrelation function
    double  aci, acj, dmi, dmj;
    double  dmlim = 100.0;                  //1/0.1^2
    double  icorrect = 0.00;
    double  scale = 0.4;
    size_t  no_pairs = 0;   //number of pairs, which is wsize! / ((wsize-2)! 2! )
    int     pribeg = 0;

    aci = acj = ac;

    if( wsize == 1 ) {
        if( GetUsePosACcorrections()) {
            aci = acj = posco[0];
        }
//         aci = GetScaleFactor() * 1.5 * ( icorrect - info[0] ) * LN2 / LOSCORES.StatisParam( Ungapped, Lambda );
//         aci = GetScaleFactor() * ( scale - ( scale * exp( scale * ( maxs[0] / GetScaleFactor()))));
        ro += ( maxs[0] + aci )*( maxs[0] + aci );
        no_pairs++;
    }
    else
        for( int priwin = wsize; 0 < priwin; priwin-- ) {
            int winpar = ( priwin  & 1 ) ^ 1;
            int midwin = ( priwin >> 1 ) - winpar;

            for( int i = pribeg; i < pribeg + priwin; i++ ) {
                int j = i + midwin;
                if( priwin <= j )
                    j -= priwin;

                if( GetUsePosACcorrections()) {
                    aci = posco[i];
                    acj = posco[j];
                }

//                 dmi = maxs[i] / GetScaleFactor(); dmj = maxs[j] / GetScaleFactor();
//                 if( 0.0 < dmi ) dmi = 1.0 / ( dmi * 0.4 ); else dmi = dmlim;
//                 if( 0.0 < dmj ) dmj = 1.0 / ( dmj * 0.4 ); else dmj = dmlim;
//                 dmi *= GetScaleFactor(); dmj *= GetScaleFactor();

//                 ro += ( maxs[i] + dmi + aci )*( maxs[j] + dmj + acj );
//                 ro += ( TIMES2( maxs[i] )+ aci )*( TIMES2( maxs[j] ) + acj );
                ro += ( maxs[i] + aci )*( maxs[j] + acj );//original
                no_pairs++;

//             for( int j = i + 1; j < wsize; j++ ) {
//                 if( !maxs[i] || !maxs[j] )
//                     continue;
//                 aci = GetScaleFactor() * 1.5 * ( icorrect - info[i] ) * LN2 / LOSCORES.StatisParam( Ungapped, Lambda );
//                 acj = GetScaleFactor() * 1.5 * ( icorrect - info[j] ) * LN2 / LOSCORES.StatisParam( Ungapped, Lambda );
//                 aci = GetScaleFactor() * ( scale - ( scale * exp( 0.2 * ( maxs[i] / GetScaleFactor()))));
//                 acj = GetScaleFactor() * ( scale - ( scale * exp( 0.2 * ( maxs[j] / GetScaleFactor()))));
//                 ro += ( maxs[i] + aci )*( maxs[j] + acj );
//                 no_pairs++;
//             }
            }
            if( priwin & 1 )
                pribeg++;
        }

    if( no_pairs )
        return ro / ( double )no_pairs;
    return ro;
}

// -------------------------------------------------------------------------
// AdjustCosts: adjusts gap costs according to the (max) scores obtained for
//     query or subject for each position of profile scoring matrix
// -------------------------------------------------------------------------

void GapScheme::AdjustCosts(
    double* maxscores,
    size_t  szvalues,
    const   LogOddsMatrix& logo,
    int     window,
    bool*   masks,
    size_t  szmasks,
    double  relent )
{
    if( szvalues == 0 )
        return;

    const double*   infovect = logo.GetInformation();
    const double*   corrects = GetPosACcorrection();

    SetThickness( logo.GetEffNoSequences());

    if( maxscores == NULL || infovect == NULL )
        throw myruntime_error( mystring( "GapScheme: Memory access error." ));

    if( szvalues < window )
        window = szvalues;

    if( window <= 0 )
        throw myruntime_error( mystring( "GapScheme: Non-positive autocorrelation window size." ));

    int     winpar = ( window  & 1 ) ^ 1;
    int     midwin = ( window >> 1 ) - winpar;

    double  extcost = GetScaledExtendCost();
    int     gapsize = GetColumns();

    int     m;
    double  value = 1.0, preval;
    double  avgmaxsc = 0.0; //average of max scores
    double  stddev = 0.0;   //standard deviation
    double  ro = 0.0;       //computed sum of autocorrelations
    double  sqrtro, roa;
    double  rolim = 100.0 * GetScaleFactor();   //1/0.1^2
    bool    winmasked;

    if( gapsize != szvalues || szvalues != szmasks ||
        gapsize != logo.GetColumns())
        throw myruntime_error( mystring( "GapScheme: Unable to compute gap costs." ));


    //window must be strictly positive !!
    for( m = 0; m < gapsize - window + 1; m++ ) {
        winmasked = false;
        for( int j = m; j < m + window; j ++ )
            if( masks[j] == true ) {
                winmasked = true;
                break;
            }

        ro = Autocorrelation( maxscores + m, window, infovect + m, corrects + m );
        //set positional gap costs for subject
        if( 0.0 < ro ) {
            preval = value;
//             value = -( int )rint( sqrt( ro ));
//             sqrtro = sqrt( ro );
//             roa = 0.0;
//             roa = sqrtro / GetScaleFactor();
//             roa = 1.0 / ( roa * 0.3 );
//             roa *= GetScaleFactor();
//             value = - TIMES2( sqrtro ) - roa;
            value = -sqrt( ro ); //original
//             if( extcost < value )
//                 value = extcost;
            //if simple derivative is 0
//             if( abs( value - preval ) <= 1 * GetScaleFactor() && 2 * extcost < value )
//                 value = extcost * 2;
            SetPosOpenAt( value, m + midwin );
        } else
//             SetPosOpenAt( preval = rolim, m + midwin );
            SetPosOpenAt( preval = 0.0/*extcost*/, m + midwin );//original
    }

    if( midwin < gapsize )
        for( m = 0; m < midwin && m < gapsize; m++ )
            SetPosOpenAt( GetPosOpenAt( midwin ), m );

    if( window - midwin < gapsize )
        for( m = gapsize - window + midwin + 1; m < gapsize; m++ )
            SetPosOpenAt( GetPosOpenAt( gapsize - window + midwin ), m );

    ComputeCosts();
    //reset context for proceeding calls of gap cost adjustment
    ResetContext();

//{{TEST***
// fprintf( stderr, "\nGap costs:\n" );
// for( m = 0; m < gapsize; m++ )
//     fprintf( stderr, "%4d  %7.2f  %7.2f  %4.2fi %4.2fd\n", m + 1,
//         GetOpenAt( m )/GetScaleFactor(), maxscores[m]/GetScaleFactor(), GetProbabilityAt( m ), GetDeleteOpenProbAt( m ));
//     fprintf( stderr, "%7.2f  %7.2f  %7.2f  %7.2f  %5.2f -- %d %5.2f\n",
//         maxscores[m], GetPosOpenAt( m ), GetOpenAt( m ), GetExtendAt( m ), GetWeightsAt( m ), masks[m], infovect[m] );
// fprintf( stderr, "\n//\n" );
//}}

}

// -------------------------------------------------------------------------
// ComputeCosts: computes position-specific gap costs given gap weights
// -------------------------------------------------------------------------

void GapScheme::ComputeCosts()
{
    double  minopen = GetScaledOpenCost();      //this is minimum score value a gap penalty can acquire
    double  maxopen = 0; //GetScaledExtendCost();//maximum score value a gap penalty can acquire
    double  minextn = GetScaledExtendCost();    //minimum score value a gap extension penalty can acquire
    double  maxextn = 0;                //maximum score value a gap extension penalty can acquire
    double  midextn = ( double )( maxextn - minextn ) / 2.0;
    double  midopen = ( double )( maxopen - minopen ) / 2.0;
    double  cost = 0;

#ifdef __DEBUG__
    if( posvec == NULL || posext == NULL ||
        vector == NULL || extvec == NULL || weights == NULL )
        throw myruntime_error( mystring( "GapScheme: Memory access error." ));
#endif

    if( maxopen < minopen )
        throw myruntime_error(
            mystring( "GapScheme: Invalid gap opening costs." ));

    for( int n = 0; n < GetColumns(); n++ ) {
        cost = minopen = GetPosOpenAt( n );
        if( !GetFixed()) {
                //do NOT uncomment: gap costs are computed by using insertion and deletion probabilities
//             cost =  minopen + ( maxopen - minopen ) * GetWeightsAt( n );
//             cost = -midopen + ( maxopen - minopen ) * SigmoidResponse( GetWeightsAt( n ));
        }
        if( 0.0 < cost )
            cost = 0.0;
        SetOpenAt( cost, n );

        cost = minextn = GetPosExtendAt( n );
        if( !GetFixed()) {
                //do NOT uncomment: gap costs are computed by using insertion and deletion probabilities
//             cost =  minextn + ( maxextn - minextn ) * GetWeightsAt( n );
//             cost = -midextn + ( maxextn - minextn ) * SigmoidResponse( GetWeightsAt( n ));
        }
        if( 0.0 < cost )
            cost = 0.0;
//         if( cost < GetOpenAt( n ))
//             cost = GetOpenAt( n );
        SetExtendAt( cost, n );
    }
}

// -------------------------------------------------------------------------
// OutputGapScheme: output full information of the gap costs
// -------------------------------------------------------------------------

void GapScheme::OutputGapScheme( const char* filename )
{
    FILE*   fp = stderr;
    size_t  l = 0;
    int     t, p;

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw myruntime_error( "Failed to open file for writing gap costs." );

    fprintf( fp, "Gap opening cost, %d\nExtending gap cost, %d\n\n", 
           ( int )rint( GetOpenCost()), ( int )rint( GetExtendCost()));
    fprintf( fp, "%9c", 32 );
    for( t = 0; t < P_NSTATES; t++ )
        fprintf( fp, " %5s", gTPTRANS_NAMES[t]);
    fprintf( fp, " %8ccosts weights deletions / int slope\n", 32 );

    for( p = -1; p < GetColumns(); p++ ) 
    {
        if( 0 <= p )
            fprintf( fp, "\n%5d %c   ", ++l, DehashCode( AAcid( p )));
        else
            fprintf( fp, "%10c", 32 );
        for( t = 0; t < P_NSTATES; t++ )
            fprintf( fp, "%6.3f", GetTransProbsAt/*GetOrgTrProbsAt*/( t, p ));
        if( 0 <= p )
            fprintf( fp, "  %4d %5.3f %6.3f - %-6.3f / %-3d %6.3f",
                ( *this )[p], GetWeightsAt( p ),
                GetDeletesBegAt( p ), GetDeletesEndAt( p ),
                GetDeletesIntervalAt( p ), GetDeletesSlopeAt( p ));
    }
    fprintf( fp, "\n" );
    if( fp != stderr )
        fclose( fp );
}

// -------------------------------------------------------------------------
// CheckIntegrity: checks integrity of the structure
// -------------------------------------------------------------------------

void GapScheme::CheckIntegrity() const
{
    int st, states;
    int pos;
    double  sum;
    const int       maxtr = gTPTRANS_NTPS;
    const double    accuracy = 1.0e-6;

    for( pos = -1; pos < GetColumns(); pos++ ) {
        if( 0 <= pos )
            if( NUMALPH <= AAcid( pos ) ||
                1.0 < GetWeightsAt( pos ) || GetWeightsAt( pos ) < 0.0,
                1.0 < GetDeletesBegAt( pos ) || GetDeletesBegAt( pos ) < 0.0 ||
                1.0 < GetDeletesEndAt( pos ) || GetDeletesEndAt( pos ) < 0.0 ||
                GetDeletesIntervalAt( pos ) < 0
            )
                throw myruntime_error( "GapScheme: Data corrupted." );

        //verify probability conservation
        for( states = 0; states < P_NSTATES; states += maxtr ) {
            sum = 0.0;
            for( st = states; st < states + maxtr && st < P_NSTATES; st++ )
                sum += GetOrgTrProbsAt( st, pos );
            if( sum < 1.0 - accuracy || sum > 1.0 + accuracy )
                throw myruntime_error( "GapScheme: Probabilities corrupted." );
        }
    }
}

// -------------------------------------------------------------------------
// Serialize: write the class data to binary file for reading them later
// -------------------------------------------------------------------------

void GapScheme::Serialize( Serializer& serializer ) const
{
    serializer.Write(( char* )&length, sizeof( length ), 1 );

    if( length <= 0 )
        return;

    int tmpopencost = ( int )rint( GetOpenCost());
    int tmpextendcost = ( int )rint( GetExtendCost());

    serializer.Write(( char* )orgtrobs_, sizeof( double ) * P_NSTATES, length + 1 );
//     serializer.Write(( char* )vector, sizeof( int ), length );
    serializer.Write(( char* )weights, sizeof( double ), length );
    serializer.Write(( char* )deletes[DBeg], sizeof( double ), length );
    serializer.Write(( char* )deletes[DEnd], sizeof( double ), length );
    serializer.Write(( char* )deleteInt, sizeof( int ), length );
    serializer.Write(( char* )aacids, sizeof( char ), length );

//     serializer.Write(( char* )&openCost, sizeof( openCost ), 1 );
//     serializer.Write(( char* )&extendCost, sizeof( extendCost ), 1 );

    serializer.Write(( char* )&tmpopencost, sizeof( tmpopencost ), 1 );
    serializer.Write(( char* )&tmpextendcost, sizeof( tmpextendcost ), 1 );
}

// -------------------------------------------------------------------------
// Deserialize: read data from binary file into the class members
// -------------------------------------------------------------------------

void GapScheme::Deserialize( Serializer& serializer )
{
    serializer.Read(( char* )&length, sizeof( length ), 1 );

    if( MAXCOLUMNS < length )
        throw myruntime_error(
            mystring( "GapScheme: Number of positions read from file is larger than maximum allowed." ));

    if( length <= 0 )
        throw myruntime_error(
            mystring( "GapScheme: Invalid number of positions read from file." ));

    int tmpopencost = SCORE_MIN;
    int tmpextendcost = SCORE_MIN;

    Reserve( length );

    serializer.Read(( char* )orgtrobs_, sizeof( double ) * P_NSTATES, length + 1 );
//     serializer.Read(( char* )vector, sizeof( int ), length );
    serializer.Read(( char* )weights, sizeof( double ), length );
    serializer.Read(( char* )deletes[DBeg], sizeof( double ), length );
    serializer.Read(( char* )deletes[DEnd], sizeof( double ), length );
    serializer.Read(( char* )deleteInt, sizeof( int ), length );
    serializer.Read(( char* )aacids, sizeof( char ), length );

//     serializer.Read(( char* )&openCost, sizeof( openCost ), 1 );
//     serializer.Read(( char* )&extendCost, sizeof( extendCost ), 1 );

    serializer.Read(( char* )&tmpopencost, sizeof( tmpopencost ), 1 );
    serializer.Read(( char* )&tmpextendcost, sizeof( tmpextendcost ), 1 );

    SetOpenCost( tmpopencost );
    SetExtendCost( tmpextendcost );
    //
    CheckIntegrity();
    Initialize();
}

