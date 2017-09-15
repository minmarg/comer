/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include <ctype.h>
#include "needconfig.h"
#include "libpro/srcaln/IMACounts.h"
#include "logitnormal.h"
#include "datapro.h"

const bool          SCORES_IN_PROFILE = true;

_TCTXVECT           CVS;

IMACounts           COUNTS;
_TLOG_FREQUENCIES   LOG_FREQUENCIES;

_TCOMPUTED_BLOSUM80 COMPUTED_BLOSUM80;
_TCOMPUTED_BLOSUM62 COMPUTED_BLOSUM62;
_TCOMPUTED_BLOSUM45 COMPUTED_BLOSUM45;
_TCOMPUTED_PSCORES_ COMPUTED_PSCORES_;
_TCOMPUTED_GONNET_  COMPUTED_GONNET;

_TRANS_PROBS        TRANSPROBS;

_LOSCORES           LOSCORES;


// =========================================================================
// _TCTXVECT: set constants related to processing context vectors
//
_TCTXVECT::_TCTXVECT()
:
    AVGLEN( 1/*7/*5*/ ),//length for averaging
    AVGCWGT( 1./*0.20/*0.40*/ ),//central weight in averaging
    AVGCOEFFS( AVGLEN, AVGCWGT ),//coefficients in averaging
    hAVGLEN( AVGLEN >> 1 ),
    //
    CTXLEN( 15/*9*/),//context length
    CWGT( 0.40/*0.30*/),//central weight of context coefficients
    MIX( true/*false*/),//mix normal vectors within the context
    AUG( false ),//augment vector 3-fold: into left, center, right partitions
    MEAN( true ),//mean vector in use
    //
    DIM( MIX? (NUMAA-1)*(AUG? 3: 1): CTXLEN*(NUMAA-1)),//dimensions of context vector
    KAPPA0( 0.1 ),//strictly positive
    NU0( double(DIM)),//must be greater than dim-1
    MU0( NUMAA-1 ),
    loKAPPA0( 1./KAPPA0 ),
    PowerNU0( -0.5*(NU0+2.)),
    CTERM(//constant term in calculating log probability of observations
      (0.5*DIM+PowerNU0)*(log(KAPPA0)-log(KAPPA0+2.)) )
{
    int ii;
    const int noeffress = NUMAA;
    double  probs[noeffress];
    double  sum = 0.0;
    for( ii = 1; ii <= DIM; ii++ )
        sum += log(NU0+1.-ii) - LN2;
    const_cast<double&>(CTERM) += sum;
    //
    if( MEAN ) {
        for( ii = 0; ii < noeffress; ii++ )
            probs[ii] = Robinson_PROBS[ii];
        //transform to normal
        ::LogitNormal2Normal( probs, noeffress, 1.e-1, false );
        for( ii = 0; ii < noeffress-1; ii++ )
            const_cast<Pslvector&>(MU0).SetValueAt( ii, probs[ii]);
    }
    //
    const_cast<CtxtCoefficients&>(AVGCOEFFS).FindCoefficients();

}

// -------------------------------------------------------------------------
// _TLOG_FREQUENCIES: precomputes log-frequency values
//
_TLOG_FREQUENCIES::_TLOG_FREQUENCIES()
{
    double  sum = 0.0;
    for( int v = 0; v < NUMFREQ; v++ )
    {
        if( 0 < v )
            sum += data[v] = log( v );
        else//factorial of zero (for integers) is equal to 1
            data[v] = 0.0;

        sums[v] = sum;
    }
}

// =========================================================================
// _TRANS_PROBS: constructor
//
_TRANS_PROBS::_TRANS_PROBS()
{
    effnos_[PS_M] = 1.0;
    effnos_[PS_I] = 0.0;
    effnos_[PS_D] = 0.0;
    InitOnce();
    Initialize();
    SetClass( Dirichlet );
}

// InitOnce: initialize data once
//
void _TRANS_PROBS::InitOnce()
{
    int st, states, n;
    double sum;

    memset( priors_, 0, sizeof( double ) * P_NSTATES );

    for( n = 0, states = 0; states < P_NSTATES; states += gTPTRANS_NTPS, n++ ) {
        sum = 0.0;
        for( st = states; st < states + gTPTRANS_NTPS && st < P_NSTATES; st++ )
            sum += gTPTRANS_ALPHAS[st];
        if( n < PS_NSTATES )
            sumalphas_[n] = sum;
    }
}

// Initialize: initialize or reset private data
//
void _TRANS_PROBS::Initialize()
{
    memset( pmestimators_, 0, sizeof( double ) * P_NSTATES );
    memset( logpmestimators_, 0, sizeof( double ) * P_NSTATES );
}

// SetClass: set class of priors
//
void _TRANS_PROBS::SetClass( TPriorClass value )
{
    double  dummy[P_NSTATES];
    class_ = value;

    memset( dummy, 0, sizeof( double ) * P_NSTATES );
    PME( &dummy );
    if( GetPMEstimators())
        memcpy( priors_, *GetPMEstimators(), sizeof( double ) * P_NSTATES );
}

// GetEffNoSequences: get effective no. observations for state `mid'
//
double _TRANS_PROBS::GetEffNoSequences( int mid ) const
{
    if( mid < 0 || PS_NSTATES <= mid )
        throw myruntime_error("_TRANS_PROBS: GetEffNoSequences: Invalid state.");
    return effnos_[mid];
}

// SetEffNoSequences: set effective no. observations for all states
//
void _TRANS_PROBS::SetEffNoSequences( double nm, double ni, double nd )
{
    effnos_[PS_M] = nm;
    effnos_[PS_I] = ni;
    effnos_[PS_D] = nd;
}

// PME: calculate posterior mean estimates of probabilities
//
void _TRANS_PROBS::PME( const double( *obsfreqs )[P_NSTATES] )
{
    Initialize();
    switch( GetClass()) {
        case Dirichlet: DirPME( obsfreqs ); break;
        default:
            throw myruntime_error( mystring( "_TRANS_PROBS: Unknown class of priors." ));
    };
    return;
}

// DirPME: apply Dirichlet priors to posterior mean estimate
// NOTE: have to have counts rather than frequencies to compute properly;
//  if not, effective number of sequences from which the observations were 
//  derived should be provided
//
void _TRANS_PROBS::DirPME( const double( *obsfreqs )[P_NSTATES] )
{
    int         n, st, states;
    const int   maxtr = gTPTRANS_NTPS;
    double      est, logest;
    double      sum, tsum;
    double      mid[PS_NSTATES];
    double      midobs, nobs;
    const double    accuracy = 1.0e-6;

    if( obsfreqs == NULL )
        throw myruntime_error( mystring( "_TRANS_PROBS: Memory access error." ));

    mid[PS_M] = GetEffNoSequences( PS_M );
    mid[PS_I] = GetEffNoSequences( PS_I );
    mid[PS_D] = GetEffNoSequences( PS_D );
    midobs = mid[PS_M] + mid[PS_I] + mid[PS_D];

    for( n = 0, states = 0; states < P_NSTATES; states += maxtr, n++ ) {
        sum = 0.0;
        for( st = states; st < states + maxtr && st < P_NSTATES; st++ ) {
//             if( st == P_ID || st == P_DI ) continue;//IGNORE states P_ID and P_DI
            sum += ( *obsfreqs )[st];
        }
        if( PS_NSTATES <= n ) {
            warning( "_TRANS_PROBS: Unexpected number of states." );
            break;
        }
        nobs = mid[n];
        tsum = nobs * sum + sumalphas_[n];
        if( 0.0 <= sum && 0.0 < tsum )
            for( st = states; st < states + maxtr && st < P_NSTATES; st++ ) {
                est = 0.0;
                logest = ( double )LOG_PROB_MIN;
//                 if( st != P_ID || st != P_DI ) {//IGNORE states P_ID and P_DI
                    est = ( nobs * ( *obsfreqs )[st] + gTPTRANS_ALPHAS[st] )/ tsum;
                    if( 0.0 < est ) logest = log( est );
//                 }
                pmestimators_[st] = est;
                logpmestimators_[st] = logest;
            }
    }

    //verify probability conservation
    for( n = 0, states = 0; states < P_NSTATES; states += maxtr, n++ ) {
        sum = 0.0;
        for( st = states; st < states + maxtr && st < P_NSTATES; st++ )
//             if( st != P_ID && st != P_DI )//omit to-be-ignored states
                sum += pmestimators_[st];
        if( sum < 1.0 - accuracy || sum > 1.0 + accuracy )
            throw myruntime_error( "_TRANS_PROBS: Transition probability estimators are not conserved." );
    }
}





// =========================================================================
//
void SetLOSCORES( const mystring& value, const char* filename )
{
    mystring msg = "Reading ";
    mystring lcval = value;
    lcval.lower();
    if( lcval == "blosum80")
        LOSCORES.SetClass( _LOSCORES::Class80 );
    else if( lcval == "blosum62")
        LOSCORES.SetClass( _LOSCORES::Class62 );
    else if( lcval == "blosum45")
        LOSCORES.SetClass( _LOSCORES::Class45 );
    else if( lcval == "gonnet")
        LOSCORES.SetClass( _LOSCORES::ClassGn );
    else if( lcval == "pscores") {
        msg += filename;
        message( msg.c_str());
        LOSCORES.SetClass( filename );
    }
    else {
        warning("Unknown substitution matrix class. Using default.");
    }
}

// =========================================================================
// _LOSCORES: constructor
//
_LOSCORES::_LOSCORES()
:   auxprobabs_( NULL )
{
    SetClass( ClassGn );
    ComputeLogProbabilities();
    InitLogProbabs( logprobs_1_ );
    InitLogProbabs( logprobs_2_ );
    InitLogProbabs( auxlogps_ );
}

// SetClass: sets specific class of scores and tries to read scores from
//     file
//
void _LOSCORES::SetClass( const char* filename )
{
    ReadScoresFile( filename );
    SetClass( ClassPS );
    ComputeLogProbabilities();
    COMPUTED_PSCORES_.Compute();
}

// ReadScoresFile: tries to read scores from a file
//
void _LOSCORES::ReadScoresFile( const char* filename )
{
    mystring    errstr1, errstr2, errstr3;
    mystring    fullname = filename;
    mystring    alt_tablename = mystring( GetFullParamDirname()) + DIRSEP + filename;
    mystring    altalt_taname = mystring( PROGDIR ) + DIRSEP +
                            UPDIR + DIRSEP +
                            GetParamDirectory() + DIRSEP +
                            filename;
    try {
        COUNTS.ReadScores( fullname.c_str());
        return;
    } catch( myexception const& ex ) {
        errstr1 = ex.what();
    }

    try {
        COUNTS.ReadScores( alt_tablename.c_str());
        return;
    } catch( myexception const& ex ) {
        errstr2 = ex.what();
    }

    try {
        COUNTS.ReadScores( altalt_taname.c_str());
        return;
    } catch( myexception const& ex ) {
        errstr3 = ex.what();
    }

    if( !errstr1.empty())   error( errstr1.c_str());
    if( !errstr2.empty())   error( errstr2.c_str());
    if( !errstr3.empty())   error( errstr3.c_str());

    throw myruntime_error( mystring( "_LOSCORES: Unable to read scores." ));
}

// StoreProbabilities_1: stores probabilities of the first profile
void _LOSCORES::StoreProbabilities_1( const double* probs )
{
    probs_1_ = probs;
    ComputeLogProbabilities( probs_1_, logprobs_1_ );
}
// RestoreProbabilities: restore probabilities
void _LOSCORES::RestoreProbabilities_1()
{
    probs_1_ = NULL;
    InitLogProbabs( logprobs_1_ );
}
// StoreProbabilities_2: stores probabilities of the second profile
void _LOSCORES::StoreProbabilities_2( const double* probs )
{
    probs_2_ = probs;
    ComputeLogProbabilities( probs_2_, logprobs_2_ );
}
// RestoreProbabilities: restore probabilities
void _LOSCORES::RestoreProbabilities_2()
{
    probs_2_ = NULL;
    InitLogProbabs( logprobs_2_ );
}
// StoreProbabilities: temporalily stores probabilities
void _LOSCORES::StoreProbabilities( const double* probs )
{
    auxprobabs_ = probs;
    ComputeLogProbabilities( auxprobabs_, auxlogps_ );
}
// RestoreProbabilities: restore probabilities
void _LOSCORES::RestoreProbabilities()
{
    auxprobabs_ = NULL;
    InitLogProbabs( auxlogps_ );
}

// InitLogProbabs: reset logs of tempralily stored probabilities
//
void _LOSCORES::InitLogProbabs( double* logprobs )
{
    if( !logprobs )
        return;
    for( unsigned char a = 0; a < NUMALPH; a++ )
        logprobs[a] = ( double ) LOG_PROB_MIN;
}

// PROBABility: background probability
//
double _LOSCORES::PROBABility( int a )
{
#ifdef __DEBUG__
    if( a < 0 || NUMALPH <= a )
        throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
#endif

    if( GetProbabilities())
        return GetProbabilities()[a];

    switch( GetClass()) {
        case Class80:
        case Class62:
        case Class45:   return Robinson_PROBS[a];
        case ClassGn:   return Robinson_PROBS[a];//Gonnet_PROBS[a];
        case ClassPS:   return COUNTS.GetBackProbability( a );
        default:
            throw myruntime_error( mystring( "_LOSCORES: Invalid class of scores." ));
    };
    return 0.0;
}

// LogPROBABility: log of background probability
//
double _LOSCORES::LogPROBABility( int a )
{
#ifdef __DEBUG__
    if( a < 0 || NUMALPH <= a )
        throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
#endif

    if( GetProbabilities())
        return auxlogps_[a];

    return logprobs_[a];
}

// PROBABILITY_1: background probability of the first profile
double _LOSCORES::PROBABILITY_1( int a )
{
    if( a < 0 || NUMALPH <= a )
        throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
    if( GetProbabilities_1())
        return probs_1_[a];
    return PROBABility( a );
}
// LogPROBABILITY_1: log of background probability of the first profile
double _LOSCORES::LogPROBABILITY_1( int a )
{
    if( a < 0 || NUMALPH <= a )
        throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
    if( GetProbabilities_1())
        return logprobs_1_[a];
    return LogPROBABility( a );
}

// PROBABILITY_2: background probability of the second profile
double _LOSCORES::PROBABILITY_2( int a )
{
    if( a < 0 || NUMALPH <= a )
        throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
    if( GetProbabilities_2())
        return probs_2_[a];
    return PROBABility( a );
}
// LogPROBABILITY_2: log of background probability of the second profile
double _LOSCORES::LogPROBABILITY_2( int a )
{
    if( a < 0 || NUMALPH <= a )
        throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
    if( GetProbabilities_2())
        return logprobs_2_[a];
    return LogPROBABility( a );
}

// ComputeLogProbabilities: recomputes logs of probabilities
//
void _LOSCORES::ComputeLogProbabilities()
{
    for( unsigned char a = 0; a < NUMALPH; a++ )
        if( 0.0 < PROBABility( a ))
            logprobs_[a] = log( PROBABility( a ));
        else
            logprobs_[a] = ( double ) LOG_PROB_MIN;
}

// ComputeLogProbabilities: compute logs of tempralily stored probabilities
//
void _LOSCORES::ComputeLogProbabilities( const double* probs, double* logprobs )
{
    if( !probs || !logprobs )
        throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));

    for( unsigned char a = 0; a < NUMALPH; a++ )
        if( 0.0 < probs[a] )
            logprobs[a] = log( probs[a] );
        else
            logprobs[a] = ( double )LOG_PROB_MIN;
}

// FreqRatio: return appropriate frequency ratio entry related to the 
//  table class specified
//
double _LOSCORES::FreqRatio( int a, int b )
{
#ifdef __DEBUG__
    if( a < 0 || NUMALPH <= a ||
        b < 0 || NUMALPH <= b )
        throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
#endif

//     if( GetProbabilities())
//         throw myruntime_error( mystring( "_LOSCORES: Restore probabilities first before calling FreqRatio." ));

    switch( GetClass()) {
        case Class80:   return BLOSUM80_FREQRATIOS[a][b];
        case Class62:   return BLOSUM62_FREQRATIOS[a][b];
        case Class45:   return BLOSUM45_FREQRATIOS[a][b];
        case ClassGn:   return GONNET_FREQRATIOS[a][b];
        case ClassPS:   return COUNTS.GetOddsRatio( a, b );//PSCORES_FREQRATIOS[a][b];
        default:
            throw myruntime_error( mystring( "_LOSCORES: Wrong class of scores to be used." ));
    };
    return -1.0;
}

// PrecomputedEntry: returns precomputed entry corresponding to the
//     appropriate table
//
double _LOSCORES::PrecomputedEntry( int a, int b )
{
#ifdef __DEBUG__
    if( a < 0 || NUMALPH <= a ||
        b < 0 || NUMALPH <= b )
        throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
#endif

//     if( GetProbabilities())
//         throw myruntime_error( mystring( "_LOSCORES: Restore probabilities first before calling PrecomputedEntry." ));

    switch( GetClass()) {
        case Class80:   return COMPUTED_BLOSUM80( a, b );
        case Class62:   return COMPUTED_BLOSUM62( a, b );
        case Class45:   return COMPUTED_BLOSUM45( a, b );
        case ClassGn:   return COMPUTED_GONNET( a, b );
        case ClassPS:   return COMPUTED_PSCORES_( a, b );
        default:
            throw myruntime_error( mystring("_LOSCORES: PrecomputedEntry: Invalid class of scores."));
    };
    return ( double )SCORE_MIN;
}

// Entry: returns entry rounded to the nearest integer
//
int _LOSCORES::Entry( int a, int b )
{
#ifdef __DEBUG__
    if( a < 0 || NUMALPH <= a ||
        b < 0 || NUMALPH <= b )
        throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
#endif

//     if( GetProbabilities())
//         throw myruntime_error( mystring( "_LOSCORES: Restore probabilities first before calling Entry." ));

    switch( GetClass()) {
        case Class80:   return BLOSUM80[a][b];
        case Class62:   return BLOSUM62[a][b];
        case Class45:   return BLOSUM45[a][b];
        case ClassPS:   return ( int )COUNTS.GetSubstScore( a, b );//PSCORES_[a][b];
        case ClassGn:   throw myruntime_error("_LOSCORES: Entry: Use precomputed values for Gonnet matrix.");
        default:
            throw myruntime_error( mystring( "_LOSCORES: Entry: Invalid class of scores." ));
    };
    return SCORE_MIN;
}

// StatisParam: returns appropriate statistical parameter as indicated by
//     field and gap scheme for the corresponding table
//
double _LOSCORES::StatisParam( int scheme, int field )
{
#ifdef __DEBUG__
    if( field < 0 || NumFields <= field )
        throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
#endif

//     if( GetProbabilities())
//         throw myruntime_error( mystring( "_LOSCORES: Restore probabilities first before calling StatisParam." ));

    switch( GetClass()) {
        case Class80:
            if( scheme < 0 || Num80Entries <= scheme )
                throw myruntime_error( mystring( "_LOSCORES: StatisParam: Memory access error." ));
            return BLOSUM80_VALUES[scheme][field];
        case Class62:
            if( scheme < 0 || NumEntries <= scheme )
                throw myruntime_error( mystring( "_LOSCORES: StatisParam: Memory access error." ));
            return BLOSUM62_VALUES[scheme][field];
        case Class45:
            if( scheme < 0 || Num45Entries <= scheme )
                throw myruntime_error( mystring( "_LOSCORES: StatisParam: Memory access error." ));
            return BLOSUM45_VALUES[scheme][field];
        case ClassGn:
            if( scheme < 0 || NumGnEntries <= scheme )
                throw myruntime_error( mystring( "_LOSCORES: StatisParam: Memory access error." ));
            return GONNET_VALUES[scheme][field];
        case ClassPS:
            break;
        default:
            throw myruntime_error( mystring( "_LOSCORES: StatisParam: Invalid class of scores." ));
    };

    if( scheme < 0 || NumPSEntries <= scheme )
        throw myruntime_error( mystring( "_LOSCORES: StatisParam: Memory access error." ));

    switch( field ) {
        case Lambda:return  COUNTS.GetStatLambda();
        case K:     return  COUNTS.GetStatK();
        case H:     return  COUNTS.GetStatH();

        case alpha: return  1.0;//return for compatibility; actually it is not used
        case beta:  return -1.0;//same as above
        case Open:
        case Extend:
        default:
            //should not be referenced
            throw myruntime_error( mystring( "_LOSCORES: StatisParam: Memory access error." ));
    }

    return -1.0;
}

// -------------------------------------------------------------------------
// _TCOMPUTED_BLOSUM80: precomputes BLOSUM80 values given frequency ratios
// -------------------------------------------------------------------------

_TCOMPUTED_BLOSUM80::_TCOMPUTED_BLOSUM80()
{
    for( unsigned char a = 0; a < NUMALPH; a++ )
        for( unsigned char b = 0; b < NUMALPH; b++ )
            if( BLOSUM80_FREQRATIOS[a][b] > 0.0 )
                data[a][b] = log( BLOSUM80_FREQRATIOS[a][b] ) / LN2 * Blosum80ScalingConstant;
            else
                data[a][b] = ( double ) SCORE_MIN;
}

// -------------------------------------------------------------------------
// _TCOMPUTED_BLOSUM62: precomputes BLOSUM62 values given frequency ratios
// -------------------------------------------------------------------------

_TCOMPUTED_BLOSUM62::_TCOMPUTED_BLOSUM62()
{
    for( unsigned char a = 0; a < NUMALPH; a++ )
        for( unsigned char b = 0; b < NUMALPH; b++ )
            if( BLOSUM62_FREQRATIOS[a][b] > 0.0 )
                data[a][b] = log( BLOSUM62_FREQRATIOS[a][b] ) / LN2 * Blosum62ScalingConstant;
            else
                data[a][b] = ( double ) SCORE_MIN;
}

// -------------------------------------------------------------------------
// _TCOMPUTED_BLOSUM45: precomputes BLOSUM45 values given frequency ratios
// -------------------------------------------------------------------------

_TCOMPUTED_BLOSUM45::_TCOMPUTED_BLOSUM45()
{
    for( unsigned char a = 0; a < NUMALPH; a++ )
        for( unsigned char b = 0; b < NUMALPH; b++ )
            if( BLOSUM45_FREQRATIOS[a][b] > 0.0 )
                data[a][b] = log( BLOSUM45_FREQRATIOS[a][b] ) / LN2 * Blosum45ScalingConstant;
            else
                data[a][b] = ( double ) SCORE_MIN;
}

// -------------------------------------------------------------------------
// _TCOMPUTED_GONNET_: precompute scaled Gonnet matrix values given
//  frequency ratios
//
_TCOMPUTED_GONNET_::_TCOMPUTED_GONNET_()
{
    unsigned char   a, b;
    for( a = 0; a < NUMALPH; a++ )
        for( b = 0; b < NUMALPH; b++ )
            if( 0.0 < GONNET_FREQRATIOS[a][b])
                data[a][b] = log( GONNET_FREQRATIOS[a][b] ) / LN2 * Gonnet_ScalingConstant;
            else
                data[a][b] = ( double )SCORE_MIN;
}

// -------------------------------------------------------------------------
// _TCOMPUTED_PSCORES_: PSCORES_ constructor: initializes scores
// -------------------------------------------------------------------------

_TCOMPUTED_PSCORES_::_TCOMPUTED_PSCORES_()
{
    for( unsigned char a = 0; a < NUMALPH; a++ )
        for( unsigned char b = 0; b < NUMALPH; b++ )
            data[a][b] = ( double ) SCORE_MIN;
}

// Compute: recomputes scores by using COUNTS information
//
void _TCOMPUTED_PSCORES_::Compute()
{
    for( unsigned char a = 0; a < NUMALPH; a++ )
        for( unsigned char b = 0; b < NUMALPH; b++ )
            if( 0.0 < COUNTS.GetOddsRatio( a, b ) && 0.0 < COUNTS.GetScoreScale())
                data[a][b] = log( COUNTS.GetOddsRatio( a, b )) / LN2 * COUNTS.GetScoreScale();
            else
                data[a][b] = ( double ) SCORE_MIN;
}

