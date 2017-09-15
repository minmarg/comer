/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include "libpro/srcsco/ProfileMatrix.h"
#include "IMACounts.h"


//minimum score: don't penalize it too much
const double    IMACounts::minscore_ = -20.0;
//scale factor of scores
double          IMACounts::scorescale_ = -1.0;
//automatically computed scale factor
bool            IMACounts::autoscale_ = true;


////////////////////////////////////////////////////////////////////////////
// CLASS IMACounts
//
// Constructor
//

IMACounts::IMACounts()
:   totresidcount_( 0.0 ),
    totpaircount_( 0.0 ),
    scoresread_( false ),
    scorentropy_( 0.0 ),
    expectedscore_( 0.0 ),
    minsubscore_( 0.0 ),
    maxsubscore_( 0.0 ),
    statlambda_( -1.0 ),
    statK_( -1.0 ),
    statH_( -1.0 )
{
    ResetResidCounts();
    ResetPairCounts();
    ResetTargetFrequencies();
    ResetMargProbabilities();
    ResetBackProbabilities();
    ResetSubstScores();
    ResetOddsRatios();
}

// -------------------------------------------------------------------------
// Destructor
// -------------------------------------------------------------------------

IMACounts::~IMACounts()
{
}

// -------------------------------------------------------------------------
// ResetResidCounts: Resets observed amino acid counts
// -------------------------------------------------------------------------

void IMACounts::ResetResidCounts()
{
    for( size_t n = 0; n < NUMALPH; n++ )
        obsresidues_[n] = 0.0;
}

// -------------------------------------------------------------------------
// ResetPairCounts: Resets observed amino acid pair counts
// -------------------------------------------------------------------------

void IMACounts::ResetPairCounts()
{
    for( size_t n = 0; n < NUMALPH; n++ )
        for( size_t m = 0; m < NUMALPH; m++ )
            obscounts_[n][m] = 0.0;
}

// -------------------------------------------------------------------------
// ResetPairCounts: Resets observed amino acid pair counts
// -------------------------------------------------------------------------

void IMACounts::ResetTargetFrequencies()
{
    for( size_t n = 0; n < NUMALPH; n++ )
        for( size_t m = 0; m < NUMALPH; m++ )
            targetfreqs_[n][m] = 0.0;
}

// -------------------------------------------------------------------------
// ResetMargProbabilities: Resets marginal amino acid probabilities
// -------------------------------------------------------------------------

void IMACounts::ResetMargProbabilities()
{
    for( size_t n = 0; n < NUMALPH; n++ )
        marginprobs_[n] = 0.0;
}

// -------------------------------------------------------------------------
// ResetBackProbabilities: Resets background amino acid probabilities
// -------------------------------------------------------------------------

void IMACounts::ResetBackProbabilities()
{
    for( size_t n = 0; n < NUMALPH; n++ )
        backgrprobs_[n] = 0.0;
}

// -------------------------------------------------------------------------
// ResetSubstScores: Resets substitution scores
// -------------------------------------------------------------------------

void IMACounts::ResetSubstScores()
{
    for( size_t n = 0; n < NUMALPH; n++ )
        for( size_t m = 0; m < NUMALPH; m++ )
            subscores_[n][m] = GetMinScore() / LN2;
}

// -------------------------------------------------------------------------
// ResetOddsRatios: Resets odds ratios
// -------------------------------------------------------------------------

void IMACounts::ResetOddsRatios()
{
    for( size_t n = 0; n < NUMALPH; n++ )
        for( size_t m = 0; m < NUMALPH; m++ )
            oddratios_[n][m] = 0.0;
}

// -------------------------------------------------------------------------
// SumUp: Sums up counts
// -------------------------------------------------------------------------

void IMACounts::SumUp( const IMACounts& right )
{
    size_t  n, m;

    for( n = 0; n < NUMALPH; n++ )
        IncResidCountBy( n, right.GetResidCount( n ));

    for( n = 0; n < NUMALPH; n++ )
        for( m = 0; m < NUMALPH; m++ )
            IncPairCountBy( n, m, right.GetPairCount( n, m ));
}

// -------------------------------------------------------------------------
// ComputeScores: Computes probabilities and derives scores
// -------------------------------------------------------------------------

void IMACounts::ComputeScores()
{
    size_t  effsize = GetEffectiveNoResids();
    double  totcount = 0.0;
    double  consv = 0.0;
    double  value;
    size_t  n, m;
    char    strbuf[BUF_MAX];

    for( n = 0; n < effsize; n++ )
        totcount += GetResidCount( n );

    if( totcount <= 0.0 )
        throw myruntime_error( mystring( "IMACounts: Total sum of counts is zero." ));

    SetTotalResidCount( totcount );
    totcount = 0.0;

    for( n = 0; n < effsize; n++ )
        for( m = 0; m < effsize; m++ )
            totcount += GetPairCount( n, m );

    if( totcount <= 0.0 )
        throw myruntime_error( mystring( "IMACounts: Total sum of substitution counts is zero." ));

    SetTotalPairCount( totcount );
    ResetMargProbabilities();

    //background probabilities...
    //
    for( n = 0; n < effsize; n++ ) {
        value = GetResidCount( n ) / GetTotalResidCount();
        SetBackProbability( n, value );
        consv += value;
    }

    if( consv < 0.9999 || consv > 1.0001 ) {
        sprintf( strbuf, "Background probabilities found to be not conserved: %f.", consv );
        warning( strbuf );
    }

    consv = 0.0;
    //target frequencies...
    //marginal probabilities...
    //
    for( n = 0; n < effsize; n++ )
        for( m = 0;  m <= n; m++ ) {
            if( n == m ) {
                value = GetPairCount( n, m ) / GetTotalPairCount();
            }
            else {
                value = ( GetPairCount( n, m ) + GetPairCount( m, n )) / TIMES2( GetTotalPairCount());
                SetTargetFrequency( m, n, value );
                IncMargProbabilityBy( m, value );
                consv += value;
            }
            SetTargetFrequency( n, m, value );
            IncMargProbabilityBy( n, value );

            consv += value;
        }

    if( consv < 0.9999 || consv > 1.0001 ) {
        sprintf( strbuf, "Target frequencies found to be not conserved: %f.", consv );
        warning( strbuf );
    }

    SetScoresRead( false );

    double  score;
    double  product;
    double  odds;
    double  value1, value2;
    double  entropy = 0.0;
    double  expected = 0.0;
    //bit scores...
    //
    for( n = 0; n < effsize; n++ )
        for( m = 0;  m <= n; m++ ) {
            odds = 0.0;
            score = GetMinScore();

            product = GetMargProbability( n ) * GetMargProbability( m );
            if( 0.0 < product )
                odds = GetTargetFrequency( n, m ) / product;
            if( 0.0 < odds )
                score = log( odds );

            score /= LN2;

            if( score < GetMinSubstScore())
                SetMinSubstScore( score );

            if( GetMaxSubstScore() < score )
                SetMaxSubstScore( score );

            value1 = score * GetTargetFrequency( n, m );
            value2 = score * product;

            entropy += value1;
            expected += value2;
            SetSubstScore( n, m, score );

            if( n != m ) {
                entropy += value1;
                expected += value2;
                SetSubstScore( m, n, score );
            }
        }

    SetEntropy( entropy );
    SetExpected( expected );
    DetermineScale();

    ComputeAbstractScores();
}

// -------------------------------------------------------------------------
// DetermineScale: Determines scale factor by using entropy value
// -------------------------------------------------------------------------

void IMACounts::DetermineScale() const
{
    if( !GetAutoScoreScale() && 0.0 < GetScoreScale())
        return;

    if( GetEntropy() <= 0.0 ) {
        SetScoreScale( 1.0 );
        return;
    }

    double value = rint( 2.0 / sqrt( GetEntropy()));

    if( value < 2.0 )
        value = 2.0;

    SetScoreScale( value );
}

// -------------------------------------------------------------------------
// DetermineScaleThroughScores: Determines scale factor by using score
//     values
// -------------------------------------------------------------------------

void IMACounts::DetermineScaleThroughScores() const
{
    size_t  res = 0;
    double  score = GetSubstScore( res, res );
    double  target = GetTargetFrequency( res, res );
    double  margin = SQUARE( GetMargProbability( res ));
    double  scale = -1.0;

    if( score < GetMinScore() || margin <= 0.0 || target <= 0.0 ) {
        scale = rint( score * LN2 / GetMinScore());
    }
    else {
        scale = rint( score * LN2 / log( target / margin ));
    }

    SetScoreScale( scale );
}

// -------------------------------------------------------------------------
// ComputeAbstractScores: Computes scores for B, Z, X
// -------------------------------------------------------------------------

void IMACounts::ComputeAbstractScores()
{
    static int  sB = HashAlphSymbol('B');
    static int  sZ = HashAlphSymbol('Z');
    static int  rN = HashAlphSymbol('N');
    static int  rD = HashAlphSymbol('D');
    static int  rQ = HashAlphSymbol('Q');
    static int  rE = HashAlphSymbol('E');

    size_t  effsize = GetEffectiveNoResids();
    size_t  n, m;
    double  value;
    double  value1, probs;
    double  lminscore = GetMinScore() / LN2;


    for( n = 0; n < effsize; n++ )
    {
        value = GetBackProbability( rN ) * GetSubstScore( rN, n ) +
                GetBackProbability( rD ) * GetSubstScore( rD, n );
        probs = GetBackProbability( rN ) + GetBackProbability( rD );
        SetSubstScore( sB, n, probs? value / probs: lminscore );
        SetSubstScore( n, sB, probs? value / probs: lminscore );

        value = GetBackProbability( rQ ) * GetSubstScore( rQ, n ) +
                GetBackProbability( rE ) * GetSubstScore( rE, n );
        probs = GetBackProbability( rQ ) + GetBackProbability( rE );
        SetSubstScore( sZ, n, probs? value / probs: lminscore );
        SetSubstScore( n, sZ, probs? value / probs: lminscore );
    }

    value  = GetBackProbability( rN ) * GetBackProbability( rN ) * GetSubstScore( rN, rN );
    value += GetBackProbability( rD ) * GetBackProbability( rD ) * GetSubstScore( rD, rD );
    value += TIMES2( GetBackProbability( rN ) * GetBackProbability( rD ) * GetSubstScore( rN, rD ));
    probs  = SQUARE( GetBackProbability( rN ) + GetBackProbability( rD ));
    SetSubstScore( sB, sB, probs? value / probs: lminscore );

    value  = GetBackProbability( rQ ) * GetBackProbability( rQ ) * GetSubstScore( rQ, rQ );
    value += GetBackProbability( rE ) * GetBackProbability( rE ) * GetSubstScore( rE, rE );
    value += TIMES2( GetBackProbability( rQ ) * GetBackProbability( rE ) * GetSubstScore( rQ, rE ));
    probs  = SQUARE( GetBackProbability( rQ ) + GetBackProbability( rE ));
    SetSubstScore( sZ, sZ, probs? value / probs: lminscore );

    value  = TIMES2( GetBackProbability( rN ) * GetBackProbability( rQ ) * GetSubstScore( rN, rQ ));
    value += TIMES2( GetBackProbability( rN ) * GetBackProbability( rE ) * GetSubstScore( rN, rE ));
    value += TIMES2( GetBackProbability( rD ) * GetBackProbability( rQ ) * GetSubstScore( rD, rQ ));
    value += TIMES2( GetBackProbability( rD ) * GetBackProbability( rE ) * GetSubstScore( rD, rE ));
    probs  = TIMES2(( GetBackProbability( rN ) + GetBackProbability( rD )) *
                    ( GetBackProbability( rQ ) + GetBackProbability( rE )));
    SetSubstScore( sB, sZ, probs? value / probs: lminscore );
    SetSubstScore( sZ, sB, probs? value / probs: lminscore );


    value1 = 0.0;
    probs = 0.0;
    for( n = 0; n < effsize; n++ )
    {
        value = 0.0;
        for( m = 0; m < effsize; m++ ) {
            value += GetBackProbability( m ) * GetSubstScore( n, m );
            value1 += GetBackProbability( n ) * GetBackProbability( m ) * GetSubstScore( n, m );
            probs += GetBackProbability( n ) * GetBackProbability( m );
        }
        SetSubstScore( n, X, value );
        SetSubstScore( X, n, value );
    }
    SetSubstScore( X, X, probs? value1 / probs: lminscore );


    value = value1 = 0.0;
    for( n = 0; n < effsize; n++ )
    {
        value += GetBackProbability( n ) * GetBackProbability( rN ) * GetSubstScore( n, rN );
        value += GetBackProbability( n ) * GetBackProbability( rD ) * GetSubstScore( n, rD );
        value1 += GetBackProbability( n ) * GetBackProbability( rQ ) * GetSubstScore( n, rQ );
        value1 += GetBackProbability( n ) * GetBackProbability( rE ) * GetSubstScore( n, rE );
    }
    probs = GetBackProbability( rN ) + GetBackProbability( rD );
    SetSubstScore( sB, X, probs? value / probs: lminscore );
    SetSubstScore( X, sB, probs? value / probs: lminscore );

    probs = GetBackProbability( rQ ) + GetBackProbability( rE );
    SetSubstScore( sZ, X, probs? value1 / probs: lminscore );
    SetSubstScore( X, sZ, probs? value1 / probs: lminscore );

    for( n = 0; n < NUMALPH - 1; n++ )
    {
        SetSubstScore( n, ASTERISK, GetMinSubstScore());
        SetSubstScore( ASTERISK, n, GetMinSubstScore());
    }
    SetSubstScore( ASTERISK, ASTERISK, 1.0 / GetScoreScale());
}

// -------------------------------------------------------------------------
// ComputeParameters: Computes statistical parameters
// -------------------------------------------------------------------------

void IMACounts::ComputeParameters()
{
    Configuration   config[NoSchemes];  //no filename, raise an error if attempting reading or writing
    SetUngappedParams( config[ProcomUngapped] );//fill parameters with values of ungapped configuration
    SetUngappedParams( config[ProcomGapped] );  //make gapped configuration identical to ungapped one

    const size_t    effsize = GetEffectiveNoResids();
    char            resids[effsize];
    double          locscaledscores[NUMALPH][NUMALPH];  //scaled scores
    size_t          n, m;

    for( n = 0; n < effsize; n++ )
        resids[n] = n;

    for( n = 0; n < NUMALPH; n++ )
        for( m = 0; m < NUMALPH; m++ )
            locscaledscores[n][m] = GetScoreScale() * GetSubstScore( n, m );

    //store probabilities needed to compute parameters
//     LOSCORES.StoreProbabilities( GetBackProbabilities());
    LOSCORES.StoreProbabilities( GetMargProbabilities());

    ProfileMatrix   natable( locscaledscores, resids, effsize, config/*, AbstractScoreMatrix::NoScaling */);
    natable.ComputeStatisticalParameters();

    //restore temporalily stored probabilities
    LOSCORES.RestoreProbabilities();

    SetStatLambda( natable.GetLambda());
    SetStatK( natable.GetK());
    SetStatH( natable.GetEntropy());
}

// =========================================================================
// Read: Reads observed amino acid count values from file
// =========================================================================

void IMACounts::Read( const char* filename )
{
    FILE*   fp = NULL;

    if( !filename || !strlen( filename ))
        throw myruntime_error( mystring( "IMACounts: Wrong filename." ));

    fp = fopen( filename, "r" );

    if( fp == NULL )
        throw myruntime_error( mystring( "IMACounts: Failed to open file for reading." ));

    ResetResidCounts();
    ResetPairCounts();

    mystring        errstr;
    size_t          len;
    const size_t    locsize = TIMES2( KBYTE );
    char            locbuffer[locsize] = {0};

    size_t          countsmax = GetEffectiveNoResids(); //number of counts per line to be read from file

    try {
        SkipComments( fp, locbuffer, locsize, &len );

        if( feof( fp ))
            throw myruntime_error( mystring( "IMACounts: No any counts read from file." ));

        ReadCounts( locbuffer, len, GetResidCountsAddr(), GetEffectiveNoResids());

        for( size_t n = 0; n < countsmax; n++ )
        {
            SkipComments( fp, locbuffer, locsize, &len );
            if( feof( fp ))
                throw myruntime_error( mystring( "IMACounts: No substitution counts read from file." ));
            ReadCounts( locbuffer, len, GetPairCountsAddrAt( n ), GetEffectiveNoResids());
        }
        if( !feof( fp ))
            ReadFooter( fp, locbuffer, locsize, &len );

    } catch( myexception const& ex )
    {
        errstr = ex.what();
    }

    fclose( fp );

    if( !errstr.empty())
        throw myruntime_error( errstr );
}

// -------------------------------------------------------------------------
// ReadFooter: Reads and checks footer
// -------------------------------------------------------------------------

void IMACounts::ReadFooter( FILE* fp, char* buffer, size_t bufsize, size_t* readlen )
{
    size_t  len;
    char*   p = buffer;
    const char* separator = "//";
    size_t      seplen = strlen( separator );

#ifdef __DEBUG__
    if( !fp || !buffer || !readlen )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif

    p = fgets( buffer, bufsize, fp );
    len = *readlen = strlen( buffer );

    if( p == NULL && feof( fp ))
        throw myruntime_error( mystring( "IMACounts: Missing end separator after data." ));

    if( p == NULL && ferror( fp ))
        throw myruntime_error( mystring( "IMACounts: Reading error." ));

    if(!( seplen <= len && strncmp( buffer, separator, seplen ) == 0 ))
        throw myruntime_error( mystring( "IMACounts: Wrong file format." ));
}

// -------------------------------------------------------------------------
// SkipComments: Skips comments from reading the file
// -------------------------------------------------------------------------

void IMACounts::SkipComments( FILE* fp, char* buffer, size_t bufsize, size_t* readlen )
{
    size_t  len;
    bool    fulline = true;
    char*   p = buffer;

#ifdef __DEBUG__
    if( !fp || !buffer || !readlen )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif

    while( !feof( fp )) {
        p = fgets( buffer, bufsize, fp );

        len = *readlen = strlen( buffer );

        if( p == NULL && feof( fp ))
            break;

        if( p == NULL && ferror( fp ))
            throw myruntime_error( mystring( "IMACounts: Reading error." ));

        if( fulline ) {
            for( p = buffer; *p == ' ' || *p == '\t'; p++ );
            if( *p != '#' && *p != '\n' && *p != '\r' && len )
                break;
        }
        if( len && ( buffer[len-1] != '\n' && buffer[len-1] != '\r' ))
            fulline = false;
    }
}

// -------------------------------------------------------------------------
// ReadCounts: Reads weighted observed counts to buffer given
// -------------------------------------------------------------------------

void IMACounts::ReadCounts(
    const char* readfrom,
    size_t      readlen,
    double*     membuffer,
    size_t      no_cnts,
    bool        allownegs )
{
//     size_t      left = readlen;
    const char* p = readfrom;
    const char* pbeg = NULL;
    const char* pend = NULL;

    mystring    number;
    double      tmpval;
    char*       paux;

    size_t      current = 0;    //serial number of counts
    size_t      countsmax = no_cnts; //GetEffectiveNoResids(); //number of counts to be read from file

#ifdef __DEBUG__
    if( !readfrom || !membuffer )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif

    while( *p && *p != '\n' && *p != '\r') {
        for( ; *p == ' ' || *p == '\t'; p++ );
        for( pbeg = p; *p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r'; p++ );

        if( pbeg == p )
            continue;

        if( countsmax <= current )
            throw myruntime_error( mystring( "IMACounts: Too many numbers in one line." ));

        number = mystring( pbeg, size_t( p - pbeg ));

        if( !number.empty()) {
            errno = 0;  //NOTE: Thread unsafe
            tmpval = strtod( number.c_str(), &paux );

            if( errno || *paux )
                throw myruntime_error( mystring( "IMACounts: Invalid number read from file." ));
            if( !allownegs && tmpval < 0.0 )
                throw myruntime_error( mystring( "IMACounts: Negative number read from file." ));

            membuffer[current++] = tmpval;
        }
    }

    if( current < countsmax )
        throw myruntime_error( mystring( "IMACounts: Not enough amount of numbers in one line." ));
}

// =========================================================================
// Write: Writes observed amino acid count values to file
// =========================================================================

void IMACounts::Write( FILE* fp )
{
    if( fp == NULL )
        return;

    size_t  effecout = GetEffectiveNoResids();

    fprintf( fp, "## Weighted observed counts\n##\n##" );

    fprintf( fp, "%13c ", DehashCode( 0 ));
    for( size_t n = 1; n < effecout; n++ )
        fprintf( fp, "%15c ", DehashCode( n ));

    fprintf( fp, "\n" );

    for( size_t n = 0; n < effecout; n++ )
        fprintf( fp, " %15.4lf", GetResidCount( n ));

    fprintf( fp, "\n## Weighted observed substitution counts\n##\n" );

    for( size_t n = 0; n < effecout; n++ ) {
        for( size_t m = 0; m < effecout; m++ )
            fprintf( fp, " %15.4lf", GetPairCount( n, m ));
        fprintf( fp, "\n" );
    }

    fprintf( fp, "//\n" );
}

// =========================================================================
// ReadScores: Reads scores and related information from file
// =========================================================================

void IMACounts::ReadScores( const char* filename )
{
    FILE*   fp = NULL;

    if( !filename || !strlen( filename ))
        throw myruntime_error( mystring( "IMACounts: Wrong filename." ));

    fp = fopen( filename, "r" );

    if( fp == NULL )
        throw myruntime_error( mystring( "IMACounts: Failed to open file " ) + filename + " for reading." );

    ResetTargetFrequencies();
    ResetMargProbabilities();
    ResetBackProbabilities();
    ResetSubstScores();
    ResetOddsRatios();

    SetScoresRead( true );

    mystring        errstr;
    size_t          len, n;
    const size_t    locsize = TIMES2( KBYTE );
    char            locbuffer[locsize] = {0};

    size_t          effecout = GetEffectiveNoResids();  //effective number of residues to be read
    size_t          printcnt = NUMALPH - 1; //number of residues to be read for scores

    try {
        SkipComments( fp, locbuffer, locsize, &len );
        if( feof( fp ))
            throw myruntime_error( mystring( "IMACounts: No any scores read from file." ));

        ReadStatisParams( locbuffer, len );

        for( n = 0; n < printcnt; n++ )
        {
            SkipComments( fp, locbuffer, locsize, &len );
            if( feof( fp ))
                throw myruntime_error( mystring( "IMACounts: No scores read from file." ));
            ReadCounts( locbuffer, len, GetSubstScoresAddrAt( n ), printcnt, true );
        }

        SkipComments( fp, locbuffer, locsize, &len );
        if( feof( fp ))
            throw myruntime_error( mystring( "IMACounts: No background probabilities read from file." ));
        ReadCounts( locbuffer, len, GetBackProbsAddr(), effecout );

        SkipComments( fp, locbuffer, locsize, &len );
        if( feof( fp ))
            throw myruntime_error( mystring( "IMACounts: No marginal probabilities read from file." ));
        ReadCounts( locbuffer, len, GetMargProbsAddr(), effecout );

        for( n = 0; n < effecout; n++ )
        {
            SkipComments( fp, locbuffer, locsize, &len );
            if( feof( fp ))
                throw myruntime_error( mystring( "IMACounts: No target frequencies read from file." ));
            ReadCounts( locbuffer, len, GetTargetFreqsAddrAt( n ), effecout );
        }

        for( n = 0; n < printcnt; n++ )
        {
            SkipComments( fp, locbuffer, locsize, &len );
            if( feof( fp ))
                throw myruntime_error( mystring( "IMACounts: No odds ratios read from file." ));
            ReadCounts( locbuffer, len, GetOddsRatiosAddrAt( n ), printcnt );
        }


        if( !feof( fp ))
            ReadFooter( fp, locbuffer, locsize, &len );

    } catch( myexception const& ex )
    {
        errstr = ex.what();
    }

    fclose( fp );

    if( !errstr.empty())
        throw myruntime_error( errstr );

    DetermineScaleThroughScores();
}

// -------------------------------------------------------------------------
// ReadStatisParams: Reads values of statistical parameters
// -------------------------------------------------------------------------

void IMACounts::ReadStatisParams( const char* readfrom, size_t readlen )
{
    const char* p = readfrom;

#ifdef __DEBUG__
    if( !readfrom )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));
#endif

    ReadAssocValue( &p, GetLambdaName(),GetLambdaAddr());
    ReadAssocValue( &p, GetKName(),     GetKAddr());
    ReadAssocValue( &p, GetHName(),     GetHAddr());

    for( ; *p == ' ' || *p == '\t'; p++ );
    if( *p && ( *p != '\n' && *p != '\r' ))
        throw myruntime_error( mystring( "IMACounts: Unrecognized parameters." ));
}

// -------------------------------------------------------------------------
// ReadStatisParams: Reads values of statistical parameters
// -------------------------------------------------------------------------

void IMACounts::ReadAssocValue( const char** p, mystring name, double* value )
{
    if( p == NULL || *p == NULL || value == NULL )
        throw myruntime_error( mystring( "IMACounts: Memory access error." ));

    const char* pbeg = NULL;
    const char* pend = NULL;
    mystring    locvalue;
    double      tmpval;
    char*       paux;

    for( ; **p == ' ' || **p == '\t'; (*p)++ );
    for( pbeg = *p; **p && **p != ' ' && **p != '\t' && **p != '\n' && **p != '\r'; (*p)++ );

    if( pbeg == *p )
        throw myruntime_error( mystring( "IMACounts: No parameter name found." ));

    locvalue = mystring( pbeg, size_t( *p - pbeg ));

    if( !name.empty()) {
        if( name != locvalue )
            throw myruntime_error( mystring( "IMACounts: Invalid parameter name." ));
        for( ; **p == ' ' || **p == '\t'; (*p)++ );
        if( **p != '=' )
            throw myruntime_error( mystring( "IMACounts: No assignment operator." ));
        (*p)++;
    }

    for( ; **p == ' ' || **p == '\t'; (*p)++ );
    for( pbeg = *p; **p && **p != ' ' && **p != '\t' && **p != ',' && **p != '\n' && **p != '\r'; (*p)++ );

    if( pbeg == *p )
        throw myruntime_error( mystring( "IMACounts: No parameter value found." ));

    locvalue = mystring( pbeg, size_t( *p - pbeg ));

    if( locvalue.empty())
        throw myruntime_error( mystring( "IMACounts: No parameter value found in file." ));

    errno = 0;  //NOTE: Thread unsafe
    tmpval = strtod( locvalue.c_str(), &paux );

    if( errno || *paux )
        throw myruntime_error( mystring( "IMACounts: Invalid parameter value read from file." ));
    if( tmpval < 0.0 )
        throw myruntime_error( mystring( "IMACounts: Negative parameter value read from file." ));

    if( **p == ',' )
        (*p)++;
    *value = tmpval;
}

// =========================================================================
// WriteScores: Writes scores and related information to file
// =========================================================================

void IMACounts::WriteScores( FILE* fp )
{
    if( fp == NULL )
        return;

    size_t  effecout = GetEffectiveNoResids();
    size_t  printcnt = NUMALPH - 1;
    size_t  n, m;

    if( 0.0 < GetScoreScale())
        fprintf( fp, "## Scores in 1/%d bit units\n", ( int )rint( GetScoreScale()));
    fprintf( fp, "## Entropy = %7.4lf, Expected = %8.4lf\n", GetEntropy(), GetExpected());
    fprintf( fp, "## Minimum = %4d, Maximum = %4d\n",
                ( int )rint( GetScoreScale() * GetMinSubstScore()),
                ( int )rint( GetScoreScale() * GetMaxSubstScore()));
    fprintf( fp, "## Reference scoring statistical parameters:\n" );
    fprintf( fp, "     Lambda = %6.4lf, K = %6.4lf, H = %6.4lf\n", GetStatLambda(), GetStatK(), GetStatH());
    fprintf( fp, "# " );

    fprintf( fp, "%c ", DehashCode( 0 ));
    for( size_t n = 1; n < printcnt; n++ )
        fprintf( fp, "%3c ", DehashCode( n ));

    fprintf( fp, "\n" );

    for( size_t n = 0; n < printcnt; n++ ) {
        for( size_t m = 0; m < printcnt; m++ ) {
            if( 0.0 < GetScoreScale() && !GetScoresRead())
                fprintf( fp, " %3d", ( int )rint( GetScoreScale() * GetSubstScore( n, m )));
            else
                fprintf( fp, " %3d", ( int )rint( GetSubstScore( n, m )));
        }
        fprintf( fp, "\n" );
    }

    fprintf( fp, "##\n## Background Probabilities:\n# " );

    fprintf( fp, "%6c ", DehashCode( 0 ));
    for( size_t n = 1; n < effecout; n++ )
        fprintf( fp, "%8c ", DehashCode( n ));

    fprintf( fp, "\n" );

    for( size_t n = 0; n < effecout; n++ )
        fprintf( fp, " %8.6lf", GetBackProbability( n ));

    fprintf( fp, "\n##\n## Marginal Probabilities:\n# \n" );

    for( size_t n = 0; n < effecout; n++ )
        fprintf( fp, " %8.6lf", GetMargProbability( n ));

    fprintf( fp, "\n" );

    fprintf( fp, "\n##\n## Target Frequencies:\n# " );

    fprintf( fp, "%6c ", DehashCode( 0 ));
    for( size_t n = 1; n < effecout; n++ )
        fprintf( fp, "%8c ", DehashCode( n ));

    fprintf( fp, "\n" );

    for( size_t n = 0; n < effecout; n++ ) {
        for( size_t m = 0; m < effecout; m++ )
            fprintf( fp, " %8.6lf", GetTargetFrequency( n, m ));
        fprintf( fp, "\n" );
    }

    fprintf( fp, "##\n## Odds Ratios:\n# " );

    fprintf( fp, "%5c ", DehashCode( 0 ));
    for( size_t n = 1; n < printcnt; n++ )
        fprintf( fp, "%7c ", DehashCode( n ));

    fprintf( fp, "\n" );

    for( size_t n = 0; n < printcnt; n++ ) {
        for( size_t m = 0; m < printcnt; m++ ) {
//             if( n < effecout && m < effecout &&
//                 0.0 < GetMargProbability( n ) &&
//                 0.0 < GetMargProbability( m ))
//                 fprintf( fp, " %7.4lf", GetTargetFrequency( n, m ) / ( GetMargProbability( n ) * GetMargProbability( m )));
//             else
            if( GetScoresRead())
                fprintf( fp, " %7.4lf", GetOddsRatio( n, m ));
            else
                fprintf( fp, " %7.4lf", exp2( GetSubstScore( n, m )));
        }
        fprintf( fp, "\n" );
    }

    fprintf( fp, "//\n" );
}

