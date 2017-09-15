/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "version.h"
#include "mystring.h"
#include "myexcept.h"

#include "ext/psl.h"
#include "ext/gamma.h"

#include "logitnormal.h"
#include "CtxtCoefficients.h"

#include "GapScheme.h"
#include "DistributionMatrix.h"
#include "Serializer.h"



// static character buffer...
char ExtendedDistributionMatrix::private_buffer[MAX_DESCRIPTION_LENGTH];


// -------------------------------------------------------------------------
// default constructor: values are initialized to zero, no memory allocation
//     is performed
// -------------------------------------------------------------------------

DistributionMatrix::DistributionMatrix()
:   values( 0 ),
    aacids( 0 ),
    columns( 0 ),
    allocated( 0 )
{
    init();
}

// -------------------------------------------------------------------------
// destructor: deallocates memory if it was allocated before
// -------------------------------------------------------------------------

DistributionMatrix::~DistributionMatrix()
{
    destroy();
}

// -------------------------------------------------------------------------
// init: initialization of members
// -------------------------------------------------------------------------

void DistributionMatrix::init()
{
    values = NULL;
    aacids = NULL;
    columns = 0;
    allocated = 0;
}

// -------------------------------------------------------------------------
// IsCompatible: verifies whether this matrix is compatible with another
//     one with respect to the amino acid vectors both matrices own
// -------------------------------------------------------------------------

bool DistributionMatrix::IsCompatible( const DistributionMatrix& one ) const
{
    if( GetColumns() != one.GetColumns())
        return false;

    for( int n = 0; n < GetColumns(); n++ )
        if( aacids[n] != one.aacids[n] )
            return false;

    return true;
}

// -------------------------------------------------------------------------
// IsCompatible: verifies whether this matrix is compatible with gap opening
//     cost vector in amino acid composition both types of structures have
// -------------------------------------------------------------------------

bool DistributionMatrix::IsCompatible( const GapScheme& goc ) const
{
    if( GetColumns() != goc.GetColumns())
        return false;

    for( int n = 0; n < GetColumns(); n++ )
        if( aacids[n] != goc.AAcid( n ))
            return false;

    return true;
}

// -------------------------------------------------------------------------
// destroy: deallocate memory and reset values
// -------------------------------------------------------------------------

void DistributionMatrix::destroy()
{
    if( values ) { free( values ); values = NULL; }
    if( aacids ) { free( aacids ); aacids = NULL; }
    columns = 0;
    allocated = 0;
}

// -------------------------------------------------------------------------
// Clear: erases all information contained in this class but leaves the
//     space allocated
// -------------------------------------------------------------------------

void DistributionMatrix::Clear()
{
    if( allocated ) {
        memset( values, 0, sizeof( double ) * NUMALPH * allocated );
        memset( aacids, 0, sizeof( char ) * allocated );
    }

    columns = 0;
}

// -------------------------------------------------------------------------
// reallocate: allocate necessary memory for the matrix 'values' and
//     vector 'aacids'
// -------------------------------------------------------------------------

void DistributionMatrix::reallocate( int howmuch )
{
    double   ( *tmp_values )[NUMALPH];
    char*       tmp_aacids;

    if( howmuch <= allocated )
        return;

    if( allocated <= 0 ) {
        tmp_values = ( double(*)[NUMALPH] )malloc( sizeof( double ) * NUMALPH * howmuch );
        tmp_aacids = ( char* )malloc( sizeof( char ) * howmuch );

    } else {
        tmp_values = ( double(*)[NUMALPH] )realloc( values, sizeof( double ) * NUMALPH * howmuch );
        tmp_aacids = ( char* )realloc( aacids, sizeof( char ) * howmuch );
    }

    if( !tmp_values || !tmp_aacids )
        throw myruntime_error( mystring( "DistributionMatrix: Not enough memory." ));

    values = tmp_values;
    aacids = tmp_aacids;

    // fill uninitialized memory with zeros
    memset( values + allocated, 0, sizeof( double ) * NUMALPH * ( howmuch - allocated ));
    memset( aacids + allocated, 0, sizeof( char ) * ( howmuch - allocated ));

    allocated = howmuch;
}

// -------------------------------------------------------------------------
// Push: pushes a vector of values to be saved in the structure
// -------------------------------------------------------------------------

void DistributionMatrix::Push( const double posvalues[NUMALPH], char aa )
{
    if( allocated < GetColumns() + 1 ) {
        reallocate(  allocated * 2 );
    }

    PushAt( posvalues, aa, GetColumns());
}

// -------------------------------------------------------------------------
// PushAt: pushes a vector of values at the given position
// -------------------------------------------------------------------------

void DistributionMatrix::PushAt( const double posvalues[NUMALPH], char aa, int position )
{
    if( allocated < position + 1 )
        throw myruntime_error(
            mystring( "DistributionMatrix: Failed to insert vector." ));

    for( int i = 0; i < NUMALPH; i++ )
        values[position][i] = posvalues[i];

    aacids[position] = aa;

    if( GetColumns() < position + 1 )
        SetColumns( position + 1 );
}

// -------------------------------------------------------------------------
// OutputMatrix: output all information this matrix possesses
// -------------------------------------------------------------------------

void DistributionMatrix::OutputMatrix( const char* filename ) const
{
    int     efective_number = NUMALPH - 1; //effective number of residues
    //
    FILE*   fp = stdout;
    size_t  l = 0;

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw myruntime_error(
            mystring( "DistributionMatrix: Failed to open file to write matrix." ));

    fprintf( fp, "%43c Weighted observed frequencies\n", 32 );
    fprintf( fp, "%9c", 32 );

    for( char r = 0; r < efective_number; r++ )
        fprintf( fp, "%4c", DehashCode( r ) );

    for( int p = 0; p < GetColumns(); p++ ) {
        fprintf( fp, "\n%5d %c   ", ++l, DehashCode(( *this )[p] ));

        for( char r = 0; r < efective_number; r++ )
            fprintf( fp, "%4d", ( int )rint( 100 * ( *this )( p, r )) );

    }
    fprintf( fp, "\n" );

    if( fp != stdout )
        fclose( fp );
}

// -------------------------------------------------------------------------
// CheckIntegrity: checks integrity of the structure; teh amino acid symbols
//     must be valid
// -------------------------------------------------------------------------

void DistributionMatrix::CheckIntegrity() const
{
    for( int n = 0; n < GetColumns(); n++ )
        if( NUMALPH <= ( *this )[n] )
            throw myruntime_error(
                mystring( "DistributionMatrix: Data corruption." ));
}

// -------------------------------------------------------------------------
// CheckForAllZeros: verifies wether there extists positions with all
//     values of zero; if so, the appropriate value is set to 1, indicating
//     100 percent amiino acid frequency in that position.
//     This is necessary for scoring matrix made of two profiles to take its
//     effect.
//     The method is applicable for frequency matrix only
// -------------------------------------------------------------------------

void DistributionMatrix::CheckForAllZeros()
{
    const double    zval = 0.001;
    int             r;

    const int       efective_number = NUMAA;    // effective number of residues
    const double    eqpart = rint(( double )FREQUENCY_SUM / efective_number ) / FREQUENCY_SUM;

    static int      symB = HashAlphSymbol('B');
    static int      symZ = HashAlphSymbol('Z');
    static int      resN = HashAlphSymbol('N');
    static int      resD = HashAlphSymbol('D');
    static int      resQ = HashAlphSymbol('Q');
    static int      resE = HashAlphSymbol('E');

    static double   Bprob = ( LOSCORES.PROBABility( resN ) + LOSCORES.PROBABility( resD ));
    static double   Zprob = ( LOSCORES.PROBABility( resQ ) + LOSCORES.PROBABility( resE ));

    for( int n = 0; n < GetColumns(); n++ ) {
        for( r = 0; r < NUMALPH; r++ )
            if( zval < ( *this )( n, r ))
                break;
        if( r == NUMALPH ) {
            // If all values are zero
//             for( r = 0; r < NUMALPH; r++ )
                //set weighted frequencies to the background frequencies
//                 ( *this )( n, r ) = LOSCORES.PROBABility( r );
            r = GetResidueAt( n );
            if( r == X ) {
                //X is at this position; set all amino acids equally probable
                for( int e = 0; e < efective_number; e++ )
                    ( *this )( n, e ) = LOSCORES.PROBABility( e );//eqpart;
            } else
            if( r == symB ) {
                //B is at the position; make appropriate amino acids available and
                //round a floating point number to precision of 2
                ( *this )( n, resN ) = rint( LOSCORES.PROBABility( resN ) * FREQUENCY_SUM / Bprob ) / FREQUENCY_SUM;
                ( *this )( n, resD ) = rint( LOSCORES.PROBABility( resD ) * FREQUENCY_SUM / Bprob ) / FREQUENCY_SUM;
            } else
            if( r == symZ ) {
                //Z is at the position; make appropriate amino acids available and
                //round a floating point number to precision of 2
                ( *this )( n, resQ ) = rint( LOSCORES.PROBABility( resQ ) * FREQUENCY_SUM / Zprob ) / FREQUENCY_SUM;
                ( *this )( n, resE ) = rint( LOSCORES.PROBABility( resE ) * FREQUENCY_SUM / Zprob ) / FREQUENCY_SUM;
            } else
                //set corresponding amino acid fully conserved else
                ( *this )( n, r ) = 1.0; //set 1 so that scores for one profile repeat the LOSCORES scores
        }
    }
}

// -------------------------------------------------------------------------
// Serialize: write the class data to file for reading them later
// -------------------------------------------------------------------------

void DistributionMatrix::Serialize( Serializer& serializer ) const
{
    serializer.Write(( char* )&columns, sizeof( columns ), 1 );

    for( int n = 0; n < columns; n++ )
        serializer.Write(( char* )values[n], sizeof( double ), NUMALPH );

    if( columns > 0 )
        serializer.Write( aacids, sizeof( char ), columns );
}

// -------------------------------------------------------------------------
// Deserialize: read data into the class members
// -------------------------------------------------------------------------

void DistributionMatrix::Deserialize( Serializer& serializer )
{
    //there's no need for destroying variables since memory allocated previously is reused
//     destroy();

    serializer.Read(( char* )&columns, sizeof( columns ), 1 );

    if( columns > MAXCOLUMNS )
        throw myruntime_error(
            mystring( "DistributionMatrix: Number of positions read from file is larger than the maximum allowed." ));
        
    if( columns <= 0 )
        throw myruntime_error(
            mystring( "DistributionMatrix: Invalid number of positions read from file." ));

    Reserve( columns ); // memory allocation

    for( int n = 0; n < columns; n++ )
        serializer.Read(( char* )values[n], sizeof( double ), NUMALPH );

    serializer.Read( aacids, sizeof( char ), columns );
    //
    CheckIntegrity();
}





////////////////////////////////////////////////////////////////////////////
// CLASS ExtendedDistributionMatrix
//
// Constructor
ExtendedDistributionMatrix::ExtendedDistributionMatrix()
{
    init();
}

// Destructor

ExtendedDistributionMatrix::~ExtendedDistributionMatrix()
{
    destroy();
}

// -------------------------------------------------------------------------
// init: initialization of members
// -------------------------------------------------------------------------

void ExtendedDistributionMatrix::init()
{
    DistributionMatrix::init();
    //
    trgfrqs_ = NULL;
    freqweights = NULL;
    information = NULL;
    expMIDs_ = NULL;
    //
    bppprob_ = NULL;
    ppprobs_ = NULL;
    pppndxs_ = NULL;
    noppps_ = NULL;
    //
    ctpsset_ = false;
    ctbppprob_ = NULL;
    ctppprobs_ = NULL;
    ctpppndxs_ = NULL;
    noctppps_ = NULL;
    //
    ctxvecset_ = false;
    lppfact_ = 0.0;
    ctxvecnorm2_ = NULL;
    ctxveclpprb_ = NULL;
    ctxvecs_ = NULL;
    ctxvecsize_ = 0;
    //
    sssset_ = false;
    sssp3_ = false;
    sss_ = NULL;
    sssprob_ = NULL;
    //
    name = NULL;
    description = NULL;

    szname = 0;
    szdescription = 0;

    nosequences = 0;
    effnosequences = 0.0;

    memset( backprobs_, 0, sizeof( double ) * NUMALPH );
    memset( postprobs_, 0, sizeof( double ) * NUMALPH );

    referenceLambda = -1.0;
    referenceK = -1.0;

    lambda = -1.0;
    entropy = -1.0;
    parameterK = -1.0;
    expscore = -1.0;

//     CalcLogPProbFactor();
}

// -------------------------------------------------------------------------
// destroy: deallocate memory and reset values
// -------------------------------------------------------------------------

void ExtendedDistributionMatrix::destroy()
{
    int n;
    if( trgfrqs_ )      { free( trgfrqs_ ); trgfrqs_ = NULL; }
    if( freqweights )   { free( freqweights ); freqweights = NULL; }
    if( information )   { free( information ); information = NULL; }
    if( expMIDs_ )      { free( expMIDs_ ); expMIDs_ = NULL; }
    if( bppprob_ )      { free( bppprob_ ); bppprob_ = NULL; }
    if( noppps_ )       { free( noppps_ ); noppps_ = NULL; }
    if( ctbppprob_ )    { free( ctbppprob_ ); ctbppprob_ = NULL; }
    if( noctppps_ )     { free( noctppps_ ); noctppps_ = NULL; }
    if( ctxvecnorm2_ )  { free( ctxvecnorm2_ ); ctxvecnorm2_ = NULL; }
    if( ctxveclpprb_ )  { free( ctxveclpprb_ ); ctxveclpprb_ = NULL; }
    if( sss_ )          { free( sss_ ); sss_ = NULL; }
    if( sssprob_ )      { free( sssprob_ ); sssprob_ = NULL; }
    if( name )          { free( name ); name = NULL; szname = 0; }
    if( description )   { free( description ); description = NULL; szdescription = 0; }
    //
    if( ppprobs_ ) {
        for( n = 0; n < columns; n++ )
            if( ppprobs_[n])
                free( ppprobs_[n]);
        free( ppprobs_ );
        ppprobs_ = NULL;
    }
    if( pppndxs_ ) {
        for( n = 0; n < columns; n++ )
            if( pppndxs_[n])
                free( pppndxs_[n]);
        free( pppndxs_ );
        pppndxs_ = NULL;
    }
    //
    if( ctppprobs_ ) {
        for( n = 0; n < columns; n++ )
            if( ctppprobs_[n])
                free( ctppprobs_[n]);
        free( ctppprobs_ );
        ctppprobs_ = NULL;
    }
    if( ctpppndxs_ ) {
        for( n = 0; n < columns; n++ )
            if( ctpppndxs_[n])
                free( ctpppndxs_[n]);
        free( ctpppndxs_ );
        ctpppndxs_ = NULL;
    }
    //
    if( ctxvecs_ ) {
        for( n = 0; n < columns; n++ )
            if( ctxvecs_[n])
                free( ctxvecs_[n]);
        free( ctxvecs_ );
        ctxvecs_ = NULL;
    }
    //
    DistributionMatrix::destroy();
}

// -------------------------------------------------------------------------
// Clear: erases all information contained in this class but leaves the
//     space allocated
// -------------------------------------------------------------------------

void ExtendedDistributionMatrix::Clear()
{
    int n;
    if( allocated ) {
        for( n = 0; n < columns; n++ ) {
            if( ppprobs_ && ppprobs_[n]) free( ppprobs_[n]);
            if( pppndxs_ && pppndxs_[n]) free( pppndxs_[n]);
            if( ctppprobs_ && ctppprobs_[n]) free( ctppprobs_[n]);
            if( ctpppndxs_ && ctpppndxs_[n]) free( ctpppndxs_[n]);
            if( ctxvecs_ && ctxvecs_[n]) free( ctxvecs_[n]);
        }
        memset( trgfrqs_, 0, sizeof( double ) * NUMAA * allocated );
        memset( freqweights, 0, sizeof( double ) * allocated );
        memset( information, 0, sizeof( double ) * allocated );
        memset( expMIDs_, 0, sizeof( double ) * PS_NSTATES * ( allocated + 1 ));
        memset( bppprob_, 0, sizeof( double ) * allocated );
        memset( ppprobs_, 0, sizeof( double*) * allocated );
        memset( pppndxs_, 0, sizeof( int*) * allocated );
        memset( noppps_,  0, sizeof( size_t ) * allocated );
        memset( ctbppprob_, 0, sizeof( double ) * allocated );
        memset( ctppprobs_, 0, sizeof( double*) * allocated );
        memset( ctpppndxs_, 0, sizeof( int*) * allocated );
        memset( noctppps_,  0, sizeof( size_t ) * allocated );
        memset( ctxvecnorm2_,  0, sizeof( double ) * allocated );
        memset( ctxveclpprb_,  0, sizeof( double ) * allocated );
        memset( ctxvecs_,  0, sizeof( double* ) * allocated );
        memset( sss_, 0, sizeof( char ) * allocated );
        memset( sssprob_, 0, sizeof( double ) * SS_NSTATES * allocated );
    }

    if( name )          { free( name ); name = NULL; szname = 0; }
    if( description )   { free( description ); description = NULL; szdescription = 0; }

    nosequences = 0;
    effnosequences = 0.0;

    //do not erase probabilities!
//     memset( backprobs_, 0, sizeof( double ) * NUMALPH );
//     memset( postprobs_, 0, sizeof( double ) * NUMALPH );

    referenceLambda = -1.0;
    referenceK = -1.0;

    lambda = -1.0;
    entropy = -1.0;
    parameterK = -1.0;
    expscore = -1.0;

    DistributionMatrix::Clear();
}

// -------------------------------------------------------------------------
// reallocate: allocate or reallocate necessary memory for the class members
// -------------------------------------------------------------------------

void ExtendedDistributionMatrix::reallocate( int howmuch )
{
    if( howmuch <= allocated )
        return;

    if( allocated <= 0 ) {
        trgfrqs_ = ( double(*)[NUMAA] )malloc( sizeof( double ) * NUMAA * howmuch );
        freqweights = ( double* )malloc( sizeof( double ) * howmuch );
        information = ( double* )malloc( sizeof( double ) * howmuch );
        expMIDs_ = ( double(*)[PS_NSTATES])malloc( sizeof(double)*PS_NSTATES*(howmuch+1));
        bppprob_ = ( double* )malloc( sizeof( double ) * howmuch );
        ppprobs_ = ( double** )malloc( sizeof( double*) * howmuch );
        pppndxs_ = ( int** )malloc( sizeof( int*) * howmuch );
        noppps_ = ( size_t* )malloc( sizeof( size_t ) * howmuch );
        ctbppprob_ = ( double* )malloc( sizeof( double ) * howmuch );
        ctppprobs_ = ( double** )malloc( sizeof( double*) * howmuch );
        ctpppndxs_ = ( int** )malloc( sizeof( int*) * howmuch );
        noctppps_ = ( size_t* )malloc( sizeof( size_t ) * howmuch );
        ctxvecnorm2_ = ( double* )malloc( sizeof( double ) * howmuch );
        ctxveclpprb_ = ( double* )malloc( sizeof( double ) * howmuch );
        ctxvecs_ = ( double** )malloc( sizeof( double*) * howmuch );
        sss_ = ( char* )malloc( sizeof( char ) * howmuch );
        sssprob_ = ( double(*)[SS_NSTATES])malloc( sizeof(double)*SS_NSTATES*howmuch );
    } else {
        trgfrqs_ = ( double(*)[NUMAA] )realloc( trgfrqs_, sizeof( double ) * NUMAA * howmuch );
        freqweights = ( double* )realloc( freqweights, sizeof( double ) * howmuch );
        information = ( double* )realloc( information, sizeof( double ) * howmuch );
        expMIDs_ = ( double(*)[PS_NSTATES])realloc( expMIDs_, sizeof(double)*PS_NSTATES*(howmuch+1));
        bppprob_ = ( double* )realloc( bppprob_, sizeof( double ) * howmuch );
        ppprobs_ = ( double** )realloc( ppprobs_, sizeof( double*) * howmuch );
        pppndxs_ = ( int** )realloc( pppndxs_, sizeof( int*) * howmuch );
        noppps_ = ( size_t* )realloc( noppps_, sizeof( size_t ) * howmuch );
        ctbppprob_ = ( double* )realloc( ctbppprob_, sizeof( double ) * howmuch );
        ctppprobs_ = ( double** )realloc( ctppprobs_, sizeof( double*) * howmuch );
        ctpppndxs_ = ( int** )realloc( ctpppndxs_, sizeof( int*) * howmuch );
        noctppps_ = ( size_t* )realloc( noctppps_, sizeof( size_t ) * howmuch );
        ctxvecnorm2_ = ( double* )realloc( ctxvecnorm2_, sizeof( double ) * howmuch );
        ctxveclpprb_ = ( double* )realloc( ctxveclpprb_, sizeof( double ) * howmuch );
        ctxvecs_ = ( double** )realloc( ctxvecs_, sizeof( double*) * howmuch );
        sss_ = ( char* )realloc( sss_, sizeof( char ) * howmuch );
        sssprob_ = ( double(*)[SS_NSTATES])realloc( sssprob_, sizeof(double)*SS_NSTATES*howmuch );
    }

    if( !trgfrqs_ || !freqweights || !information || !expMIDs_ )
        throw myruntime_error("ExtendedDistributionMatrix: Not enough memory.");
    if( !bppprob_ || !ppprobs_ || !pppndxs_ || !noppps_ )
        throw myruntime_error("ExtendedDistributionMatrix: Not enough memory.");
    if( !ctbppprob_ || !ctppprobs_ || !ctpppndxs_ || !noctppps_ )
        throw myruntime_error("ExtendedDistributionMatrix: Not enough memory.");
    if( !ctxvecnorm2_ || !ctxveclpprb_ || !ctxvecs_ )
        throw myruntime_error("ExtendedDistributionMatrix: Not enough memory.");
    if( !sss_ || !sssprob_ )
        throw myruntime_error("ExtendedDistributionMatrix: Not enough memory.");

    double    (*tmp_trgfrqs)[NUMAA] = trgfrqs_;
    double*     tmp_freqweights = freqweights;
    double*     tmp_information = information;
    double    (*tmp_expMIDs)[PS_NSTATES] = expMIDs_ + 1;
    double*     tmp_bppprob = bppprob_;
    double**    tmp_ppprobs = ppprobs_;
    int**       tmp_pppndxs = pppndxs_;
    size_t*     tmp_noppps = noppps_;
    double*     tmp_ctbppprob = ctbppprob_;
    double**    tmp_ctppprobs = ctppprobs_;
    int**       tmp_ctpppndxs = ctpppndxs_;
    size_t*     tmp_noctppps = noctppps_;
    double*     tmp_ctxvecnorm2 = ctxvecnorm2_;
    double*     tmp_ctxveclpprb = ctxveclpprb_;
    double**    tmp_ctxvecs = ctxvecs_;
    char*       tmp_sss = sss_;
    double    (*tmp_sssprob)[SS_NSTATES] = sssprob_;

    if( 0 < allocated ) {
        tmp_trgfrqs += allocated;
        tmp_freqweights += allocated;
        tmp_information += allocated;
        tmp_expMIDs += allocated;
        tmp_bppprob += allocated;
        tmp_ppprobs += allocated;
        tmp_pppndxs += allocated;
        tmp_noppps += allocated;
        tmp_ctbppprob += allocated;
        tmp_ctppprobs += allocated;
        tmp_ctpppndxs += allocated;
        tmp_noctppps += allocated;
        tmp_ctxvecnorm2 += allocated;
        tmp_ctxveclpprb += allocated;
        tmp_ctxvecs += allocated;
        tmp_sss += allocated;
        tmp_sssprob += allocated;
    } else {
        memset( expMIDs_, 0, sizeof( double ) * PS_NSTATES );
    }

    // fill uninitialized memory with zeros
    memset( tmp_trgfrqs, 0, sizeof( double ) * NUMAA * ( howmuch - allocated ));
    memset( tmp_freqweights, 0, sizeof( double ) * ( howmuch - allocated ));
    memset( tmp_information, 0, sizeof( double ) * ( howmuch - allocated ));
    memset( tmp_expMIDs, 0, sizeof( double ) * PS_NSTATES * ( howmuch - allocated ));
    memset( tmp_bppprob, 0, sizeof( double ) * ( howmuch - allocated ));
    memset( tmp_ppprobs, 0, sizeof( double*) * ( howmuch - allocated ));
    memset( tmp_pppndxs, 0, sizeof( int*) * ( howmuch - allocated ));
    memset( tmp_noppps, 0, sizeof( size_t ) * ( howmuch - allocated ));
    memset( tmp_ctbppprob, 0, sizeof( double ) * ( howmuch - allocated ));
    memset( tmp_ctppprobs, 0, sizeof( double*) * ( howmuch - allocated ));
    memset( tmp_ctpppndxs, 0, sizeof( int*) * ( howmuch - allocated ));
    memset( tmp_noctppps, 0, sizeof( size_t ) * ( howmuch - allocated ));
    memset( tmp_ctxvecnorm2, 0, sizeof( double ) * ( howmuch - allocated ));
    memset( tmp_ctxveclpprb, 0, sizeof( double ) * ( howmuch - allocated ));
    memset( tmp_ctxvecs, 0, sizeof( double*) * ( howmuch - allocated ));
    memset( tmp_sss, 0, sizeof( char ) * ( howmuch - allocated ));
    memset( tmp_sssprob, 0, sizeof( double ) * SS_NSTATES * ( howmuch - allocated ));
    //
    DistributionMatrix::reallocate( howmuch );
}

// -------------------------------------------------------------------------
// Push: save a vector of values;
//     weight -- weight of frequencies
//     info -- information content
//
void ExtendedDistributionMatrix::Push( const double posvalues[NUMALPH], char aa, 
                                       double weight, double info, double expnobs[PS_NSTATES] )
{
    if( allocated < GetColumns() + 1 ) {
        reallocate( TIMES2( allocated +1 ));
    }

    PushAt( posvalues, aa, weight, info, expnobs, GetColumns());
}

// -------------------------------------------------------------------------
// PushAt: save a vector of values at the given position;
//     weight -- weight of frequencies
//     info -- information content
//
void ExtendedDistributionMatrix::PushAt( const double posvalues[NUMALPH], char aa, 
                                         double weight, double info, double expnobs[PS_NSTATES], 
                                         int position )
{
    if( allocated < position + 1 )
        throw myruntime_error("ExtendedDistributionMatrix: Failed to insert vector of values.");

    DistributionMatrix::PushAt( posvalues, aa, position );
    freqweights[position] = weight;
    information[position] = info;
    SetMIDExpNoObservationsAt( position, expnobs );
    //
}

// -------------------------------------------------------------------------
// PushAt: pushes vector of values at the given position
//
void ExtendedDistributionMatrix::PushAt( const double posvalues[NUMALPH], char aa, int position )
{
    if( allocated < position + 1 )
        throw myruntime_error("ExtendedDistributionMatrix: Memory access error.");
    DistributionMatrix::PushAt( posvalues, aa, position );
}

// -------------------------------------------------------------------------
// Finalize: finalize data; calculate target frequencies from scores
//
void ExtendedDistributionMatrix::Finalize()
{
    CalcTrgFrequencies();
}





// -------------------------------------------------------------------------
// CalcLogPProbFactor: calculate log of prior probability factor used in 
//  processing context vectors
//
void ExtendedDistributionMatrix::CalcLogPProbFactor()
{
    double  lpfact = 0.0;
    double  term1, term2;
    double  operr;
    int     err;
    const double    do2 = (double)CVS.DIM * 0.5;
    const double    nu0p1o2 = (CVS.NU0+1.) * 0.5;
    const double    kp0 = CVS.KAPPA0;
    const double    kp0p1 = kp0+1.;

    const double    nu_t_o2 = nu0p1o2;//nu over 2

    const double    dofdom_o2 = nu_t_o2 - do2;
    const double    dofdom_pdim_o2 = dofdom_o2 + do2;//(deg.of freedom + dim.) over 2

    if(( err = psl_lngamma_e( dofdom_pdim_o2, &term1, &operr )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));
    if(( err = psl_lngamma_e( dofdom_o2, &term2, &operr )) != 0 )
        throw myruntime_error( TranslatePSLError( err ));

    lpfact = term1 - term2;
//     lpfact -= do2 * SLC_LNPI;//NOTE:do not include pi^{-1/2A}
    lpfact += do2 * (log(kp0)-log(kp0p1));

    lppfact_ = lpfact;
}





// -------------------------------------------------------------------------
// CalcTgtFrequencies: calculate target frequencies from scores
//
void ExtendedDistributionMatrix::CalcTrgFrequencies()
{
    int     effnoress = NUMAA;
    double  (*trgfs )[NUMAA] = GetTrgFreqs();
    char    res;
    int     m, r;

    static const double lsclambda = LOSCORES.StatisParam( Ungapped, Lambda );

    if( trgfs == NULL )
        throw myruntime_error("ExtendedDistributionMatrix: CalcTrgFrequencies: Null target frequencies.");

    for( m = 0; m < GetColumns(); m++, trgfs++ ) {
        res = GetResidueAt( m );
        // omit unused positions and gaps in query
        if( res == GAP )
            continue;
        for( r = 0; r < effnoress; r++ )
            (*trgfs)[r] = exp( GetValueAt(m,r)*lsclambda ) * LOSCORES.PROBABility(r);
    }
}

// -------------------------------------------------------------------------
// CalcScores: calculate scores from target frequencies
//
void ExtendedDistributionMatrix::CalcScores()
{
    int     noress = NUMALPH;
    int     effnoress = NUMAA;
    double  (*trgfs)[NUMAA] = GetTrgFreqs();
    double  (*scos)[NUMALPH] = GetVector();
    double  bprob, ratv;
    char    res;
    int     m, r;

    static const double lsclambda = LOSCORES.StatisParam( Ungapped, Lambda );

    if( scos == NULL )
        throw myruntime_error("ExtendedDistributionMatrix: CalcScores: Null scores.");
    if( trgfs == NULL )
        throw myruntime_error("ExtendedDistributionMatrix: CalcScores: Null target frequencies.");
    if( lsclambda <= 0.0 )
        throw myruntime_error("ExtendedDistributionMatrix: CalcScores: Invalid value of lambda.");

    for( m = 0; m < GetColumns(); m++, scos++, trgfs++ ) {
        res = GetResidueAt( m );
        // omit unused positions and gaps in query
        if( res == GAP )
            continue;
        for( r = 0; r < effnoress; r++ ) {
            bprob = LOSCORES.PROBABility(r);
            if((*trgfs)[r] <= 0.0 || 1.0 <= (*trgfs)[r])
                throw myruntime_error("ExtendedDistributionMatrix: CalcScores: Invalid target frequencies.");
            if( bprob <= 0.0 || 1.0 <= bprob )
                throw myruntime_error("ExtendedDistributionMatrix: CalcScores: Invalid back. probabilities.");
            (*scos)[r] = log((*trgfs)[r] / bprob ) / lsclambda;
        }
        for( ; r < noress; r++ ) {
            ratv = LOSCORES.FreqRatio( res, r );
            if( ratv <= 0.0 )
                (*scos)[r] = SCORE_MIN;
            else
                (*scos)[r] = LOSCORES.PrecomputedEntry( res, r );
        }
    }
}

// =========================================================================
// MixTrgFrequencies: mix target frequencies under HDP framework
//
void ExtendedDistributionMatrix::MixTrgFrequencies( const HDPbase* hdpbase )
{
    mystring errstr;
    mystring preamb = "MixTrgFrequencies: ";

    const int       noress = NUMALPH;
    const int       noeffress = NUMAA;
    double          pme = 0.0;//estimated probability
    double          tfr = 0.0;//target frequency
    double*         infvv = GetInformation();
    double  (*trgfs)[NUMAA] = GetTrgFreqs();

    if( infvv == NULL )
        throw myruntime_error( preamb + "Null information values.");
    if( trgfs == NULL )
        throw myruntime_error( preamb + "Null target frequencies.");
    if( hdpbase == NULL )
        throw myruntime_error( preamb + "Null HDP structure.");

    int     ctxtsz = hdpbase->GetCtxtSize();
    int     prolen = GetColumns();
    int     p;
    int     left, right, r;
    int     parity = ( ctxtsz & 1 ) ^ 1;
    int     hlf = ctxtsz >> 1;
    int     mid = hlf - parity;

    if( ctxtsz < 1 )
        throw myruntime_error( preamb + "Wrong HDP context size.");
    if( prolen < ctxtsz ) {
        warning("Profile length is less than context size. Not mixing.");
        return;
    }
    Pslmatrix   promtx( prolen, noeffress );//profile matrix
    Pslmatrix   ctxmtx;//context matrix
    Pslvector   ctxsck;//context stack
    Pslvector   ctxnrm( noeffress*ctxtsz );//normal transform
    Pslvector   mixed( noeffress );//mixed vector
    double      infrm;

    //make matrix representation of profile
    for( p = 0; p < prolen; p++, trgfs++ ) {
        if( !IsValidResSym( GetResidueAt(p)))
            throw myruntime_error( preamb + "Not a residue at profile position.");

        for( r = 0; r < noeffress; r++ ) {
            if((*trgfs)[r] <= 0.0 || 1.0 <= (*trgfs)[r])
                throw myruntime_error( preamb + "Invalid target frequencies.");
            promtx.SetValueAt( p, r, (*trgfs)[r]);
        }
    }

    //iterate over all profile positions
    trgfs = GetTrgFreqs();
    for( p = 0; p < prolen; p++, trgfs++, infvv++ )
    {
        right = ( int )SLC_MIN( prolen-1, p+hlf );
        left = SLC_MAX( 0, right-ctxtsz+1 );

        ctxmtx = promtx.SubMatrix( left, 0, ctxtsz, noeffress );
        ctxsck = ctxmtx.Stack();
        ctxnrm.Copy( ctxsck );

        //mix central vector of context;
        //input vector ctxnrm will be transformed;
        //mixed will be returned in logit-normal space
        hdpbase->MixCLNVar( ctxnrm, &mixed );

        //write PME back to target frequencies
        infrm = 0.0;
        for( r = 0; r < noeffress; r++ ) {
            pme = mixed.GetValueAt(r);
            (*trgfs)[r] = pme;
            //calculate relative entropy
            if( LOSCORES.PROBABility(r) <= 0.0 )
                continue;
            if( 0.0 < pme )
                infrm += pme * log( pme / LOSCORES.PROBABility(r));
        }
        //save relative entropy
        infrm /= LN2;
        infrm = SLC_MAX( 0.0, infrm );
        *infvv = infrm;
    }

    //recalculate scores
    CalcScores();
}

// =========================================================================
// CalcTFPostPredictives: calculate posterior predictives probabilities of 
//  target frequencies under HDP framework
//
void ExtendedDistributionMatrix::CalcTFPostPredictives( const HDPbase* hdpbase )
{
    mystring errstr;
    mystring preamb = "CalcTFPostPredictives: ";

    const int       noress = NUMALPH;
    const int       noeffress = NUMAA;
    int             ndx = 0;//index of cluster
    double          ppr = 0.0;//posterior probability
    double          tfr = 0.0;//target frequency
    double  (*trgfs)[NUMAA] = GetTrgFreqs();

    if( trgfs == NULL )
        throw myruntime_error( preamb + "Null target frequencies.");
    if( hdpbase == NULL )
        throw myruntime_error( preamb + "Null HDP structure.");

    int     nosupcls = hdpbase->GetNoSupClusters();
    int     ctxtsz = hdpbase->GetCtxtSize();
    int     prolen = GetColumns();
    int     p;
    int     left, right, r, n;
    int     parity = ( ctxtsz & 1 ) ^ 1;
    int     hlf = ctxtsz >> 1;
    int     mid = hlf - parity;

    if( nosupcls < 1 && nosupcls != -1 )
        throw myruntime_error( preamb + "Wrong number of HDP support clusters.");
    if( ctxtsz < 1 )
        throw myruntime_error( preamb + "Wrong HDP context size.");
    if( prolen < ctxtsz ) {
        warning("Profile length is less than context size. Not using posteriors.");
        return;
    }
    Pslmatrix   promtx( prolen, noeffress );//profile matrix
    Pslmatrix   ctxmtx;//context matrix
    Pslvector   ctxsck;//context stack
    Pslvector   ctxnrm( noeffress*ctxtsz );//normal transform
    Pslvector   ppprobs;//vector of posterior probabilities
    Ivector     cindcs;//indices of clusters

    //make matrix representation of profile
    for( p = 0; p < prolen; p++, trgfs++ ) {
        if( !IsValidResSym( GetResidueAt(p)))
            throw myruntime_error( preamb + "Not a residue at profile position.");

        for( r = 0; r < noeffress; r++ ) {
            if((*trgfs)[r] <= 0.0 || 1.0 <= (*trgfs)[r])
                throw myruntime_error( preamb + "Invalid target frequencies.");
            promtx.SetValueAt( p, r, (*trgfs)[r]);
        }
    }

    //iterate over all profile positions
    for( p = 0; p < prolen; p++ )
    {
        right = ( int )SLC_MIN( prolen-1, p+hlf );
        left = SLC_MAX( 0, right-ctxtsz+1 );

        ctxmtx = promtx.SubMatrix( left, 0, ctxtsz, noeffress );
        ctxsck = ctxmtx.Stack();
        ctxnrm.Copy( ctxsck );

        //calculate posterior predictive probabilities;
        //input vector ctxnrm will be transformed;
        hdpbase->CalcPPProbs( ctxnrm, &ppr, &ppprobs, &cindcs, 0.02/*lpfact*/);

        if( /*ppprobs.GetSize() < 1 || */ppprobs.GetSize() != cindcs.GetSize())
            throw myruntime_error( preamb + "Invalid number of posterior predictive probabilities.");

        //save posteriors
        SetBckPPProbAt( ppr, p );
        SetPPProbsAt( p, ppprobs.GetVector(), cindcs.GetVector(), ppprobs.GetSize());
    }
}

// -------------------------------------------------------------------------
// Calc_ctTFPostPredictives: calculate posterior predictives 
//  probabilities of target frequencies under HDP framework; environmental 
//  case
//
void ExtendedDistributionMatrix::Calc_ctTFPostPredictives( const HDPbase* hdpctbase )
{
    mystring errstr;
    mystring preamb = "Calc_ctTFPostPredictives: ";

    const int       noress = NUMALPH;
    const int       noeffress = NUMAA;
    int             ndx = 0;//index of cluster
    double          ppr = 0.0;//posterior probability
    double          tfr = 0.0;//target frequency
    double          cfc;
    double  (*trgfs)[NUMAA] = GetTrgFreqs();
    const size_t    cnCTXLEN = 21;//TODO: read ctx length from parameter file
    const double    cdCTCWGT = 0.5;//TODO: read ctx central weight from parameter file
    CtxtCoefficients coeffs( cnCTXLEN, cdCTCWGT );
    coeffs.FindCoefficients();

    if( trgfs == NULL )
        throw myruntime_error( preamb + "Null target frequencies.");
    if( hdpctbase == NULL )
        throw myruntime_error( preamb + "Null HDP structure.");

    int     nosupcls = hdpctbase->GetNoSupClusters();
    int     ctxtsz = hdpctbase->GetCtxtSize();
    int     prolen = GetColumns();
    int     p, pn, c;
    int     left, right, r, n;
    int     parity = ( cnCTXLEN & 1 ) ^ 1;
    int     hlf = cnCTXLEN >> 1;
    int     mid = hlf - parity;

    if( nosupcls < 1 && nosupcls != -1 )
        throw myruntime_error( preamb + "Wrong number of HDP support clusters.");
    if( ctxtsz < 1 || 1 < ctxtsz )
        throw myruntime_error( preamb + "Wrong HDP context size.");
    if( prolen < ctxtsz || prolen < (cnCTXLEN-1)>>1 + 1 ) {
        warning("Profile length is less than context size. Not using ct. posteriors.");
        return;
    }
    Pslmatrix   promtx( prolen, noeffress-1 );//transformed profile matrix
    Pslvector   tfrqs( noeffress );//target frequencies
    Pslvector   ctxnrm( noeffress-1 );//positional representative
    Pslvector   ppprobs;//vector of posterior probabilities
    Ivector     cindcs;//indices of clusters

    //make matrix representation of profile
    for( p = 0; p < prolen; p++, trgfs++ ) {
        if( !IsValidResSym( GetResidueAt(p)))
            throw myruntime_error( preamb + "Not a residue at profile position.");

        for( r = 0; r < noeffress; r++ ) {
            if((*trgfs)[r] <= 0.0 || 1.0 <= (*trgfs)[r])
                throw myruntime_error( preamb + "Invalid target frequencies.");
            tfrqs.SetValueAt( r, (*trgfs)[r]);
        }
        //transform to normal
        ::LogitNormal2Normal( tfrqs.GetVector(), noeffress, 1.e-1, false );
        for( r = 0; r < noeffress-1; r++ ) {
            promtx.SetValueAt( p, r, tfrqs.GetValueAt(r));
        }
    }

    SetctPsSet( true );
    //iterate over all profile positions
    for( p = 0; p < prolen; p++ )
    {
        right = p + hlf;//( int )SLC_MIN( prolen-1, p+hlf );
        left = p - hlf + parity;//SLC_MAX( 0, right-ctxtsz+1 );

        ctxnrm.Zero();
        for( r = left, c = 0; r <= right && c < coeffs.GetLength(); r++, c++ ) {
            pn = r;
            if( prolen <= pn || pn < 0 )
                pn = p + p - pn;
            cfc = coeffs.GetCoefficientAt(c);
            ctxnrm.Superposition( cfc, promtx.RowVector(pn));
        }

        //calculate posterior predictive probabilities;
        hdpctbase->CalcPPProbs( ctxnrm, &ppr, &ppprobs, &cindcs,
                                0.02/*lpfact*/, true/*usepriors*/, false/*tonormal*/ );

        if( /*ppprobs.GetSize() < 1 || */ppprobs.GetSize() != cindcs.GetSize())
            throw myruntime_error( preamb + "Invalid number of posterior predictive probabilities.");

        //save posteriors
        SetctBckPPProbAt( ppr, p );
        SetctPPProbsAt( p, ppprobs.GetVector(), cindcs.GetVector(), ppprobs.GetSize());
    }
}





// =========================================================================
// CalcCtxVector: calculate context vector, its prior probability, and its 
//  squared norm for each position
//
void ExtendedDistributionMatrix::CalcCtxVector()
{
    mystring errstr;
    mystring preamb = "CalcCtxVector: ";

    const int       noress = NUMALPH;
    const int       noeffress = NUMAA;
    double          cfc;
    double  (*trgfs)[NUMAA] = GetTrgFreqs();
    CtxtCoefficients coeffs( CVS.CTXLEN, CVS.CWGT );
    coeffs.FindCoefficients();

    if( CVS.CTXLEN < 1 )
        throw myruntime_error( preamb + "Invalid context length.");
    if( CVS.CWGT < 0.0 || 1.0 < CVS.CWGT )
        throw myruntime_error( preamb + "Invalid context central weight.");

    CalcLogPProbFactor();

    if( trgfs == NULL )
        throw myruntime_error( preamb + "Null target frequencies.");

    int     prolen = GetColumns();
    int     p, pn, c, cc, d;
    int     left, right, r, n, err;
    int     parity = ( CVS.CTXLEN & 1 ) ^ 1;
    int     hlf = CVS.CTXLEN >> 1;
    int     mid = hlf - parity;

    if( prolen < CVS.CTXLEN || prolen < (CVS.CTXLEN-1)>>1 + 1 ) {
        warning("Profile length is less than context size; not calculating vectors.");
        return;
    }
    Pslmatrix   promtx( prolen, noeffress-1 );//transformed profile matrix
    Pslvector   tfrqs( noeffress );//target frequencies
    Pslvector   ctxnrm( CVS.DIM );//positional representative
    double      norm2;//squared norm of vector
    double      pprob;//prior probability of vector

    //make matrix representation of profile
    for( p = 0; p < prolen; p++, trgfs++ ) {
        if( !IsValidResSym( GetResidueAt(p)))
            throw myruntime_error( preamb + "Not a residue at profile position.");

        for( r = 0; r < noeffress; r++ ) {
            if((*trgfs)[r] <= 0.0 || 1.0 <= (*trgfs)[r])
                throw myruntime_error( preamb + "Invalid target frequencies.");
            tfrqs.SetValueAt( r, (*trgfs)[r]);
        }
        //transform to normal
        ::LogitNormal2Normal( tfrqs.GetVector(), noeffress, 1.e-1, false );
        for( r = 0; r < noeffress-1; r++ ) {
            promtx.SetValueAt( p, r, tfrqs.GetValueAt(r));
        }
    }

    SetCtxVecSet( true );
    //iterate over all profile positions
    for( p = 0; p < prolen; p++ )
    {
        right = p + hlf;
        left = p - hlf + parity;

        ctxnrm.Zero(); cc = 0;
        for( r = left, c = 0; r <= right && c < coeffs.GetLength(); r++, c++ ) {
            pn = r;
            if( prolen <= pn || pn < 0 )
                pn = p + p - pn;
            cfc = coeffs.GetCoefficientAt(c);
            if( CVS.MIX ) {
                if( CVS.AUG ) {
                    cc = 0;
                    if( mid == c ) cc = promtx.GetNoCols();
                    else if( mid < c ) cc = TIMES2( promtx.GetNoCols());
                    //
                    for( d = 0; d < promtx.GetNoCols(); d++ )
                        ctxnrm.AddValueAt( cc+d, cfc*(promtx.GetValueAt(pn,d)-CVS.MU0.GetValueAt(d)) );
                }
                else
                    if(( err = ctxnrm.Superposition( cfc, promtx.RowVector(pn))) != 0 )
                        throw myruntime_error( preamb + TranslatePSLError( err ));
            }
            else {
                for( d = 0; d < promtx.GetNoCols(); d++ )
                    ctxnrm.SetValueAt( cc++, cfc*(promtx.GetValueAt(pn,d)-CVS.MU0.GetValueAt(d)) );
            }
        }

        if( CVS.MIX && !CVS.AUG && CVS.MEAN ) {
            if(( err = ctxnrm.Superposition( -1., CVS.MU0 )) != 0 )
                throw myruntime_error( preamb + TranslatePSLError( err ));
        }

        //vector's squared norm
        norm2 = ctxnrm.Norm2();
        norm2 *= norm2;

        //vector's prior probability scaled by pi^{1/2dim}, where dim=noeffress-1
        pprob = GetLogPProbFactor();
        pprob -= 0.5*(CVS.NU0+1.) * log(1.+CVS.KAPPA0/(CVS.KAPPA0+1.)*norm2);

        //save vector and calculated magnitudes
        SetCtxVecPlusAt( p, norm2, pprob, ctxnrm.GetVector(), ctxnrm.GetSize());
    }
}





// -------------------------------------------------------------------------
// SetName: sets name of the multiple alignment the matrix was
//     constructed by
//
void ExtendedDistributionMatrix::SetName( const char* newname )
{
    if( name ) {
        free( name );
        name = NULL;
        SetNameSize( 0 );
    }

    if( !newname || !*newname )
        return;

    size_t  newlength = strlen( newname );
    if( !newlength )
        throw myruntime_error( mystring( "ExtendedDistributionMatrix: Invalid argument of filename." ));

    name = ( char* )malloc( newlength + 1 );
    if( !name )
        throw myruntime_error( mystring( "ExtendedDistributionMatrix: Not enough memory." ));

    strncpy( name, newname, newlength + 1 ); // include the terminating null symbol
    SetNameSize( newlength );
}

// -------------------------------------------------------------------------
// SetTitle: replaces text of the title to the new one
//
void ExtendedDistributionMatrix::SetDescription( const char* newdesc )
{
    if( description ) {
        free( description );
        description = NULL;
        SetDescriptionSize( 0 );
    }

    if( !newdesc )
        return;

    size_t  newlength = strlen( newdesc );
    if( !newlength )
        throw myruntime_error( mystring( "ExtendedDistributionMatrix: Wrong title argument." ));

    description = ( char* )malloc( newlength + 1 );
    if( !description )
        throw myruntime_error( mystring( "ExtendedDistributionMatrix: Not enough memory." ));

    strncpy( description, newdesc, newlength + 1 ); // include the terminating null symbol
    SetDescriptionSize( newlength );
}

// -------------------------------------------------------------------------
// OutputMatrix: output all information to file
//
void ExtendedDistributionMatrix::OutputMatrix( const char* filename ) const
{
    int     efective_number = NUMALPH - 1; //effective number of residues
    //
    FILE*   fp = stdout;
    size_t  l = 0;

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw myruntime_error(
            mystring( "ExtendedDistributionMatrix: Failed to open file to write matrix." ));

    fprintf( fp, "Multiple alignment (no. sequences %d, eff. no. sequences %.3f)", GetNoSequences(), GetEffNoSequences());

    if( GetName())
        fprintf( fp, ": %s", GetName());

    fprintf( fp, "\n" );

    if( GetDescription())
        fprintf( fp, "First sequence description: %s\n\n", GetDescription());

    fprintf( fp, "%28c Position-specific scoring matrix %22c Freq weights %c Information %c Exp.Obs.\n", 32, 32, 32, 32 );
    fprintf( fp, "%9c", 32 );

    for( char r = 0; r < efective_number; r++ )
        fprintf( fp, "%3c", DehashCode( r ) );

    for( int p = 0; p < GetColumns(); p++ ) {
        fprintf( fp, "\n%5d %c   ", ++l, DehashCode(( *this )[p] ));

        for( char r = 0; r < efective_number; r++ )
            fprintf( fp, "%3d", ( int )rint(( *this )( p, r )) );

        fprintf( fp, " %11d", ( int )rint( GetFrequencyWeightAt( p )));
        fprintf( fp, " %13.2f", GetInformationAt( p ));
        fprintf( fp, " %7.2f", GetMIDExpNoObservationsAt( p, PS_M ));
        fprintf( fp, " %7.2f", GetMIDExpNoObservationsAt( p, PS_I ));
        fprintf( fp, " %7.2f", GetMIDExpNoObservationsAt( p, PS_D ));
    }
    fprintf( fp, "\n" );

    if( fp != stdout )
        fclose( fp );
}

// -------------------------------------------------------------------------
// OutputProfile: output profile information in the text format
//
void OutputProfile( const char* filename,
        const FrequencyMatrix& freq,
        const LogOddsMatrix& odds,
        const GapScheme& gaps )
{
    FILE*           fp = stdout;
    size_t          l = 0;
    const size_t    res_count = NUMALPH - 1; //exclude gap symbol

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw myruntime_error(
            mystring( "Failed to open file for writing." ));

    fprintf( fp, "Multiple alignment (no. sequences %d, eff. no. sequences %.3f)",
            odds.GetNoSequences(), odds.GetEffNoSequences());

    if( odds.GetName())
        fprintf( fp, ": %s", odds.GetName());

    fprintf( fp, "\n" );

    if( odds.GetDescription())
        fprintf( fp, "First sequence description: %s\n\n", odds.GetDescription());


    fprintf( fp,"%28c Position-specific scoring matrix "
                "%53c Weighted observed frequencies %30c Gap weights %c Freq weights %c Information\n",
                32, 32, 32, 32, 32 );

    fprintf( fp, "%9c", 32 );

    for( unsigned char r = 0; r < res_count; r++ )
        fprintf( fp, "%3c", DehashCode( r ) );

    for( unsigned char r = 0; r < res_count; r++ )
        fprintf( fp, "%4c", DehashCode( r ) );

    for( int p = 0; p < odds.GetColumns(); p++ ) {
        // omit unused positions and gaps in query
        if( odds.GetResidueAt( p ) == GAP )
            continue;

        fprintf( fp, "\n%5d %c   ", ++l, DehashCode( odds.GetResidueAt( p )));

        for( unsigned char r = 0; r < res_count; r++ )
            fprintf( fp, "%2d ", ( int )rint( odds.GetValueAt( p, r )));

        for( unsigned char r = 0; r < res_count; r++ )
            fprintf( fp, "%4d", ( int )rint( 100 * freq.GetValueAt( p, r )));

        fprintf( fp, " %6d", ( int )rint( 100 * gaps.GetWeightsAt( p )));
        fprintf( fp, " %14d", ( int )rint( odds.GetFrequencyWeightAt( p )));

        fprintf( fp, " %13.2f", odds.GetInformationAt( p ));
    }

    fprintf( fp, "\n\n" );

// #ifdef SCALE_PROFILES
    odds.PrintParameters( fp );
// #endif

    if( fp != stdout )
        fclose( fp );
}

// -------------------------------------------------------------------------
// GetMaxAnnotationWidth: get maximum width of annotation 
//
size_t ExtendedDistributionMatrix::GetMaxAnnotationWidth() const
{
    static size_t   width = 2 * OUTPUTINDENT + OUTPUTWIDTH + 2;
    return width;
}

// -------------------------------------------------------------------------
// PrintAnnotation: print short annotation of the name and description to
//     string stream; space of stream must be pre-allocated!!
//
void ExtendedDistributionMatrix::PrintAnnotation( char* sp ) const
{
    if( !sp )
        return;
    *sp = 0;//to ensure printing to the end of the stream
    PrintAnnotation( &string_print, sp );
}

// -------------------------------------------------------------------------
// PrintAnnotation: print short annotation to file
//
void ExtendedDistributionMatrix::PrintAnnotation( FILE* fp ) const
{
    PrintAnnotation( &file_print, fp );
}

// -------------------------------------------------------------------------
// PrintAnnotation: format and print short annotation of the name and
//     description
//
void ExtendedDistributionMatrix::PrintAnnotation( TPrintFunction print_func, void* vpn ) const
{
    static size_t   preamble = 0;
    static size_t   textwidth = OUTPUTINDENT + OUTPUTWIDTH + 2;
    static size_t   width = textwidth + preamble;
    static size_t   max_length = width;
    static size_t   max_rows = 1;

    PrintDescriptionHelper( print_func, vpn, preamble, textwidth, width, max_rows, max_length, true );
}

// -------------------------------------------------------------------------
// GetMinimumRequiredSizeForDescription: get minimum size required to
//     contain full description of the profile
//
size_t ExtendedDistributionMatrix::GetMinimumRequiredSizeForDescription() const
{
    static size_t   indent = OUTPUTINDENT;              //indentation
    static size_t   textwidth = OUTPUTWIDTH;            //output width

#ifdef __DEBUG__
    if( !textwidth )
        return 0;
#endif

    size_t  name_size = GetNameSize();              //actual size of name
    size_t  desc_size = GetDescriptionSize();       //actual size of description
    size_t  max_desc_size = GetPrivateBufferSize(); //maximum size of description

    size_t  title_size  = ( name_size + desc_size ) * ( 1 + ( indent + 2 ) / textwidth ) + 1;

    if( title_size >= max_desc_size )
        title_size  = max_desc_size;

    return title_size;
}

// -------------------------------------------------------------------------
// PrintDescriptionFirst: format and print name and description to file
//
void ExtendedDistributionMatrix::PrintDescriptionFirst( FILE* fp ) const
{
    PrintDescriptionFirst( &file_print, fp );
}

// -------------------------------------------------------------------------
// PrintDescriptionFirst: format and print name and description of the
//     matrix
//
void ExtendedDistributionMatrix::PrintDescriptionFirst( TPrintFunction print_func, void* vpn ) const
{
    static size_t   preamble = 0;
    static size_t   textwidth = OUTPUTINDENT + OUTPUTWIDTH + 1;
    static size_t   width = textwidth + preamble;
    static size_t   max_length = GetPrivateBufferSize();
    static size_t   max_rows = max_length / textwidth + 2;

    PrintDescriptionHelper( print_func, vpn, preamble, textwidth, width, max_rows, max_length, true );
}

// -------------------------------------------------------------------------
// PrintDescription: format and print name and description to string stream;
//     space of stream must be pre-allocated!!
//
void ExtendedDistributionMatrix::PrintDescription( char* sp ) const
{
    if( !sp )
        return;
    *sp = 0;//to ensure printing to the end of the stream
    PrintDescription( &string_print, sp );
}

// -------------------------------------------------------------------------
// PrintDescription: format and print name and description to file
//
void ExtendedDistributionMatrix::PrintDescription( FILE* fp ) const
{
    PrintDescription( &file_print, fp );
}

// -------------------------------------------------------------------------
// PrintDescription: format and print name and description
//
void ExtendedDistributionMatrix::PrintDescription( TPrintFunction print_func, void* vpn ) const
{
    static size_t   preamble = OUTPUTINDENT;
    static size_t   textwidth = OUTPUTWIDTH + 1;
    static size_t   width = textwidth + preamble;
    static size_t   max_length = GetPrivateBufferSize();
    static size_t   max_rows = max_length / textwidth + 2;

    PrintDescriptionHelper( print_func, vpn, preamble, textwidth, width, max_rows, max_length, false );
}

// -------------------------------------------------------------------------
// PrintDescription: format and print method helper
//
void ExtendedDistributionMatrix::PrintDescriptionHelper(
        TPrintFunction print_func, void* vpn,
        size_t preamble, size_t textwidth, size_t width, size_t max_rows, size_t max_length,
        bool annotation ) const
{
    if( vpn == NULL )
        return;

    bool            printname = false;

    char*           buffer = GetPrivateBuffer();
    const char*     name = printname? GetName(): NULL;
    const char*     description = GetDescription();
    const char*     ending = "...";

    size_t  sz_name = printname? GetNameSize(): 0;
    size_t  sz_description = GetDescriptionSize();
    size_t  szending = strlen( ending );
    size_t  size = 0;
    size_t  pos = 0;    //current position in line
    size_t  tot = 0;    //position in the buffer
    bool    term = false;

    size  = sz_name + sz_description + ( name ? 3: 0 );     //+parentheses and space
    size += ( size / textwidth + 1 ) * ( preamble + 1 );    //#lines times preamble

    if( buffer == NULL || !max_length || GetPrivateBufferSize() < max_length )
        return;

    if( size >= max_length - max_rows * ( preamble + 1 )) {
        size  = max_length - max_rows * ( preamble + 1 ) - szending - 1;
        term = true;
    }

    if( !annotation ) {
        *buffer++ = '>'; pos++; tot++;
    }
    if( name ) {
        *buffer++ = '('; pos++; tot++;
        FormatBuffer( buffer, name, tot, pos, size, preamble, width );
        if( tot < size ) { *buffer++ = ')'; pos++; tot++; }
        if( tot < size ) { *buffer++ = 32;  pos++; tot++; }
    }
    if( description ) {
        FormatBuffer( buffer, description, tot, pos, size, preamble, width );
    }

    if( term ) {
        if( width <= pos + szending ) {
            *buffer++ = '\n';
            for( size_t k = 0; k < preamble; k++ ) *buffer++ = 32;
        }
        for( const char* p = ending; *p; *buffer++ = *p++ );
    }

    if( !annotation )
        if( *( buffer - 1 ) != '\n' )
            *buffer++ = '\n';

    *buffer++ = 0;

    print_func( vpn, "%s", GetPrivateBuffer());
}

// -------------------------------------------------------------------------
// FormatBuffer: auxiliary method to format character buffer; note that
//     space required for the first argument must be pre-allocated
//
void ExtendedDistributionMatrix::FormatBuffer( char*& format, const char* desc,
        size_t& tot, size_t& pos,
        const size_t size,
        const size_t indent, const size_t width ) const
{
    size_t  loc = 0;
    size_t  lto = 0;
    const char* beg = ( tot <= 2 )? NULL: desc;//if the very beginning, don't check for line feed

    while( *desc && tot < size ) {
        if( *desc == 32 || *desc == 9 || desc == beg ) {
            loc = pos + 1;
            lto = tot + 1;
            for( const char* p = desc + 1;
                    *p && *p != 32 && *p != 9 && loc < width && lto < size;
                     p++, loc++, lto++ );
            if( width <= loc && indent < pos ) {
                *format++ = '\n';
                for( size_t k = 0; k < indent; k++ ) *format++ = 32;
                tot += indent;
                pos  = indent;
                if( *desc == 32 || *desc == 9 )
                    desc++;
                if( beg )
                    beg = NULL;
                continue;
            }
        }
        *format++ = *desc++;
        pos++; tot++;
    }
}

// -------------------------------------------------------------------------
// Serialize: write the class data to file
//
void ExtendedDistributionMatrix::Serialize( Serializer& serializer ) const
{
    size_t  sz_name = 0;
    size_t  sz_description = 0;
    int n;
    //
    DistributionMatrix::Serialize( serializer );
    //
    if( columns > 0 ) {
        serializer.Write(( char* )freqweights, sizeof( double ), columns );
        serializer.Write(( char* )information, sizeof( double ), columns );
        for( n = 0; n < columns; n++ )
            serializer.Write(( char* )expMIDs_[n], sizeof( double ), PS_NSTATES );
    }
    //
    if( GetName())
        sz_name = strlen( GetName()) + 1; // to include null symbol

    serializer.Write(( char* )&sz_name, sizeof( sz_name ), 1 );

    if( GetName())
        serializer.Write(( char* )GetName(), 1, sz_name );
    //
    if( GetDescription())
        sz_description = strlen( GetDescription()) + 1; // to include null symbol

    serializer.Write(( char* )&sz_description, sizeof( sz_description ), 1 );

    if( GetDescription())
        serializer.Write(( char* )GetDescription(), 1, sz_description );


    serializer.Write(( char* )&nosequences,         sizeof( nosequences ), 1 );
    serializer.Write(( char* )&effnosequences,      sizeof( effnosequences ), 1 );

    serializer.Write(( char* )&referenceLambda,     sizeof( referenceLambda ), 1 );
    serializer.Write(( char* )&referenceK,          sizeof( referenceK ), 1 );
    serializer.Write(( char* )&lambda,              sizeof( lambda ), 1 );
    serializer.Write(( char* )&entropy,             sizeof( entropy ), 1 );
    serializer.Write(( char* )&parameterK,          sizeof( parameterK ), 1 );
    serializer.Write(( char* )&expscore,            sizeof( expscore ), 1 );

}

// -------------------------------------------------------------------------
// Deserialize: read data into the class members
// -------------------------------------------------------------------------

void ExtendedDistributionMatrix::Deserialize( Serializer& serializer )
{
    size_t  sz_name;
    size_t  sz_description;
    int n;
    //
    //memory allocation is there!
    DistributionMatrix::Deserialize( serializer );
    //
    if( columns > 0 ) {
        serializer.Read(( char* )freqweights, sizeof( double ), columns );
        serializer.Read(( char* )information, sizeof( double ), columns );
        for( n = 0; n < columns; n++ )
            serializer.Read(( char* )expMIDs_[n], sizeof( double ), PS_NSTATES );
    }
    //
    serializer.Read(( char* )&sz_name, sizeof( sz_name ), 1 );

    if( sz_name ) {
        char*   newname = ( char* )malloc( sz_name );
        if( !newname )
            throw myruntime_error( mystring( "ExtendedDistributionMatrix: Not enough memory." ));

        serializer.Read( newname, 1, sz_name );
        newname[ sz_name - 1 ] = 0;  // to avoid crashes

        SetName( newname );
        free( newname );
    } else
        SetName( NULL );
    //
    serializer.Read(( char* )&sz_description, sizeof( sz_description ), 1 );

    if( sz_description ) {
        char*   newdesc = ( char* )malloc( sz_description );
        if( !newdesc )
            throw myruntime_error( mystring( "ExtendedDistributionMatrix: Not enough memory." ));

        serializer.Read( newdesc, 1, sz_description );
        newdesc[ sz_description - 1 ] = 0;  // to avoid crashes

        SetDescription( newdesc );
        free( newdesc );
    } else
        SetDescription( NULL );

    serializer.Read(( char* )&nosequences,          sizeof( nosequences ), 1 );
    serializer.Read(( char* )&effnosequences,       sizeof( effnosequences ), 1 );

    serializer.Read(( char* )&referenceLambda,      sizeof( referenceLambda ), 1 );
    serializer.Read(( char* )&referenceK,           sizeof( referenceK ), 1 );
    serializer.Read(( char* )&lambda,               sizeof( lambda ), 1 );
    serializer.Read(( char* )&entropy,              sizeof( entropy ), 1 );
    serializer.Read(( char* )&parameterK,           sizeof( parameterK ), 1 );
    serializer.Read(( char* )&expscore,             sizeof( expscore ), 1 );
}





// =========================================================================

const char* patstrDATVER = "COMER profile v";
const char* patstrDESC = "DESC:";
const char* patstrFILE = "FILE:";
const char* patstrCMD = "CMD:";
const char* patstrLEN = "LEN:";
const char* patstrNOS = "NOS:";
const char* patstrEFFNOS = "EFFNOS:";
const char* patstrSCALE = "SCALE:";
const char* patstrNULL = "NULL";
const char* patstrPOST = "POST";
const char* patstrCT = "CT:";
const char* patstrCV = "CV:";
const char* patstrSS = "SS:";
const char* patstrSPCOMP = "Computed  ungapped,";
const char* patstrSPREFR = "Reference ungapped,";
const char* patstrENTROPY = "Entropy,";
const char* patstrEXPECTED = "Expected,";
const char* patstrEXPNN = "Expected score per pos. is non-negative,";
const char* patstrEND = "*";
const int   lenstrEND = strlen( patstrEND );
const int   INTSCALE = 1000;
const int   SSSCALE = 10000;

// =========================================================================
// PrintParameters: print statistical parameters of profile
//

void ExtendedDistributionMatrix::PrintParameters( FILE* fp ) const
{
    if( fp == NULL )
        return;

    if( 0.0 <= GetExpectedScore()) {
        fprintf( fp, "%s %.4f!\n", patstrEXPNN, GetExpectedScore());
        return;
    }

    fprintf( fp, "%-25s  %-6s   %-6s\n", " ", "K", "Lambda" );
    fprintf( fp, "%-25s  %6.4f   %6.4f\n", patstrSPCOMP, GetK(),    GetLambda());
    fprintf( fp, "%-25s  %6.4f   %6.4f\n", patstrSPREFR, GetRefK(), GetRefLambda());
    fprintf( fp, "%s %6.4f; %s %6.4f\n", patstrENTROPY, GetEntropy(), patstrEXPECTED, GetExpectedScore());
}

// =========================================================================
// TextWriteProfile: write profile data to file
//
void TextWriteProfile( 
    FILE* fp,
    const FrequencyMatrix& frqs, const LogOddsMatrix& pssm, const GapScheme& gaps, 
    int scale )
{
    if( !fp )
        return;

    if( frqs.GetColumns() != pssm.GetColumns() ||
        frqs.GetColumns() != gaps.GetColumns())
        throw myruntime_error( "TextWriteProfile: Inconsistent Profile data." );

    if( scale < 1 )
        throw myruntime_error( "TextWriteProfile: Invalid scale factor." );

    const int   noress = NUMALPH; //gap symbol included
    double          bppprob;//background posterior predicitive
    const double*   ppprobs;//posterior predictive probabilities
    const int*      pppndxs;//indices of posteriors
    size_t          noppps;//number of p.p.probability values
    double          norm2;//context vector's squared norm
    double          lpprb;//log of context vector's prior probability
    const double*   ctxvec;//context vector
    int             vsize = pssm.GetCtxVecSize();//size of vector
    char        res;
    int         m, r, t, n = 0;

    fprintf( fp, "%s%s\n", patstrDATVER, dataversion );
    fprintf( fp, "%-9s%s\n", patstrDESC, pssm.GetDescription()? pssm.GetDescription(): "" );
    fprintf( fp, "%-9s%s\n", patstrFILE, pssm.GetName()? pssm.GetName(): "" );

    fprintf( fp, "%-9s", patstrCMD );   print_cmdline( &file_print, fp );
    fprintf( fp, "%-9s%d\n", patstrLEN, pssm.GetColumns());
    fprintf( fp, "%-9s%d\n", patstrNOS, pssm.GetNoSequences());
    fprintf( fp, "%-9s%.1f\n", patstrEFFNOS, pssm.GetEffNoSequences());
    fprintf( fp, "%-9s%d\n", patstrSCALE, scale );

    fprintf( fp, "# Scores / Frequencies / ");
    for( t = 0; t < P_NSTATES; t++ ) fprintf( fp, "%s ", gTPTRANS_NAMES[t]);
    fprintf( fp, "/ NexpM NexpI NexpD; Weight; Information /\n");
    fprintf( fp, "# BPProb NoPosts / PProbs / PPNdxs / CT:NoPosts BPProb PProbs PPNdxs ");
    fprintf( fp, "/ CV: Prior Norm2 N Vector / SS:Pred Prob (C,E,H)\n");
//     fprintf( fp, "/ NexpM NexpI NexpD; Weight; Information; Gap weight, deletion at beg., end, interval\n");
//     fprintf( fp, "/ Weight, information, exp. no. observations; Gap weight, deletion at beg., end, interval\n");

    fprintf( fp, "%9c", 32 );

    for( r = 0; r < noress; r++ )
        fprintf( fp, " %7c", DehashCode( r ) );


    fprintf( fp, "\n%7s   ", patstrNULL );
    for( r = 0; r < noress; r++ )
        fprintf( fp, "%7d ", ( int )rint( scale * pssm.GetBackProbsAt( r )));

    fprintf( fp, "\n%7s   ", patstrPOST );
    for( r = 0; r < noress; r++ )
        fprintf( fp, "%7d ", ( int )rint( scale * pssm.GetPostProbsAt( r )));

    if( 0 < pssm.GetColumns()) {
        fprintf( fp, "\n%10c", 32 );
        for( t = 0; t < P_NSTATES; t++ )
            fprintf( fp, "%7d ", ( int )rint( scale * gaps.GetOrgTrProbsAt( t, -1 )));

        fprintf( fp, "\n%10c", 32 );
        for( t = 0; t < PS_NSTATES; t++ )
            fprintf( fp, "%7d ", ( int )rint( INTSCALE * pssm.GetMIDExpNoObservationsBeg( t )));
    }

    for( m = 0; m < pssm.GetColumns(); m++ ) {
        res = pssm.GetResidueAt( m );
        if( res != frqs.GetResidueAt( m ) || res != gaps.AAcid( m ))
            throw myruntime_error( "TextWriteProfile: Inconsistent Profile data." );
        // omit unused positions and gaps in query
        if( res == GAP )
            continue;

        fprintf( fp, "\n%5d %c   ", ++n, DehashCode( pssm.GetResidueAt( m )));

        for( r = 0; r < noress; r++ )
            fprintf( fp, "%7d ", ( int )rint( scale * pssm.GetValueAt( m, r )));

        fprintf( fp, "\n%10c", 32 );

        for( r = 0; r < noress; r++ )
            fprintf( fp, "%7d ", ( int )rint( scale * frqs.GetValueAt( m, r )));

        fprintf( fp, "\n%10c", 32 );

        for( t = 0; t < P_NSTATES; t++ )
            fprintf( fp, "%7d ", ( int )rint( scale * gaps.GetOrgTrProbsAt( t, m )));

        fprintf( fp, "\n%10c", 32 );

        for( t = 0; t < PS_NSTATES; t++ )
            fprintf( fp, "%7d ", ( int )rint( INTSCALE * pssm.GetMIDExpNoObservationsAt( m, t )));

        fprintf( fp, "%7d %7d ",
                ( int )rint( scale * pssm.GetFrequencyWeightAt( m )),
                ( int )rint( scale * pssm.GetInformationAt( m )));

//         fprintf( fp, "%7d %7d %7d %7d ",
//                 ( int )rint( scale * gaps.GetWeightsAt( m )),
//                 ( int )rint( scale * gaps.GetDeletesBegAt( m )),
//                 ( int )rint( scale * gaps.GetDeletesEndAt( m )),
//                 gaps.GetDeletesIntervalAt( m ));

        //{{HDP1
        bppprob = pssm.GetBckPPProbAt(m);
        ppprobs = pssm.GetPPProbsAt(m);
        pppndxs = pssm.GetPPPIndsAt(m);
        noppps = pssm.GetNoPPProbsAt(m);

        fprintf( fp, "\n%10c", 32 );
        fprintf( fp, "%7d %7d ", (int)rint( scale*bppprob ), noppps );

        if( 0 < (int)noppps ) {
            if( ppprobs == NULL || pppndxs == NULL )
                throw myruntime_error("TextWriteProfile: Null posterior probabilities.");

            fprintf( fp, "\n%10c", 32 );

            for( t = 0; t < noppps; t++ )
                fprintf( fp, "%7d ", (int)rint( scale*ppprobs[t]));

            fprintf( fp, "\n%10c", 32 );

            for( t = 0; t < noppps; t++ )
                fprintf( fp, "%7d ", pppndxs[t]);
        }
        //}}

        if( pssm.GetctPsSet()) {
            //HDP ctx
            bppprob = pssm.GetctBckPPProbAt(m);
            ppprobs = pssm.GetctPPProbsAt(m);
            pppndxs = pssm.GetctPPPIndsAt(m);
            noppps = pssm.GetNoctPPProbsAt(m);

            fprintf( fp, "\n%13c%s%d %7d ", 32, patstrCT, noppps, (int)rint( scale*bppprob ));

            if( 0 < (int)noppps ) {
                if( ppprobs == NULL || pppndxs == NULL )
                    throw myruntime_error("TextWriteProfile: Null ct posterior probabilities.");

                for( t = 0; t < noppps; t++ )
                    fprintf( fp, "%7d ", (int)rint( scale*ppprobs[t]));

                for( t = 0; t < noppps; t++ )
                    fprintf( fp, "%7d ", pppndxs[t]);
            }
        }

        //{{ context vector
        if( pssm.GetCtxVecSet()) {
            norm2 = pssm.GetCtxVecNorm2At(m);
            lpprb = pssm.GetCtxVecLpprobAt(m);
            ctxvec = pssm.GetCtxVecAt(m);

            fprintf( fp, "\n%14c%s %7d %7d %7d ", 32, patstrCV,
                   (int)rint( scale*lpprb ), (int)rint( scale*norm2 ), vsize );

            if( 0 < vsize ) {
                if( ctxvec == NULL )
                    throw myruntime_error("TextWriteProfile: Null context vector.");

                for( t = 0; t < vsize; t++ )
                    fprintf( fp, "%7d ", (int)rint( scale*ctxvec[t]));
            }
        }
        //}}

        //{{
        if( pssm.GetSSSSet()) {
            //use floor in rounding SS state probability
            fprintf( fp, "\n%13c%s%c", 32, patstrSS, DehashSSCode( pssm.GetSSStateAt(m)));
            if( pssm.GetSSSP3()) {
                for( t = 0; t < SS_NSTATES; t++ )
                    fprintf( fp, " %7d", (int)( SSSCALE*pssm.GetSSStateProbAt(m,t)));
            }
            else
                fprintf( fp, " %7d", (int)( SSSCALE*pssm.GetSSStateProbAt(m)));
        }
        //}}
    }

    fprintf( fp, "\n%10c%7d %7d", 32,
                ( int )rint( scale * gaps.GetOpenCost()),
                ( int )rint( scale * gaps.GetExtendCost()));
    fprintf( fp, "\n" );

// #ifdef SCALE_PROFILES
    pssm.PrintParameters( fp );
// #endif
    fprintf( fp, "%s\n", patstrEND );
}

// -------------------------------------------------------------------------
// TextReadProfile: read profile data from file
//
void TextReadProfile( FILE* fp, FrequencyMatrix& frqs, LogOddsMatrix& pssm, GapScheme& gaps )
{
    if( !fp )
        return;

    const double    accuracy = 1.0e-3;
    const int       maxnocls = 10000;
    const int       noress = NUMALPH;
    size_t          length, rbts;
    bool            lineread = false;
    const size_t    locsize = 10*KBYTE;
    char            locbuffer[locsize+1] = {0};
    char*           p;
    int             emsg;
    mystring    desc, file;
    const char* descp, *filep;
    int         prolen = 0;
    int         nos = 0;
    double      effnos = 0;
    int         scale = 0;
    char        res, sss;
    int         intval;
    double      value, consv;
    double      scores[NUMALPH];
    double      freqns[NUMALPH];
    double      trnpro[P_NSTATES];
    double      midobs[PS_NSTATES];
    double      weight, inform;
    double      gapwgt, delbeg, delend, open, extend;
    double      expnobs, sssprob[SS_NSTATES];
    int         delint;
    double      bppprob, ctbppprob;//background posterior predicitive
    Pslvector   ppprobs, ctppprobs;//posterior predictive probabilities
    Ivector     pppndxs, ctpppndxs;//indices of posteriors
    size_t      noppps, noctppps;//number of p.p.probability values
    double      norm2, lpprb;//squared norm and log of prior probability for context vector
    Pslvector   ctxvec;//context vector
    int         vsize = 0;//vector size
    double      refunlambda, refunK;
    double      cmpunlambda, cmpunK;
    double      proentropy,  proexpected;
    int         m, r, t, t0, tn, n = 0;

    memset( scores, 0, NUMALPH * sizeof( double ));
    memset( freqns, 0, NUMALPH * sizeof( double ));

    //read version number
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong profile format." );

    if(( p = strstr( locbuffer, patstrDATVER )) == NULL )
        throw myruntime_error( "Wrong profile format." );

    p += strlen( patstrDATVER );

    if( length <= size_t( p - locbuffer ))
        throw myruntime_error( mystring( "Wrong profile format." ));

    if(( p = strstr( p, dataversion )) == NULL )
        throw myruntime_error( "Wrong data version number." );


    //read description line
    if(( emsg = skip_comments( fp, desc )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( feof( fp ) || desc.empty())
        throw myruntime_error( "Wrong profile format." );

    if(( p = ( char* )strstr( desc.c_str(), patstrDESC )) == NULL )
        throw myruntime_error( "Wrong profile format." );

    descp = p + strlen( patstrDESC );
    for( ; *descp == ' ' || *descp == '\t'; descp++ );
    for( n = ( int )desc.length() - 1; 0 <= n && ( desc[n] == '\n' || desc[n] == '\r' ); n-- ) desc[n] = 0;


    //read filename
    if(( emsg = skip_comments( fp, file )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( feof( fp ) || file.empty())
        throw myruntime_error( "Wrong profile format." );

    if(( p = ( char* )strstr( file.c_str(), patstrFILE )) == NULL )
        throw myruntime_error( "Wrong profile format." );

    filep = p + strlen( patstrFILE );
    for( ; *filep == ' ' || *filep == '\t' ; filep++ );
    for( n = ( int )file.length() - 1; 0 <= n && ( file[n] == '\n' || file[n] == '\r' ); n-- ) file[n] = 0;


    //read command line
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong profile format." );

    if(( p = strstr( locbuffer, patstrCMD )) == NULL )
        throw myruntime_error( "Wrong profile format." );


    //read profile length
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong profile format." );

    if(( p = strstr( locbuffer, patstrLEN )) == NULL )
        throw myruntime_error( "Wrong profile format." );

    p += strlen( patstrLEN );

    if( length <= size_t( p - locbuffer ))
        throw myruntime_error( "Wrong profile format." );

    if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &prolen, &rbts )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( prolen < 1 )
        throw myruntime_error( "Wrong profile format: Invalid profile length." );


    //read number of sequences
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong profile format." );

    if(( p = strstr( locbuffer, patstrNOS )) == NULL )
        throw myruntime_error( "Wrong profile format." );

    p += strlen( patstrNOS );

    if( length <= size_t( p - locbuffer ))
        throw myruntime_error( "Wrong profile format." );

    if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &nos, &rbts )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( nos < 1 )
        throw myruntime_error( "Wrong profile format: Invalid value of no. sequences." );


    //read effective number of sequences
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong profile format." );

    if(( p = strstr( locbuffer, patstrEFFNOS )) == NULL )
        throw myruntime_error( "Wrong profile format." );

    p += strlen( patstrEFFNOS );

    if( length <= size_t( p - locbuffer ))
        throw myruntime_error( "Wrong profile format." );

    if(( emsg = read_double( p, length - size_t( p - locbuffer ), &effnos, &rbts )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( effnos <= 0.0 )
        throw myruntime_error( "Wrong profile format: Invalid value of eff. no. sequences." );


    //read scale factor
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong profile format." );

    if(( p = strstr( locbuffer, patstrSCALE )) == NULL )
        throw myruntime_error( "Wrong profile format." );

    p += strlen( patstrSCALE );

    if( length <= size_t( p - locbuffer ))
        throw myruntime_error( "Wrong profile format." );

    if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &scale, &rbts )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( scale < 1 )
        throw myruntime_error( "Wrong profile format: Invalid scale factor." );


    //read amino acid symbols; omitting check of ordering
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));


    //read background probabilities
    if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong profile format." );

    if(( p = strstr( locbuffer, patstrNULL )) == NULL )
        throw myruntime_error( "Wrong profile format: Background probabilities missing." );

    p += strlen( patstrNULL );

    consv = 0.0;
    for( r = 0; r < noress; r++ ) {
        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No background probabilities." );

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        p += rbts;

        if( intval < 0 || scale < intval )
            throw myruntime_error( "Wrong profile format: Invalid probability value." );

        value = ( double )intval / ( double )scale;
        consv += value;
        pssm.SetBackProbsAt( r, value );
    }

    if( consv < 1.0 - accuracy || consv > 1.0 + accuracy )
        throw myruntime_error( mystring( "Invalid profile data: Probabilities are not conserved." ));
    //


    //read posterior (generalized target) probabilities
    if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong profile format." );

    if(( p = strstr( locbuffer, patstrPOST )) == NULL )
        throw myruntime_error( "Wrong profile format: Posterior probabilities missing." );

    p += strlen( patstrPOST );

    consv = 0.0;
    for( r = 0; r < noress; r++ ) {
        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No posterior probabilities." );

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        p += rbts;

        if( intval < 0 || scale < intval )
            throw myruntime_error( "Wrong profile format: Invalid post. probability value." );

        value = ( double )intval / ( double )scale;
        consv += value;
        pssm.SetPostProbsAt( r, value );
    }

    if( consv < 1.0 - accuracy || consv > 1.0 + accuracy )
        throw myruntime_error( mystring( "Invalid profile data: Post. probabilities are not conserved." ));
    //


    //beginning transition probabilities
    if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong profile format." );

    for( t = 0; t < P_NSTATES; t++ )
    {
        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No transition probabilities." );

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        p += rbts;

        trnpro[t] = ( double )intval / ( double )scale;
    }
    //


    //beginning MID state expected observations
    if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong profile format." );

    for( t = 0; t < PS_NSTATES; t++ )
    {
        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No MID state obervations." );

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        p += rbts;

        midobs[t] = ( double )intval / ( double )INTSCALE;
    }
    //

    if( MAXCOLUMNS < prolen )
        throw myruntime_error( "Wrong profile format: Profile length is too large." );

    pssm.Clear();
    frqs.Clear();
    gaps.Clear();

    pssm.Reserve( prolen );
    frqs.Reserve( prolen );
    gaps.Reserve( prolen );

    if( *descp )
        pssm.SetDescription( descp );
    else
        pssm.SetDescription( NULL );

    if( *filep )
        pssm.SetName( filep );
    else
        pssm.SetName( NULL );

    pssm.SetNoSequences( nos );
    pssm.SetEffNoSequences( effnos );
    pssm.SetMIDExpNoObservationsBeg( midobs );

    gaps.SetOrgTrProbsBeg( &trnpro );//set beginning transition probabilities here

    lineread = false;

    for( m = 0; m < prolen; m++ )
    {
        //scores
        if( !lineread ) {
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw myruntime_error( "Wrong profile format." );
        }

        lineread = false;

        if(( emsg = read_integer( p = locbuffer, length, &n, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( n != m + 1 )
            throw myruntime_error( "Wrong profile format: Odering." );

        p += rbts;

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No residue." );

        if(( emsg = read_symbol( p, length - size_t( p - locbuffer ), &res, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        p += rbts;

        res = HashAlphSymbol( res );

        for( r = 0; r < noress; r++ )
        {
            if( length <= size_t( p - locbuffer ))
                throw myruntime_error( "Wrong profile format: No scores." );

            if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            p += rbts;

            scores[r] = ( double )intval / ( double )scale;
        }


        //frequencies
        if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error( "Wrong profile format." );

        for( r = 0; r < noress; r++ )
        {
            if( length <= size_t( p - locbuffer ))
                throw myruntime_error( "Wrong profile format: No frequencies." );

            if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            p += rbts;

            freqns[r] = ( double )intval / ( double )scale;
        }


        //transition probabilities
        if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error( "Wrong profile format." );

        for( t = 0; t < P_NSTATES; t++ )
        {
            if( length <= size_t( p - locbuffer ))
                throw myruntime_error( "Wrong profile format: No transition probabilities." );

            if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            p += rbts;

            trnpro[t] = ( double )intval / ( double )scale;
        }


        //MID state observations
        if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error( "Wrong profile format." );

        for( t = 0; t < PS_NSTATES; t++ )
        {
            if( length <= size_t( p - locbuffer ))
                throw myruntime_error( "Wrong profile format: No MID state observations." );

            if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( intval < 0 )
                throw myruntime_error( "Wrong profile format: Invalid value of MID state observations." );

            p += rbts;

            midobs[t] = ( double )intval / ( double )INTSCALE;
        }


        //weight, information
        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No weight value." );

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( intval < 0 )
            throw myruntime_error( "Wrong profile format: Negative weight." );

        weight = ( double )intval / ( double )scale;

        p += rbts;

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No information values." );

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( intval < 0 )
            throw myruntime_error( "Wrong profile format: Negative information value." );

        inform = ( double )intval / ( double )scale;

        p += rbts;


//         //gap weights, deletes at beg., end, its intervals
//         if( length <= size_t( p - locbuffer ))
//             throw myruntime_error( "Wrong profile format: No gap weights." );
// 
//         if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
//             throw myruntime_error( TranslateReadError( emsg ));
// 
//         if( intval < 0 || scale < intval )
//             throw myruntime_error( "Wrong profile format: Negative gap weights." );
// 
//         gapwgt = ( double )intval / scale;
// 
//         p += rbts;
// 
//         if( length <= size_t( p - locbuffer ))
//             throw myruntime_error( "Wrong profile format: No deletion values." );
// 
//         if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
//             throw myruntime_error( TranslateReadError( emsg ));
// 
//         if( intval < 0 || scale < intval )
//             throw myruntime_error( "Wrong profile format: Invalid deletion values." );
// 
//         delbeg = ( double )intval / scale;
// 
//         p += rbts;
// 
//         if( length <= size_t( p - locbuffer ))
//             throw myruntime_error( "Wrong profile format: No deletion values." );
// 
//         if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
//             throw myruntime_error( TranslateReadError( emsg ));
// 
//         if( intval < 0 || scale < intval )
//             throw myruntime_error( "Wrong profile format: Invalid deletion values." );
// 
//         delend = ( double )intval / scale;
// 
//         p += rbts;
// 
//         if( length <= size_t( p - locbuffer ))
//             throw myruntime_error( "Wrong profile format: No deletion interval." );
// 
//         if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
//             throw myruntime_error( TranslateReadError( emsg ));
// 
//         if( intval < 0 )
//             throw myruntime_error( "Wrong profile format: Invalid deletion interval." );
// 
//         delint = intval;


        //CLUSTERS: HDP1
        //bck posterior probability, no posterior probability values
        if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong profile format.");

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error("Wrong profile format: No HDP1 back. posterior probability.");

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( intval < 0 )
            throw myruntime_error("Wrong profile format: Negative HDP1 back. posterior probability.");

        bppprob = (double)intval / (double)scale;

        p += rbts;

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error("Wrong profile format: No no. HDP1 post. probabilities.");

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( intval < 0 )
            throw myruntime_error("Wrong profile format: Negative no. HDP1 post. probabilities.");
        if( maxnocls < intval )
            throw myruntime_error("Wrong profile format: Too large no. HDP1 post. probabilities.");

        noppps = (size_t)intval;
        if( noppps ) {
            ppprobs.Allocate((int)noppps );
            pppndxs.Allocate((int)noppps );
            ppprobs.Clear();
            pppndxs.Clear();

            //posterior predictive probabilities
            if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw myruntime_error("Wrong profile format.");

            for( t = 0; t < (int)noppps; t++ )
            {
                if( length <= size_t( p - locbuffer ))
                    throw myruntime_error("Wrong profile format: No HDP1 post. probabilities.");

                if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));

                if( intval < 0 )
                    throw myruntime_error("Wrong profile format: Negative HDP1 post. probability value.");

                p += rbts;

                ppprobs.Push(( double)intval / (double)scale );
            }

            //indices of posteriors
            if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw myruntime_error("Wrong profile format.");

            for( t = 0; t < (int)noppps; t++ )
            {
                if( length <= size_t( p - locbuffer ))
                    throw myruntime_error("Wrong profile format: No HDP1 probability indices.");

                if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));

                if( intval != -1 && intval < 0 )
                    throw myruntime_error("Wrong profile format: Negative HDP1 post. probability index.");

                p += rbts;

                pppndxs.Push( intval );
            }
        }//if( noppps )


        //HDP ctx
        if( m < 1 )
            pssm.SetctPsSet( true );//probe for HDP ctx information
        if( pssm.GetctPsSet()) {
            //bck posterior probability, no posterior probability values
            if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw myruntime_error("Wrong profile format.");

            lineread = true;

            if(( p = strstr( locbuffer, patstrCT )) == NULL ) {
                if( m < 1 )
                    pssm.SetctPsSet( false );//no HDP ctx information
                else
                    throw myruntime_error("Wrong profile format: No HDP ctx probabilities.");
            }
            else {
                lineread = false;
                p += strlen( patstrCT );

                if( length <= size_t( p - locbuffer ))
                    throw myruntime_error("Wrong profile format: No no. HDP ctx prob. values.");

                if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));

                if( intval < 0 )
                    throw myruntime_error("Wrong profile format: Negative no. HDP ctx prob. values.");
                if( maxnocls < intval )
                    throw myruntime_error("Wrong profile format: Too large no. HDP ctx prob. values.");

                noctppps = (size_t)intval;

                p += rbts;

                if( length <= size_t( p - locbuffer ))
                    throw myruntime_error("Wrong profile format: No HDP ctx bck. probability.");

                if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));

                if( intval < 0 )
                    throw myruntime_error("Wrong profile format: Negative HDP ctx bck. probability.");

                ctbppprob = (double)intval / (double)scale;

                p += rbts;

                if( noctppps ) {
                    ctppprobs.Allocate((int)noctppps );
                    ctpppndxs.Allocate((int)noctppps );
                    ctppprobs.Clear();
                    ctpppndxs.Clear();

                    //posterior predictive probabilities
                    for( t = 0; t < (int)noctppps; t++ )
                    {
                        if( length <= size_t( p - locbuffer ))
                            throw myruntime_error("Wrong profile format: No HDP ctx probabilities.");

                        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
                            throw myruntime_error( TranslateReadError( emsg ));

                        if( intval < 0 )
                            throw myruntime_error("Wrong profile format: Negative HDP ctx probability value.");

                        p += rbts;

                        ctppprobs.Push(( double)intval / (double)scale );
                    }

                    //indices of posteriors
                    for( t = 0; t < (int)noctppps; t++ )
                    {
                        if( length <= size_t( p - locbuffer ))
                            throw myruntime_error("Wrong profile format: No HDP ctx probability indices.");

                        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
                            throw myruntime_error( TranslateReadError( emsg ));

                        if( intval != -1 && intval < 0 )
                            throw myruntime_error("Wrong profile format: Negative HDP ctx probability index.");

                        p += rbts;

                        ctpppndxs.Push( intval );
                    }
                }//if( noppps )
            }//else
        }//if( pssm.GetctPsSet())


        //context vector
        if( m < 1 )
            pssm.SetCtxVecSet( true );//probe for context vector data
        if( pssm.GetCtxVecSet()) {
            if( !lineread )
                if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw myruntime_error("Wrong profile format.");

            lineread = true;

            if(( p = strstr( locbuffer, patstrCV )) == NULL ) {
                if( m < 1 )
                    pssm.SetCtxVecSet( false );//no context vector data
                else
                    throw myruntime_error("Wrong profile format: No context vector data.");
            }
            else {
                lineread = false;
                p += strlen( patstrCV );

                if( length <= size_t( p - locbuffer ))
                    throw myruntime_error("Wrong profile format: No prior for context vector.");

                if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));

                lpprb = (double)intval / (double)scale;

                p += rbts;

                if( length <= size_t( p - locbuffer ))
                    throw myruntime_error("Wrong profile format: No context vector norm.");

                if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));

                if( intval < 0 )
                    throw myruntime_error("Wrong profile format: Negative context vector norm.");

                norm2 = (double)intval / (double)scale;

                p += rbts;

                if( length <= size_t( p - locbuffer ))
                    throw myruntime_error("Wrong profile format: No context vector size.");

                if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));

                if( intval < 0 )
                    throw myruntime_error("Wrong profile format: Negative context vector size.");
                if( 1000 < intval )
                    throw myruntime_error("Wrong profile format: Too large context vector size.");
                if( vsize && vsize != intval )
                    throw myruntime_error("Wrong profile format: Inconsistent context vector size.");

                vsize = intval;

                p += rbts;

                if( vsize ) {
                    ctxvec.Allocate( vsize );
                    ctxvec.Clear();

                    //vector elements
                    for( t = 0; t < vsize; t++ )
                    {
                        if( length <= size_t( p - locbuffer ))
                            throw myruntime_error("Wrong profile format: No context vector entries.");

                        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
                            throw myruntime_error( TranslateReadError( emsg ));

                        p += rbts;

                        ctxvec.Push(( double)intval / (double)scale );
                    }
                }//if( vsize )
            }//else
        }//if( pssm.GetCtxVecSet())


        //SS state
        if( m < 1 ) {
            pssm.SetSSSSet( true );//probe SS information existance
            pssm.SetSSSP3( true );
            t0 = 0; tn = SS_NSTATES-1;
        }
        if( pssm.GetSSSSet()) {
            if( !lineread )
                if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw myruntime_error( "Wrong profile format." );

            lineread = true;

            if(( p = strstr( locbuffer, patstrSS )) == NULL ) {
                if( m < 1 )
                    pssm.SetSSSSet( false );//no SS information
                else
                    throw myruntime_error( "Wrong profile format: No SS state." );
            }
            else {
                memset( sssprob, 0, SS_NSTATES*sizeof(double));

                lineread = false;
                p += strlen( patstrSS );

                if(( emsg = read_symbol( p, length - size_t( p - locbuffer ), &sss, &rbts )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));

                p += rbts;

                sss = HashSSState( sss );

                if( !pssm.GetSSSP3()) {
                    t0 = sss; tn = sss;
                }

                for( t = t0; t <= tn; t++ )
                {
                    if( length <= size_t( p - locbuffer )) {
                        if( m < 1 && t == t0+1 ) {
                            pssm.SetSSSP3( false );
                            sssprob[sss] = ( double )intval /( double )SSSCALE;
                            break;
                        }
                        else
                            throw myruntime_error( "Wrong profile format: No SS state probability." );
                    }

                    if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
                        throw myruntime_error( TranslateReadError( emsg ));

                    p += rbts;
                    for( ; *p == ' ' || *p == '\t' || *p == '\r' || *p == '\n'; p++ );

                    sssprob[t] = ( double )intval /( double )SSSCALE;
                }
            }
        }//if( pssm.GetSSSSet())


        //SAVE profile position
        frqs.PushAt( freqns, res, m );
        pssm.PushAt( scores, res, weight, inform, midobs, m );
        pssm.SetBckPPProbAt( bppprob, m );
        pssm.SetPPProbsAt( m, ppprobs.GetVector(), pppndxs.GetVector(), (int)noppps );
        if( pssm.GetctPsSet()) {
            pssm.SetctBckPPProbAt( ctbppprob, m );
            pssm.SetctPPProbsAt( m, ctppprobs.GetVector(), ctpppndxs.GetVector(), (int)noctppps );
        }
        if( pssm.GetCtxVecSet()) pssm.SetCtxVecPlusAt( m, norm2, lpprb, ctxvec.GetVector(), vsize );
        if( pssm.GetSSSSet()) pssm.SetSSStateAt( m, sss, sssprob );
        gaps.PushAt(&trnpro, gapwgt, delbeg, delend, delint, res, m );
    }

    //open, extend
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong profile format." );

    if(( emsg = read_integer( p = locbuffer, length, &intval, &rbts )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( 0 < intval )
        throw myruntime_error( "Wrong profile format." );

    open = ( double )intval / scale;

    p += rbts;

    if( length <= size_t( p - locbuffer ))
        throw myruntime_error( "Wrong profile format." );

    if(( emsg = read_integer( p, length, &intval, &rbts )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( 0 < intval )
        throw myruntime_error( "Wrong profile format." );

    extend = ( double )intval / scale;

    gaps.SetOpenCost( open );
    gaps.SetExtendCost( extend );
    gaps.Initialize();

// #ifdef SCALE_PROFILES
    //statistical parameters
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong profile format: No statistical parameters." );

    if(( p = strstr( locbuffer, patstrEXPNN )) != NULL ) {
        p += strlen( patstrEXPNN );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No expected value." );

        if(( emsg = read_double( p, length - size_t( p - locbuffer ), &proexpected, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( proexpected < 0.0 )
            throw myruntime_error( "Wrong profile format: Wrong expected value." );

        pssm.SetExpectedScore( proexpected );
    }
    else {
        //computed
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error( "Wrong profile format: No statistical parameters." );

        if(( p = strstr( locbuffer, patstrSPCOMP )) == NULL )
            throw myruntime_error( "Wrong profile format: No statistical parameters." );

        p += strlen( patstrSPCOMP );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No statistical parameters." );

        if(( emsg = read_double( p, length - size_t( p - locbuffer ), &cmpunK, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

//         if( cmpunK < 0.0 )
//             throw myruntime_error( "Wrong profile format: Invalid statistical parameters." );

        p += rbts;

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No statistical parameters." );

        if(( emsg = read_double( p, length - size_t( p - locbuffer ), &cmpunlambda, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

//         if( cmpunlambda < 0.0 )
//             throw myruntime_error( "Wrong profile format: Invalid statistical parameters." );


        //reference
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error( "Wrong profile format: No reference parameters." );

        if(( p = strstr( locbuffer, patstrSPREFR )) == NULL )
            throw myruntime_error( "Wrong profile format: No reference parameters." );

        p += strlen( patstrSPREFR );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No reference parameters." );

        if(( emsg = read_double( p, length - size_t( p - locbuffer ), &refunK, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( refunK < 0.0 )
            throw myruntime_error( "Wrong profile format: Invalid reference parameters." );

        p += rbts;

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No reference parameters." );

        if(( emsg = read_double( p, length - size_t( p - locbuffer ), &refunlambda, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( refunlambda < 0.0 )
            throw myruntime_error( "Wrong profile format: Invalid reference parameters." );


        //entropy, expected
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error( "Wrong profile format: No entropy." );

        if(( p = strstr( locbuffer, patstrENTROPY )) == NULL )
            throw myruntime_error( "Wrong profile format: No entropy." );

        p += strlen( patstrENTROPY );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No entropy." );

        if(( emsg = read_double( p, length - size_t( p - locbuffer ), &proentropy, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

//         if( proentropy < 0.0 )
//             throw myruntime_error( "Wrong profile format: Invalid reference parameters." );

        if(( p = strstr( locbuffer, patstrEXPECTED )) == NULL )
            throw myruntime_error( "Wrong profile format: No expected value." );

        p += strlen( patstrEXPECTED );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No expected value." );

        if(( emsg = read_double( p, length - size_t( p - locbuffer ), &proexpected, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));


        pssm.SetRefLambda( refunlambda );
        pssm.SetRefK( refunK );
        pssm.SetLambda( cmpunlambda );
        pssm.SetEntropy( proentropy );
        pssm.SetK( cmpunK );
        pssm.SetExpectedScore( proexpected );
    }
// #endif

    pssm.Finalize();


    //footer
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( length < lenstrEND )
        throw myruntime_error("Wrong profile format: No end symbol.");

    for( n = 0; n < lenstrEND; n++ )
        if( locbuffer[n] != patstrEND[n] )
            throw myruntime_error("Wrong profile format: Invalid end symbol.");
}

