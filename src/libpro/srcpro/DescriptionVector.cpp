/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#include <math.h>
#include <stdlib.h>


#include "rc.h"
#include "DescriptionVector.h"




////////////////////////////////////////////////////////////////////////////
// CLASS PosDescriptionVector
//
// -------------------------------------------------------------------------
// Default constructor
// -------------------------------------------------------------------------

PosDescriptionVector::PosDescriptionVector()
{
    Init();
}

// Initialization constructor
//

PosDescriptionVector::PosDescriptionVector( unsigned reservation )
{
    Init();
    Realloc( reservation );
}

// Copy constructor
//
PosDescriptionVector::PosDescriptionVector( const PosDescriptionVector& one )
{
    Init();
    if( !one.capacity())
        throw myruntime_error( mystring( "Initialization error: wrong argument's been passed." ));

    Realloc( one.capacity());
    *this = one;
}

// -------------------------------------------------------------------------
// Destructor
// -------------------------------------------------------------------------

PosDescriptionVector::~PosDescriptionVector()
{
    if( residues ) free( residues );
    if( flags ) free( flags );
    if( weights ) free( weights );
}

// -------------------------------------------------------------------------
// Init: initialization of the class members
// -------------------------------------------------------------------------

void PosDescriptionVector::Init()
{
    residues = NULL;
    flags = NULL;
    weights = NULL;
    glbwght_ = 0.0;
    length_ = 0;
    efflength_ = 0;
    capacity_ = 0;
    used = true;
    firstused_ = SIZE_MAX;
    lastused_ = SIZE_MAX;
    counter_ = SIZE_MAX;
    ResetCluster();
}

// -------------------------------------------------------------------------
// Assignment
// -------------------------------------------------------------------------

PosDescriptionVector& PosDescriptionVector::operator=( const PosDescriptionVector& one )
{
    if( capacity_ < one.capacity_ )
        Realloc( one.capacity_ );

    memcpy( residues,   one.residues,   sizeof( unsigned char ) * one.length_ );
    memcpy( flags,      one.flags,      sizeof( unsigned char ) * one.length_ );
    memcpy( weights,    one.weights,    sizeof( double ) * one.length_ );

    glbwght_ = one.glbwght_;
    used = one.used;
    length_ = one.length_;
    efflength_ = one.efflength_;
    firstused_ = one.firstused_;
    lastused_ = one.lastused_;
    cluster_ = one.cluster_;

    memset( residues + length_,   0,  sizeof( unsigned char ) * ( capacity_ - length_ ));
    memset( flags + length_,      0,  sizeof( unsigned char ) * ( capacity_ - length_ ));
    memset( weights + length_,    0,  sizeof( double ) * ( capacity_ - length_ ));

    return *this;
}

// -------------------------------------------------------------------------
// Realloc: tries to reallocate necessary memory
// -------------------------------------------------------------------------

void PosDescriptionVector::Realloc( int newcap )
{
    if( capacity_ == 0 ) {
        residues = ( unsigned char* )malloc( sizeof( unsigned char ) * newcap );
        flags = ( unsigned char* )malloc( sizeof( unsigned char ) * newcap );
        weights = ( double* )malloc( sizeof( double ) * newcap );
    } else {
        residues = ( unsigned char* )realloc( residues, sizeof( unsigned char ) * newcap );
        flags = ( unsigned char* )realloc( flags, sizeof( unsigned char ) * newcap );
        weights = ( double* )realloc( weights, sizeof( double ) * newcap );
    }

    if( !residues || !flags || !weights )
        throw myruntime_error( mystring( "Not enough memory." ));

    unsigned char*      tress = residues;
    unsigned char*      tflgs = flags;
    double*             tweis = weights;


    if( capacity_ != 0 ) {
        tress = residues + capacity_;
        tflgs = flags + capacity_;
        tweis = weights + capacity_;
    }

    memset( tress, 0, sizeof( unsigned char ) * ( newcap - capacity_ ));
    memset( tflgs, 0, sizeof( unsigned char ) * ( newcap - capacity_ ));
    memset( tweis, 0, sizeof( double ) * ( newcap - capacity_ ));

    capacity_ = newcap;
}

// -------------------------------------------------------------------------
// push: pushes all information describing one position of the sequence
// -------------------------------------------------------------------------

void PosDescriptionVector::push( unsigned char r, unsigned char f, double w )
{
    if( capacity_ < length_ + 1 ) {
        Realloc(  capacity_ * 2 );
    }

    residues[length_] = r;
    flags[length_] = f;
    weights[length_] = w;

    length_++;

    if( r != GAP )
        efflength_++;
}

// -------------------------------------------------------------------------
// clear: clears all buffers allocated for the sequence
// -------------------------------------------------------------------------

void PosDescriptionVector::clear()
{
    memset( residues, 0, sizeof( unsigned char ) * capacity_ );
    memset( flags, 0, sizeof( unsigned char ) * capacity_ );
    memset( weights, 0, sizeof( double ) * capacity_ );

    length_ = 0;
    efflength_ = 0;
}





////////////////////////////////////////////////////////////////////////////
// CLASS ExtendedDescriptionVector
//
double* ExtendedDescriptionVector::prexpnores_ = NULL;

// -------------------------------------------------------------------------
// Constructor
//
ExtendedDescriptionVector::ExtendedDescriptionVector( unsigned reservation )
{
    Init();
    InitPrexpNoDistinctRes();
    Realloc( reservation );
}

// Constructor
//
ExtendedDescriptionVector::ExtendedDescriptionVector( const PosDescriptionVector& sequence )
{
    Init();
    InitPrexpNoDistinctRes();
    if( !sequence.capacity())
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Initialization: Wrong argument." ));

    Realloc( sequence.capacity());
    PosDescriptionVector::operator=( sequence );
}

// Default constructor
//
ExtendedDescriptionVector::ExtendedDescriptionVector()
{
    throw myruntime_error( mystring( "ExtendedDescriptionVector: Default initialization is impossible." ));
}

// -------------------------------------------------------------------------
// Destructor
//
ExtendedDescriptionVector::~ExtendedDescriptionVector()
{
    size_t  n;
    DestroyPrexpNoDistinctRes();

    if( indices ) {
        for( n = 0; n < length_; n++ )
            if( indices[n] )
                free( indices[n] );
        free( indices );
    }
    if( indicesLengths_ )       free( indicesLengths_ );
    if( indicesCapacities_ )    free( indicesCapacities_ );

    if( states )        free( states );
    if( extents )       free( extents );
    if( nodifrs_ )      free( nodifrs_ );
    if( MSextents )     free( MSextents );
    if( counts )        free( counts );

    FreeSqnWeightsAt();

    if( distribution )  free( distribution );
    if( matchWeights )  free( matchWeights );
    if( transWeights )  free( transWeights );
    if( gapWeights   )  free( gapWeights   );
    if( distinctHist )  free( distinctHist );
    if( MSdistinctHist )free( MSdistinctHist );
    if( targetFreqns )  free( targetFreqns );
    if( targetTranst )  free( targetTranst );
    if( rawPSSMatrix )  free( rawPSSMatrix );
    if( bppprob_ )      free( bppprob_ );
    if( noppps_ )       free( noppps_ );
    if( information  )  free( information  );
    if( noseqs_ )       free( noseqs_ );
    if( expnoseqs_ )    free( expnoseqs_ );
    if( expnosmid_ )    free( expnosmid_ );

    if( ppprobs_ ) {
        for( n = 0; n < length_; n++ )
            if( ppprobs_[n] )
                free( ppprobs_[n] );
        free( ppprobs_ );
    }
    if( pppndxs_ ) {
        for( n = 0; n < length_; n++ )
            if( pppndxs_[n] )
                free( pppndxs_[n] );
        free( pppndxs_ );
    }
}

// -------------------------------------------------------------------------
// FreeSqnWeightsAt: free sequence weights of all states
//
void ExtendedDescriptionVector::FreeSqnWeightsAt()
{
    int     st;
    size_t  n;
    for( st = 0; st < PS_NSTATES; st++ ) {
        if( sqnweights_[st]) {
            for( n = 0; n < length_; n++ )
                FreeSqnWeightsAt( n, st );
        }
        if( sqnweights_[st]) {
            free( sqnweights_[st]);
            sqnweights_[st] = NULL;
        }
        if( sqnweights_[st+PS_NSTATES]) {
            free( sqnweights_[st+PS_NSTATES]);
            sqnweights_[st+PS_NSTATES] = NULL;
        }
    }
}





// =========================================================================
// InitPrexpNoObservations: initialize expected number of distinct residues:
//     SUM ( 1 - ( 1 - backp[i]^n )), where n number of observations
//
void ExtendedDescriptionVector::InitPrexpNoDistinctRes( const double* backprobs )
{
    int     noelems = GetSizeOfExpNoDistinctRes();
    int     noeffres = NUMAA;

    if( noelems < 1 )
        return;

    DestroyPrexpNoDistinctRes();

    double  pterm;  //partial sum of probabilities
    int     n, r;

    prexpnores_ = ( double* )malloc( noelems * sizeof( double ));

    if( prexpnores_ == NULL )
        throw myruntime_error( "InitPrexpNoDistinctRes: Not enough memory." );

    memset( prexpnores_, 0, noelems * sizeof( double ));

    prexpnores_[0] = 0.0;

    for( n = 1; n < noelems; n++ ) {
        pterm = 0.0;
        for( r = 0; r < noeffres; r++ )
            pterm += exp( n * log( 1.0 - ( backprobs? backprobs[r]: LOSCORES.PROBABility( r )) ));
        prexpnores_[n] = noeffres - pterm;
    }
}

// DestroyPrexpNoObservations: destroy expected number of distinct residues
//
void ExtendedDescriptionVector::DestroyPrexpNoDistinctRes()
{
    if( prexpnores_ ) {
        free( prexpnores_ );
        prexpnores_ = NULL;
    }
}

// GetExpNoObservations: get expected number of observations
//     corresponding to average number of distinct observed residues
//
double ExtendedDescriptionVector::GetExpNoObservations( double avgnodistres )
{
    int     noelems = GetSizeOfExpNoDistinctRes();
    double  expnobs = 0.0;
    int     n;

    if( prexpnores_ == NULL )
        throw myruntime_error( "GetExpNoObservations: Memory access error." );

    for( n = 1; n < noelems && prexpnores_[n] <= avgnodistres; n++ );
    expnobs = ( noelems <= n )? n: n - ( prexpnores_[n] - avgnodistres )/( prexpnores_[n] - prexpnores_[n-1]);
    return expnobs;
}





// =========================================================================
// Init: reset all values
//
void ExtendedDescriptionVector::Init()
{
    int st;
    PosDescriptionVector::Init();
    //
    indices = NULL;
    indicesLengths_ = NULL;
    indicesCapacities_= NULL;
    reservation_ = ALLOCSEQ;
    //
    memset( backprobs_, 0, sizeof( double ) * NUMALPH );
    memset( postprobs_, 0, sizeof( double ) * NUMALPH );
    states = NULL;
    extents = NULL;
    nodifrs_ = NULL;
    MSextents = NULL;
    counts = NULL;
    distribution = NULL;
    for( st = 0; st < PS_NSTATES; st++ ) {
        sqnweights_[st] = NULL;
        sqnweights_[st+PS_NSTATES] = NULL;
    }
    matchWeights = NULL;
    transWeights = NULL;
    gapWeights   = NULL;
    distinctHist = NULL;
    MSdistinctHist = NULL;
    targetFreqns = NULL;
    targetTranst = NULL;
    rawPSSMatrix = NULL;
    information  = NULL;
    bppprob_ = NULL;
    ppprobs_ = NULL;
    pppndxs_ = NULL;
    noppps_ = NULL;
    noseqs_      = NULL;
    expnoseqs_   = NULL;
    expnosmid_   = NULL;
}

// -------------------------------------------------------------------------
// InitRightExtents: initialize right extent values
//
void ExtendedDescriptionVector::InitRightExtents( size_t from, size_t to )
{
    size_t  n, loclength = capacity();
    int st;

    if( extents == NULL || MSextents == NULL )
        throw myruntime_error( "ExtendedDescriptionVector: InitRightExtents: Memory access error." );
    if( to )
        loclength = to;

    for( n = from; n < loclength; n++ )
        for( st = 0; st < PS_NSTATES; st++ )
            extents[n][st][xRight] = ( size_t )SIZE_MAX;
    for( n = from; n < loclength; n++ )
        MSextents[n][xRight] = ( size_t )SIZE_MAX;
}

// -------------------------------------------------------------------------
// assignment operator
//
ExtendedDescriptionVector& ExtendedDescriptionVector::operator=( const PosDescriptionVector& sequence )
{
    throw myruntime_error("ExtendedDescriptionVector: Object assignment is not allowed.");
    return *this;
}

// -------------------------------------------------------------------------
// Realloc: reallocate memory
//
void ExtendedDescriptionVector::Realloc( int newcap )
{
    int st;
    if( newcap < capacity_ )
        return;

    if( capacity_ == 0 ) {
        indices             = ( size_t** )malloc( sizeof( void* ) * newcap );
        indicesLengths_     = ( size_t* )malloc( sizeof( size_t ) * newcap );
        indicesCapacities_  = ( size_t* )malloc( sizeof( size_t ) * newcap );
        //
        states = ( int* )malloc( sizeof( int ) * newcap );
        extents = ( size_t(*)[PS_NSTATES][xCount])malloc( sizeof( size_t ) * PS_NSTATES * xCount * newcap );
        nodifrs_ = ( double(*)[PS_NSTATES])malloc( sizeof( double ) * PS_NSTATES * newcap );
        MSextents = ( size_t(*)[xCount] )malloc( sizeof( size_t ) * xCount * newcap );
        counts  = ( size_t* )malloc( sizeof( size_t ) * newcap );

        for( st = 0; st < PS_NSTATES; st++ ) {
            sqnweights_[st] = ( double** )malloc( sizeof( double* ) * newcap );
            sqnweights_[st+PS_NSTATES] = ( double** )malloc( sizeof( double* ) * newcap );
        }

        distribution = ( size_t(*)[NUMALPH] )malloc( sizeof( size_t ) * NUMALPH * newcap );
        matchWeights = ( double(*)[NUMALPH] )malloc( sizeof( double ) * NUMALPH * newcap );
        transWeights = ( double(*)[P_NSTATES] )malloc( sizeof( double ) * P_NSTATES * ( newcap + 1 ));
        gapWeights   = ( double*            )malloc( sizeof( double ) * newcap );
        distinctHist = ( size_t(*)[PS_NSTATES][NUMALPH])malloc( sizeof(size_t)*PS_NSTATES*NUMALPH *( newcap+1 ));
        MSdistinctHist = ( size_t(*)[NUMALPH] )malloc( sizeof( size_t ) * NUMALPH * newcap );
        targetFreqns = ( double(*)[NUMALPH] )malloc( sizeof( double ) * NUMALPH * newcap );
        targetTranst = ( double(*)[P_NSTATES] )malloc( sizeof( double ) * P_NSTATES * ( newcap + 1 ));
        rawPSSMatrix = ( double(*)[NUMALPH] )malloc( sizeof( double ) * NUMALPH * newcap );
        information  = ( double*            )malloc( sizeof( double ) * newcap );
        bppprob_     = ( double*            )malloc( sizeof( double ) * newcap );
        ppprobs_     = ( double**           )malloc( sizeof( double*) * newcap );
        pppndxs_     = ( int**              )malloc( sizeof( int*) * newcap );
        noppps_      = ( size_t*            )malloc( sizeof( size_t ) * newcap );
        noseqs_      = ( size_t*            )malloc( sizeof( size_t ) * newcap );
        expnoseqs_   = ( double*            )malloc( sizeof( double ) * newcap );
        expnosmid_ = ( double(*)[PS_NSTATES])malloc( sizeof( double ) * PS_NSTATES *( newcap + 1 ));

    } else {
        indices             = ( size_t** )realloc( indices, sizeof( void* ) * newcap );
        indicesLengths_     = ( size_t* )realloc( indicesLengths_, sizeof( size_t ) * newcap );
        indicesCapacities_  = ( size_t* )realloc( indicesCapacities_, sizeof( size_t ) * newcap );
        //
        states = ( int* )realloc( states, sizeof( int ) * newcap );
        extents = ( size_t(*)[PS_NSTATES][xCount])realloc( extents, sizeof(size_t)*PS_NSTATES*xCount * newcap );
        nodifrs_ = ( double(*)[PS_NSTATES])realloc( nodifrs_, sizeof( double ) * PS_NSTATES * newcap );
        MSextents = ( size_t(*)[xCount] )realloc( MSextents, sizeof( size_t ) * xCount * newcap );
        counts  = ( size_t* )realloc( counts, sizeof( size_t ) * newcap );

        for( st = 0; st < PS_NSTATES; st++ ) {
            sqnweights_[st] = (double**)realloc( sqnweights_[st], sizeof(double*)* newcap );
            sqnweights_[st+PS_NSTATES] = (double**)realloc( sqnweights_[st+PS_NSTATES], sizeof(double*)*newcap );
        }

        distribution = ( size_t(*)[NUMALPH] )realloc( distribution, sizeof( size_t ) * NUMALPH * newcap );
        matchWeights = ( double(*)[NUMALPH] )realloc( matchWeights, sizeof( double ) * NUMALPH * newcap );
        transWeights = ( double(*)[P_NSTATES] )realloc( transWeights, sizeof( double ) * P_NSTATES * ( newcap+1));
        gapWeights   = ( double*            )realloc( gapWeights,   sizeof( double ) * newcap );
        distinctHist = ( size_t(*)[PS_NSTATES][NUMALPH])realloc( distinctHist, sizeof(size_t)*PS_NSTATES*NUMALPH *(newcap+1));
        MSdistinctHist = ( size_t(*)[NUMALPH] )realloc( MSdistinctHist, sizeof( size_t ) * NUMALPH * newcap );
        targetFreqns = ( double(*)[NUMALPH] )realloc( targetFreqns, sizeof( double ) * NUMALPH * newcap );
        targetTranst = ( double(*)[P_NSTATES] )realloc( targetTranst, sizeof( double ) * P_NSTATES * ( newcap+1));
        rawPSSMatrix = ( double(*)[NUMALPH] )realloc( rawPSSMatrix, sizeof( double ) * NUMALPH * newcap );
        information  = ( double*            )realloc( information,  sizeof( double ) * newcap );
        bppprob_     = ( double*            )realloc( bppprob_, sizeof( double ) * newcap );
        ppprobs_     = ( double**           )realloc( ppprobs_, sizeof( double*) * newcap );
        pppndxs_     = ( int**              )realloc( pppndxs_, sizeof( int*) * newcap );
        noppps_      = ( size_t*            )realloc( noppps_, sizeof( size_t ) * newcap );
        noseqs_      = ( size_t*            )realloc( noseqs_,      sizeof( size_t ) * newcap );
        expnoseqs_   = ( double*            )realloc( expnoseqs_,   sizeof( double ) * newcap );
        expnosmid_ = ( double(*)[PS_NSTATES])realloc( expnosmid_, sizeof(double)*PS_NSTATES*(newcap+1));
    }

    if( !indices || !indicesLengths_ || !indicesCapacities_ )
        throw myruntime_error( mystring( "Not enough memory." ));

    if( !states || !extents || !nodifrs_ || !MSextents || !counts )
        throw myruntime_error( mystring( "Not enough memory." ));

    if( !distribution || !matchWeights || !transWeights || !gapWeights || 
        !distinctHist || !MSdistinctHist || !targetFreqns || !targetTranst )
        throw myruntime_error( mystring( "Not enough memory." ));

    if( !bppprob_ || !ppprobs_ || !pppndxs_ || !noppps_ )
        throw myruntime_error( mystring( "Not enough memory." ));

    if( !information || !rawPSSMatrix || !noseqs_ || !expnoseqs_ || !expnosmid_ )
        throw myruntime_error( mystring( "Not enough memory." ));

    for( st = 0; st < PS_NSTATES; st++ )
        if( !sqnweights_[st] || !sqnweights_[st+PS_NSTATES])
            throw myruntime_error( mystring( "Not enough memory." ));

    size_t**            locindices = indices;
    size_t*             locindicesLengths_ = indicesLengths_;
    size_t*             locindicesCapacities_ = indicesCapacities_;
    //
    int*                locstates = states;
    size_t           ( *locextents )[PS_NSTATES][xCount] = extents;
    double           ( *locnodifrs )[PS_NSTATES] = nodifrs_;
    size_t           ( *locMSextents )[xCount] = MSextents;
    size_t*             loccounts = counts;
    double**            locsqnweights[TIMES2(PS_NSTATES)];
    size_t           ( *locdistribution )[NUMALPH] = distribution;
    double           ( *locmatchWeights )[NUMALPH] = matchWeights;
    double           ( *loctransWeights )[P_NSTATES] = transWeights + 1;
    double*             locgapWeights              = gapWeights;
    size_t           ( *locdistinctHist )[PS_NSTATES][NUMALPH] = distinctHist + 1;
    size_t           ( *locMSdistinctHist )[NUMALPH] = MSdistinctHist;
    double           ( *loctargetFreqns )[NUMALPH] = targetFreqns;
    double           ( *loctargetTranst )[P_NSTATES] = targetTranst + 1;
    double           ( *locrawPSSMatrix )[NUMALPH] = rawPSSMatrix;
    double*             locinformation             = information;
    double*             locbppprob                 = bppprob_;
    double**            locppprobs                 = ppprobs_;
    int**               locpppndxs                 = pppndxs_;
    size_t*             locnoppps                  = noppps_;
    size_t*             locnoseqs                  = noseqs_;
    double*             locexpnoseqs               = expnoseqs_;
    double           ( *locexpnosmid )[PS_NSTATES] = expnosmid_ + 1;

    for( st = 0; st < PS_NSTATES; st++ ) {
        locsqnweights[st] = sqnweights_[st];
        locsqnweights[st+PS_NSTATES] = sqnweights_[st+PS_NSTATES];
    }

    if( capacity_ != 0 ) {
        locindices += capacity_;
        locindicesLengths_ += capacity_;
        locindicesCapacities_ += capacity_;
        //
        locstates += capacity_;
        locextents += capacity_;        //!
        locnodifrs += capacity_;
        locMSextents += capacity_;
        loccounts += capacity_;

        for( st = 0; st < PS_NSTATES; st++ ) {
            locsqnweights[st] += capacity_;
            locsqnweights[st+PS_NSTATES] += capacity_;
        }
        locdistribution += capacity_;   //!
        locmatchWeights += capacity_;
        loctransWeights += capacity_;
        locgapWeights   += capacity_;
        locdistinctHist += capacity_;
        locMSdistinctHist += capacity_;
        loctargetFreqns += capacity_;
        loctargetTranst += capacity_;
        locrawPSSMatrix += capacity_;
        locinformation  += capacity_;
        locbppprob      += capacity_;
        locppprobs      += capacity_;
        locpppndxs      += capacity_;
        locnoppps       += capacity_;
        locnoseqs       += capacity_;
        locexpnoseqs    += capacity_;
        locexpnosmid    += capacity_;
    }
    else {
        //beginning transition weights and frequencies
        memset( transWeights, 0, sizeof( double ) * P_NSTATES );
        memset( targetTranst, 0, sizeof( double ) * P_NSTATES );
        memset( distinctHist, 0, sizeof( size_t ) * PS_NSTATES * NUMALPH );
        memset( expnosmid_, 0, sizeof( double ) * PS_NSTATES );
    }


    memset( locindices,             0, sizeof( size_t* ) * ( newcap - capacity_ ));
    memset( locindicesLengths_,     0, sizeof( size_t  ) * ( newcap - capacity_ ));
    memset( locindicesCapacities_,  0, sizeof( size_t  ) * ( newcap - capacity_ ));
    //
    memset( locstates,  0, sizeof( int ) * ( newcap - capacity_ ));
    memset( locextents, 0, sizeof( size_t ) * ( newcap - capacity_ ) * PS_NSTATES * xCount );
    memset( locnodifrs, 0, sizeof( double ) * ( newcap - capacity_ ) * PS_NSTATES );
    memset( locMSextents, 0, sizeof( size_t ) * ( newcap - capacity_ ) * xCount );
    memset( loccounts,  0, sizeof( size_t ) * ( newcap - capacity_ ));
    InitRightExtents( capacity_, newcap );
    // 
    for( st = 0; st < PS_NSTATES; st++ ) {
        memset( locsqnweights[st], 0, sizeof( double* ) * ( newcap - capacity_ ));
        memset( locsqnweights[st+PS_NSTATES], 0, sizeof( double* ) * ( newcap - capacity_ ));
    }
    //
    memset( locdistribution, 0, sizeof( size_t ) * ( newcap - capacity_ ) * NUMALPH );
    memset( locmatchWeights, 0, sizeof( double ) * ( newcap - capacity_ ) * NUMALPH );
    memset( loctransWeights, 0, sizeof( double ) * ( newcap - capacity_ ) * P_NSTATES );
    memset( locgapWeights,   0, sizeof( double ) * ( newcap - capacity_ ));

    memset( locdistinctHist, 0, sizeof( size_t ) * ( newcap - capacity_ ) * PS_NSTATES * NUMALPH );
    memset( locMSdistinctHist, 0, sizeof( size_t ) * ( newcap - capacity_ ) * NUMALPH );

    memset( loctargetFreqns, 0, sizeof( double ) * ( newcap - capacity_ ) * NUMALPH );
    memset( loctargetTranst, 0, sizeof( double ) * ( newcap - capacity_ ) * P_NSTATES );
    memset( locrawPSSMatrix, 0, sizeof( double ) * ( newcap - capacity_ ) * NUMALPH );
    //
    memset( locinformation,  0, sizeof( double ) * ( newcap - capacity_ ));
    memset( locbppprob,      0, sizeof( double ) * ( newcap - capacity_ ));
    memset( locppprobs,      0, sizeof( double*) * ( newcap - capacity_ ));
    memset( locpppndxs,      0, sizeof( int*) * ( newcap - capacity_ ));
    memset( locnoppps,       0, sizeof( size_t ) * ( newcap - capacity_ ));
    memset( locnoseqs,       0, sizeof( size_t ) * ( newcap - capacity_ ));
    memset( locexpnoseqs,    0, sizeof( double ) * ( newcap - capacity_ ));
    memset( locexpnosmid,    0, sizeof( double ) * ( newcap - capacity_ ) * PS_NSTATES );

    PosDescriptionVector::Realloc( newcap );
}

// -------------------------------------------------------------------------
// PushIndexAt: insert indices of sequences not participating in
//     extent at the given position
//
void ExtendedDescriptionVector::PushIndexAt( size_t value, size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif

    if( indicesCapacities_[n] < indicesLengths_[n] + 1 ) {
        ReallocIndices( indicesCapacities_[n]? indicesCapacities_[n] * 2: reservation_, n );
    }

    indices[n][ indicesLengths_[n] ] = value;

    indicesLengths_[n]++;
}

// -------------------------------------------------------------------------
// ReallocIndices: reallocate memory associated with vector of indices of
//     sequences
//
void ExtendedDescriptionVector::ReallocIndices( size_t newcap, size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif

    if( indicesCapacities_[n] == 0 ) {
        indices[n] = ( size_t* )malloc( sizeof( size_t ) * newcap );
    } else {
        indices[n] = ( size_t* )realloc( indices[n], sizeof( size_t ) * newcap );
    }

    if( !indices[n] )
        throw myruntime_error( mystring( "Not enough memory." ));

    size_t*     indarr = indices[n];

    if( indicesCapacities_[n] != 0 ) {
        indarr += indicesCapacities_[n];
    }

    memset( indarr, 0, sizeof( size_t ) * ( newcap - indicesCapacities_[n] ));

    indicesCapacities_[n] = newcap;
}

// -------------------------------------------------------------------------
// clear: clear all buffers associated with the extended description vector
//
void ExtendedDescriptionVector::clear()
{
    size_t  n;
    int     st;

    for( n = 0; n < length_; n++ )
        if( indices[n] )
            free( indices[n] );

    memset( indices,             0, sizeof( size_t* ) * capacity_ );
    memset( indicesLengths_,     0, sizeof( size_t  ) * capacity_ );
    memset( indicesCapacities_,  0, sizeof( size_t  ) * capacity_ );
    //
    memset( backprobs_, 0, sizeof( double ) * NUMALPH );
    memset( postprobs_, 0, sizeof( double ) * NUMALPH );
    memset( states,  0, sizeof( int ) * capacity_ );
    memset( extents, 0, sizeof( size_t ) * capacity_ * PS_NSTATES * xCount );
    memset( nodifrs_, 0, sizeof( double ) * capacity_ * PS_NSTATES );
    memset( MSextents, 0, sizeof( size_t ) * capacity_ * xCount );
    memset( counts,  0, sizeof( size_t ) * capacity_ );
    InitRightExtents();

    if( capacity_ ) {
        for( st = 0; st < PS_NSTATES; st++ ) {
            for( n = 0; n < length_; n++ )
                FreeSqnWeightsAt( n, st );
            memset( sqnweights_[st], 0, sizeof( double* ) * capacity_ );
            memset( sqnweights_[st+PS_NSTATES], 0, sizeof( double* ) * capacity_ );
        }
        for( n = 0; n < length_; n++ ) {
            if( ppprobs_ && ppprobs_[n]) free( ppprobs_[n]);
            if( pppndxs_ && pppndxs_[n]) free( pppndxs_[n]);
        }
        memset( distribution, 0, sizeof( size_t ) * capacity_ * NUMALPH );
        memset( matchWeights, 0, sizeof( double ) * capacity_ * NUMALPH );
        memset( transWeights, 0, sizeof( double ) * ( capacity_+1 ) * P_NSTATES );
        memset( gapWeights,   0, sizeof( double ) * capacity_           );
        memset( distinctHist, 0, sizeof( size_t ) * ( capacity_+1 )* PS_NSTATES*NUMALPH );
        memset( MSdistinctHist, 0, sizeof( size_t ) * capacity_ * NUMALPH );
        memset( targetFreqns, 0, sizeof( double ) * capacity_ * NUMALPH );
        memset( targetTranst, 0, sizeof( double ) * ( capacity_+1 ) * P_NSTATES );
        memset( rawPSSMatrix, 0, sizeof( double ) * capacity_ * NUMALPH );
        memset( information,  0, sizeof( double ) * capacity_           );
        memset( bppprob_,     0, sizeof( double ) * capacity_           );
        memset( ppprobs_,     0, sizeof( double*) * capacity_           );
        memset( pppndxs_,     0, sizeof( int*) * capacity_              );
        memset( noppps_,      0, sizeof( size_t ) * capacity_           );
        memset( noseqs_,      0, sizeof( size_t ) * capacity_           );
        memset( expnoseqs_,   0, sizeof( double ) * capacity_           );
        memset( expnosmid_,   0, sizeof( double ) * ( capacity_+1 ) * PS_NSTATES );
    }

    PosDescriptionVector::clear();
}

// -------------------------------------------------------------------------
// PrintMatchWeights: outputs match weights which correspond to the observed
//     weighted frequencies
// -------------------------------------------------------------------------

void ExtendedDescriptionVector::PrintMatchWeights( FILE* fp )
{
    size_t  l = 0;

    if( fp == NULL )
        return;

    fprintf( fp, "%9c", 32 );

    for( unsigned char r = 0; r < NUMALPH; r++ )
        fprintf( fp, "%3c%c", 32, DehashCode( r ) );

    for( size_t p = 0; p < size(); p++ ) {
        // omit unused positions and gaps in query
        if( /* !IsUsedAt( p ) || */ResidueAt( p ) == GAP )
            continue;

        fprintf( fp, "\n%5d %c   ", ++l, DehashCode( ResidueAt( p )));

        for( unsigned char r = 0; r < NUMALPH; r++ )
            fprintf( fp, "%4d", ( int )rint( 100 * GetMatchWeightsAt( r, p )));
    }
    fprintf( fp, "\n" );
}

// -------------------------------------------------------------------------
// PrintTransWeights: output weighted transition weights
//
void ExtendedDescriptionVector::PrintTransWeights( FILE* fp )
{
    int     n;
    size_t  l = 0;

    if( fp == NULL )
        return;

    fprintf( fp, "%9c", 32 );

    for( n = 0; n < P_NSTATES; n++ )
        fprintf( fp, " %3s", gTPTRANS_NAMES[n]);

    fprintf( fp, "\n%5d %c   ", l, 32 );

    for( n = 0; n < P_NSTATES; n++ )
        fprintf( fp, "%4d", ( int )rint( 100 * GetTransWeightsBeg( n )));

    for( size_t p = 0; p < size(); p++ ) {
        // omit unused positions and gaps in query
        if( /* !IsUsedAt( p ) || */ResidueAt( p ) == GAP )
            continue;

        fprintf( fp, "\n%5d %c   ", ++l, DehashCode( ResidueAt( p )));

        for( n = 0; n < P_NSTATES; n++ )
            fprintf( fp, "%4d", ( int )rint( 100 * GetTransWeightsAt( n, p )));
    }
    fprintf( fp, "\n" );
}

// -------------------------------------------------------------------------
// PrintTransWeights: output target transition weights
//
void ExtendedDescriptionVector::PrintTargetTranst( FILE* fp )
{
    int     n;
    size_t  l = 0;

    if( fp == NULL )
        return;

    fprintf( fp, "%9c", 32 );

    for( n = 0; n < P_NSTATES; n++ )
        fprintf( fp, " %3s", gTPTRANS_NAMES[n]);

    fprintf( fp, "\n%5d %c   ", l, 32 );

    for( n = 0; n < P_NSTATES; n++ )
        fprintf( fp, "%4d", ( int )rint( 100 * GetTargetTranstBeg( n )));

    for( size_t p = 0; p < size(); p++ ) {
        // omit unused positions and gaps in query
        if( /* !IsUsedAt( p ) || */ResidueAt( p ) == GAP )
            continue;

        fprintf( fp, "\n%5d %c   ", ++l, DehashCode( ResidueAt( p )));

        for( n = 0; n < P_NSTATES; n++ )
            fprintf( fp, "%4d", ( int )rint( 100 * GetTargetTranstAt( n, p )));
    }
    fprintf( fp, "\n" );
}

// -------------------------------------------------------------------------
// PrintMIDExpNoObservations: output weighted MID state observations
//
void ExtendedDescriptionVector::PrintMIDExpNoObservations( FILE* fp )
{
    int     n;
    size_t  l = 0;

    if( fp == NULL )
        return;

    fprintf( fp, "%9c", 32 );

    for( n = 0; n < PS_NSTATES; n++ )
        fprintf( fp, " %11s", gTPSTATES_NAMES[n]);

    fprintf( fp, "\n%5d %c   ", l, 32 );

    for( n = 0; n < PS_NSTATES; n++ )
        fprintf( fp, "%12g", GetMIDExpNoObservationsAt( -1, n ));

    for( size_t p = 0; p < size(); p++ ) {
        // omit unused positions and gaps in query
        if( /* !IsUsedAt( p ) || */ResidueAt( p ) == GAP )
            continue;

        fprintf( fp, "\n%5d %c   ", ++l, DehashCode( ResidueAt( p )));

        for( n = 0; n < PS_NSTATES; n++ )
            fprintf( fp, "%12g", GetMIDExpNoObservationsAt( p, n ));
    }
    fprintf( fp, "\n" );
}

// -------------------------------------------------------------------------
// PrintPSSMatrix: outputs raw PSSM matrix
// -------------------------------------------------------------------------

void ExtendedDescriptionVector::PrintPSSMatrix( FILE* fp )
{
    size_t  l = 0;

    if( fp == NULL )
        return;

    fprintf( fp, "%9c", 32 );

    for( unsigned char r = 0; r < NUMALPH; r++ )
        fprintf( fp, "%3c", DehashCode( r ) );

    for( size_t p = 0; p < size(); p++ ) {
        // omit unused positions and gaps in query
        if( /* !IsUsedAt( p ) || */ResidueAt( p ) == GAP )
            continue;

        fprintf( fp, "\n%5d %c  ", ++l, DehashCode( ResidueAt( p )));

        for( unsigned char r = 0; r < NUMALPH; r++ )
            fprintf( fp, "%3d", ( int )rint( GetPSSMEntryAt( r, p )));
//             fprintf( fp, "%f ", GetTargetFreqnsAt( r, p ));
    }
    fprintf( fp, "\n" );
}

// -------------------------------------------------------------------------
// PrintPSSMandWeights: outputs PSSM matrix together with match weights
// -------------------------------------------------------------------------

void ExtendedDescriptionVector::PrintSuppressedPSSMandWeights( FILE* fp )
{
    int             n;
    size_t          p, l = 0;
    const size_t    effective_nr = NUMAA;
    unsigned char   r = 0;

    if( fp == NULL )
        return;

    fprintf( fp,"%22c Position-specific scoring matrix "
                "%39c Weighted observed frequencies "
                "%18c Target transition frequencies %7c Information\n", 32, 32, 32, 32 );

    fprintf( fp, "%9c", 32 );

    for( r = 0; r < effective_nr; r++ )
        fprintf( fp, "%3c", DehashCode( r ) );

    for( r = 0; r < effective_nr; r++ )
        fprintf( fp, "%4c", DehashCode( r ) );

    for( n = 0; n < P_NSTATES; n++ )
        fprintf( fp, "%4s", gTPTRANS_NAMES[n]);

    for( p = 0; p < size(); p++ ) {
        // omit unused positions and gaps in query
        if( /* !IsUsedAt( p ) || */ResidueAt( p ) == GAP )
            continue;

        fprintf( fp, "\n%5d %c   ", ++l, DehashCode( ResidueAt( p )));

        for( r = 0; r < effective_nr; r++ )
            fprintf( fp, "%2d ", ( int )rint( GetPSSMEntryAt( r, p )));

        for( r = 0; r < effective_nr; r++ )
            fprintf( fp, "%4d", ( int )rint( 100 * GetMatchWeightsAt( r, p )));

        for( n = 0; n < P_NSTATES; n++ )
            fprintf( fp, "%4d", ( int )rint( 100 * GetTargetTranstAt( n, p )));

        fprintf( fp, " %5.2f", GetInformationAt( p ));
    }
    fprintf( fp, "\n" );
}

// -------------------------------------------------------------------------
// PrintProfile: outputs all needed information about the profile
// -------------------------------------------------------------------------

void ExtendedDescriptionVector::PrintProfile( FILE* fp )
{
    int             n;
    size_t          p, l = 0;
    const size_t    res_count = NUMALPH - 1;
    unsigned char   r = 0;

    if( fp == NULL )
        return;

    fprintf( fp,"%22c Position-specific scoring matrix "
                "%39c Weighted observed frequencies "
                "%18c Target transition frequencies %7c Gap weights %c Information\n", 32, 32, 32, 32, 32 );

    fprintf( fp, "%9c", 32 );

    for( r = 0; r < res_count; r++ )
        fprintf( fp, "%3c", DehashCode( r ) );

    for( r = 0; r < res_count; r++ )
        fprintf( fp, "%4c", DehashCode( r ) );

    for( n = 0; n < P_NSTATES; n++ )
        fprintf( fp, "%4s", gTPTRANS_NAMES[n]);

    for( p = 0; p < size(); p++ ) {
        // omit unused positions and gaps in query
        if( /* !IsUsedAt( p ) || */ResidueAt( p ) == GAP )
            continue;

        fprintf( fp, "\n%5d %c   ", ++l, DehashCode( ResidueAt( p )));

        for( r = 0; r < res_count; r++ )
            fprintf( fp, "%2d ", ( int )rint( GetPSSMEntryAt( r, p )));

        for( r = 0; r < res_count; r++ )
            fprintf( fp, "%4d", ( int )rint( 100 * GetMatchWeightsAt( r, p )));

        for( n = 0; n < P_NSTATES; n++ )
            fprintf( fp, "%4d", ( int )rint( 100 * GetTargetTranstAt( n, p )));

        fprintf( fp, " %6d", ( int )rint( 100 * GetGapWeightsAt( p )));

        fprintf( fp, " %13.2f", GetInformationAt( p ));
    }
    fprintf( fp, "\n" );
}

