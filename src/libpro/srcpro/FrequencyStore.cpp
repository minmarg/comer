/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "mystring.h"
#include "myexcept.h"
#include "Serializer.h"
#include "FrequencyStore.h"


TFVectorProbabilities   FrequencyStore::probtype = DISCRETE;
THashFunction FrequencyStore::func1stlevel = FrequencyVector::RJHashing;
THashFunction FrequencyStore::func2ndlevel = FrequencyVector::SBoxHashing;//;

THashFunction FrequencyVector::hash_FUNCTIONS[] = {
    RJHashing,
    SBoxHashing,
    OATHashing,
    MD5Hashing,
    CRCHashing
};

const size_t FrequencyVector::no_hash_FUNCTIONS = sizeof( hash_FUNCTIONS ) / sizeof( THashFunction );

// -------------------------------------------------------------------------
// constructor: allocates memory required by vector and initializes it
// -------------------------------------------------------------------------

FrequencyVector::FrequencyVector(
    const double    freq[NUMALPH],
    double          weight,
    size_t          thickn,
    const double    scores[NUMALPH],
    double          info )
:
    vector( NULL )
{
    Init( freq, weight, thickn, scores, info );
}

// -------------------------------------------------------------------------
// constructor: this constructor for reading purposes; no memory is
//     allocated within
// -------------------------------------------------------------------------

FrequencyVector::FrequencyVector( const char* vect )
:   vector( const_cast<char*>( vect ))
{
}

// -------------------------------------------------------------------------
// constructor: default
// -------------------------------------------------------------------------

FrequencyVector::FrequencyVector()
:   vector( NULL )
{
//     throw myruntime_error(
//         mystring( "FrequencyVector: Default construction is not allowed." ));
}

// -------------------------------------------------------------------------
// destructor: it does nothing for one; destruction should be done
//     MANUALLY!! by calling method Destroy()
// -------------------------------------------------------------------------

FrequencyVector::~FrequencyVector()
{
}

// -------------------------------------------------------------------------
// Init: initialization without copying of values
// -------------------------------------------------------------------------

void FrequencyVector::Init()
{
    if( vector == NULL ) {
        vector = ( char* )malloc( Overall );
        if( vector == NULL )
            throw myruntime_error(
                mystring( "FrequencyVector: Not enough memory." ));
    }

    memset( vector, 0, Overall );
    *( unsigned int* )vector = constant;
}

// -------------------------------------------------------------------------
// Init: initialization and copy of values are performed
// -------------------------------------------------------------------------

void FrequencyVector::Init(
    const double    freq[NUMALPH],
    double          weight,
    size_t          thickn,
    const double    scores[NUMALPH],
    double          info )
{
    Init();

    size_t          thickness = thickn;
    unsigned int    infcontent = ( unsigned int )rint( info * INFO_SCALE_CONSTANT );

    if( MAXINFVALUE < infcontent )
        infcontent = MAXINFVALUE;

    if( MAXNOSVALUE < thickness )
        thickness = MAXNOSVALUE;

    SetWeight(( WGHTTYPE )rint( weight ));
    SetThickness(( NOSTYPE ) thickness );
    SetInfContent(( INFMTYPE ) infcontent );

    for( size_t n = 0; n < no_Elems; n++ ) {
#ifdef __DEBUG__
        if( 1.0 < freq[n] )
            throw myruntime_error( mystring( "FrequencyVector: Values of frequencies corrupted." ));
#endif
        SetValueAt( n, ( FREQTYPE )/*( int )*/rint( FREQUENCY_SUM * freq[n] ));
        SetScoreAt( n, ( short )rint( scores[n] * SCALE_CONSTANT ));
    }

#if 0   //if use interference of adjacent bytes
    size_t  k = 3, r, l;

    if( no_Elems < k )
        k = no_Elems;

    size_t  left = k >> 1;
    size_t  rght = k - left;
    size_t  sum;

    for( size_t n = 0; n < no_Elems; n++ ) {
        l = ( left < n )? n - left: 0;
        r = ( n + rght < no_Elems )? n + rght: no_Elems;
        sum = 0;

        for( ; l < r; l++ )
            sum += GetValueAt( l );

        SetValueAt( n, sum );
    }
#endif
}

// -------------------------------------------------------------------------
// SetVector: changes contents of the vector with the values of frequencies,
//     scores specified
// -------------------------------------------------------------------------

void FrequencyVector::SetVector(
    const double    freq[NUMALPH],
    double          weight,
    size_t          thickn,
    const double    scores[NUMALPH],
    double          info )
{
    Init( freq, weight, thickn, scores, info );
}

// -------------------------------------------------------------------------
// Destroy: destroys allocated memory; object destruction must go through
//     this method
// -------------------------------------------------------------------------

void FrequencyVector::Destroy()
{
    if( vector ) {
        free( vector );
        vector = NULL;
    }
}

// -------------------------------------------------------------------------
// ComputeProVectorProbability: compute probability of occurence of the
//     frequency vector; frequencies must be written in the vector before
//     calling this method; probability is computed according to profile
//     vector theory
//
//    __    (Fja)
//   |  | Pa
//    a
//
//  where Pa are background probabilities, Fja are observed frequencies from
//  vector j
// -------------------------------------------------------------------------

double FrequencyVector::ComputeProVectorProbability()
{
    double      probv = 0.0;
    double      freqv;

#ifdef __DEBUG__
    if( NUMALPH < no_Elems )
        throw myruntime_error( mystring( "FrequencyVector: Wrong size of vector." ));
#endif

    for( int r = 0; r < no_Elems; r++ ) {
        if( LOSCORES.PROBABility( r ) <= 0.0 )
            //consider only residues which have background probabilities
            continue;

        freqv  = ( double )( unsigned )GetValueAt( r );
        probv += freqv / FREQUENCY_SUM * LOSCORES.LogPROBABility( r );
    }

#ifdef __DEBUG__
    if( 0.0 < probv )
        throw myruntime_error( mystring( "FrequencyVector: Error in computing probability." ));
#endif

    probv = exp( probv );
    return probv;
}

// -------------------------------------------------------------------------
// ComputeMultinomialProbability: compute probability of occurence of the
//     frequency vector; frequencies must be written in the vector before
//     calling this method; probability is computed using MULTINOMIAL
//     distribution:
//       _
//    ( \  Fja ) !
//      /_ a       __    (Fja)
//    ----------- |  | Pa
//      __         a
//     |  | Fja!
//      a
//
//  where Pa are background probabilities, Fja are observed frequencies from
//  vector j
// NOTE: such scheme of computing probabilities is not very suitable
//     because of differences between vectors of different compositions; for
//     example vectors of uniformly distributed frequencies dominate
//     against the others
// -------------------------------------------------------------------------

double FrequencyVector::ComputeMultinomialProbability()
{
    double      probv = 0.0;
    unsigned    freqv;

#ifdef __DEBUG__
    if( NUMALPH < no_Elems )
        throw myruntime_error( mystring( "FrequencyVector: Wrong size of vector." ));
#endif

    for( int r = 0; r < no_Elems; r++ ) {
        if( LOSCORES.PROBABility( r ) <= 0.0 )
            //consider only residues which have background probabilities
            continue;

        freqv  = ( unsigned )GetValueAt( r );
        probv += freqv * LOSCORES.LogPROBABility( r ) - LOG_FREQUENCIES.SumOf( freqv );
    }

    probv += LOG_FREQUENCIES.Total();//log gamma(sum freqv)

#ifdef __DEBUG__
    if( 0.0 < probv )
        throw myruntime_error( mystring( "FrequencyVector: Error in computing probability." ));
#endif

    probv = exp( probv );
    return probv;
}

// -------------------------------------------------------------------------
// ComputeProbability: compute probability of occurence of the frequency
//     vector; frequencies must be written in the vector before calling this
//     method;
// -------------------------------------------------------------------------

double FrequencyVector::ComputeProbability()
{
    double      probv = 0.0;

    switch( FrequencyStore::GetDistributionType()) {
        case MULTINOMIAL:
            probv = ComputeMultinomialProbability();
            break;

        case PROVECTOR:
            probv = ComputeProVectorProbability();
            break;

        case DISCRETE:
            //do nothing
            return probv;

        default:
            throw myruntime_error( mystring( "FrequencyVector: Unknown vector distribution type." ));
    }

    SetProbability( probv );
    return probv;
}

// -------------------------------------------------------------------------
// NormalizeProbability: normalize probability
// -------------------------------------------------------------------------

void FrequencyVector::NormalizeProbability( double probnorm )
{
    if( FrequencyStore::GetDistributionType() != MULTINOMIAL &&
        FrequencyStore::GetDistributionType() != PROVECTOR )
        return;

#ifdef __DEBUG__
    if( probnorm < 0.0 )
        throw myruntime_error( mystring( "FrequencyVector: Wrong normalization term." ));
#endif

    double  probv = GetProbability();
    if( probnorm )
        SetProbability( probv / probnorm );
}

// -------------------------------------------------------------------------
// UpdateProbability: increase probability for the frequency vector to
//     occur; this method is an alternative to the ComputeProbability method
// -------------------------------------------------------------------------

void FrequencyVector::UpdateProbability()
{
    //the method is for uniform distribution only
    if( FrequencyStore::GetDistributionType() != DISCRETE )
        return;

    double  probv = GetProbability() + 1.0;
    SetProbability( probv );
}

// -------------------------------------------------------------------------
// FinalProbability: compute finally probability given total number of
//     frequency vectors; this method is to be used as an alternative to the
//     multinomial computations
// -------------------------------------------------------------------------

void FrequencyVector::FinalProbability( size_t total )
{
    //the method is for uniform distribution only
    if( FrequencyStore::GetDistributionType() != DISCRETE )
        return;

#ifdef __DEBUG__
    if( total == 0 )
        throw myruntime_error( mystring( "FrequencyVector: Illegal operation." ));
#endif

    double  probv = GetProbability() + 1.0;
    SetProbability( probv / total );
}



// -------------------------------------------------------------------------
// Serialize: writes frequencies for effective number of residues to file
// -------------------------------------------------------------------------

void FrequencyVector::Serialize( Serializer& serializer ) const
{
    if( !vector )
        return;

    //omit preamble, write alternate information
    serializer.Write( vector + preamble, sizeof( char ), Overall - preamble );
}

// -------------------------------------------------------------------------
// Deserialize: read data into the class member, vector of frequencies
// -------------------------------------------------------------------------

void FrequencyVector::Deserialize( Serializer& serializer )
{
    if( vector )
        //vector supposed to be not initialized yet
        throw myruntime_error( mystring( "FrequencyVector: Deserialize: Vector is non-empty." ));

    Init();
    serializer.Read( vector + preamble, sizeof( char ), Overall - preamble );
}



// -------------------------------------------------------------------------

const int FrequencyVector_Prob_scaling = 65536;

// -------------------------------------------------------------------------
// TextWriteVector: writes full vector data to file
//
void FrequencyVector::TextWriteVector( FILE* fp ) const
{
    if( !vector || !fp )
        return;

    const size_t    locsize = KBYTE;
    char            locbuffer[locsize+1] = {0};
    char*           p = NULL;
    size_t          n;

    //preamble omitted
    p = locbuffer;
    for( n = 0; n < no_Elems; n++ ) {
        sprintf( p, "%d ", GetValueAt( n ));
        p += strlen( p );
    }
    sprintf( p, "%d ", GetWeight()); p += strlen( p );
    sprintf( p, "%d ", GetThickness()); p += strlen( p );
    for( n = 0; n < no_Elems; n++ ) {
        sprintf( p, "%d ", GetScoreAt( n ));
        p += strlen( p );
    }
    sprintf( p, "%d ", GetInfContent()); p += strlen( p );
    sprintf( p, "%.8g", GetProbability());
//     sprintf( p, "%d", ( int )rint( FrequencyVector_Prob_scaling * GetProbability()));
    fprintf( fp, "%s\n", locbuffer );
}

// -------------------------------------------------------------------------
// TextReadVector: read vector data from file
//
void FrequencyVector::TextReadVector( FILE* fp )
{
    if( !fp )
        return;
    if( vector )
        //vector supposed to be not initialized yet
        throw myruntime_error( mystring( "FrequencyVector: TextReadVector: Vector is non-empty." ));

    size_t          length, rbts;
    const size_t    locsize = KBYTE;
    char            locbuffer[locsize+1] = {0};
    char*           p;
    size_t          n;
    int             emsg;

    double          prob;
    int             intval;
//     int             freqs[no_Elems];
//     int             scores[no_Elems];
//     int             weight, thickn, info;

    Init();

    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong format of data vectors." );

    p = locbuffer;
    //frequencies
    for( n = 0; n < no_Elems; n++ ) {
        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong format of data vectors: No frequencies." );

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( FREQUENCY_SUM < intval || intval < 0 )
            throw myruntime_error( "Wrong format of data vectors: Invalid frequencies." );

        SetValueAt( n,( FREQTYPE )intval );
        p += rbts;
    }

    //weight
    if( length <= size_t( p - locbuffer ))
        throw myruntime_error( "Wrong format of data vectors: No frequency weight." );

    if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if(( int )( unsigned WGHTTYPE )( -1 ) < intval || intval < 0 )
        throw myruntime_error( "Wrong format of data vectors: Invalid frequency weight." );

    SetWeight(( WGHTTYPE )intval );
    p += rbts;

    //thickness
    if( length <= size_t( p - locbuffer ))
        throw myruntime_error( "Wrong format of data vectors: No eff. thickness value." );

    if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( MAXNOSVALUE < intval || intval < 1 )
        throw myruntime_error( "Wrong format of data vectors: Invalid eff. thickness value." );

    SetThickness(( NOSTYPE )intval );
    p += rbts;

    //scores
    for( n = 0; n < no_Elems; n++ ) {
        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong format of data vectors: No scores." );

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        SetScoreAt( n, ( short )intval );
        p += rbts;
    }

    //information
    if( length <= size_t( p - locbuffer ))
        throw myruntime_error( "Wrong format of data vectors: No information content." );

    if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( MAXINFVALUE < intval || intval < 0 )
        throw myruntime_error( "Wrong format of data vectors: Invalid information content value." );

    SetInfContent(( INFMTYPE )intval );
    p += rbts;

    //log-probability
    if( length <= size_t( p - locbuffer ))
        throw myruntime_error( "Wrong format of data vectors: No log probability." );

    if(( emsg = read_double( p, length - size_t( p - locbuffer ), &prob, &rbts )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

//     prob = ( double )intval / FrequencyVector_Prob_scaling;
    SetProbability( prob );
    p += rbts;
}



// -------------------------------------------------------------------------
// Print: print frequencies in human readable format
// -------------------------------------------------------------------------

void FrequencyVector::Print( FILE* fp ) const
{
    if( !fp )
        return;

    fprintf( fp, "freq: " );

    for( size_t n = 0; n < no_Elems; n++ )
        fprintf( fp, "%4d", ( int )GetValueAt( n ));

    fprintf( fp, ", alpha: %2d, thickn:%4d, scores: ", ( int )GetWeight(), ( int )GetThickness());

    for( size_t n = 0; n < no_Elems; n++ )
        fprintf( fp, "%3d", ( int )rint(( double ) GetScoreAt( n ) / SCALE_CONSTANT ));

    fprintf( fp, ", info: %.2f, prob: %.4f", ( double ) GetInfContent() / INFO_SCALE_CONSTANT, GetProbability());

    fprintf( fp, "\n" );
}

// -------------------------------------------------------------------------
// STATIC members
//
// FreqENOSComparer: compare effective number of sequences of vectors
//
int FrequencyVector::FreqENOSComparer( const void* vect1, const void* vect2 )
{
    if( vect1 == NULL || vect2 == NULL )
        throw myruntime_error("FrequencyVector:FreqENOSComparer: Null vectors.");
    NOSTYPE enos1 = *( NOSTYPE* )(( const char* )vect1 + preamble + vectorSz + alphaSz );
    NOSTYPE enos2 = *( NOSTYPE* )(( const char* )vect2 + preamble + vectorSz + alphaSz );
    return enos1 - enos2;
}

// -------------------------------------------------------------------------
// FrequencyComparer: Comparison function for two frequency vectors
//
int FrequencyVector::FrequencyComparer( const void* vect1, const void* vect2 )
{
    return memcmp( ExtractKey( vect1 ), ExtractKey( vect2 ), SizeofKey());
}

// -------------------------------------------------------------------------
// FrequencyValueComparer: compare values of two frequency vectors
//
int FrequencyVector::FrequencyValueComparer( const void* vect1, const void* vect2 )
{
    return memcmp( ExtractValue( vect1 ), ExtractValue( vect2 ), SizeofValue());
}

// -------------------------------------------------------------------------
// RJHashing: Jenkins' hashing
// -------------------------------------------------------------------------

size_t FrequencyVector::RJHashing( const void* vect )
{
    static  Uint32      initval = 0;//zero works best
    size_t  address = hashlittle( ExtractKey( vect ), SizeofKey(), initval );
    return  address;
}

// -------------------------------------------------------------------------
// SBoxHashing: S-Box hashing
// -------------------------------------------------------------------------

size_t FrequencyVector::SBoxHashing( const void* vect )
{
    size_t  address = sboxhash( ExtractKey( vect ), SizeofKey());
    return  address;
}

// -------------------------------------------------------------------------
// OATHashing: one at a time hashing
// -------------------------------------------------------------------------

size_t FrequencyVector::OATHashing( const void* vect )
{
    size_t  address = oneatatime( ExtractKey( vect ), SizeofKey());
    return  address;
}

// -------------------------------------------------------------------------
// MD5Hashing: MD5 cryptographic hashing
// -------------------------------------------------------------------------

size_t FrequencyVector::MD5Hashing( const void* vect )
{
    Uint32  cipher[SZCIPHER];  //the algorithm produces a 128-bit cipher
    md5hashing( cipher, ExtractKey( vect ), SizeofKey());

    size_t  which = 0;
#if 1   //if use selection of double word
    size_t  max = 0;
    size_t  rgh = SizeofKey() / SZCIPHER;
    size_t  sum = 0;
    for( size_t n = 0, cnt = 1, part = 0; n < SizeofKey(); n++, cnt++ )
        if( cnt == rgh ) {
            if( max < sum ) { max = sum; which = part; }
            sum = 0;
            cnt = 1;
            part++;
        }
#endif
    return cipher[which];//return the first byte of the obtained cipher
}

// -------------------------------------------------------------------------
// CRCHashing: hashing by CRC
// -------------------------------------------------------------------------

size_t FrequencyVector::CRCHashing( const void* vect )
{
    size_t  address = crchash( ExtractKey( vect ), SizeofKey());
    return  address;
}

// -------------------------------------------------------------------------
// MixHashing: double hashing algorithm
// -------------------------------------------------------------------------

size_t FrequencyVector::MixHashing( const void* vect )
{
    size_t  code1 = crchash( ExtractKey( vect ), SizeofKey());
//     size_t  highest = (( code1 >> 27 ) & 1 ) | (( code1 >> 29 ) & 2 ) | (( code1 >> 31 ) & 4 );
    size_t  highest = ( code1 >> 22 ) & 7;
    size_t  address = ( *hash_FUNCTIONS[ highest % no_hash_FUNCTIONS ] )( vect );
    return  address;
}





// /////////////////////////////////////////////////////////////////////////
// CLASS FrequencyStore
//
// constructor: default
//
FrequencyStore::FrequencyStore( int fstype )
:   primaryHash( NULL ),
    frequencies( NULL ),
    no_collisions1( 0 ),
    no_collisions2( 0 ),
    freqstructtype_( fstype )
{
    Init();
}

// -------------------------------------------------------------------------
// destructor:
// -------------------------------------------------------------------------

FrequencyStore::~FrequencyStore()
{
    Destroy();
    DestroyFrequencies();
}

// -------------------------------------------------------------------------
// Init: initialization
// -------------------------------------------------------------------------

void FrequencyStore::Init()
{
    //create primary hash strucutre using Jenkins' hash function;
    //hash size is defined as 2 to FirstLevelHashBits
    primaryHash = new HashTable( FirstLevelHashing, hashsize( FirstLevelHashBits ));

    if( primaryHash == NULL )
        throw myruntime_error(
            mystring( "FrequencyStore: Not enough memory." ));
}

// -------------------------------------------------------------------------
// SetNoFrequencyVectors: resets number that affects size of vector
//     frequencies; all information in frequencies will be destroyed
// -------------------------------------------------------------------------

void FrequencyStore::SetNoFrequencyVectors( size_t newsize )
{
    DestroyFrequencies();
    NewFrequencies( newsize );
}

// -------------------------------------------------------------------------
// Destroy: destroys all elements contained in the two-level hash system
// -------------------------------------------------------------------------

void FrequencyStore::Destroy()
{
    if( !primaryHash )
        return;

    for( size_t n = 0; n < primaryHash->GetHashSize(); n++ )
    {   //iterate over all hash table...
        const void* value = primaryHash->GetValueAt( n );
        if( !value )
            continue;
        FrequencyVector frequencies(( const char* )value );

        if( frequencies.IsVector())
            frequencies.Destroy();
        else
            DestroySecondaryHash( value );
    }

    delete primaryHash;
}

// -------------------------------------------------------------------------
// DestroySecondaryHash: destroys secondary hash table at the address
//     specified
// -------------------------------------------------------------------------

void FrequencyStore::DestroySecondaryHash( const void* value )
{
    const HashTable*    hash = ( const HashTable* )value;

    if( !hash )
        return;

    for( size_t n = 0; n < hash->GetHashSize(); n++ )
    {   //iterate over all hash table...
        const void* value = hash->GetValueAt( n );
        if( !value )
            continue;

        FrequencyVector frequencies(( const char* )value );

        if( frequencies.IsVector())
            frequencies.Destroy();
        else
            DestroyBSearchStructure( value );
    }

    delete hash;
}

// -------------------------------------------------------------------------
// DestroyBSearchStructure: destroys binary search structure at the address
//     specified
// -------------------------------------------------------------------------

void FrequencyStore::DestroyBSearchStructure( const void* value )
{
    const BinarySearchStructure*
        bins = ( const BinarySearchStructure* )value;

    if( !bins )
        return;

    for( size_t n = 0; n < bins->GetSize(); n++ )
    {   //iterate over all hash table...
        const void* value = bins->GetValueAt( n );
        if( !value )
            continue;

        FrequencyVector frequencies(( const char* )value );

        if( frequencies.IsVector())
            frequencies.Destroy();
#ifdef __DEBUG__
        else
            throw myruntime_error( mystring( "FrequencyStore: Internal error." ));
#endif
    }

    delete bins;
}

// -------------------------------------------------------------------------
// Store: stores frequency vector (key) into the hash system; at first the
//     key is probed to put at the computed address in the primary hash
//     table; if the position is occupied, secondary hash table is created
//     if not exists and the occupation along with the key are put
//     in the hash. If secondary hash table exists, then binary search
//     structure is created and the occupation along with the key are put in
//     there;
//  returns address of the vector to be stored if frequencies have been
//  stored in the system, otherwise returns address of the duplicate if
//  there is one;
//  returns NULL on error
// -------------------------------------------------------------------------

const void* FrequencyStore::Store( const FrequencyVector& freq )
{
    if( !primaryHash )
        return NULL;

    size_t      address = primaryHash->GetAddress( freq.GetVector());
    const void* value = primaryHash->GetValueAt( address );


    if( value == NULL ) {
        primaryHash->SetValueAt( address, freq.GetVector());
        PushInFrequencies( freq.GetVector());
        return freq.GetVector();
    }

    HashTable*              secondHash = ( HashTable* )value;
    const FrequencyVector   freq_quest (( const char* )value );

    if( freq_quest.IsVector()) {
        //check if it is not a duplicate first
        if( VectorsAreEqual( freq, freq_quest ))
            //it is a duplicate, omit further processing
            return freq_quest.GetVector();

        //if there is only one frequency vector at the moment,
        //  we need to create a secondary hash table next
        secondHash = new HashTable( SecndLevelHashing, hashsize( SecndLevelHashBits ));

        if( secondHash == NULL )
            throw myruntime_error( mystring( "FrequencyStore: Not enough memory." ));

        //put the existing key in the secondary hash first
        if( ! StoreInSecondaryHash( secondHash, freq_quest ))
            throw myruntime_error( mystring( "FrequencyStore: Store failed." ));

        //now this position of the primary hash will point to the secondary hash
        primaryHash->SetValueAt( address, secondHash );
    }

    //put the original key in the secondary hash
    const void* inserted = StoreInSecondaryHash( secondHash, freq );

    //if the vector has been inserted
    if( inserted == freq.GetVector()) {
        no_collisions1++;
        PushInFrequencies( freq.GetVector());
    }

    return  inserted;
}

// -------------------------------------------------------------------------
// StoreInSecondaryHashAt: stores frequency vector (key) into the secondary
//     hash table
//  returns address of the vector to be stored if frequencies have been
//  stored in the hash, otherwise returns address of the duplicate if
//  there is one;
//  returns NULL on error
// -------------------------------------------------------------------------

const void* FrequencyStore::StoreInSecondaryHash( HashTable* hash, const FrequencyVector& freq )
{
    if( !hash )
        return NULL;

    //compute address for the key
    size_t  address = hash->GetAddress( freq.GetVector());
    const void* value = hash->GetValueAt( address );

    if( value == NULL ) {
        hash->SetValueAt( address, freq.GetVector());
        return freq.GetVector();
    }

    BinarySearchStructure*  bins = ( BinarySearchStructure* )value;
    const FrequencyVector   freq_quest (( const char* )value );

    if( freq_quest.IsVector()) {
        //check if it is not a duplicate first
        if( VectorsAreEqual( freq, freq_quest ))
            //it is a duplicate, omit further processing
            return freq_quest.GetVector();

        //if there is only one frequency vector at the moment,
        //  we need to create a binary search structure next
        bins = new BinarySearchStructure( FrequencyVector::FrequencyComparer, binSize );

        if( bins == NULL )
            throw myruntime_error( mystring( "FrequencyStore: Not enough memory." ));

        //put the existing key in the binary search structure first
        if( ! StoreInBSearchStructure( bins, freq_quest ))
            throw myruntime_error( mystring( "FrequencyStore: Store failed." ));

        //now this position of the secondary hash will point to the binary search structure
        hash->SetValueAt( address, bins );
    }

    //put the original key in the secondary hash
    const void* inserted = StoreInBSearchStructure( bins, freq );

    //if the vector has been inserted
    if( inserted == freq.GetVector())
        no_collisions2++;

    return  inserted;
}

// -------------------------------------------------------------------------
// BinarySearchStructure: stores frequency vector (key) into the binary
//     search structure
//  returns address of the vector to be stored if frequencies have been
//  stored in the binary search structure, otherwise returns address of the
//  duplicate if there is one;
//  returns NULL on error
// -------------------------------------------------------------------------

const void* FrequencyStore::StoreInBSearchStructure( BinarySearchStructure* bins, const FrequencyVector& freq )
{
    if( !bins )
        return NULL;

    int     location = -1;
    //push the key into the appropriate place; the key
    //  may not be inserted if it is a duplicate of
    //  existing one
    bool    inserted = bins->Push( freq.GetVector(), &location );
    if( !inserted && 0 <= location ) {
        //if not inserted, construct a vector object
        const FrequencyVector   freq_quest(( const char* )bins->GetValueAt( location ));
        //check whether it is a duplicate and the values (scores) are equal if so
        /*inserted = */VectorsAreEqual( freq, freq_quest );
        return freq_quest.GetVector();
    }
    return  inserted? freq.GetVector(): NULL;
}

// -------------------------------------------------------------------------
// Find: finds frequency vector in the two-level hash system
// -------------------------------------------------------------------------

const void* FrequencyStore::Find( const FrequencyVector& freq ) const
{
    if( !primaryHash )
        return NULL;

    size_t      address = primaryHash->GetAddress( freq.GetVector());
    const void* value = primaryHash->GetValueAt( address );

    if( value == NULL ) {
        return NULL;
    }

    const HashTable*        secondHash = ( const HashTable* )value;
    const FrequencyVector   freq_quest (( const char* )value );

    if( freq_quest.IsVector())
        return value;

    //must search in the secondary hash table
    return FindInSecondaryHash( secondHash, freq );
}

// -------------------------------------------------------------------------
// FindInSecondaryHash: finds frequency vector in the secondary hash table
// -------------------------------------------------------------------------

const void* FrequencyStore::FindInSecondaryHash( const HashTable* hash, const FrequencyVector& freq ) const
{
    if( !hash )
        return NULL;

    //compute address for the key
    size_t  address = hash->GetAddress( freq.GetVector());
    const void* value = hash->GetValueAt( address );

    if( value == NULL ) {
        return NULL;
    }

    const BinarySearchStructure*    bins = ( const BinarySearchStructure* )value;
    const FrequencyVector           freq_quest (( const char* )value );

    if( freq_quest.IsVector())
        return value;

    //must search in the binary search structure
    return FindInBSearchStructure( bins, freq );
}

// -------------------------------------------------------------------------
// FindInBSearchStructure: finds frequency vector (key) in the binary
//     search structure
// -------------------------------------------------------------------------

const void* FrequencyStore::FindInBSearchStructure( const BinarySearchStructure* bins, const FrequencyVector& freq ) const
{
    if( !bins )
        return NULL;

    int location = -1;

    if( bins->Find( freq.GetVector(), &location ))
        return bins->GetValueAt( location );

    return NULL;
}

