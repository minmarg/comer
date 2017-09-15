/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __FrequencyStore__
#define __FrequencyStore__

#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "rc.h"
#include "mystring.h"
#include "myexcept.h"
#include "BinarySearchStructure.h"
#include "libhsh/HashTable.h"
#include "libhsh/RJHashing.h"
#include "libhsh/SBoxHashing.h"
#include "libhsh/HashFunctions.h"
#include "libhsh/MD5Hashing.h"
#include "libhsh/CRCHashing.h"

#include "datapro.h"


class Serializer;

//to avoid using templates
#define FREQTYPE    char
#define SIZETYPE  ( sizeof( FREQTYPE ))

#define WGHTTYPE    short
#define WGHTSIZE  ( sizeof( WGHTTYPE ))

#define INFMTYPE    unsigned char
#define INFMSIZE  ( sizeof( INFMTYPE ))
#define MAXINFVALUE 255

#define NOSTYPE     unsigned short
#define NOSSIZE   ( sizeof( NOSTYPE ))
#define MAXNOSVALUE 65535


// -------------------------------------------------------------------------
// Type determines how frequency vector probabilities shoud be computed
//

#define DTEXT_DISCRETE      ( "Discrete" )
#define DTEXT_PROVECTOR     ( "Profile" )
#define DTEXT_MULTINOMIAL   ( "Multinomial" )

enum TFVectorProbabilities {
    DISCRETE,
    PROVECTOR,
    MULTINOMIAL,
    DTypeUNKNOWN
};

// _________________________________________________________________________
// Class FrequencyVector
//
class FrequencyVector
{
    enum {  constant = 0xffffffff };                    //flag value
    enum {  preamble = sizeof( unsigned int ),          //dummy field for indication of frequency vector existance in memory
            no_Elems = NUMAA,                           //number of elements in a vector
            vectorSz = no_Elems * SIZETYPE,             //size of vector = frequencies
            alphaSz  = WGHTSIZE,                        //positional weight for frequencies, alpha
            thickSz  = NOSSIZE,                         //effective thickness
            scoreSz  = no_Elems * sizeof( short ),      //scores for the frequencies
            inforSz  = INFMSIZE,                        //information content, relative entropy
            lprobSz  = sizeof( double ),                //log-probability to occur for this vector of frequencies
            Overall  = preamble + vectorSz + alphaSz + thickSz + scoreSz + inforSz + lprobSz
    };

public:                                                 //for writing
    FrequencyVector( const double freq[NUMALPH], double weight, size_t thickn, const double scores[NUMALPH], double info );
    FrequencyVector( const char*  vect );               //for reading
    explicit FrequencyVector();
    ~FrequencyVector();

    void        Destroy();
                                                        //change contents of the vector
    void        SetVector( const double freq[NUMALPH], double weight, size_t thickn, const double scores[NUMALPH], double info );

    INFMTYPE    GetInfContent() const;                  //get information content
    NOSTYPE     GetThickness() const;                   //get effective thickness
    WGHTTYPE    GetWeight() const;                      //get weight (known as alpha)
    FREQTYPE    GetValueAt( size_t n ) const;
    short       GetScoreAt( size_t n ) const;
    double      GetProbability() const;
    double      ComputeProVectorProbability();
    double      ComputeMultinomialProbability();
    double      ComputeProbability();                   //compute log-probability of vector to occur
    void        UpdateProbability();                    //increase probability of vector to occur (non-multinomial case)
    void        NormalizeProbability( double );         //normalize probability
    void        FinalProbability( size_t );             //set final probability value (non-multinomial case)
                                                        //number of logical components in the vector
    static int  GetNoElems()            { return no_Elems; }
    const char* GetVector() const       { return vector; }

    bool        IsVector() const;

    void        Serialize( Serializer& ) const;         //data serialization
    void        Deserialize( Serializer& );             //data deserialization

    void        TextWriteVector( FILE* fp ) const;
    void        TextReadVector( FILE* fp );


    static int      FreqENOSComparer( const void* vect1, const void* vect2 );
    static int      FrequencyComparer( const void* vect1, const void* vect2 );
    static int      FrequencyValueComparer( const void* vect1, const void* vect2 );

    static size_t   RJHashing( const void* key );
    static size_t   SBoxHashing( const void* key );
    static size_t   OATHashing( const void* key );
    static size_t   MD5Hashing( const void* key );
    static size_t   CRCHashing( const void* key );
    static size_t   MixHashing( const void* key );

    void        Print( FILE* ) const;

protected:
    void        Init();
    void        Init( const double[NUMALPH], double, size_t, const double[NUMALPH], double );

    void        SetValueAt( size_t n, FREQTYPE value );
    void        SetScoreAt( size_t n, short );

    void        SetInfContent( INFMTYPE value );
    void        SetThickness( NOSTYPE value );
    void        SetWeight( WGHTTYPE value );
    void        SetProbability( double value );         //set log-probability value

    static const char*  ExtractKey( const void* vect )  { return ( const char* )vect + preamble; }
    static size_t       SizeofPrivateKey()              { return vectorSz + alphaSz + thickSz; }
    static size_t       SizeofKey()                     { return vectorSz + alphaSz + thickSz + scoreSz; }
    //scoreSz above in SizeofKey is added artificially and is necessary in case where scaling of individual profiles is used

    static const char*  ExtractValue( const void* vect )  { return ( const char* )vect + preamble + SizeofPrivateKey(); }
    static size_t       SizeofValue()                     { return scoreSz; }

private:
    char*                   vector;                     //frequencies for effective number of residues
    static THashFunction    hash_FUNCTIONS[];           //distinct hash functions
    static const size_t     no_hash_FUNCTIONS;          //number of hash functions
};

// _________________________________________________________________________
// Class FrequencyStore
//
//  This class is for storage of frequency vectors. Frequency vectors are
//  hashed in a two-level hashing structure with binary search structures
//  in case of collisions. The top-level hash function gives address of the
//  secondary hash structure designed for the case of collisions. If
//  collisions appear in the secondary hash structures anyway, dynamic
//  binary search structures are used.
//
class FrequencyStore
{
public:
    enum {
        TFreqSimpleVector,
        TFreqBSEnosStructure
    };
    enum {
        FirstLevelHashBits = 22,    //number of significant bits used for the first-level hash size
        SecndLevelHashBits = 7,     //number of significant bits used for the second-level hash size
        binSize = 20,               //initial size of binary search structure
        FrequenciesBits = 21        //number of significant bits used for the size of vector frequencies;
                                    //  make size two times smaller than one of the primary hash
    };

public:
    FrequencyStore( int fstype = TFreqSimpleVector );
    ~FrequencyStore();

    const void* Store( const FrequencyVector& freq );       //store frequency vector in the two-level hash system
    const void* Find( const FrequencyVector& freq ) const;  //find frequency vector

    size_t  GetCollisions1() const          { return no_collisions1; }
    size_t  GetCollisions2() const          { return no_collisions2; }


    const SimpleVector* GetFrequencies() const { return frequencies; }
    void                SetNoFrequencyVectors( size_t );    //this number affects size of vector frequencies

    bool    IsConsistent( double count ) const;             //verify consistency of frequency vectors
    bool    IsConsistent( double count, size_t ) const;     //verify consistency of frequency vectors

    static  mystring                GetDistributionText( int type );
    static  TFVectorProbabilities   GetDistributionType( const mystring& distrstr );
    static  TFVectorProbabilities   GetDistributionType()           { return probtype; }
    static  void SetDistributionType( TFVectorProbabilities type )  { probtype = type; }

protected:
    void    Init();
    void    Destroy();                                      //destroy whole two-level hash system
    void    DestroySecondaryHash( const void* );            //destroy secondary hash structure at the address
    void    DestroyBSearchStructure( const void* );         //destroy binary search structure at the address

    const void* StoreInSecondaryHash( HashTable*, const FrequencyVector& );
    const void* StoreInBSearchStructure( BinarySearchStructure*, const FrequencyVector& );

    const void* FindInSecondaryHash( const HashTable*, const FrequencyVector& ) const;
    const void* FindInBSearchStructure( const BinarySearchStructure*, const FrequencyVector& ) const;

    void    DestroyFrequencies();                           //destroy vector frequencies
    void    NewFrequencies( size_t );                       //allocate space for frequencies
    void    PushInFrequencies( const void* );               //store address of frequency vector

                                                            //whether two vectors are equal
    bool    VectorsAreEqual( const FrequencyVector&, const FrequencyVector& );

    static size_t   FirstLevelHashing( const void* key );
    static size_t   SecndLevelHashing( const void* key );

    int     GetFreqStructType() const { return freqstructtype_; }

private:
    static TFVectorProbabilities    probtype;       //probability distribution type
    static THashFunction    func1stlevel;   //first-level hash function
    static THashFunction    func2ndlevel;   //second level hash function
    HashTable*      primaryHash;            //hash to store frequency vectors for instant access
    SimpleVector*   frequencies;            //contains pointers to frequency vectors
    size_t      no_collisions1;             //1st-level collisions
    size_t      no_collisions2;             //2nd-level collisions
    const int   freqstructtype_;            //type of stucture to hold frequencies
};

////////////////////////////////////////////////////////////////////////////
// Class FrequencyVector inlines
//
// GetScore: get appropriate score in the vector

inline
short FrequencyVector::GetScoreAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !vector || no_Elems <= n )
        throw myruntime_error(
            mystring( "FrequencyVector: Memory access error." ));
#endif
    static int  beg = preamble + SizeofPrivateKey();
    return *( short* )( vector + beg + n * sizeof( short ));
}

// -------------------------------------------------------------------------
// SetScore: set score at the position

inline
void FrequencyVector::SetScoreAt( size_t n, short value )
{
#ifdef __DEBUG__
    if( !vector || no_Elems <= n )
        throw myruntime_error(
            mystring( "FrequencyVector: Memory access error." ));
#endif
    static int  beg = preamble + SizeofPrivateKey();
    *( short* )( vector + beg + n * sizeof( short )) = value;
}

// -------------------------------------------------------------------------
// IsVector: return flag indicating that this is a vector

inline
bool FrequencyVector::IsVector() const
{
#ifdef __DEBUG__
    if( !vector )
        throw myruntime_error(
            mystring( "FrequencyVector: Memory access error." ));
#endif

    return *( unsigned* )vector == constant;
}

// -------------------------------------------------------------------------
// GetValueAt: get frequency value at the appropriate position

inline
FREQTYPE FrequencyVector::GetValueAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !vector || no_Elems <= n )
        throw myruntime_error(
            mystring( "FrequencyVector: Memory access error." ));
#endif

    return *( FREQTYPE* )( vector + preamble + n * SIZETYPE );
}

// -------------------------------------------------------------------------
// SetValueAt: set frequency value at the appropriate position

inline
void FrequencyVector::SetValueAt( size_t n, FREQTYPE value )
{
#ifdef __DEBUG__
    if( !vector || no_Elems <= n )
        throw myruntime_error(
            mystring( "FrequencyVector: Memory access error." ));
#endif

    *( FREQTYPE* )( vector + preamble + n * SIZETYPE ) = value;
}

// -------------------------------------------------------------------------
// GetWeight: get weight value

inline
WGHTTYPE FrequencyVector::GetWeight() const
{
#ifdef __DEBUG__
    if( !vector )
        throw myruntime_error(
            mystring( "FrequencyVector: Memory access error." ));
#endif
    return *( WGHTTYPE* )( vector + preamble + vectorSz );
}

// -------------------------------------------------------------------------
// SetWeight: set weight value

inline
void FrequencyVector::SetWeight( WGHTTYPE value )
{
#ifdef __DEBUG__
    if( !vector )
        throw myruntime_error(
            mystring( "FrequencyVector: Memory access error." ));
#endif
    *( WGHTTYPE* )( vector + preamble + vectorSz ) = value;
}

// -------------------------------------------------------------------------
// GetThickness: get value of effective thickness

inline
NOSTYPE FrequencyVector::GetThickness() const
{
#ifdef __DEBUG__
    if( !vector )
        throw myruntime_error(
            mystring( "FrequencyVector: Memory access error." ));
#endif
    return *( NOSTYPE* )( vector + preamble + vectorSz + alphaSz );
}

// -------------------------------------------------------------------------
// SetThickness: set value for effective thickness

inline
void FrequencyVector::SetThickness( NOSTYPE value )
{
#ifdef __DEBUG__
    if( !vector )
        throw myruntime_error(
            mystring( "FrequencyVector: Memory access error." ));
#endif
    *( NOSTYPE* )( vector + preamble + vectorSz + alphaSz ) = value;
}

// -------------------------------------------------------------------------
// GetInfContent: extract information content

inline
INFMTYPE FrequencyVector::GetInfContent() const
{
#ifdef __DEBUG__
    if( !vector )
        throw myruntime_error(
            mystring( "FrequencyVector: Memory access error." ));
#endif
    return *( INFMTYPE* )( vector + preamble + SizeofKey());
}

// -------------------------------------------------------------------------
// SetInfContent: set information content

inline
void FrequencyVector::SetInfContent( INFMTYPE value )
{
#ifdef __DEBUG__
    if( !vector )
        throw myruntime_error(
            mystring( "FrequencyVector: Memory access error." ));
#endif
    *( INFMTYPE* )( vector + preamble + SizeofKey()) = value;
}

// -------------------------------------------------------------------------
// GetProbability: get log-probability value

inline
double FrequencyVector::GetProbability() const
{
#ifdef __DEBUG__
    if( !vector )
        throw myruntime_error(
            mystring( "FrequencyVector: Memory access error." ));
#endif
    return *( double* )( vector + preamble + SizeofKey() + inforSz );
}

// -------------------------------------------------------------------------
// SetProbability: set log-probability value

inline
void FrequencyVector::SetProbability( double value )
{
#ifdef __DEBUG__
    if( !vector )
        throw myruntime_error(
            mystring( "FrequencyVector: Memory access error." ));
#endif
    *( double* )( vector + preamble + SizeofKey() + inforSz ) = value;
}


// =========================================================================
// CLASS FrequencyStore inlines
//
// IsConsistent: verify probability of frequency vectors for consistency
// -------------------------------------------------------------------------

inline
bool FrequencyStore::IsConsistent( double count ) const
{
    return IsConsistent( count, GetFrequencies()->GetSize());
}

inline
bool FrequencyStore::IsConsistent( double count, size_t no_vectors ) const
{
    //the method is for uniform distribution only
    if( FrequencyStore::GetDistributionType() != DISCRETE )
        return true;

    if( !GetFrequencies())
        return false;

    return ( size_t )count == no_vectors;
}

// -------------------------------------------------------------------------
// DestroyFrequencies: destroy vector frequencies
// -------------------------------------------------------------------------

inline
void FrequencyStore::DestroyFrequencies()
{
    if( frequencies )
        delete frequencies;
}

// -------------------------------------------------------------------------
// NewFrequencies: allocates space for vector frequencies
// -------------------------------------------------------------------------

inline
void FrequencyStore::NewFrequencies( size_t size )
{
    if( GetFreqStructType() == TFreqBSEnosStructure )
        frequencies = new BinarySearchStructure( &FrequencyVector::FreqENOSComparer, size, true );
    else
        frequencies = new SimpleVector( size );

    if( frequencies == NULL )
        throw myruntime_error("FrequencyStore: Not enough memory.");
}

// -------------------------------------------------------------------------
// PushInFrequencies: push one frequency vector (distribution) into their
//     vector
// -------------------------------------------------------------------------

inline
void FrequencyStore::PushInFrequencies( const void* vect )
{
    if( !frequencies )
        NewFrequencies( hashsize( FrequenciesBits ));
    frequencies->Push( vect );
}

// -------------------------------------------------------------------------
// FirstLevelHashing: perform hash computing for the first level hash table
// -------------------------------------------------------------------------

inline
size_t FrequencyStore::FirstLevelHashing( const void* key )
{
#ifdef __DEBUG__
    if( !func1stlevel )
        throw myruntime_error(
            mystring( "FrequencyStore: Unable to perform hash computations." ));
#endif
    size_t  address = ( *func1stlevel )( key );
//     return  ( address >> (( sizeof( address ) << 3 ) - FirstLevelHashBits ))& hashmask( FirstLevelHashBits );
    return  address & hashmask( FirstLevelHashBits );
}

// -------------------------------------------------------------------------
// PushInFrequencies: perform hash computing for the second level hash table
// -------------------------------------------------------------------------

inline
size_t FrequencyStore::SecndLevelHashing( const void* key )
{
#ifdef __DEBUG__
    if( !func2ndlevel )
        throw myruntime_error(
            mystring( "FrequencyStore: Unable to perform hash computations." ));
#endif
    size_t  address = ( *func2ndlevel )( key );
    return  address & hashmask( SecndLevelHashBits );
}

// -------------------------------------------------------------------------
// VectorsAreEqual: verify if two vectors are equal; if so check that
//     values of the keys are equal as well
// -------------------------------------------------------------------------

inline
bool FrequencyStore::VectorsAreEqual( const FrequencyVector& one, const FrequencyVector& another )
{
    if( FrequencyVector::FrequencyComparer( one.GetVector(), another.GetVector()) == 0 ) {
#ifdef __DEBUG__
        //if scaling of profiles has been performed, there may be cases when scores are not equal
        //  even the keys (frequency vectors) are equal but scores then differ marginally and have
        //  no significant effect on tendency of score distribution
        if( FrequencyVector::FrequencyValueComparer( one.GetVector(), another.GetVector()) != 0 ) {
// one.Print( stderr );
// another.Print( stderr );
            error( "FrequencyStore: Scores of two equal frequency vectors are not equal." );
        }
#endif
        return true;
    }
    return false;
}

// -------------------------------------------------------------------------
// GetDistributionText: return distribution description given type
//
inline
mystring FrequencyStore::GetDistributionText( int distr )
{
    switch( distr ) {
        case DISCRETE:      return DTEXT_DISCRETE;
        case PROVECTOR:     return DTEXT_PROVECTOR;
        case MULTINOMIAL:   return DTEXT_MULTINOMIAL;
    }
    return mystring();
}

// GetDistributionType: return distribution type given description
//
inline
TFVectorProbabilities FrequencyStore::GetDistributionType( const mystring& distrstr )
{
    if( distrstr.empty())
        return DTypeUNKNOWN;

    if( !distrstr.ncompare( DTEXT_DISCRETE ))       return DISCRETE;
    else if( !distrstr.ncompare( DTEXT_PROVECTOR )) return PROVECTOR;
    else if( !distrstr.ncompare( DTEXT_MULTINOMIAL))return MULTINOMIAL;

    return DTypeUNKNOWN;
}


#endif//__FrequencyStore__
