/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __AbstractUniversalScoreMatrix__
#define __AbstractUniversalScoreMatrix__

#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "mystring.h"
#include "myexcept.h"
#include "libpro/srcpro/FrequencyStore.h"
#include "AbstractScoreMatrix.h"


extern double rint( double x );

class Serializer;

class FrequencyStore;
class BinarySearchStructure;

// _________________________________________________________________________
// Class ScoreProbability
//
class ScoreProbability
{
public:
    ScoreProbability( double, double );
    explicit ScoreProbability();
    ~ScoreProbability();

    int             GetScore() const                { return score; }
    double          GetDoubleScore() const          { return double_score; }
    double          GetProbability() const          { return probb; }

    void            SetScores( double value )       { double_score = value; score = ( int )rint( double_score ); }
    void            SetProbability( double value )  { probb = value; }

    void            IncProbability( double );

    static  int                 ScoreComparer( const void* key1, const void* key2 );
    static ScoreProbability*    NewScoreProbability();
    static void                 Destroy( const ScoreProbability* );

protected:

private:
    int     score;          //score
    double  double_score;   //same score of type double
    double  probb;          //probability
};


// _________________________________________________________________________
// Class AbstractUniversalScoreMatrix
//
// This class contains data and methods for manipulation of frequency pool
// in which distinct vectors represent distribution of possible positional
// combinations; the vectors are used to comprise scoring system for
// profile-profile comparisons
//
class AbstractUniversalScoreMatrix: public AbstractScoreMatrix
{
public:
    enum {
        ProbsHeader = 0xffffffff,       //header of the formatted probabilities message
        ScaleHeader = 0xfffffff7        //header of the formatted scale factor message
    };

public:
    AbstractUniversalScoreMatrix(
            const FrequencyStore*   store,
            TType           type,
            Configuration   config[NoSchemes],
            TBehaviour      beh,
            TScaling        scl,
            TMask           msk,
            TFVectorProbabilities,
            bool cpu = true
    );
    explicit AbstractUniversalScoreMatrix();
    virtual ~AbstractUniversalScoreMatrix();

    static bool         IsFormattedDataValid( const char* );    //verifies the formatted message for validness

    bool                PreferCPU() const                   { return cpu_preference; }
                                                                //output parameters in final print
    void                PrintFinal( FILE* fp ) const        { AbstractScoreMatrix::PrintFinal( fp ); }

    static size_t       GetMaxElemsOfProbCalculator();          //maximum number of elements in the probability structure
    static size_t       GetMaxSizeOfProbCalculator();           //maximum size required to contain the probability structure
    static size_t       GetMinSizeOfProbCalculator();           //minimum size of message can be valid

    static size_t       GetSizeOfScaleMessage();                //size of formatted message for scale factor

    const TFVectorProbabilities GetDistributionType() const { return distribution; }

    virtual void        Serialize( Serializer& ) const;         //serialization of the parameters
    virtual void        Deserialize( Serializer& );             //deserialization of the parameters

    virtual void        ComputeProfileScoringMatrix( bool final = false ) = 0;
                                                                //compute probabilities to observe scores at each position
    virtual void        ComputeScoreProbabilities( AttributableScores* );
                                                                //compute positional probabilities of scores
    virtual void        ComputePositionalScoreProbs( AttributableScores* );

protected:
    virtual double      GetScoreBase( int m, int n ) const = 0;

    virtual bool        ComputeScoreProbabilitiesCPU( AttributableScores* ) = 0;    //abstract method of ComputeScoreProbabilities
    virtual bool        ComputeScoreProbabilitiesMem( AttributableScores* ) = 0;    //abstract method of ComputeScoreProbabilities
                                                                                    //process scores and probabilities
    void                ProcessScoreProbabilities( AttributableScores*, double = 0.0 );

    const FrequencyStore*           GetStore() const            { return freqstore;         }
    BinarySearchStructure*          GetProbCalculator()         { return prob_calculator;   }
    const BinarySearchStructure*    GetProbCalculator() const   { return prob_calculator;   }
    static void                     ClearProbCalculator( BinarySearchStructure* );  //destroy members of probCalculator and clear it
    static void                     DestroyProbCalculator( BinarySearchStructure* );//destroy totally probCalculator

    bool                    Push( const ScoreProbability* );        //push score-probability strucutre into a private class member

    size_t                  FormatTermMessage( char* ) const;       //commpose termination message using the scale factor format
    bool                    IsThatTermMessage( const char* ) const; //whether the message is termination message

    size_t                  FormatScaleFactor( char*, double value ) const;
    size_t                  FormatScaleFactor( char* ) const;       //format message with a scale facor in it
    void                    DeformatScaleFactor( const char* );     //deformat message to extract scale factor written in it

    size_t                  FormatProbCalculator( char* ) const;    //format score probabilities and write them to the string stream
    BinarySearchStructure*  DeformatProbCalculator( const char* );  //deformat score probabilities and re-create binary search structure
                                                                    //infuse probabilities with new ones given with the argument
    bool                    InfuseProbCalculator( const BinarySearchStructure*, double*, bool delete_duplicates = true );

    void                    PrintProbCalculator( FILE* );                       //print information of scores and corresponding probabilities
    virtual void            USM_THROW( const char*, int = NOCLASS ) const = 0;  //throw method private for classes of this type

private:
                                                                    //for testing purposes
    void                    PrintStringinHex( const char*, size_t ) const;

private:
    const FrequencyStore*   freqstore;          //big array of frequency vectors
    BinarySearchStructure*  prob_calculator;    //vector used to calculate probabilities each time the scaling is performed
    const TFVectorProbabilities distribution;   //frequency vector distribution type
    bool                    cpu_preference;     //applying universal score system requires a large amount of memory;
                                                //  an alternative is to use more intense computations instead, i.e.
                                                //  CPU preference
    static const double     fake_scale_factor;  //fake scale factor used to indicate the termination
};


// INLINES ...
// CLASS ScoreProbability
//
// -------------------------------------------------------------------------
// ScoreComparer: the method is for comparison of scores; returns 0 if the
// keys are equal, value >0 if key1 > key2, and value <0 otherwise
//
inline
int ScoreProbability::ScoreComparer( const void* key1, const void* key2 )
{
    const ScoreProbability* scp1 = ( const ScoreProbability* )key1;
    const ScoreProbability* scp2 = ( const ScoreProbability* )key2;

    if( !scp1 || !scp2 )
        throw myruntime_error( mystring( "ScoreProbability: Unable to compare objects." ));

    return scp1->GetScore() - scp2->GetScore();
}

// -------------------------------------------------------------------------
// NewScoreProbability: creates new default object of class ScoreProbability
//
inline
ScoreProbability* ScoreProbability::NewScoreProbability()
{
    ScoreProbability*   scp = new ScoreProbability;
    if( scp == NULL )
        throw myruntime_error( mystring( "ScoreProbability: Not enough memory." ));
    return scp;
}

// -------------------------------------------------------------------------
// Destroy: destroys the object of class ScoreProbability
//
inline
void ScoreProbability::Destroy( const ScoreProbability* scp )
{
    if( scp )
        delete scp;
}

// -------------------------------------------------------------------------
// IncProbability: increases probability with the given value
//
inline
void ScoreProbability::IncProbability( double value )
{
    probb += value;
}

// =========================================================================
// INLINES of CLASS AbstractUniversalScoreMatrix 
//
// GetMaxElemsOfProbCalculator: returns maximum number of elements
//     allowed in the structure
//
inline
size_t AbstractUniversalScoreMatrix::GetMaxElemsOfProbCalculator()
{
    return MAX_RANGE_OF_SCORES;
}

// GetMaxElemsOfProbCalculator: returns maximum number of elements
//     allowed in the structure
//
inline
size_t AbstractUniversalScoreMatrix::GetMaxSizeOfProbCalculator()
{
    return  ( sizeof( double )  + sizeof( double )) * GetMaxElemsOfProbCalculator() +
            ( sizeof( int ))    +   //header of the message
            ( sizeof( size_t )) +   //space required to contain number of elements
            ( sizeof( int ));       //space required to contain CRC
}

// GetMinSizeOfProbCalculator: returns minimum size of message could be
//     valid in receipt
//
inline
size_t AbstractUniversalScoreMatrix::GetMinSizeOfProbCalculator()
{
    return  ( sizeof( int )) +      //header of the message
            ( sizeof( size_t )) +   //space required to contain number of elements
            ( sizeof( int ));       //space required to contain CRC
}

// GetSizeOfScaleMessage: size of formatted message for scale factor
//
inline
size_t AbstractUniversalScoreMatrix::GetSizeOfScaleMessage()
{
    return  ( sizeof( int )) +      //header of the message
            ( sizeof( double ));    //size of scale factor itself
}

// -------------------------------------------------------------------------
// IsFormattedDataValid: verifies wthether the message formatted is valid
//
inline
bool AbstractUniversalScoreMatrix::IsFormattedDataValid( const char* strstream )
{
    if( !strstream )
        return false;

    if( *( int* )strstream != ProbsHeader )
        return false;

    return true;
}


#endif//__AbstractUniversalScoreMatrix__
