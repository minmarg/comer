/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/



#ifndef __AttributableScores__
#define __AttributableScores__

#include <stdlib.h>

#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "data.h"
#include "mystring.h"
#include "myexcept.h"
#include "libpro/srcpro/GapScheme.h"

#define ACWINDOW    4

// extern double rint( double x );

class AbstractScoreMatrix;


////////////////////////////////////////////////////////////////////////////
// CLASS AttributableScores
// This class is to contain a particular score table and to perform
// computations to derive statistical parameters subject to the table
//
class AttributableScores
{
public:
    typedef void    ( AbstractScoreMatrix::*TrogressFunction )( AttributableScores* );
    typedef bool    ( AbstractScoreMatrix::*TScoreReductFunction )( AttributableScores* );
    typedef void    ( AbstractScoreMatrix::*TProbabilityFunction )( AttributableScores* );
    typedef void    ( AttributableScores::*TConservationFunction )( double, double*, double*, int, void* );
    typedef double  ( AttributableScores::*TGetProbFunction )( TScore, int n ) const;
    typedef double  ( AttributableScores::*TGetProbNormFunction )( int n ) const;

    AttributableScores(
            AbstractScoreMatrix* prnt,
            TProbabilityFunction probfun,
            TScoreReductFunction redtfun,
            TrogressFunction progressfun,
            double reflambda,
            double refH,
            double refK,
            const double** im,
            const TMask** msk,
            bool keep,
            int as_factor = 1 );
    virtual ~AttributableScores();

    bool                IsValid() const;                                //whether score matrix is valid

    double              GetImageScore( int m, int n ) const;            //returns image score at specified profile positions

    bool                GetMaskedUnmasked( int m, int n  ) const    { return GetMasked( m, n  ) == Unmasked; }
    bool                GetMaskedToIgnore( int m, int n  ) const    { return GetMasked( m, n  ) == MaskToIgnore; }
    bool                GetMaskedToConsider( int m, int n  ) const  { return GetMasked( m, n  ) == MaskToConsider; }

    TMask               GetMasked( int m, int n ) const;                //score mask at the position

    virtual double              GetScore( int m, int n ) const = 0;     //returns score at specified profile positions
    virtual void                SetScore( double, int m, int n ) = 0;   //set score at specified positions

    const double* const         GetQueryInfContent() const             { return queryinfcontent; }
    const double* const         GetSbjctInfContent() const             { return sbjctinfcontent; }

    int                 GetQuerySize() const                        { return querySize;   }
    int                 GetSubjectSize() const                      { return subjectSize; }

    double              GetRefLambda() const                        { return referenceLambda; }
    double              GetRefH() const                             { return referenceH; }
    double              GetRefK() const                             { return referenceK; }

    double              GetLambda() const                           { return lambda; }      //computed ungapped Lambda
    double              GetH() const                                { return entropy; }     //computed relative entropy
    double              GetK() const                                { return parameterK; }  //computed ungapped parameter K
    double              GetExpectedScore() const                    { return expscore; }    //expected score
    void                SetExpectedScore( double E )                { expscore = E; }

    bool                GetAllNegatives() const                     { return allnegatives; }
    void                SetAllNegatives( bool value )               { allnegatives = value; }

    bool                GetKeepInMemory() const                     { return keepinmem; }

    double              GetConstantForHAdjustment() const           { return constant_c; }

    virtual void                Init( int sz_query, int sz_sbjct ) = 0;

    bool                SearchForHSPs( double minhspscore, int hsplen, int nohsps, int maxdist, int* = NULL, int* = NULL );

    void                AdjustGaps(
                                const LogOddsMatrix& qlogo, const LogOddsMatrix& slogo,
                                GapScheme& qg, GapScheme& sg,
                                bool autoc, int acwindow ) const;

    void                ScaleScoringMatrix();
    void                ComputeStatisticalParameters( bool = true, bool wrn = true);//calc. statistical parameters
    void                ComputeProbabilitiesAndLambda( bool wrn = true);//compute score probabilities and statistical parameter lambda

    virtual void                ComputeHConstant();                     //compute constant c used in the adjustment computations of H

    int                 GetAutoScalingFactor() const            { return auto_scaling_factor; }

    double              GetPrivateMultiplier() const            { return private_multiplier; }
    void                SetPrivateMultiplier( double value )    { private_multiplier = value; }

    void                MultiplyScoresByPrivateMultiplier();            //multiply scores by the private multiplier

    TScore              GetMinScore() const { return min_score; }       //obtain min score
    TScore              GetMaxScore() const { return max_score; }       //obtain max score

    void                SetMinMaxScores( double min, double max )       { SetMinMaxScores(( TScore )rint( min ), ( TScore )rint( max )); }
    void                SetProbabilityOf( double value, double sc )     { SetProbabilityOf( value, ( TScore )rint( sc )); }
    void                IncProbabilityOf( double value, double sc )     { IncProbabilityOf( value, ( TScore )rint( sc )); }
    void                DivideProbabilityOf( double value, double sc )  { DivideProbabilityOf( value, ( TScore )rint( sc )); }

    void                SetMinMaxScores( TScore min, TScore max );      //set min max scores and allocates space for prob vector if needed
    double              GetProbabilityOf( TScore ) const;               //get probability given score
    void                SetProbabilityOf( double value, TScore );       //set probability value
    void                IncProbabilityOf( double value, TScore );       //increment probability value
    void                DivideProbabilityOf( double value, TScore );    //divide probability by value


    void                ComputePositionalInfContents();                 //compute entropies at each position of score system

    TScore              GetQueryMinScoreAt( int n ) const;              //get min score at query position n
    TScore              GetQueryMaxScoreAt( int n ) const;              //get max score at query position n
    TScore              GetSbjctMinScoreAt( int m ) const;              //get min score at subject position n
    TScore              GetSbjctMaxScoreAt( int m ) const;              //get max score at subject position n

    void                SetQueryInfMinMaxScoresAt( TScore min, TScore max, int n ); //set min max scores at pos. n
    void                SetQueryInfProbNorm( double value, int n );         //set norm. term of probability at pos. n
    double              GetQueryInfProbNorm( int n ) const;                 //get norm. term of probability at pos. n
    double              GetQueryInfProbOf( TScore, int n ) const;           //get probability at pos. n given score
    double              IncQueryInfProbOf( double value, TScore, int n );   //increment probability value for score at pos. n
    double              DivQueryInfProbOf( double value, TScore, int n );   //divide probability of score by value at pos. n

    void                SetSbjctInfMinMaxScoresAt( TScore min, TScore max, int m ); //set min max scores at pos. m
    void                SetSbjctInfProbNorm( double value, int m );         //set norm. term of probability at pos. m
    double              GetSbjctInfProbNorm( int m ) const;                 //get norm. term of probability at pos. m
    double              GetSbjctInfProbOf( TScore, int m ) const;           //get probability at pos. m given score
    double              IncSbjctInfProbOf( double value, TScore, int m );   //increment probability value for score at pos. m
    double              DivSbjctInfProbOf( double value, TScore, int m );   //divide probability of score by value at pos. m

    void                CheckforExpectedScore();//indirectly adjust score by expected value

protected:
    enum TTable {       //type indicating the presence of table type: image, scores
        ImageTable,     //  use double**image
        ScoreTable      //  use TScore**scores
    };

    explicit AttributableScores();

    void                SetQuerySize( int value )                   { querySize = value;   }
    void                SetSubjectSize( int value )                 { subjectSize = value; }

    void                Init();                             //IMPORTANT: query and subject sizes are supposed to be initialized

    double              GetQueryInfContentAt( int n ) const;        //information content at specified position
    double              GetSbjctInfContentAt( int m ) const;        //information content at specified position

    void                SetQueryInfContentAt( int n, double );      //set information content at specified position
    void                SetSbjctInfContentAt( int m, double );      //set information content at specified position

    void                SetConstantForHAdjustment( double value )   { constant_c = value; }

    void                SetRefLambda( double value )                { referenceLambda = value; }
    void                SetRefH( double value )                     { referenceH = value; }
    void                SetRefK( double value )                     { referenceK = value; }

    void                SetLambda( double newlambda )       { lambda = newlambda; }
    void                SetH( double H )                    { entropy = H; }
    void                SetK( double K )                    { parameterK = K; }

    void                Autocorrelation( double ( *maxs )[ACWINDOW], int length, double* output, int outlen ) const;

    void                AdjustGapsByAutocorr(
                                const LogOddsMatrix& qlogo, const LogOddsMatrix& slogo,
                                GapScheme& qg, GapScheme& sg, int acwindow ) const;
    void                AdjustGapsByAutocorr2(
                                const LogOddsMatrix& qlogo, const LogOddsMatrix& slogo,
                                GapScheme& qg, GapScheme& sg, int acwindow ) const;

    void                ComputeKHGivenLambda();                     //compute entropy and Karlin's parameter K

    virtual double              FindLambdaRootAtPosition(           //find lambda at specified position given probs.
        TConservationFunction,
        int gcd,
        double x1,
        double x2,
        double xacc,
        int maxit,
        int minscore,
        int maxscore,
        int pos,
        TGetProbFunction ) = 0;

    virtual double              findLambdaRoot( TConservationFunction, int, double, double, double, int, void* = NULL ) = 0;
    virtual double              findConstantRoot( TConservationFunction, int, double, double, double, int, void* = NULL ) = 0;
    double              rootByNRandBisection( TConservationFunction, int, double, double, double, int, void*, bool = true );
    void                testRoot();

    //functions describing conservation equations
    virtual void                lambdaconspos( double x, double* f, double* df, int step, void* params ) = 0;
    virtual void                conservation( double, double*, double*, int step = 1, void* = NULL ) = 0;
    virtual void                relentropy_conservation( double x, double* f, double* df, int step = 1, void* = NULL ) = 0;
    //

    virtual void                ScaleToAttainLambda() = 0;          //perform iterative scaling until required lambda is attained

    void                ComputeGCD();                               //computation of the greatest common divisor

    void                ComputeConstantToAdjustH();                 //compute constant c to be able to adjust relative entropy of the scores

    void                ComputePositionalScoreProbs();              //compute positional probabilities of scores
    virtual void        ComputeScoreProbabilities();                //compute probabilities to observe scores at each position
    void                ComputeScalingLambda( bool wrn = true );    //compute scaling parameter lambda
    virtual void                ComputeEntropyGivenLambda() = 0;    //compute relative entropy given previously computed lambda
    virtual double              ComputePositionalEntropy(           //compute entropy at specified position of score system
                            double lambda, int minscore, int maxscore, int pos,
                            TGetProbFunction, TGetProbNormFunction ) const = 0;
    virtual void        ComputeKarlinsK();                          //compute parameter K

    void                AdjustScoresInMatrix();                     //adjust scores in the matrix given scaling parameter lambda
    void                AdjustScoresInMatrix( double, TTable, bool direct = true ); //adjust scores in the matrix overloaded

    void                NewProbabilities( size_t );                         //allocate memory for probabilities
    void                DestroyProbabilities();                             //deallocate memory for probabilities


    TScore              GetGCD() const      { return score_gcd; }           //get the greatest common divisor

    virtual void                MultiplyScoreBy( double, int m, int n ) = 0;//multiply score at the specified position by a factor

    void                SetMinScore( TScore value ) { min_score = value; }  //set min score
    void                SetMaxScore( TScore value ) { max_score = value; }  //set max score
    void                SetGCD( TScore value )      { score_gcd = value; }  //set the greatest common divisor

    void                PrintProbabilities( FILE* );                        //print probability values

    //Private routines for manipulation of probability vector used to compute parameter K
    void                NewPrivateProbVector( size_t maxit, size_t range );
    void                DestroyPrivateProbVector();
    void                ResetPrivateProbVector( size_t = 0 );
    double*             GetPrivateProbVector() const        { return priv_prob_vector; }
    size_t              GetPrivateProbVectorSize() const    { return prob_vector_size; }
    //End of provate routines.


protected:
    virtual void** const    GetScores() const = 0;

    AbstractScoreMatrix*        GetParent()                 { return parent; }
    const AbstractScoreMatrix*  GetParent() const           { return parent; }
    TProbabilityFunction    GetParentProbFunction()         { return parent_probfunction; }
    TScoreReductFunction    GetParentRedtFunction()         { return parent_redtfunction; }
    TrogressFunction        GetParentProgressFunction()     { return parent_progressfunction; }

    const double** const    GetImage() const                { return image; }
    const TMask** const     GetMask() const                 { return mask; }

    void                NewQueryInfContent( int size );
    void                NewSbjctInfContent( int size );
    void                DestroyQueryInfContent();
    void                DestroySbjctInfContent();

    void                NewQueryInfProbsAt( int size, int n );      //allocate memory of probabilities for inf. content at pos. n
    void                NewQueryInfProbs();                         //allocate memory of probabilities for inf. content
    void                DestroyQueryInfProbs();                     //deallocate memory of probabilities
    void                NewSbjctInfProbsAt( int size, int n );      //allocate memory of probabilities for inf. content at pos. m
    void                NewSbjctInfProbs();                         //allocate memory of probabilities for inf. content
    void                DestroySbjctInfProbs();                     //deallocate memory of probabilities

    void                NewQueryMinMaxScores();                     //allocate memory to keep min/max scores for query positions
    void                NewSbjctMinMaxScores();                     //allocate memory to keep min/max scores for subject positions
    void                DestroyQueryMinMaxScores();                 //deallocate memory of min/max scores for query positions
    void                DestroySbjctMinMaxScores();                 //deallocate memory of min/max scores for subject positions

    void                MakeQuerySbjctMasks( bool* queryvect, int qlen, bool* sbjctvect, int slen ) const;

private:
    double          private_multiplier;  
    const int       auto_scaling_factor;

private:
    AbstractScoreMatrix*    parent;             //Parent score matrix
    TProbabilityFunction    parent_probfunction;//Parent probability function
    TScoreReductFunction    parent_redtfunction;//parent score reduction function
    TrogressFunction        parent_progressfunction;//Parent function to show progress of the values of score multiplier
    const double**          image;              //private matrix of double scores
//  void**                  scores;             //type of matrix of scores depends on implementation of this class
    const TMask**           mask;               //matrix of score masks

    double*                 queryinfcontent;    //information content along the query positions 10.21
    double*                 sbjctinfcontent;    //information content along the subject positions 10.21
    double**                queryinfprobabs;    //probabilities of information content for query positions
    double**                sbjctinfprobabs;    //probabilities of information content for subject positions
    int*                    queryminmaxscores;  //min/max scores along the query positions
    int*                    sbjctminmaxscores;  //min/max scores along the subject positions

    double*                 probabilities;      //score probabilities
    TScore                  min_score;          //minimum score value
    TScore                  max_score;          //maximum score value
    TScore                  score_gcd;          //the greatest common divisor of scores

    double*                 priv_prob_vector;   //private probability vector
    size_t                  prob_vector_size;   //size of private probability vector

    double                  expscore;           //expected score per column pair
    bool                    allnegatives;       //whether scores are all negative

    double                  referenceLambda;    //reference lambda parameter
    double                  referenceH;         //reference parameter H
    double                  referenceK;         //reference parameter K

    double                  lambda;             //scaling parameter, lambda
    double                  entropy;            //entropy describing information per aligned pair of positions (parameter H)
    double                  parameterK;         //Karlin's parameter K

    double                  constant_c;         //constant c used to adjust relative entropy of scores

    int                     querySize;          //length of query sequence (profile)
    int                     subjectSize;        //length of subject sequence (profile)

    bool                    keepinmem;          //whether to keep storage in memory
};

// INLINES ...

inline bool AttributableScores::IsValid() const
{
    if( GetKeepInMemory())
        return GetScores() != NULL && 0 < GetQuerySize() && 0 < GetSubjectSize();
    return true;
}

// -------------------------------------------------------------------------
// GetImageScore: image score at the profile positions
// -------------------------------------------------------------------------

inline
double AttributableScores::GetImageScore( int m, int n ) const
{
#ifdef __DEBUG__
    if( !image || subjectSize <= m || m < 0 || querySize <= n || n < 0 )
        throw myruntime_error(
            mystring( "AttributableScores: Memory access error." ));
#endif

    return image[m][n];
}

// -------------------------------------------------------------------------
// GetMasked: score mask at the position
// -------------------------------------------------------------------------

inline
TMask AttributableScores::GetMasked( int m, int n ) const
{
#ifdef __DEBUG__
    if( !mask || subjectSize <= m || m < 0 || querySize <= n || n < 0 )
        throw myruntime_error(
            mystring( "AttributableScores: Memory access error." ));
#endif

    return mask[m][n];
}

// -------------------------------------------------------------------------
// GetQueryInfContentAt: information content at specified query position
// -------------------------------------------------------------------------
inline
double AttributableScores::GetQueryInfContentAt( int n ) const
{
#ifdef __DEBUG__
    if( !queryinfcontent || querySize <= n || n < 0 )
        throw myruntime_error(
            mystring( "AttributableScores: Memory access error." ));
#endif

    return queryinfcontent[n];
}

// -------------------------------------------------------------------------
// GetSbjctInfContentAt: information content at specified subject position
// -------------------------------------------------------------------------
inline
double AttributableScores::GetSbjctInfContentAt( int m ) const
{
#ifdef __DEBUG__
    if( !sbjctinfcontent || subjectSize <= m || m < 0 )
        throw myruntime_error(
            mystring( "AttributableScores: Memory access error." ));
#endif

    return sbjctinfcontent[m];
}

// -------------------------------------------------------------------------
// SetQueryInfContentAt: set information content at query position
// -------------------------------------------------------------------------
inline
void AttributableScores::SetQueryInfContentAt( int n, double value )
{
#ifdef __DEBUG__
    if( !queryinfcontent || querySize <= n || n < 0 )
        throw myruntime_error(
            mystring( "AttributableScores: Memory access error." ));
#endif

    queryinfcontent[n] = value;
}

// -------------------------------------------------------------------------
// SetSbjctInfContentAt: set information content at specified position
// -------------------------------------------------------------------------
inline
void AttributableScores::SetSbjctInfContentAt( int m, double value )
{
#ifdef __DEBUG__
    if( !sbjctinfcontent || subjectSize <= m || m < 0 )
        throw myruntime_error(
            mystring( "AttributableScores: Memory access error." ));
#endif

    sbjctinfcontent[m] = value;
}

// -------------------------------------------------------------------------
// MultiplyScoresByPrivateMultiplier: multiply scores by the private
//     multiplier
// -------------------------------------------------------------------------

inline
void AttributableScores::MultiplyScoresByPrivateMultiplier()
{
    AdjustScoresInMatrix( GetPrivateMultiplier(), ImageTable, false/*direct*/ );
}

// -------------------------------------------------------------------------
// SetMinMaxScores(): set min max scores and allocates space for probability
//     vector if needed
// -------------------------------------------------------------------------

inline
void AttributableScores::SetMinMaxScores( TScore min, TScore max )
{
    if( max < min )
        throw myruntime_error( mystring( "AttributableScores: Max score is less than min score." ));

    NewProbabilities( max - min + 1 );

    SetMinScore( min );
    SetMaxScore( max );
}

// -------------------------------------------------------------------------
// NewQueryInfContent(): allocate memory of information content for query
//     positions
// -------------------------------------------------------------------------
inline
void AttributableScores::NewQueryInfContent( int size )
{
    if( size <= 0 ) {
        warning( "AttributableScores: Memory allocation size <=0 requested." );
        return;
    }
    DestroyQueryInfContent();

    queryinfcontent = ( double* )malloc( sizeof( double ) * size );

    if( !queryinfcontent ) 
        throw myruntime_error( mystring( "AttributableScores: Not enough memory." ));

    memset( queryinfcontent, 0, sizeof( double ) * size );
}

// -------------------------------------------------------------------------
// DestroyQueryInfContent(): free memory of information content for query pos.
// -------------------------------------------------------------------------

inline
void AttributableScores::DestroyQueryInfContent()
{
    if( queryinfcontent )
        free( queryinfcontent );
    queryinfcontent = NULL;
}

// -------------------------------------------------------------------------
// NewSbjctInfContent(): allocate memory of information content for subject
//     positions
// -------------------------------------------------------------------------
inline
void AttributableScores::NewSbjctInfContent( int size )
{
    if( size <= 0 ) {
        warning( "AttributableScores: Memory allocation size <=0 requested." );
        return;
    }
    DestroySbjctInfContent();

    sbjctinfcontent = ( double* )malloc( sizeof( double ) * size );

    if( !sbjctinfcontent )
        throw myruntime_error( mystring( "AttributableScores: Not enough memory." ));

    memset( sbjctinfcontent, 0, sizeof( double ) * size );
}

// -------------------------------------------------------------------------
// DestroySbjctInfContent(): free memory of information content for subject pos.
// -------------------------------------------------------------------------

inline
void AttributableScores::DestroySbjctInfContent()
{
    if( sbjctinfcontent )
        free( sbjctinfcontent );
    sbjctinfcontent = NULL;
}

// -------------------------------------------------------------------------
// GetQueryMinScoreAt/GetQueryMaxScoreAt: returns min/max scores at query
//     position n
//
inline
TScore AttributableScores::GetQueryMinScoreAt( int n ) const
{
    if( !queryminmaxscores || GetQuerySize() <= n || n < 0 )
        throw myruntime_error(
            mystring( "AttributableScores: GetQueryMinScoreAt: Memory access error." ));

    return queryminmaxscores[n+n];
}
inline
TScore AttributableScores::GetQueryMaxScoreAt( int n ) const
{
    if( !queryminmaxscores || GetQuerySize() <= n || n < 0 )
        throw myruntime_error(
            mystring( "AttributableScores: GetQueryMaxScoreAt: Memory access error." ));

    return queryminmaxscores[n+n+1];
}

// -------------------------------------------------------------------------
// GetSbjctMinScoreAt/GetSbjctMaxScoreAt: returns min/max scores at subject
//     position m
//
inline
TScore AttributableScores::GetSbjctMinScoreAt( int m ) const
{
    if( !sbjctminmaxscores || GetSubjectSize() <= m || m < 0 )
        throw myruntime_error(
            mystring( "AttributableScores: GetSbjctMinScoreAt: Memory access error." ));

    return sbjctminmaxscores[m+m];
}
inline
TScore AttributableScores::GetSbjctMaxScoreAt( int m ) const
{
    if( !sbjctminmaxscores || GetSubjectSize() <= m || m < 0 )
        throw myruntime_error(
            mystring( "AttributableScores: GetSbjctMaxScoreAt: Memory access error." ));

    return sbjctminmaxscores[m+m+1];
}

// -------------------------------------------------------------------------
// SetQueryInfMinMaxScoresAt: set min max scores at query position n
// -------------------------------------------------------------------------
inline
void AttributableScores::SetQueryInfMinMaxScoresAt( TScore min, TScore max, int n )
{
    if( max < min )
        throw myruntime_error( mystring( "AttributableScores: Wrong min/max scores." ));

    if( !queryminmaxscores || GetQuerySize() <= n || n < 0 )
        throw myruntime_error(
            mystring( "AttributableScores: SetQueryInfMinMaxScoresAt: Memory access error." ));

    NewQueryInfProbsAt( max - min + 1, n );

    queryminmaxscores[n+n] = min;
    queryminmaxscores[n+n+1] = max;
}

// -------------------------------------------------------------------------
// SetSbjctInfMinMaxScoresAt: set min max scores at subject position m
// -------------------------------------------------------------------------
inline
void AttributableScores::SetSbjctInfMinMaxScoresAt( TScore min, TScore max, int m )
{
    if( max < min )
        throw myruntime_error( mystring( "AttributableScores: Wrong min/max scores." ));

    if( !sbjctminmaxscores || GetSubjectSize() <= m || m < 0 )
        throw myruntime_error(
            mystring( "AttributableScores: SetSbjctInfMinMaxScoresAt: Memory access error." ));

    NewSbjctInfProbsAt( max - min + 1, m );

    sbjctminmaxscores[m+m] = min;
    sbjctminmaxscores[m+m+1] = max;
}

// -------------------------------------------------------------------------
// SetQueryInfProbNorm/SetSbjctInfProbNorm: set normalizing term of
//     probability at query/subject position n/m
// -------------------------------------------------------------------------
inline
void AttributableScores::SetQueryInfProbNorm( double value, int n )
{
#ifdef __DEBUG__
    if( !queryinfprobabs || GetQuerySize() <= n || n < 0 )
        throw myruntime_error( mystring( "AttributableScores: SetQueryInfProbNorm: Memory access error." ));
#endif
    if( !queryinfprobabs[n] )
        throw myruntime_error( mystring( "AttributableScores: SetQueryInfProbNorm: Memory access error." ));

    queryinfprobabs[n][0] = value;
}

inline
void AttributableScores::SetSbjctInfProbNorm( double value, int m )
{
#ifdef __DEBUG__
    if( !sbjctinfprobabs || GetSubjectSize() <= m || m < 0 )
        throw myruntime_error( mystring( "AttributableScores: SetSbjctInfProbNorm: Memory access error." ));
#endif
    if( !sbjctinfprobabs[m] )
        throw myruntime_error( mystring( "AttributableScores: SetSbjctInfProbNorm: Memory access error." ));

    sbjctinfprobabs[m][0] = value;
}

// -------------------------------------------------------------------------
// GetQueryInfProbNorm/GetSbjctInfProbNorm: get normalizing term of
//     probability at query/subject position n/m
// -------------------------------------------------------------------------
inline
double AttributableScores::GetQueryInfProbNorm( int n ) const
{
#ifdef __DEBUG__
    if( !queryinfprobabs || GetQuerySize() <= n || n < 0 )
        throw myruntime_error( mystring( "AttributableScores: GetQueryInfProbNorm: Memory access error." ));
#endif
    if( !queryinfprobabs[n] )
        throw myruntime_error( mystring( "AttributableScores: GetQueryInfProbNorm: Memory access error." ));

    return queryinfprobabs[n][0];
}

inline
double AttributableScores::GetSbjctInfProbNorm( int m ) const
{
#ifdef __DEBUG__
    if( !sbjctinfprobabs || GetSubjectSize() <= m || m < 0 )
        throw myruntime_error( mystring( "AttributableScores: GetSbjctInfProbNorm: Memory access error." ));
#endif
    if( !sbjctinfprobabs[m] )
        throw myruntime_error( mystring( "AttributableScores: GetSbjctInfProbNorm: Memory access error." ));

    return sbjctinfprobabs[m][0];
}


// -------------------------------------------------------------------------
// GetQueryInfProbOf: get probability of score at query position n
// -------------------------------------------------------------------------
inline
double AttributableScores::GetQueryInfProbOf( TScore score, int n ) const
{
#ifdef __DEBUG__
    if( !queryinfprobabs || GetQuerySize() <= n || n < 0 )
        throw myruntime_error( mystring( "AttributableScores: GetQueryInfProbOf: Memory access error." ));
#endif
    TScore  minscore = GetQueryMinScoreAt( n );

    if( !queryinfprobabs[n] || score - minscore < 0 )
        throw myruntime_error( mystring( "AttributableScores: GetQueryInfProbOf: Memory access error." ));

    return queryinfprobabs[n][ score - minscore + 1 ];//NOTE: position 0 is reserved for norm. term
}

// -------------------------------------------------------------------------
// GetSbjctInfProbOf: get probability of score at subject position m
// -------------------------------------------------------------------------
inline
double AttributableScores::GetSbjctInfProbOf( TScore score, int m ) const
{
#ifdef __DEBUG__
    if( !sbjctinfprobabs || GetSubjectSize() <= m || m < 0 )
        throw myruntime_error( mystring( "AttributableScores: GetSbjctInfProbOf: Memory access error." ));
#endif
    TScore  minscore = GetSbjctMinScoreAt( m );

    if( !sbjctinfprobabs[m] || score - minscore < 0 )
        throw myruntime_error( mystring( "AttributableScores: GetSbjctInfProbOf: Memory access error." ));

    return sbjctinfprobabs[m][ score - minscore + 1 ];
}

// -------------------------------------------------------------------------
// IncQueryInfProbOf: increment probability of score by value at query
//     position n
// -------------------------------------------------------------------------
inline
double AttributableScores::IncQueryInfProbOf( double value, TScore score, int n )
{
#ifdef __DEBUG__
    if( !queryinfprobabs || GetQuerySize() <= n || n < 0 )
        throw myruntime_error( mystring( "AttributableScores: IncQueryInfProbOf: Memory access error." ));
#endif
    TScore  minscore = GetQueryMinScoreAt( n );

    if( !queryinfprobabs[n] || score - minscore < 0 )
        throw myruntime_error( mystring( "AttributableScores: IncQueryInfProbOf: Memory access error." ));

    queryinfprobabs[n][ score - minscore + 1 ] += value;
    return queryinfprobabs[n][ score - minscore + 1 ];
}

// -------------------------------------------------------------------------
// IncSbjctInfProbOf: increment probability of score by value at subject
//     position m
// -------------------------------------------------------------------------
inline
double AttributableScores::IncSbjctInfProbOf( double value, TScore score, int m )
{
#ifdef __DEBUG__
    if( !sbjctinfprobabs || GetSubjectSize() <= m || m < 0 )
        throw myruntime_error( mystring( "AttributableScores: IncSbjctInfProbOf: Memory access error." ));
#endif
    TScore  minscore = GetSbjctMinScoreAt( m );

    if( !sbjctinfprobabs[m] || score - minscore < 0 )
        throw myruntime_error( mystring( "AttributableScores: IncSbjctInfProbOf: Memory access error." ));

    sbjctinfprobabs[m][ score - minscore + 1 ] += value;
    return sbjctinfprobabs[m][ score - minscore + 1 ];
}

// -------------------------------------------------------------------------
// DivQueryInfProbOf: divide probability of score at query position n
//     by value
// -------------------------------------------------------------------------

inline
double AttributableScores::DivQueryInfProbOf( double value, TScore score, int n )
{
#ifdef __DEBUG__
    if( !queryinfprobabs || GetQuerySize() <= n || n < 0 )
        throw myruntime_error( mystring( "AttributableScores: DivQueryInfProbOf: Memory access error." ));
#endif
    TScore  minscore = GetQueryMinScoreAt( n );

    if( !queryinfprobabs[n] || score - minscore < 0 )
        throw myruntime_error( mystring( "AttributableScores: DivQueryInfProbOf: Memory access error." ));

    if( value == 0.0 )
        throw myruntime_error( mystring( "AttributableScores: DivQueryInfProbOf: Illegal operation." ));

    queryinfprobabs[n][ score - minscore + 1 ] /= value;
    return queryinfprobabs[n][ score - minscore + 1 ];
}

// -------------------------------------------------------------------------
// DivSbjctInfProbOf: divide probability of score at subject position m
//     by value
// -------------------------------------------------------------------------

inline
double AttributableScores::DivSbjctInfProbOf( double value, TScore score, int m )
{
#ifdef __DEBUG__
    if( !sbjctinfprobabs || GetSubjectSize() <= m || m < 0 )
        throw myruntime_error( mystring( "AttributableScores: DivSbjctInfProbOf: Memory access error." ));
#endif
    TScore  minscore = GetSbjctMinScoreAt( m );

    if( !sbjctinfprobabs[m] || score - minscore < 0 )
        throw myruntime_error( mystring( "AttributableScores: DivSbjctInfProbOf: Memory access error." ));

    if( value == 0.0 )
        throw myruntime_error( mystring( "AttributableScores: DivSbjctInfProbOf: Illegal operation." ));

    sbjctinfprobabs[m][ score - minscore + 1 ] /= value;
    return sbjctinfprobabs[m][ score - minscore + 1 ];
}



// -------------------------------------------------------------------------
// NewQueryMinMaxScores: allocate memory for min/max scores of query
//     positions
// -------------------------------------------------------------------------
inline
void AttributableScores::NewQueryMinMaxScores()
{
    if( GetQuerySize() <= 0 ) {
        warning( "AttributableScores: NewQueryMinMaxScores: Wrong size." );
        return;
    }
    DestroyQueryMinMaxScores();

    queryminmaxscores = ( int* )malloc( sizeof( int ) * GetQuerySize()* 2 );

    if( !queryminmaxscores ) 
        throw myruntime_error( mystring( "AttributableScores: Not enough memory." ));

    memset( queryminmaxscores, 0, sizeof( int ) * GetQuerySize()* 2 );
}

// -------------------------------------------------------------------------
// DestroyQueryMinMaxScores: free memory of min/max scores for query pos.
// -------------------------------------------------------------------------
inline
void AttributableScores::DestroyQueryMinMaxScores()
{
    if( queryminmaxscores )
        free( queryminmaxscores );
    queryminmaxscores = NULL;
}

// -------------------------------------------------------------------------
// NewSbjctMinMaxScores: allocate memory for min/max scores of subject
//     positions
// -------------------------------------------------------------------------
inline
void AttributableScores::NewSbjctMinMaxScores()
{
    if( GetSubjectSize() <= 0 ) {
        warning( "AttributableScores: NewSbjctMinMaxScores: Wrong size." );
        return;
    }
    DestroySbjctMinMaxScores();

    sbjctminmaxscores = ( int* )malloc( sizeof( int ) * GetSubjectSize()* 2 );

    if( !sbjctminmaxscores ) 
        throw myruntime_error( mystring( "AttributableScores: Not enough memory." ));

    memset( sbjctminmaxscores, 0, sizeof( int ) * GetSubjectSize()* 2 );
}

// -------------------------------------------------------------------------
// DestroySbjctMinMaxScores: free memory of min/max scores for subject pos.
// -------------------------------------------------------------------------
inline
void AttributableScores::DestroySbjctMinMaxScores()
{
    if( sbjctminmaxscores )
        free( sbjctminmaxscores );
    sbjctminmaxscores = NULL;
}




// -------------------------------------------------------------------------
// NewProbabilities(): allocate memory for probabilities
// -------------------------------------------------------------------------

inline
void AttributableScores::NewProbabilities( size_t count )
{
    if( probabilities )
        DestroyProbabilities();

    probabilities = ( double* )malloc( sizeof( double ) * count );

    if( !probabilities ) 
        throw myruntime_error( mystring( "AttributableScores: Not enough memory." ));

    memset( probabilities, 0, sizeof( double ) * count );
}

// -------------------------------------------------------------------------
// DestroyProbabilities(): deallocate memory previously allocated for
//     probabilities
// -------------------------------------------------------------------------

inline
void AttributableScores::DestroyProbabilities()
{
    if( probabilities )
        free( probabilities );
    probabilities = NULL;
}

// -------------------------------------------------------------------------
// GetProbabilityOf(): get score probability given score value
// -------------------------------------------------------------------------

inline
double AttributableScores::GetProbabilityOf( TScore score ) const
{
#ifdef __DEBUG__
    if( !probabilities )
        throw myruntime_error( mystring( "AttributableScores: Memory access error." ));

    if( score - min_score < 0 )
        throw myruntime_error( mystring( "AttributableScores: Memory access error." ));
#endif
    return probabilities[score-min_score];
}

// -------------------------------------------------------------------------
// SetProbabilityAt(): set probability of score
// -------------------------------------------------------------------------

inline
void AttributableScores::SetProbabilityOf( double value, TScore score )
{
#ifdef __DEBUG__
    if( !probabilities )
        throw myruntime_error( mystring( "AttributableScores: Memory access error." ));

    if( score - min_score < 0 )
        throw myruntime_error( mystring( "AttributableScores: Memory access error." ));
#endif
    probabilities[score-min_score] = value;
}

// -------------------------------------------------------------------------
// IncProbabilityAt(): increment score probability with the value given
// -------------------------------------------------------------------------

inline
void AttributableScores::IncProbabilityOf( double value, TScore score )
{
#ifdef __DEBUG__
    if( !probabilities )
        throw myruntime_error( mystring( "AttributableScores: Memory access error." ));

    if( score - min_score < 0 )
        throw myruntime_error( mystring( "AttributableScores: Memory access error." ));
#endif
    probabilities[score-min_score] += value;
}

// -------------------------------------------------------------------------
// DivideProbabilityAt(): divide score probability by the value given
// -------------------------------------------------------------------------

inline
void AttributableScores::DivideProbabilityOf( double value, TScore score )
{
#ifdef __DEBUG__
    if( !probabilities )
        throw myruntime_error( mystring( "AttributableScores: Memory access error." ));

    if( score - min_score < 0 )
        throw myruntime_error( mystring( "AttributableScores: Memory access error." ));
#endif
    if( value == 0.0 )
        throw myruntime_error( mystring( "AttributableScores: Illegal operation." ));

    if( probabilities[score-min_score])
        probabilities[score-min_score] /= value;
}

// -------------------------------------------------------------------------
// InitializePrivateProbVector(): initialize private probability vector
// -------------------------------------------------------------------------

inline
void AttributableScores::NewPrivateProbVector( size_t maxit, size_t range )
{
    if( priv_prob_vector )
        DestroyPrivateProbVector();

    size_t  size = maxit * range + 1;   //plus to include the zero value

    priv_prob_vector = ( double* )malloc( sizeof( double ) * size );

    if( !priv_prob_vector )
        throw myruntime_error( mystring( "AttributableScores: Not enough memory." ));

    prob_vector_size = size;
    ResetPrivateProbVector();
}

// -------------------------------------------------------------------------
// DestroyPrivateProbVector(): deallocate memory previously allocated for
//     private probability vector
// -------------------------------------------------------------------------

inline
void AttributableScores::DestroyPrivateProbVector()
{
    if( priv_prob_vector )
        free( priv_prob_vector );
    priv_prob_vector = NULL;
    prob_vector_size = 0;
}

// -------------------------------------------------------------------------
// ResetPrivateProbVector(): reset values of private probability vector
// -------------------------------------------------------------------------

inline
void AttributableScores::ResetPrivateProbVector( size_t count )
{
#ifdef __DEBUG__
    if( !priv_prob_vector )
        throw myruntime_error( mystring( "AttributableScores: Memory access error." ));

    if( prob_vector_size < count ) {
        warning( "AttributableScores: Wrong argument." );
        count = prob_vector_size;
    }
#endif
    if( count == 0 )
        count = prob_vector_size;
    memset( priv_prob_vector, 0, sizeof( double ) * count );
}


#endif//__AttributableScores__
