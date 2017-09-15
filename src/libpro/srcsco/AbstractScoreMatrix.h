/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/



#ifndef __AbstractScoreMatrix__
#define __AbstractScoreMatrix__

#include <stdlib.h>
#include <math.h>

#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "data.h"
#include "mystring.h"
#include "myexcept.h"

#include "libpro/srcpro/Configuration.h"
#include "libpro/srcpro/GapScheme.h"
#include "libpro/srcpro/DistributionMatrix.h"

//Types for masking of scores
enum TMask {
        Unmasked = 0,       //unmasked score
        MaskToIgnore = 1,   //masked score with flag of excluding it from the comp. statistics
        MaskToConsider = 2  //masked score with having no effect on comp. statistics but the alignment procedure
};

#include "AttributableScores.h"


// extern double rint( double x );

////////////////////////////////////////////////////////////////////////////
// CLASS AbstractScoreMatrix
// This is a base class for various implementations of Score Matrices
//
class AbstractScoreMatrix
{
public:
    enum TType {
            PositionSpecific,
            ProfileSpecific,
            ProfSpecLSO,
            HDPProfileSpecific,
            HDP0CtxProfileSpecific,
            Universal,
            ParallelUniversal,
            ParallelLSOUniversal,
            ParallelHDPUniversal,
            AdjustedProfileSpecific,
            Dummy
    };

    enum TBehaviour {
            StatisticsGiven,
            ComputeStatistics
    };

    enum TScaling {
            NoScaling = 0,
            AutoScalling = 1,
            FPScaling = 2
    };

public:
    AbstractScoreMatrix( TType, Configuration config[NoSchemes], TBehaviour, TScaling a_scaling, TMask maskapp );
    virtual ~AbstractScoreMatrix();

    TType               GetType() const                 { return matrixType; }
    TBehaviour          GetBehaviour() const            { return behaviour; }
    TMask               GetMaskingApproach() const      { return commonmasking; }

    virtual const char* GetMethodName() const = 0;

    const char*         GetName() const { return name_; }
    void                SetName( const char* value ) { name_ = value; }

    double              GetDeletionCoefficient() const          { return deletioncoeff; }
    void                SetDeletionCoefficient( double value )  { deletioncoeff = value; }

    double              GetInformationThreshold() const { return information_content; }
    void                SetInfoThresholdByEval( double eval );

    bool                GetAutocorrectionPositional() const             { return  positcorrects; }
    void                SetAutocorrectionPositional( bool value )       { positcorrects = value; }


    double              GetInfoCorrectionUpperBound2nd() const          { return ainfoupper2nd; }
    void                SetInfoCorrectionUpperBound2nd( double value )  { ainfoupper2nd = value; }

    double              GetInfoCorrectionNumerator2nd() const           { return ainfonumerator2nd; }
    void                SetInfoCorrectionNumerator2nd( double value )   { ainfonumerator2nd = value; }

    double              GetInfoCorrectionScale2nd() const               { return ainfoscale2nd; }
    void                SetInfoCorrectionScale2nd( double value )       { ainfoscale2nd = value; }

    double              GetInfoCorrectionNumeratorAlt() const           { return ainfonumeratoralt; }
    void                SetInfoCorrectionNumeratorAlt( double value )   { ainfonumeratoralt = value; }

    double              GetInfoCorrectionScaleAlt() const               { return ainfoscalealt; }
    void                SetInfoCorrectionScaleAlt( double value )       { ainfoscalealt = value; }

    bool                GetSupportOptimFreq() const                     { return supportoptfreq_; }
    virtual void        OptimizeTargetFrequencies()                     { return; }

    bool                GetScoreAdjRelEnt() const { return scadjrelent_; }
    void                SetScoreAdjRelEnt( bool value ) { scadjrelent_ = value; }

    bool                GetScoreAdjCorrel() const { return scadjcorrsc_; }
    void                SetScoreAdjCorrel( bool value ) { scadjcorrsc_ = value; }


//     double              operator()( int m, int n ) const{ return GetScore( m, n ); }
    operator            bool() const                    { return IsValid();  }

    bool                GetFPScaling() const            { return auto_scaling == FPScaling; }
    bool                GetAutoScaling() const          { return auto_scaling == AutoScalling; }
    int                 GetAutoScalingFactor() const    { return auto_scaling_constant; }

    bool                IsValid() const;                            //whether score matrix is valid
                                                                    //returns score at specified profile positions
    virtual double      GetScore( int m, int n ) const  { return GetScoreBase( m, n ); }
    double              GetImageScore( int m, int n ) const;        //returns image score at specified profile positions
    double              GetModImageScore( int m, int n ) const;
    double              GetModScore( int m, int n ) const;

    bool                GetUseModImage() const { return usemodimage_; }
    void                SetUseModImage( bool value ) { usemodimage_ = value; }

    TMask               GetMasked( int m, int n  ) const;           //get mask value at

    virtual int         GetQuerySize() const            { return GetPrivateQuerySize(); }
    virtual int         GetSubjectSize() const          { return GetPrivateSecndSize(); }
                                                                    //perform search of HSPs
    bool                ScanForHSPs( double minhspscore, int hsplen, int nohsps, int maxdist, int* = NULL, int* = NULL );

    double              GetRefLambda() const            { return referenceLambda; }
    double              GetRefH() const                 { return referenceH; }
    double              GetRefK() const                 { return referenceK; }
    double              GetExpGappedLambda() const      { return experimentalGappedLambda; }
    double              GetExpGappedK() const           { return experimentalGappedK; }
    double              GetDerivedGappedLambda() const  { return derivedGappedLambda; }
    double              GetDerivedGappedK() const       { return derivedGappedK; }

    double              GetLambda() const               { return lambda; }      //computed ungapped Lambda
    double              GetEntropy() const              { return entropy; }     //computed relative entropy
    double              GetK() const                    { return parameterK; }  //computed ungapped parameter K
    double              GetExpectedScore() const        { return expscore; }    //expected score

    double              GetPrelimScore() const          { return prescore; }    //preliminary score
    double              GetPrelimExpect() const         { return prexpect; }    //expectation value of preliminary score

    TScore              GetMinScore() const             { return min_score; }   //get min score
    TScore              GetMaxScore() const             { return max_score; }   //get max score

    double              GetMultiplier() const           { return score_multiplier; }

    virtual void        ComputeProfileScoringMatrix( bool final = false ) = 0;

    void                PostScalingProc(
                                const LogOddsMatrix& qlogo, const LogOddsMatrix& slogo,
                                GapScheme& qg, GapScheme& sg,
                                bool autoc, int acwindow, bool firstpass );
    int                 AdjustGapCost( int cost ) const;
    double              GetFinalScore( double value ) const;

    virtual void        ScaleScoringMatrix();
    void                ComputeStatisticalParameters();             //compute all statistical parameters that are required
    virtual double      ComputeExpectation( double, double* = NULL, double* = NULL, double* = NULL, double* = NULL ) const;
    double              GetMeanLengthGivenExpectation( double, double* ) const;

    void                PrintParameterTable( char* ) const;         //print table of computed statistical parameter values
    void                PrintParameterTable( FILE* ) const;         //print table of computed statistical parameter values
    virtual void        PrintParameterTable( TPrintFunction, void* vpn ) const = 0;
    void                PrintFinal( char* ) const;                  //final print if needed
    void                PrintFinal( FILE* ) const;                  //final print if needed
    virtual void        PrintFinal( TPrintFunction, void* vpn ) const = 0;

    virtual void        PrintScoringMatrix( FILE* ) = 0;      //print the computed scoring system

    //STATIC public methods...
    static  bool        ComputeLengthAdjustment(                    //compute length adjustment to correct edge-effects
                const Configuration& , size_t query_len, Uint64 db_len, size_t no_sequences );

    static  size_t      GetDeltaLength() { return deltaLength; }    //get the length adjustment value
    static  Uint64      GetSearchSpace() { return searchSpace; }    //get the search space computed
    static  Uint64      GetRawSearchSpace() { return raw_search_space; }

                                                                    //print table of statistical parameter values...
    void                PrintReferenceParameterTable( char* ) const;
    void                PrintReferenceParameterTable( FILE* ) const;
    void                PrintReferenceParameterTable( TPrintFunction, void* vpn ) const;


                                                        //compute positional probabilities of scores
    virtual void        ComputePositionalScoreProbs( AttributableScores* ) = 0;

    const double* const         GetQueryInfContent() const;         //information contents corresponding query positions
    const double* const         GetSbjctInfContent() const;         //information contents corresponding subject positions

protected:
    explicit AbstractScoreMatrix();

    void                Init( int sz_query, int sz_sbjct );
    void                InitializeSSParameters();

    virtual double      GetScoreBase( int m, int n ) const;         //returns score at specified profile positions

    void                SetScore( double, int m, int n );
    void                SetModScore( double, int m, int n );

    int                 GetPrivateQuerySize() const                 { return querySize;   }
    int                 GetPrivateSecndSize() const                 { return subjectSize; }

    void                SetInformationThreshold( double value )     { information_content = value; }

    bool                GetMaskedUnmasked( int m, int n  ) const    { return GetMasked( m, n  ) == Unmasked; }
    bool                GetMaskedToIgnore( int m, int n  ) const    { return GetMasked( m, n  ) == MaskToIgnore; }
    bool                GetMaskedToConsider( int m, int n  ) const  { return GetMasked( m, n  ) == MaskToConsider; }

    void                SetMasked( TMask, int m, int n  );          //set mask value at

    const Configuration&    GetConfiguration( TProcomSchemes ps ) const;

    virtual void        MultiplierProgress( AttributableScores* );  //method used to show progress of the values of score_multiplier

    bool                GetAllNegatives() const                     { return allnegatives; }
    void                SetAllNegatives( bool value );

    void                SetPrelimScore( double value )              { prescore = value; }   //set preliminary score
    void                SetPrelimExpect( double value )             { prexpect = value; }   //set expect of prelimianry score

    void                SetMinScore( TScore value )                 { min_score = value; }  //set min score
    void                SetMaxScore( TScore value )                 { max_score = value; }  //set max score

    void                SetMultiplier( double value )               { score_multiplier = value; }

                                                        //compute probabilities of scores at each position
    virtual void        ComputeScoreProbabilities( AttributableScores* ) = 0;
    virtual bool        ReduceScores( AttributableScores* ) { return false; }

    void                SetSupportOptimFreq( bool value )           { supportoptfreq_ = value; }

private:
    void                SetRefLambda( double value )                { referenceLambda = value; }
    void                SetRefH( double value )                     { referenceH = value; }
    void                SetRefK( double value )                     { referenceK = value; }
public:
    void                SetExpGappedLambda( double value )          { experimentalGappedLambda = value; }
    void                SetExpGappedK( double value )               { experimentalGappedK = value; }
private:
    void                SetDerivedGappedLambda( double value )      { derivedGappedLambda = value; }
    void                SetDerivedGappedK( double value )           { derivedGappedK = value; }

    void                SetLambda( double newlambda )               { lambda = newlambda; }
    void                SetEntropy( double H )                      { entropy = H; }
    void                SetK( double K )                            { parameterK = K; }
    void                SetExpectedScore( double E )                { expscore = E; }

    void                TransformScoresToAchieveH();                //transform scores to acquire desired relative entropy
    void                TransformScoresInMatrix( double c );        //transform scores in the matrix given constant c for relative entropy

    void                DeriveGappedParameters();                   //evaluate estimated gapped parameters Lambda and K


    //STATIC protected methods...
                                                //obtain statistical reference parameters given gap cost scheme
    static  void        ObtainReferenceParameters( int ss,
            double* l, double* k, double* h, double* a, double* b,
            double* L = NULL, double* K = NULL, double* H = NULL, double* A = NULL, double* B = NULL );

    static  bool        ComputeLengthAdjustment(    //compute length adjustment to correct edge-effects
                double lambda, double K, double alpha, double beta,
                size_t query_len, Uint64 db_len, size_t no_sequences );

    static  void        SetDeltaLength( size_t value ) { deltaLength = value; } //save the length adjustment value
    static  void        SetSearchSpace( Uint64 value ) { searchSpace = value; } //set the search space value
    static  void        SetRawSearchSpace( Uint64 value )   { raw_search_space = value; }

protected:
    const AttributableScores*   GetCorresScores() const             { return corres_scores; }
    AttributableScores*         GetCorresScores()                   { return corres_scores; }
    const AttributableScores*   GetScaledScores() const             { return scaled_scores; }
    AttributableScores*         GetScaledScores()                   { return scaled_scores; }

private:
    const double**              GetImage() const                    { return ( const double** )image; }
    const TMask**               GetMask() const                     { return ( const TMask** )mask; }

    void                        SetQuerySize( int qsize )           { querySize = qsize;   }
    void                        SetSubjectSize( int rsize )         { subjectSize = rsize; }

private:
    const TType         matrixType;                 //type of scoring method and matrix
    const TBehaviour    behaviour;                  //behaviour about computations to perform
    const TMask         commonmasking;              //common masking approach
    const char*         name_;

    double**            image;                      //private matrix of double scores
    double**            modimage_;                  //modified matrix of double scores
    bool                usemodimage_;               //whether using modified matrix
    TMask**             mask;           //masking of scores; some of masked scores should not be included in commp. statistics
    AttributableScores* corres_scores;  //matrix of scores corresponding to image (computed by processing two profiles) (representative)
    AttributableScores* scaled_scores;  //scores multiplied by a factor to increase precision of calculations (for computations)

    double              score_multiplier;           //the last value of multiplier used in scaling of scores

    double              deletioncoeff;              //deletion probability weight

    double              information_content;        //information content threshold
    double              ainfoupper2nd;              //upper bound of information content threshold used in 2nd-pass computations
    double              ainfonumerator2nd;          //numerator of expression to compute 2nd-pass information content threshold
    double              ainfoscale2nd;              //logarithmic scale to compute 2nd-pass inf. content threshold
    double              ainfonumeratoralt;          //numerator of alternative expression to compute inf. content threshold
    double              ainfoscalealt;              //logarithmic scale to alternatively compute inf. content threshold
    bool                positcorrects;              //whether corrections computed positionally by entropies are in use

    bool                scadjrelent_;               //relative entropy score adjustment
    bool                scadjcorrsc_;               //correlated score adjustment

    double              prescore;                   //preliminary obtained score
    double              prexpect;                   //expectation value corresponding preliminary obtained score

    TScore              min_score;                  //minimum score value; just for information
    TScore              max_score;                  //maximum score value; just for information

    TScaling            auto_scaling;               //whether to use auto scaling to increase precision of scores
    const int           auto_scaling_constant;      //constant factor used to scale values in order to gain precision

    double              expscore;                   //expected score per column pair
    bool                allnegatives;               //whether scores are all negative

    double              referenceLambda;            //reference lambda parameter
    double              referenceH;                 //reference parameter H
    double              referenceK;                 //reference parameter K
    double              experimentalGappedLambda;   //experimental gapped lambda
    double              experimentalGappedK;        //experimental gapped K
    double              derivedGappedLambda;        //gapped lambda derived for this scoring matrix
    double              derivedGappedK;             //gapped K derived for this scoring matrix

    double              lambda;                     //analitically computed ungapped Lambda
    double              entropy;                    //computed relative entropy
    double              parameterK;                 //analitically computed ungapped parameter K

    int                 querySize;                  //length of query sequence (profile)
    int                 subjectSize;                //length of subject sequence (profile)

    bool                supportoptfreq_;            //supports optimization of target frequencies

    const Configuration*    configuration;          //parameter configuration

    static  size_t      deltaLength;                //length adjustment for edge-effect correction
    static  Uint64      searchSpace;                //search space computed taking into account deltaLength
    static  Uint64      raw_search_space;           //raw search space without any corrections
};

// INLINES ...

// -------------------------------------------------------------------------
// IsValid: whether the class object is valid
// -------------------------------------------------------------------------

inline bool AbstractScoreMatrix::IsValid() const
{
    return  GetCorresScores() != NULL && (
            GetScaledScores() != NULL || ! GetAutoScaling()
        );
}

// -------------------------------------------------------------------------
// GetScoreBase: used to access score at the profile positions
// -------------------------------------------------------------------------

inline
double AbstractScoreMatrix::GetScoreBase( int m, int n ) const
{
#ifdef __DEBUG__
    if( GetCorresScores() == NULL || ( GetAutoScaling() && GetScaledScores() == NULL ))
        throw myruntime_error(
            mystring( "AbstractScoreMatrix: Memory access error." ));
#endif

    if( GetAutoScaling())
        return GetScaledScores()->GetScore( m, n );

    return GetCorresScores()->GetScore( m, n );
}

// -------------------------------------------------------------------------
// GetImageScore: image score at the profile positions
// -------------------------------------------------------------------------

inline
double AbstractScoreMatrix::GetImageScore( int m, int n ) const
{
#ifdef __DEBUG__
    if( !image || subjectSize <= m || m < 0 || querySize <= n || n < 0 )
        throw myruntime_error(
            mystring( "AbstractScoreMatrix: Memory access error." ));
#endif

    return image[m][n];
}

// -------------------------------------------------------------------------
// GetModImageScore: modified image score at the profile positions
//
inline
double AbstractScoreMatrix::GetModImageScore( int m, int n ) const
{
#ifdef __DEBUG__
    if( !modimage_ || subjectSize <= m || m < 0 || querySize <= n || n < 0 )
        throw myruntime_error("AbstractScoreMatrix: Memory access error.");
#endif
    return modimage_[m][n];
}

// -------------------------------------------------------------------------
// GetModScore: get scaled modified score at the profile positions
//
inline
double AbstractScoreMatrix::GetModScore( int m, int n ) const
{
    int sf = 1;
    double  score = GetModImageScore( m, n );
#ifdef __DEBUG__
    if( GetCorresScores() == NULL || ( GetAutoScaling() && GetScaledScores() == NULL ))
        throw myruntime_error("AbstractScoreMatrix: Memory access error.");
#endif
    if( GetAutoScaling())
        sf = GetScaledScores()->GetAutoScalingFactor();
    else
        sf = GetCorresScores()->GetAutoScalingFactor();

    if( SCORE_MIN < score && sf != 1 )
        score *= ( double )sf;
    return score;
}

// -------------------------------------------------------------------------
// GetMasked: get score mask at the position of the image
// -------------------------------------------------------------------------

inline
TMask AbstractScoreMatrix::GetMasked( int m, int n  ) const
{
#ifdef __DEBUG__
    if( !mask || subjectSize <= m || m < 0 || querySize <= n || n < 0 )
        throw myruntime_error(
            mystring( "AbstractScoreMatrix: Memory access error." ));
#endif

    return mask[m][n];
}

// -------------------------------------------------------------------------
// GetQueryInfContent/GetSbjctInfContent: returns information content
//     vector of query/subject positions
//
inline
const double* const AbstractScoreMatrix::GetQueryInfContent() const
{
#ifdef __DEBUG__
    if( GetCorresScores() == NULL || ( GetAutoScaling() && GetScaledScores() == NULL ))
        throw myruntime_error(
            mystring( "AbstractScoreMatrix: GetQueryInfContent: Memory access error." ));
#endif
    if( GetAutoScaling())
        return GetScaledScores()->GetQueryInfContent();
    return GetCorresScores()->GetQueryInfContent();
}

inline
const double* const AbstractScoreMatrix::GetSbjctInfContent() const
{
#ifdef __DEBUG__
    if( GetCorresScores() == NULL || ( GetAutoScaling() && GetScaledScores() == NULL ))
        throw myruntime_error(
            mystring( "AbstractScoreMatrix: GetQueryInfContent: Memory access error." ));
#endif
    if( GetAutoScaling())
        return GetScaledScores()->GetSbjctInfContent();
    return GetCorresScores()->GetSbjctInfContent();
}

// -------------------------------------------------------------------------
// PostScalingProc: performs post scaling procedure
// -------------------------------------------------------------------------
inline
void AbstractScoreMatrix::PostScalingProc(
    const LogOddsMatrix&    querylogo,
    const LogOddsMatrix&    sbjctlogo,
    GapScheme&              querygaps,
    GapScheme&              sbjctgaps,
    bool    autoc,
    int     acwindow,
    bool    firstpass )
{
#ifdef __DEBUG__
    if( GetCorresScores() == NULL || ( GetAutoScaling() && GetScaledScores() == NULL ))
        throw myruntime_error(
            mystring( "AbstractScoreMatrix: Memory access error." ));
#endif

    if( firstpass ) {
        if( GetAutocorrectionPositional()) {
            if( GetAutoScaling())
                GetScaledScores()->ComputePositionalInfContents();
            else
                if( autoc )
                    GetCorresScores()->ComputePositionalInfContents();
        }
    } else {
        ComputeProfileScoringMatrix( true /*final*/ );
    }

    if( GetAutoScaling())
        GetScaledScores()->AdjustGaps( querylogo, sbjctlogo, querygaps, sbjctgaps, autoc, acwindow );
    else
        if( autoc )
            GetCorresScores()->AdjustGaps( querylogo, sbjctlogo, querygaps, sbjctgaps, autoc, acwindow );
}

// -------------------------------------------------------------------------
// AdjustGapCost: adjusts gap cost so that it is mulitplied by a factor when
//     scaling is used to gain precision
//
inline
int AbstractScoreMatrix::AdjustGapCost( int cost ) const
{
    int cost__ = cost;

    if( GetAutoScaling())
        cost__ *= GetAutoScalingFactor();

    return cost__;
}
// -------------------------------------------------------------------------
// GetFinalScore: adjusts final alignment score so that it is divided by a
//     factor when scaling is used to gain precision
//
inline
double AbstractScoreMatrix::GetFinalScore( double value ) const
{
    double  val = value;

    if( val && GetAutoScaling())
        val /= GetAutoScalingFactor();

    return val;
}

// -------------------------------------------------------------------------
// SetImageScore: change the image at the specified location of the matrix
// -------------------------------------------------------------------------

inline
void AbstractScoreMatrix::SetScore( double value, int m, int n )
{
#ifdef __DEBUG__
    if( !image || subjectSize <= m || m < 0 || querySize <= n || n < 0 )
        throw myruntime_error( mystring( "AbstractScoreMatrix: Memory access error." ));
#endif
    image[m][n] = value;

    if( GetCorresScores()) GetCorresScores()->SetScore( value, m, n );
    if( GetScaledScores()) GetScaledScores()->SetScore( value, m, n );
}

// -------------------------------------------------------------------------
// SetModScore: set modified score at the specified position
//
inline
void AbstractScoreMatrix::SetModScore( double value, int m, int n )
{
#ifdef __DEBUG__
    if( !modimage_ || subjectSize <= m || m < 0 || querySize <= n || n < 0 )
        throw myruntime_error( mystring( "AbstractScoreMatrix: Memory access error." ));
#endif
    modimage_[m][n] = value;
}

// -------------------------------------------------------------------------
// SetMasked: sets score mask at the position of the image
// -------------------------------------------------------------------------

inline
void AbstractScoreMatrix::SetMasked( TMask value, int m, int n  )
{
#ifdef __DEBUG__
    if( !mask || subjectSize <= m || m < 0 || querySize <= n || n < 0 )
        throw myruntime_error( mystring( "AbstractScoreMatrix: Memory access error." ));
#endif

    mask[m][n] = value;
}

// -------------------------------------------------------------------------
// SetAllNegatives: set flag indicating the presence of all negative scores
// -------------------------------------------------------------------------

inline
void AbstractScoreMatrix::SetAllNegatives( bool value )
{
    if( GetCorresScores()) GetCorresScores()->SetAllNegatives( value );
    if( GetScaledScores()) GetScaledScores()->SetAllNegatives( value );
    allnegatives = value;
}

// -------------------------------------------------------------------------
// GetConfiguration: returns reference to parameter configuration object
// -------------------------------------------------------------------------

inline
const Configuration& AbstractScoreMatrix::GetConfiguration( TProcomSchemes ps ) const
{
    if( !configuration )
        throw myruntime_error( mystring( "AbstractScoreMatrix: Configuration is null." ));

    if( NoSchemes <= ps )
        throw myruntime_error( mystring( "AbstractScoreMatrix: Wrong argument while accessing configuration." ));

    return configuration[ps];
}


#endif//__AbstractScoreMatrix__
