/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __ProfileSearching__
#define __ProfileSearching__

#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "rc.h"
#include "data.h"
#include "mystring.h"
#include "myexcept.h"

#include "libHDP/HDPbase.h"
#include "libpro/srcpro/Configuration.h"
#include "libpro/srcpro/Database.h"
#include "libpro/srcpro/GapScheme.h"
#include "libpro/srcpro/DistributionMatrix.h"
#include "libpro/srcpro/SEGProfile.h"
#include "libpro/srcsco/ScoringMatrix.h"
#include "libpro/srcsco/LSOSMatrix.h"
#include "libpro/srcsco/HDP0ScoreMatrix.h"
#include "libpro/srcsco/HDP0CtxScoreMatrix.h"
#include "libpro/srcsco/AdjustedScoreMatrix.h"
#include "libpro/srcsco/UniversalScoreMatrix.h"
#include "libpro/srcaln/ProfileAlignment.h"
#include "libpro/srcaln/HybridAlignment.h"

// _________________________________________________________________________
// Class HitInformation
//

class HitInformation {
public:
    HitInformation( double sc, double eval, double refeval, char* annot, char* fullinfo, bool destroy = true );
    ~HitInformation();

    bool        operator<( const HitInformation& ) const;

    double      GetScore() const        { return score; }
    double      GetEvalue() const       { return evalue; }
    double      GetRefEvalue() const    { return ref_evalue; }

    const char* GetAnnotation() const   { return annotation; }
    const char* GetFullAlignment() const{ return fullalignment; }

//     void    SetScore( double value )    { score = value; }
//     void    SetEvalue( double value )   { evalue = value; }

protected:
    explicit HitInformation();

private:
    double  score;          //raw alignment score
    double  evalue;         //e-value of the score
    double  ref_evalue;     //reference e-value with respect to which comparisons are accomplished
    char*   annotation;     //short (up to 80 ch.) annotation of the hit
    char*   fullalignment;  //alignment and additional information
    bool    autodestroy;    //whether to deallocate memory used by the class members
};


// _________________________________________________________________________
// Class ProfileSearching
//

class ProfileSearching
{
public:
    ProfileSearching(
            const char* paramconfigfile,
            const char* input,
            const char* database,
            const char* output,
            double eval_threshold,
            int no_hits,
            int no_alns,
            double ident_threshold,
            double infrm_threshold,
            double mask_percnt,
            int gapopen,
            int gapextend,
            bool showpars,
            bool fixedcosts,
            AbstractScoreMatrix::TType,
            AbstractScoreMatrix::TBehaviour,
            AbstractScoreMatrix::TScaling,
            TMask
    );
    ~ProfileSearching();

    void            Run();

    const char*     GetParamConfigFile() const  { return paramConfigFile; }
    const char*     GetInput() const            { return input_name; }
    const char*     GetDatabase() const         { return database_name; }
    const char*     GetOutput() const           { return output_name; }

    double          GetEvalueThreshold() const  { return max_evalue; }
    int             GetNoHitsThreshold() const  { return max_no_hits; }
    int             GetNoAlnsThreshold() const  { return max_no_alns; }
    double          GetIdentityThreshold()const { return seqidentity; }
    double          GetInformationThreshold() const { return information_content; }
    double          GetMaskscalePercents() const{ return maskscale_percnt; }
    int             GetGapOpenCost() const      { return gapopencost; }
    int             GetGapExtnCost() const      { return gapextncost; }
    bool            ToShowPars() const          { return show_pars; }
    bool            GetFixedCosts() const       { return fixed_costs; }

    bool            UngappedAlignments() const  { return configuration[ProcomGapped].IsUngapped(); }

    int             SetTarFrMix( int value ) { tfrmix_ = value; }
    bool            GetTarFrMixHDPCtx() const { return tfrmix_ == tfrmixHDPCtx; }

    int             SetScoAdjment( int value ) { scoadj_ = value; }
    bool            GetScoAdjmentHDPCtx() const { return scoadj_ == scoadjHDPCtx || scoadj_ == scoadjHDPsco; }

    void            SetHDPbase( const HDPbase* value ) { HDPbase_ = value; }
    const HDPbase*  GetHDPbase() const { return HDPbase_; }

    void            SetHDPctbase( const HDPbase* value ) { HDPctbase_ = value; }
    const HDPbase*  GetHDPctbase() const { return HDPctbase_; }

    void            SetSSEModel( int modndx ) { ssemodel_ = modndx; }
    int             GetSSEModel() const { return ssemodel_; }

    bool            GetHCSeg() const            { return hcseg; }
    size_t          GetHCWinLength() const      { return hcwinlen; }
    double          GetHCLowEntropy() const     { return hclowentropy; }
    double          GetHCHighEntropy() const    { return hchighentropy; }

    bool            GetUsingSeg() const         { return usingseg; }
    size_t          GetSegWinLength() const     { return segwinlen; }
    double          GetSegLowEntropy() const    { return seglowentropy; }
    double          GetSegHighEntropy() const   { return seghighentropy; }
    double          GetSegDistance() const      { return segdistance; }

    bool            GetUsingSeqSeg() const          { return usingseqseg; }
    size_t          GetSeqSegWinLength() const      { return seqsegwinlen; }
    double          GetSeqSegLowEntropy() const     { return seqseglowent; }
    double          GetSeqSegHighEntropy() const    { return seqseghighent; }

    bool            GetAutoGapCosts()           { return autogapcosts; }
    int             GetAutocorrWinsize()        { return autocorrwinsize; }
    void            SetAutoGapCosts( bool agc, int ac_winsize ) {
                        autogapcosts = agc;
                        autocorrwinsize = ac_winsize;
                    }


    void            SetHCParameters( size_t winlen, double lowent, double highent )  {
                        SetHCSeg( true );
                        SetHCWinLength( winlen );
                        SetHCLowEntropy( lowent );
                        SetHCHighEntropy( highent );
                    }

    void            SetSegParameters( size_t winlen, double lowent, double highent, double distance )  {
                        SetUsingSeg( true );
                        SetSegWinLength( winlen );
                        SetSegLowEntropy( lowent );
                        SetSegHighEntropy( highent );
                        SetSegDistance( distance );
                        ProfileAlignment::SetSEGdistanceThreshold( distance );
                    }
    void            SetSeqSegParameters( size_t winlen, double lowent, double highent )  {
                        SetUsingSeqSeg( true );
                        SetSeqSegWinLength( winlen );
                        SetSeqSegLowEntropy( lowent );
                        SetSeqSegHighEntropy( highent );
                    }


    bool            GetComputeDELETEstates() const          { return deletestates_; }
    void            SetComputeDELETEstates( bool value )    { deletestates_ = value; }

    double          GetPseudoCountWeight() const            { return pseudocntweight_; }
    void            SetPseudoCountWeight( double value )    { pseudocntweight_ = value; }

    double          GetExtentMinSeqPercentage() const       { return extminseqperc_; }
    void            SetExtentMinSeqPercentage( double value ) { extminseqperc_ = value; }

    size_t          GetExtentMinWindow() const              { return extminwindow_; }
    void            SetExtentMinWindow( size_t value )      { extminwindow_ = value; }


    int             GetThicknessNumber() const              { return thickness_number; }
    void            SetThicknessNumber( int value )         { thickness_number = value; }

    double          GetThicknessPercents() const            { return thickness_percnt; }
    void            SetThicknessPercents( double value )    { thickness_percnt = value; }

    bool            GetAutocorrectionPositional() const     { return  positcorrects; }
    void            SetAutocorrectionPositional( bool value){ positcorrects = value; }

    bool            GetAutoACcorrection() const             { return autoaccorrection; }
    void            SetAutoACcorrection( bool value )       { autoaccorrection = value; }

    double          GetDeletionCoefficient() const          { return deletioncoeff; }
    void            SetDeletionCoefficient( double value )  { deletioncoeff = value; }


    double          GetGapProbabFactorEvalue() const            { return gapprobfactevalue; }
    void            SetGapProbabFactorEvalue( double value )    { gapprobfactevalue = value; }

    double          GetGapProbabFactorWeight() const            { return gapprobfactweight; }
    void            SetGapProbabFactorWeight( double value )    { gapprobfactweight = value; }

    double          GetGapProbabFactorShift() const             { return gapprobfactshift; }
    void            SetGapProbabFactorShift( double value )     { gapprobfactshift = value; }


    double          GetAutocorrectionNumerator1st() const           { return acorrnumerator1st; }
    void            SetAutocorrectionNumerator1st( double value )   { acorrnumerator1st = value; }

    double          GetAutocorrectionNumerator2nd() const           { return acorrnumerator2nd; }
    void            SetAutocorrectionNumerator2nd( double value )   { acorrnumerator2nd = value; }

    double          GetAutocorrectionLogScale() const           { return acorrlogscale; }
    void            SetAutocorrectionLogScale( double value )   { acorrlogscale = value; }

    double          GetAutocorrectionDenomScale() const         { return acorrdenomscale; }
    void            SetAutocorrectionDenomScale( double value ) { acorrdenomscale = value; }


    double          GetInfoCorrectionUpperBound2nd() const          { return ainfoupper2nd; }
    void            SetInfoCorrectionUpperBound2nd( double value )  { ainfoupper2nd = value; }

    double          GetInfoCorrectionNumerator2nd() const           { return ainfonumerator2nd; }
    void            SetInfoCorrectionNumerator2nd( double value )   { ainfonumerator2nd = value; }

    double          GetInfoCorrectionScale2nd() const           { return ainfoscale2nd; }
    void            SetInfoCorrectionScale2nd( double value )   { ainfoscale2nd = value; }

    double          GetInfoCorrectionNumeratorAlt() const           { return ainfonumeratoralt; }
    void            SetInfoCorrectionNumeratorAlt( double value )   { ainfonumeratoralt = value; }

    double          GetInfoCorrectionScaleAlt() const           { return ainfoscalealt; }
    void            SetInfoCorrectionScaleAlt( double value )   { ainfoscalealt = value; }


    int             GetHSPLength() const                        { return hsplength_; }
    void            SetHSPLength( int value )                   { hsplength_ = value; }

    int             GetHSPScore() const                         { return hspscore_; }
    void            SetHSPScore( int value )                    { hspscore_ = value; }

    int             GetHSPDistance() const                      { return hspdistance_; }
    void            SetHSPDistance( int value )                 { hspdistance_ = value; }

    int             GetHSPNoHSPs() const                        { return hspnohsps_; }
    void            SetHSPNoHSPs( int value )                   { hspnohsps_ = value; }

    ProfileAlignment::TAlgorithm    GetAlnAlgorithm() const                 { return alnalgo_; }
    void            SetAlnAlgorithm( ProfileAlignment::TAlgorithm value )   { alnalgo_ = value; }


    double          GetExpectForAlnLength() const           { return evalue_alnlength; }
    void            SetExpectForAlnLength( double value )   { evalue_alnlength = value; }

    void            PrintMethodName( FILE* fp ) const;          //printing of the method name used in scoring alignments
    void            PrintParameterTable( FILE* ) const;         //printing of parameter table

protected:
    explicit ProfileSearching();
                                                                //create alternative score system given subject profile
    void                        CreateScoreSystem( const FrequencyMatrix&, const LogOddsMatrix& );
    void                        CreateScoreSystem();            //create member score system
    void                        DestroyScoreSystem();           //destroy score system
    void                        ComputeScoreSystem();           //compute score system if needed
    void                        ScaleScoreSystem();             //scale score system
    bool    ScanForHSPs( double minhspscore, int hsplen, int nohsps, int mindist, int* possbjct = NULL, int* posquery = NULL );
    bool                        PreprocessSubject(              //preprocessing of subject profile
        const FrequencyMatrix&, const LogOddsMatrix&, GapScheme&, bool firstpass );

    const AbstractScoreMatrix*  GetScoreSystem() const  { return scoreSystem; }
    AbstractScoreMatrix*        GetScoreSystem()        { return scoreSystem; }

    void                        Realloc( size_t newcap );       //memory reallocation
    void                        Push( HitInformation* hit );    //push hit into the list

    bool                        IsCompatible( GapScheme&, GapScheme& ) const;
    void                        ComputationLogicWithProfiles( FrequencyMatrix&, LogOddsMatrix&, GapScheme& );
    void                        PostComputationLogic();


    void            SetHCSeg( bool value )              { hcseg = value; }
    void            SetHCWinLength( size_t value )      { hcwinlen = value; }
    void            SetHCLowEntropy( double value )     { hclowentropy = value; }
    void            SetHCHighEntropy( double value )    { hchighentropy = value; }

    void            SetUsingSeg( bool value )           { usingseg = value; }
    void            SetSegWinLength( size_t value )     { segwinlen = value; }
    void            SetSegLowEntropy( double value )    { seglowentropy = value; }
    void            SetSegHighEntropy( double value )   { seghighentropy = value; }
    void            SetSegDistance( double value )      { segdistance = value; }

    void            SetUsingSeqSeg( bool value )            { usingseqseg = value; }
    void            SetSeqSegWinLength( size_t value )      { seqsegwinlen = value; }
    void            SetSeqSegLowEntropy( double value )     { seqseglowent = value; }
    void            SetSeqSegHighEntropy( double value )    { seqseghighent = value; }

    void            PrintSearchingHeader( FILE* );  //print general information before searching against the database
    void            PrintHits( FILE* );             //print searching results
    void            PrintSearchingFooter( FILE* );  //print detail information of the searching

    AbstractScoreMatrix::TType      GetMethod() const               { return method; }
    AbstractScoreMatrix::TBehaviour GetBehaviour() const            { return statbehaviour; }
    AbstractScoreMatrix::TScaling   GetScaling() const              { return scaling; }
    TMask                           GetMasking() const              { return masking; }

                                                    //configuration access routines
    const Configuration&        GetConfiguration( TProcomSchemes ) const;
    Configuration&              GetConfiguration( TProcomSchemes );
    Configuration* const        GetConfiguration()              { return configuration; }

    void            SetGapOpenCost( int value )         { gapopencost = value; }
    void            SetGapExtnCost( int value )         { gapextncost = value; }

    void                        ProcessQuery(); //read query information from the file
    void                        FillMatrices(   //fill profile matrices by reading and processing information from file
                                    FrequencyMatrix&, LogOddsMatrix&, GapScheme&,
                                    const char* filename );


    size_t          GetNoSequences() const      { return no_sequences; }            //number of sequences in database
    Uint64          GetDbSize() const           { return db_size; }                 //size of database

    void                        SetNoSequences( size_t value )  { no_sequences = value; }
    void                        SetDbSize( Uint64 value )       { db_size = value; }


    HitInformation*&            GetHitAt( size_t );                                 //hit at the specified location
    const HitInformation*       GetHitAt( size_t ) const;                           //hit at the specified location
    size_t                      GetHitlistCapacity() const  { return capacity; }    //capacity of hitListing
    size_t                      GetHitlistSize() const      { return size; }        //size of hitListing

    void                        QSortHits();                                        //quick sort of hits

private:
    const char*             paramConfigFile;//full pathname to parameter configuration file
    const char*             input_name;     //input profile's name
    const char*             database_name;  //profile database name
    const char*             output_name;    //output file name, null if standard output
    AbstractScoreMatrix*    scoreSystem;    //score system used to align profiles
    Configuration           configuration[NoSchemes];   //parameter configuration

    double                  max_evalue;     //e-value threshold used for outputing of the alignments
    int                     max_no_hits;    //maximum number of hits to show in the result list
    int                     max_no_alns;    //maximum number of alignments to show in the output
    double                  seqidentity;    //sequence identity used to make profiles
    const double            information_content;//information content threshold
    bool                    deletestates_;  //flag of computing delete states

    double                  pseudocntweight_;   //weight for pseudo count frequencies
    double                  extminseqperc_;     //minimum required sequence percentage an extent must cover
    size_t                  extminwindow_;      //minimum required extent window length in positions

    int                     thickness_number;   //thickness in number of sequences a position must have (rarely used)
    double                  thickness_percnt;   //thickness in percentage (rarely used)

    const double            maskscale_percnt;   //scaling of masked positions in percentage
    int                     gapopencost;    //gap open cost
    int                     gapextncost;    //gap extend cost
    bool                    show_pars;      //show statistical significance parameters below the alignments
    bool                    fixed_costs;    //use fixed (not position-specific) gap cost scheme

    HitInformation**        hitListing;     //the resulting hit listing
    size_t                  capacity;       //capacity of hitListing
    size_t                  size;           //current size of hitListing

    FrequencyMatrix         query_freq;     //frequency matrix for the query
    LogOddsMatrix           query_pssm;     //PSSM matrix for the query
    GapScheme               query_gaps;     //position-specific gap costs for the query

    Database                profile_db;     //profile database
    size_t                  no_sequences;   //number of sequences within the database
    Uint64                  db_size;        //size of the database

    AbstractScoreMatrix::TType      method;         //scoring method
    AbstractScoreMatrix::TBehaviour statbehaviour;  //statistical behaviour
    AbstractScoreMatrix::TScaling   scaling;        //strategy of precision
    TMask                           masking;        //masking approach
    ProfileAlignment::TAlgorithm    alnalgo_;       //alignment algorithm

    int                     tfrmix_;        //mixing of target frequencies
    int                     scoadj_;        ////score adjustment
    const HDPbase*          HDPbase_;       //HDP base structure with data
    const HDPbase*          HDPctbase_;     //HDP ctx base structure

    int                     ssemodel_;      //index of a model for the estimation of stat. sign.

    bool                    autogapcosts;   //automatically computed gap opening costs
    int                     autocorrwinsize;//autocorrelation window size

    bool                    hcseg;          //HC seg in use
    size_t                  hcwinlen;       //HC seg window length
    double                  hclowentropy;   //HC seg low entropy threshold
    double                  hchighentropy;  //HC seg  high entropy threshold

    bool                    usingseg;       //seg in use
    size_t                  segwinlen;      //seg window length
    double                  seglowentropy;  //seg low entropy threshold
    double                  seghighentropy; //seg  high entropy threshold
    double                  segdistance;    //seg  distance of equivalence

    bool                    usingseqseg;    //whether to use seg for sequences in multiple alignment
    size_t                  seqsegwinlen;   //window length
    double                  seqseglowent;   //low entropy threshold
    double                  seqseghighent;  //high entropy threshold

    double                  deletioncoeff;          //deletion probability weight

    double                  gapprobfactevalue;      //evalue threshold for computing of gap probability factor
    double                  gapprobfactweight;      //argument weight in expression of gap probability factor
    double                  gapprobfactshift;       //argument shift in expression of gap probability factor

    double                  acorrnumerator1st;      //numerator of expression to compute 1st-pass autocorrection
    double                  acorrnumerator2nd;      //numerator of expression to compute 2nd-pass upper bound for autocorrection
    double                  acorrlogscale;          //logarithmic scale to compute 2nd-pass autocorrection
    double                  acorrdenomscale;        //denominator scale to compute 2nd-pass autocorrection

    double                  ainfoupper2nd;          //upper bound of information content threshold used in 2nd-pass computations
    double                  ainfonumerator2nd;      //numerator of expression to compute 2nd-pass information content threshold
    double                  ainfoscale2nd;          //logarithmic scale to compute 2nd-pass inf. content threshold
    double                  ainfonumeratoralt;      //numerator of alternative expression to compute inf. content threshold
    double                  ainfoscalealt;          //logarithmic scale to alternatively compute inf. content threshold

    int                     hsplength_;             //lengths of HSPs
    int                     hspscore_;              //minimum HSP score
    int                     hspdistance_;           //distance between the HSPs
    int                     hspnohsps_;             //number of HSPs in the diagonal

    double                  evalue_alnlength;       //evalue to compute expected mean alignment length
    bool                    positcorrects;          //whether corrections computed positionally by entropies are to be used
    bool                    autoaccorrection;       //if auto correction for autocorrelation gap cost function is in effect
};


////////////////////////////////////////////////////////////////////////////
// INLINES
//
// -------------------------------------------------------------------------
// comparison operator for the HitInformation ojects
//
inline
bool HitInformation::operator<( const HitInformation& right ) const
{
    return GetRefEvalue() < right.GetRefEvalue();
}

// -------------------------------------------------------------------------
// CLASS ProfileSearching
//
// IsCompatible: whether the profiles are mutually compatible
// -------------------------------------------------------------------------

inline
bool ProfileSearching::IsCompatible( GapScheme& first, GapScheme& secnd ) const
{
    return first.IsCompatible( secnd );
}

// -------------------------------------------------------------------------
// GetConfiguration: returns reference to parameter configuration object
// -------------------------------------------------------------------------

inline
const Configuration& ProfileSearching::GetConfiguration( TProcomSchemes ps ) const
{
    if( NoSchemes <= ps )
        throw myruntime_error( mystring( "ProfileSearching: Wrong argument while accessing configuration." ));

    return configuration[ps];
}

// GetConfiguration: overloaded
//

inline
Configuration& ProfileSearching::GetConfiguration( TProcomSchemes ps )
{
    if( NoSchemes <= ps )
        throw myruntime_error( mystring( "ProfileSearching: Wrong argument while accessing configuration." ));

    return configuration[ps];
}

// -------------------------------------------------------------------------
// obtain hit at the specified location

inline
const HitInformation* ProfileSearching::GetHitAt( size_t n ) const
{
#ifdef __DEBUG__
    if( n < size )
#endif
        return hitListing[n];

    throw myruntime_error( mystring( "ProfileSearching: Memory access error." ));
}

// -------------------------------------------------------------------------
// returns address of hit at the specified location to be modified

inline
HitInformation*& ProfileSearching::GetHitAt( size_t n )
{
#ifdef __DEBUG__
    if( n < size )
#endif
        return hitListing[n];

    throw myruntime_error( mystring( "ProfileSearching: Memory access error." ));
}

// -------------------------------------------------------------------------
// CreateScoreSystem: creates score system used to align profiles
//
inline
void ProfileSearching::CreateScoreSystem()
{
    if( GetMethod() != AbstractScoreMatrix::Universal )
        return;

    scoreSystem = new UniversalScoreMatrix(
                query_freq,
                query_pssm,
                profile_db.GetStore(),
                GetInformationThreshold(),
                GetThicknessNumber(),
                GetThicknessPercents(),
                GetConfiguration(),
                GetBehaviour(), //ComputeStatistics or StatisticsGiven
                GetScaling(),   //AutoScalling or NoScaling
                GetMasking(),
                Database::GetDistributionType()
    );

    if( scoreSystem == NULL )
        throw myruntime_error( mystring( "ProfileSearching: Not enough memory." ));

    scoreSystem->SetDeletionCoefficient( GetDeletionCoefficient());
    scoreSystem->SetInfoCorrectionUpperBound2nd( GetInfoCorrectionUpperBound2nd());
    scoreSystem->SetInfoCorrectionNumerator2nd( GetInfoCorrectionNumerator2nd());
    scoreSystem->SetInfoCorrectionScale2nd( GetInfoCorrectionScale2nd());
    scoreSystem->SetInfoCorrectionNumeratorAlt( GetInfoCorrectionNumeratorAlt());
    scoreSystem->SetInfoCorrectionScaleAlt( GetInfoCorrectionScaleAlt());
    scoreSystem->SetAutocorrectionPositional( GetAutocorrectionPositional());

    ScaleScoreSystem();
}

// CreateScoreSystem: creates score system of alternative type
//
inline
void ProfileSearching::CreateScoreSystem(
        const FrequencyMatrix& sbjctfreq,
        const LogOddsMatrix& sbjctpssm )
{
    switch( GetMethod()) {
        case AbstractScoreMatrix::ProfileSpecific:
                scoreSystem = new ScoringMatrix(
                    query_freq,
                    query_pssm,
                    sbjctfreq,
                    sbjctpssm,
                    GetInformationThreshold(),
                    GetThicknessNumber(),
                    GetThicknessPercents(),
                    GetMaskscalePercents(),
                    GetConfiguration(),
                    GetBehaviour(), //ComputeStatistics or StatisticsGiven
                    GetScaling(),   //FPScaling, AutoScalling, or NoScaling
                    GetMasking()
                );
        break;
        case AbstractScoreMatrix::ProfSpecLSO:
                scoreSystem = new LSOSMatrix(
                    query_freq,
                    query_pssm,
                    sbjctfreq,
                    sbjctpssm,
                    GetScoAdjmentHDPCtx()? GetHDPbase(): NULL,
                    GetScoAdjmentHDPCtx()? GetHDPctbase(): NULL,
                    GetInformationThreshold(),
                    GetThicknessNumber(),
                    GetThicknessPercents(),
                    GetMaskscalePercents(),
                    GetConfiguration(),
                    GetBehaviour(), //ComputeStatistics or StatisticsGiven
                    GetScaling(),   //FPScaling, AutoScalling, or NoScaling
                    GetMasking()
                );
        break;
        case AbstractScoreMatrix::AdjustedProfileSpecific:
                scoreSystem = new AdjustedScoreMatrix(
                    query_freq,
                    query_pssm,
                    sbjctfreq,
                    sbjctpssm,
                    profile_db.GetStore(),
                    GetInformationThreshold(),
                    GetThicknessNumber(),
                    GetThicknessPercents(),
                    GetMaskscalePercents(),
                    GetConfiguration(),
                    GetBehaviour(), //ComputeStatistics or StatisticsGiven
                    GetScaling(),   //FPScaling, AutoScalling, or NoScaling
                    GetMasking()
                );
        break;
        case AbstractScoreMatrix::HDPProfileSpecific:
                scoreSystem = new HDP0ScoreMatrix(
                    GetHDPbase(),
                    query_freq,
                    query_pssm,
                    sbjctfreq,
                    sbjctpssm,
                    GetInformationThreshold(),
                    GetThicknessNumber(),
                    GetThicknessPercents(),
                    GetMaskscalePercents(),
                    GetConfiguration(),
                    GetBehaviour(), //ComputeStatistics or StatisticsGiven
                    GetScaling(),   //FPScaling, AutoScalling, or NoScaling
                    GetMasking()
                );
        break;
        case AbstractScoreMatrix::HDP0CtxProfileSpecific:
                scoreSystem = new HDP0CtxScoreMatrix(
                    GetHDPbase(),
                    query_freq,
                    query_pssm,
                    sbjctfreq,
                    sbjctpssm,
                    GetInformationThreshold(),
                    GetThicknessNumber(),
                    GetThicknessPercents(),
                    GetMaskscalePercents(),
                    GetConfiguration(),
                    GetBehaviour(), //ComputeStatistics or StatisticsGiven
                    GetScaling(),   //FPScaling, AutoScalling, or NoScaling
                    GetMasking()
                );
        break;
        default:    return;
    };

    if( scoreSystem == NULL )
        throw myruntime_error( mystring( "ProfileSearching: Not enough memory." ));

    scoreSystem->SetDeletionCoefficient( GetDeletionCoefficient());
    scoreSystem->SetInfoCorrectionNumerator2nd( GetInfoCorrectionNumerator2nd());
    scoreSystem->SetInfoCorrectionUpperBound2nd( GetInfoCorrectionUpperBound2nd());
    scoreSystem->SetInfoCorrectionScale2nd( GetInfoCorrectionScale2nd());
    scoreSystem->SetInfoCorrectionNumeratorAlt( GetInfoCorrectionNumeratorAlt());
    scoreSystem->SetInfoCorrectionScaleAlt( GetInfoCorrectionScaleAlt());
    scoreSystem->SetAutocorrectionPositional( GetAutocorrectionPositional());
    scoreSystem->SetName( sbjctpssm.GetName());
}

// -------------------------------------------------------------------------
// DestroyScoreSystem: destroys score system
//
inline
void ProfileSearching::DestroyScoreSystem()
{
    if( scoreSystem )
        delete scoreSystem;
}

// -------------------------------------------------------------------------
// computeProfileScoringMatrix: calls appropriate method of ScoresMatrix object
// -------------------------------------------------------------------------

inline
void ProfileSearching::ComputeScoreSystem()
{
#ifdef __DEBUG__
    if( !scoreSystem )
        throw myruntime_error( mystring( "ProfileSearching: Unable to compute score matrix." ));
#endif
    scoreSystem->ComputeProfileScoringMatrix();
}

// -------------------------------------------------------------------------
// ScanForHSPs: run algorithm of hsps
//
inline
bool ProfileSearching::ScanForHSPs(
    double minhspscore, int hsplen, int nohsps, int mindist,
    int* possbjct, int* posquery )
{
#ifdef __DEBUG__
    if( !scoreSystem )
        throw myruntime_error( mystring( "ProfileSearching: Unable to scan score matrix for HSPs." ));
#endif
    if( minhspscore <= 0.0 )
        return true;

    switch( scoreSystem->GetType()) {
        case AbstractScoreMatrix::ProfileSpecific:
        case AbstractScoreMatrix::ProfSpecLSO:
        case AbstractScoreMatrix::AdjustedProfileSpecific:
        case AbstractScoreMatrix::HDPProfileSpecific:
        case AbstractScoreMatrix::HDP0CtxProfileSpecific:
            return scoreSystem->ScanForHSPs( minhspscore, hsplen, nohsps, mindist, possbjct, posquery );
//             break;

        case AbstractScoreMatrix::Universal:
            break;

        default:
            break;
    };
    return true;
}

// -------------------------------------------------------------------------
// ScaleScoreSystem: scale score system
//
inline
void ProfileSearching::ScaleScoreSystem()
{
#ifdef __DEBUG__
    if( !scoreSystem )
        throw myruntime_error( mystring( "ProfileSearching: Unable to scale score matrix." ));
#endif
    switch( scoreSystem->GetType()) {
        case AbstractScoreMatrix::ProfileSpecific:
            scoreSystem->ScaleScoringMatrix();
            break;

        case AbstractScoreMatrix::ProfSpecLSO:
            scoreSystem->ScaleScoringMatrix();
            break;

        case AbstractScoreMatrix::AdjustedProfileSpecific:
            scoreSystem->ScaleScoringMatrix();
            break;

        case AbstractScoreMatrix::HDPProfileSpecific:
        case AbstractScoreMatrix::HDP0CtxProfileSpecific:
            scoreSystem->ScaleScoringMatrix();
            break;

        case AbstractScoreMatrix::Universal:
//             scoreSystem->ComputeStatisticalParameters();
            scoreSystem->ScaleScoringMatrix();
            break;

        default:
            break;
    };
}

// -------------------------------------------------------------------------
// PreprocessSubject: prepares scoring system for alignment with subject
//     profile
//
inline
bool ProfileSearching::PreprocessSubject(
    const FrequencyMatrix&  freq,
    const LogOddsMatrix&    pssm,
    GapScheme&              gaps,
    bool                    firstpass )
{
    if( firstpass ) {
        if( scoreSystem == NULL ||
            scoreSystem->GetType() == AbstractScoreMatrix::ProfileSpecific ||
            scoreSystem->GetType() == AbstractScoreMatrix::ProfSpecLSO ||
            scoreSystem->GetType() == AbstractScoreMatrix::AdjustedProfileSpecific ||
            scoreSystem->GetType() == AbstractScoreMatrix::HDPProfileSpecific ||
            scoreSystem->GetType() == AbstractScoreMatrix::HDP0CtxProfileSpecific )
        {
            DestroyScoreSystem();
            CreateScoreSystem( freq, pssm );
            ComputeScoreSystem();
            if( !ScanForHSPs( GetHSPScore(), GetHSPLength(), GetHSPNoHSPs(), GetHSPDistance()))
                return false;
            if( scoreSystem->GetSupportOptimFreq())///originally commented
                scoreSystem->OptimizeTargetFrequencies();
            ScaleScoreSystem();
        }
#ifdef __DEBUG__
        if( !scoreSystem )
            throw myruntime_error( mystring( "ProfileSearching: Unable to preprocess subject profile." ));
#endif

        if( scoreSystem->GetType() == AbstractScoreMatrix::Universal )
            dynamic_cast<UniversalScoreMatrix*>( scoreSystem )->PreserveSubject( freq, pssm );
    }

    gaps.SetUsePosACcorrections( GetAutocorrectionPositional());
    query_gaps.SetUsePosACcorrections( GetAutocorrectionPositional());

    scoreSystem->PostScalingProc(
        query_pssm, pssm,
        query_gaps, gaps,
        GetAutoGapCosts(),
        GetAutocorrWinsize(),
        firstpass
    );
    return true;
}

// -------------------------------------------------------------------------
// PrintMethodName: prints name of method used to score profile positions
//
inline
void ProfileSearching::PrintMethodName( FILE* fp ) const
{
    const size_t    padding = OUTPUTINDENT;
    char            blanks[padding+1];
    size_t          n = 0;

    for( n = 0; n < padding; blanks[n++] = 32 );
    blanks[n] = 0;

    if( scoreSystem ) {
        fprintf( fp, "Scoring method: %s\n", scoreSystem->GetMethodName());
        if( scoreSystem->GetAutoScaling())
            fprintf( fp, "%s(Increased precision)\n", blanks );
        else if( scoreSystem->GetFPScaling())
                fprintf( fp, "%s(Floating-point precision)\n", blanks );
        fprintf( fp, "\n");
        return;
    }
}

// -------------------------------------------------------------------------
// PrintParameterTable: prints parameter table
//
inline
void ProfileSearching::PrintParameterTable( FILE* fp ) const
{
    if( !scoreSystem )
        return;

    scoreSystem->PrintFinal( fp );
}


#endif//__ProfileSearching__
