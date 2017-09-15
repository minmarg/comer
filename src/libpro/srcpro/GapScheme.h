/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __GapScheme__
#define __GapScheme__

#include <math.h>

#include "debug.h"
#include "rc.h"
#include "mystring.h"
#include "myexcept.h"

#include "datapro.h"
#include "Configuration.h"
#include "DistributionMatrix.h"

#define DEFAULT_AC_CORRECTION   ( 1.2 )

class Serializer;

// enumeration for delete state boundaries
enum TDeleteBounds {
    DBeg,   //beginning
    DEnd,   //end
    DSlope, //slope of probability
    DCnt    //number of terms
};

////////////////////////////////////////////////////////////////////////////
// CLASS GapScheme
// Information of and manipulation with positional gaps
//
class GapScheme
{
public:
    GapScheme( double gap_open = DEFAULTGAPOPENCOST, double gap_extend = DEFAULTEXTENDCOST );
    ~GapScheme();

    const Configuration*    GetConfiguration() const { return configuration; }
    void        SetConfiguration( const Configuration* config );

    int         GetColumns() const          { return length; }

    bool        GetFixed() const            { return fixed; }
    void        SetFixed( bool value )      { fixed = value; }

    double      GetOpenCost() const         { return openCost; }
    double      GetExtendCost() const       { return extendCost; }

    bool        GetUsePosACcorrections() const                  { return useposaccorrect; }
    void        SetUsePosACcorrections( bool value )            { useposaccorrect = value; }


    double          GetGapProbabFactorEvalue() const            { return gapprobfactevalue; }
    void            SetGapProbabFactorEvalue( double value )    { gapprobfactevalue = value; }

    double          GetGapProbabFactorWeight() const            { return gapprobfactweight; }
    void            SetGapProbabFactorWeight( double value )    { gapprobfactweight = value; }

    double          GetGapProbabFactorShift() const             { return gapprobfactshift; }
    void            SetGapProbabFactorShift( double value )     { gapprobfactshift = value; }


    double      GetAutocorrectionNumerator1st() const           { return acorrnumerator1st; }
    void        SetAutocorrectionNumerator1st( double value )   { acorrnumerator1st = value; }

    double      GetAutocorrectionNumerator2nd() const           { return acorrnumerator2nd; }
    void        SetAutocorrectionNumerator2nd( double value )   { acorrnumerator2nd = value; }

    double      GetAutocorrectionLogScale() const               { return acorrlogscale; }
    void        SetAutocorrectionLogScale( double value )       { acorrlogscale = value; }

    double      GetAutocorrectionDenomScale() const             { return acorrdenomscale; }
    void        SetAutocorrectionDenomScale( double value )     { acorrdenomscale = value; }

    void        AdjustContextByEval( double evalue, double relent, const double* = NULL );
    void        ResetContext();

    static double   GetTransPowIndex( int tr );
    static void     SetTransPowIndex( int tr, double value );

    double          GetOrgTrProbsAt( int trans, int pos ) const;
    const double ( *GetOrgTrProbsAt( int pos ) const )[P_NSTATES];

    double          GetTransProbsAt( int trans, int pos ) const;
    const double ( *GetTransProbsAt( int pos ) const )[P_NSTATES];

    double          GetLogTransProbsAt( int trans, int pos ) const;
    const double ( *GetLogTransProbsAt( int pos ) const )[P_NSTATES];

    double          GetNSLogTransProbsAt( int trans, int pos ) const;
    const double ( *GetNSLogTransProbsAt( int pos ) const )[P_NSTATES];

    double      GetDeleteOpenWeightAt( int pos ) const  { return GetDeletesBegAt( pos ); }
    double      GetDeleteOpenProbAt( int pos ) const;
    double      GetDeleteExtendWeightAt( int pos, int time ) const;
    double      GetDeleteExtensionProbAt( int pos, int time ) const;

    void        SetOpenCost( double value )                     { openCost = value; }
    void        SetExtendCost( double value )                   { extendCost = value; }

    double      GetPosOpenAt( int pos ) const;
    double      GetPosExtendAt( int pos ) const;

    void        SetPosOpenAt( double value, int pos );
    void        SetPosExtendAt( double value, int pos );
    void        Initialize();                       //initialization of data given partial information

    double      GetOpenAt( int ) const;             //gap opening cost at the position
    double      GetExtendAt( int ) const;           //gap extend cost at the position
    double      GetWeightsAt( int ) const;          //gap weight at the position
    double      GetProbabilityAt( int ) const;      //gap probability at the position

    double      GetDeletesBegAt( int ) const;       //starting deletion weight at the position
    double      GetDeletesEndAt( int ) const;       //final deletion weight at the position
    double      GetDeletesSlopeAt( int ) const;     //slope of deletion variability at the position

    int         GetDeletesIntervalAt( int ) const;  //interval of deletion variability at the position

    double      operator[]( int ind ) const;

    void        Clear();
    void        Push( const double (*)[P_NSTATES], double w, double delbeg, double delend, int delta, char );
    void        PushAt( const double (*)[P_NSTATES], double w, double delbeg, double delend, int delta, char, int pos );
    void        SetOrgTrProbsBeg( const double (*)[P_NSTATES]);

    double      GetScoresMultiplier()               { return scoremult_; }
    void        SetScoresMultiplier( double value ) { scoremult_ = value; }

    void        Prepare( int );                                     //prepare gap costs for using for alignments
    void        Prepare( int, const double* posinfcontents );       //prepare augmented with init of positional values of correct.
    void        Prepare( double opencost, double extncost, bool );  //prepare: overloaded
                                                                //adjust gap costs with autocorrelation function
    void        AdjustCosts( double*, size_t, const LogOddsMatrix&, int window, bool* masks, size_t, double );

    void        Serialize( Serializer& ) const;
    void        Deserialize( Serializer& );

    bool        IsCompatible( const GapScheme& one ) const;         //is it compatible with other object?
    bool        IsComposIdentical( const GapScheme& one ) const;    //whether it is compositionally identical

    void        OutputGapScheme( const char* = NULL );

    char        AAcid( int ind ) const;

    void        Reserve( int amount ) { reallocate( amount ); }

protected:
    void        destroy();                          //memory deallocation and reset of values
    void        reallocate( int howmuch );          //memory allocation
    void        SetColumns( int col )               { length = col; }

    void        ReevaluateTransProbs();
    void        InitializePosGapCosts();            //initialize positional gap costs (upper bound)
    void        InitializePosGapOpenCosts();
    void        InitializePosGapExtendCosts();
                                                    //initialize positional values of correction
    void        InitializePosACcorrection( const double* posinfcontents, double = -1.0 );

    void        CheckIntegrity() const;             //check integrity of the structure

    void        ResetExtendCost();
    void        SetExtendCostByEval( double evalue, double relent );

    double      GetACcorrection() const             { return autcorrection; }
    void        ResetACcorrection();                //reset alternative correction term for autocorrelation function
    void        ResetPosACcorrection();             //reset positional values of correction for autocorrelation function
    void        SetACcorrection( double value );
    void        SetACcorrectionByEval( double evalue, double relent, const double* = NULL );

    double      GetScaledACcorrection() const           { return ac_correction; }
    void        SetScaledACcorrection( double value )   { ac_correction = value; }

    int         GetScaleFactor() const              { return scalefactor; }
    void        SetScaleFactor( int value )         { scalefactor = value; }

    double      GetScaledOpenCost() const           { return scaledopenCost; }
    void        SetScaledOpenCost( double value )   { scaledopenCost = value; }

    double      GetScaledExtendCost() const         { return scaledextendCost; }
    void        SetScaledExtendCost( double value ) { scaledextendCost = value; }

    void        SetOrgTrProbsAt( double value, int trans, int pos );
    void        SetOrgTrProbsAt( const double (*)[P_NSTATES], int pos );

    void        SetTransProbsAt( double value, int trans, int pos );
    void        SetTransProbsAt( const double (*)[P_NSTATES], int pos );

    void        SetLogTransProbsAt( double value, int trans, int pos );
    void        SetLogTransProbsAt( const double (*)[P_NSTATES], int pos );

    void        SetNSLogTransProbsAt( double value, int trans, int pos );
    void        SetNSLogTransProbsAt( const double (*)[P_NSTATES], int pos );

    void        SetOpenAt( double value, int pos );
    void        SetExtendAt( double value, int pos );
    void        SetWeightsAt( double value, int pos );

    const double*   GetPosACcorrection() const      { return positionalaccorr; }
    double          GetPosACcorrectionAt( int pos ) const;          //get correction for autocorrelation function at position pos
    void            SetPosACcorrectionAt( double value, int pos );  //set correction for autocorrelation function at position pos

    void        ComputeDeleteSlopes();                      //compute delete-state slopes for all positions

    void        SetDeletesBegAt( double value, int );       //set starting deletion weight at the position
    void        SetDeletesEndAt( double value, int );       //set final deletion weight at the position
    void        ComputeDeletesSlopeAt( int );               //compute slope of deletion variability at the position

    void        SetDeletesIntervalAt( int value, int );     //set interval of deletion variability at the position
                                                            //set delete-state information at the position
    void        SetDeleteStateAt( double beg, double end, int delta, int );

                                                    //autocorrelation function
    double      Autocorrelation( double*, int window, const double* info, const double* posco );
    void        ComputeCosts();                     //compute position-specific gap costs

    size_t      GetThickness() const                { return thickness; }
    void        SetThickness( size_t value )        { thickness = value; }
    void        ResetThickness()                    { thickness = 0; }

    bool        GetContextOn() const                { return contextadjusted; }
    void        SetContextOn()                      { contextadjusted = true; }
    void        SetContextOff()                     { contextadjusted = false; }

    double      GetContextEvalue() const            { return contextevalue; }
    void        SetContextEvalue( double evalue )   { contextevalue = evalue; }
    void        ResetContextEvalue()                { contextevalue = -1.0; }

    double      GetProbabilityFactor() const;
    static double   SigmoidResponse( double );

private:
    bool    fixed;                  //whether gap costs are fixed at the positions
    double* posvec;                 //vector of positional gap opening costs (upper bound)
    double* posext;                 //vector of positional gap extension costs (upper bound)

    static double s_trnpowind_[P_NSTATES];//transition power indices
    double (*orgtrobs_)[P_NSTATES]; //estimated (target) transition probabilities
    double (*trnprobs_)[P_NSTATES]; //recalculated transition probabilities
    double (*logtrnps_)[P_NSTATES]; //logs of recalculated transition probabilities
    double (*nsltrnps_)[P_NSTATES]; //not scaled logs of recalculated transition probabilities
    double* vector;                 //vector of gap costs for each position of sequence
    double* extvec;                 //vector of gap extend costs for each position of sequence
    double* weights;                //vector of position-specific gap weights
    double* deletes[DCnt];          //position-specific deletion weights
    int*    deleteInt;              //intervals of deletions
    int     length;                 //length of gapCostVector
    double  openCost;               //default gap opening cost
    double  extendCost;             //default cost for extending a gap
    double  scaledopenCost;         //scaled gap opening cost
    double  scaledextendCost;       //scaled cost for extending a gap
    char*   aacids;                 //sequence of amino acids
    int     allocated;              //how many positions allocated
    int     scalefactor;            //factor to scale gap-cost-related values
    double  scoremult_;             //multiplier of scores (from scaling of score system)
    const Configuration*
            configuration;          //system configuration

    double  autcorrection;          //logical value of correction used in autocorrelation function
    double  ac_correction;          //scaled correction used in autocorrelation function
    double* positionalaccorr;       //positional values of correction used in autocorrelation function

    double  gapprobfactevalue;      //evalue threshold for computing of gap probability factor
    double  gapprobfactweight;      //argument weight in expression of gap probability factor
    double  gapprobfactshift;       //argument shift in expression of gap probability factor

    double  acorrnumerator1st;      //numerator of expression to compute 1st-pass autocorrection
    double  acorrnumerator2nd;      //numerator of expression to compute 2nd-pass upper bound for autocorrection
    double  acorrlogscale;          //logarithmic scale to compute 2nd-pass autocorrection
    double  acorrdenomscale;        //denominator scale to compute 2nd-pass autocorrection

    size_t  thickness;              //effective thickness in sequences used to construct corresponding profile
    double  contextevalue;          //context evalue; -1 in case of reset context
    bool    contextadjusted;        //whether a context is adjusted
    bool    useposaccorrect;        //whether to use positional values of correction to autocorrelation function
};

// =========================================================================
// INLINES...
//
// SetConfiguration: sets the configuration and resets parameters
//     depending on it
//
inline
void GapScheme::SetConfiguration( const Configuration* config )
{
    configuration = config;
    ResetContext();
}

// -------------------------------------------------------------------------
// SigmoidResponse: gives sigmoid response corresponding to weight which is
//     supposed to be in [0,1]. Range of the returned value is (-.5,.5)
// -------------------------------------------------------------------------

inline
double GapScheme::SigmoidResponse( double weight )
{
    static double   wshift = 0.5;
    static double   rshift = 0.5;
    double          expont = 10.0;

    return 1.0 / ( 1.0 + exp( -expont * ( weight - wshift ))) - rshift;
}

// -------------------------------------------------------------------------
// GetProbabilityFactor: returns factor of probability that depends on
//     effective thickness of corresponding profile; return value is
//     always in [0,1]
//
inline
double GapScheme::GetProbabilityFactor() const
{
// return 1.0;
    double  evalthresh = GetGapProbabFactorEvalue();
    double  weightfact = GetGapProbabFactorWeight();
    double  argshift = GetGapProbabFactorShift();

    double  effthick = GetThickness();
    double  evalue = GetContextEvalue();
#ifdef __DEBUG__
    if( effthick < 0 )
        throw myruntime_error( mystring( "GapScheme: GetProbabilityFactor: Negative effective thickness." ));
#endif
    if( 0.0 <= evalue && evalue < evalthresh )
        return 1.0;

    double  factor01 = 1.0 / ( 1.0 + exp( -weightfact * effthick + argshift ));
    return  factor01;
}

// -------------------------------------------------------------------------
// operator[]: used to access a gap opening cost at the position
// -------------------------------------------------------------------------

inline
double GapScheme::operator[]( int n ) const
{
    return GetOpenAt( n );
}

// -------------------------------------------------------------------------
// AAcid: used to access amino acid at the position
// -------------------------------------------------------------------------

inline
char GapScheme::AAcid( int n ) const
{
#ifdef __DEBUG__
    if( length <= n || n < 0 || aacids == NULL )
        throw myruntime_error(
            mystring( "GapScheme: Memory access error." ));
#endif
    return aacids[n];
}

// -------------------------------------------------------------------------
// GetTransFactor: get transition probability power index
//
inline
double GapScheme::GetTransPowIndex( int tr )
{
#ifdef __DEBUG__
    if( tr < 0 || P_NSTATES <= tr )
        throw myruntime_error("GapScheme: GetTransPowIndex: Memory access error.");
#endif
    return s_trnpowind_[tr];
}

inline
void GapScheme::SetTransPowIndex( int tr, double value )
{
#ifdef __DEBUG__
    if( tr < 0 || P_NSTATES <= tr )
        throw myruntime_error("GapScheme: SetTransPowIndex: Memory access error.");
#endif
    s_trnpowind_[tr] = value;
}

// -------------------------------------------------------------------------
// GetOrgTrProbsAt: return specified target transition probability at the 
//     position
//
inline
double GapScheme::GetOrgTrProbsAt( int trans, int pos ) const
{
    if( orgtrobs_ == NULL ||
        length <= pos || pos < -1 || P_NSTATES <= trans || trans < 0 )
        throw myruntime_error( "GapScheme: Memory access error." );
    return orgtrobs_[pos+1][trans];
}

// -------------------------------------------------------------------------
// GetOrgTrProbsAt: return target transition probabilities at the position
//
inline
const double ( *GapScheme::GetOrgTrProbsAt( int pos ) const )[P_NSTATES]
{
    if( orgtrobs_ == NULL || length <= pos || pos < -1 )
        throw myruntime_error( "GapScheme: Memory access error." );
    return orgtrobs_ + ( pos+1 );
}


// -------------------------------------------------------------------------
// GetTransProbsAt: return specified recalc. target transition 
//     probability at the position
//
inline
double GapScheme::GetTransProbsAt( int trans, int pos ) const
{
    if( trnprobs_ == NULL ||
        length <= pos || pos < -1 || P_NSTATES <= trans || trans < 0 )
        throw myruntime_error( "GapScheme: Memory access error." );
    return trnprobs_[pos+1][trans];
}

// -------------------------------------------------------------------------
// GetTransProbsAt: return recalc. target transition probabilities at the 
//     position
//
inline
const double ( *GapScheme::GetTransProbsAt( int pos ) const )[P_NSTATES]
{
    if( trnprobs_ == NULL || length <= pos || pos < -1 )
        throw myruntime_error( "GapScheme: Memory access error." );
    return trnprobs_ + ( pos+1 );
}


// -------------------------------------------------------------------------
// GetLogTransProbsAt: return log of specified recalc. target transition 
//     probability at the position
//
inline
double GapScheme::GetLogTransProbsAt( int trans, int pos ) const
{
    if( logtrnps_ == NULL ||
        length <= pos || pos < -1 || P_NSTATES <= trans || trans < 0 )
        throw myruntime_error( "GapScheme: Memory access error." );
    return logtrnps_[pos+1][trans];
}

// -------------------------------------------------------------------------
// GetLogTransProbsAt: return log of recalc. target transition 
//     probabilities at the position
//
inline
const double ( *GapScheme::GetLogTransProbsAt( int pos ) const )[P_NSTATES]
{
    if( logtrnps_ == NULL || length <= pos || pos < -1 )
        throw myruntime_error( "GapScheme: Memory access error." );
    return logtrnps_ + ( pos+1 );
}


// -------------------------------------------------------------------------
// GetNSLogTransProbsAt: return non-scaled log of specified recalc. target 
//     transition probability at the position
//
inline
double GapScheme::GetNSLogTransProbsAt( int trans, int pos ) const
{
    if( nsltrnps_ == NULL ||
        length <= pos || pos < -1 || P_NSTATES <= trans || trans < 0 )
        throw myruntime_error( "GapScheme: Memory access error." );
    return nsltrnps_[pos+1][trans];
}

// -------------------------------------------------------------------------
// GetNSLogTransProbsAt: return non-scaled log of recalc. target transition 
//     probabilities at the position
//
inline
const double ( *GapScheme::GetNSLogTransProbsAt( int pos ) const )[P_NSTATES]
{
    if( nsltrnps_ == NULL || length <= pos || pos < -1 )
        throw myruntime_error( "GapScheme: Memory access error." );
    return nsltrnps_ + ( pos+1 );
}


// -------------------------------------------------------------------------
// GetPosOpenAt: returns positional gap opening cost at the position
//     (upper bound)
// -------------------------------------------------------------------------

inline
double GapScheme::GetPosOpenAt( int n ) const
{
#ifdef __DEBUG__
    if( length <= n || n < 0 || posvec == NULL )
        throw myruntime_error(
            mystring( "GapScheme: Memory access error." ));
#endif
    return posvec[n];
}

// -------------------------------------------------------------------------
// GetPosExtendAt: returns positional gap extension cost at the position
//     (upper bound)
// -------------------------------------------------------------------------

inline
double GapScheme::GetPosExtendAt( int n ) const
{
#ifdef __DEBUG__
    if( length <= n || n < 0 || posext == NULL )
        throw myruntime_error(
            mystring( "GapScheme: Memory access error." ));
#endif
    return posext[n];
}

// -------------------------------------------------------------------------
// GetOpenAt: accesses gap opening cost at the position
// -------------------------------------------------------------------------

inline
double GapScheme::GetOpenAt( int n ) const
{
#ifdef __DEBUG__
    if( length <= n || n < 0 || vector == NULL )
        throw myruntime_error(
            mystring( "GapScheme: Memory access error." ));
#endif
    return vector[n];
}

// -------------------------------------------------------------------------
// GetExtendAt: accesses gap extend cost at the position
// -------------------------------------------------------------------------

inline
double GapScheme::GetExtendAt( int n ) const
{
#ifdef __DEBUG__
    if( length <= n || n < 0 || extvec == NULL )
        throw myruntime_error(
            mystring( "GapScheme: Memory access error." ));
#endif
    return extvec[n];
}

// -------------------------------------------------------------------------
// GetWeightsAt: returns weight value at the specified position
// -------------------------------------------------------------------------

inline
double GapScheme::GetWeightsAt( int n ) const
{
#ifdef __DEBUG__
    if( length <= n || n < 0 || weights == NULL )
        throw myruntime_error(
            mystring( "GapScheme: Memory access error." ));
#endif
    return weights[n];
}

// -------------------------------------------------------------------------
// GetDeletesBegAt: returns starting deletion weight at the position
// -------------------------------------------------------------------------

inline
double GapScheme::GetDeletesBegAt( int n ) const
{
#ifdef __DEBUG__
    if( length <= n || n < 0 || deletes[DBeg] == NULL )
        throw myruntime_error(
            mystring( "GapScheme: Memory access error." ));
#endif
    return deletes[DBeg][n];
}

// -------------------------------------------------------------------------
// GetDeletesEndAt: returns final deletion weight at the position
// -------------------------------------------------------------------------

inline
double GapScheme::GetDeletesEndAt( int n ) const
{
#ifdef __DEBUG__
    if( length <= n || n < 0 || deletes[DEnd] == NULL )
        throw myruntime_error(
            mystring( "GapScheme: Memory access error." ));
#endif
    return deletes[DEnd][n];
}

// -------------------------------------------------------------------------
// GetDeletesSlopeAt: returns slope of deletion variability at the position
// -------------------------------------------------------------------------

inline
double GapScheme::GetDeletesSlopeAt( int n ) const
{
#ifdef __DEBUG__
    if( length <= n || n < 0 || deletes[DSlope] == NULL )
        throw myruntime_error(
            mystring( "GapScheme: Memory access error." ));
#endif
    return deletes[DSlope][n];
}

// -------------------------------------------------------------------------
// GetDeletesIntervalAt: returns interval of deletion variability at the
//     position
// -------------------------------------------------------------------------

inline
int GapScheme::GetDeletesIntervalAt( int n ) const
{
#ifdef __DEBUG__
    if( length <= n || n < 0 || deleteInt == NULL )
        throw myruntime_error(
            mystring( "GapScheme: Memory access error." ));
#endif
    return deleteInt[n];
}

// -------------------------------------------------------------------------
// GetDeleteExtendWeightAt: returns extend cost for deletion at the given
//     position and time moment of extend
// -------------------------------------------------------------------------

inline
double GapScheme::GetDeleteExtendWeightAt( int n, int time ) const
{
#ifdef __DEBUG__
    if( length <= n || n < 0 ||
        deletes[DBeg] == NULL ||
        deletes[DEnd] == NULL ||
        deletes[DSlope] == NULL ||
        deleteInt == NULL
    )
        throw myruntime_error(
            mystring( "GapScheme: Memory access error." ));
#endif
    double  beg = deletes[DBeg][n];
    double  end = deletes[DEnd][n];
    double  slope = deletes[DSlope][n];
    int     delta = deleteInt[n];

    if( delta <= 0 )
        return 0.0;

    if( time <= 0 )
        return beg;

    if( delta <= time )
        return end;

    return slope * ( double ) time + beg;
}

// -------------------------------------------------------------------------
// GetProbabilityAt: returns gap probability at the given position
// -------------------------------------------------------------------------
inline
double GapScheme::GetProbabilityAt( int pos ) const
{
    double  factor01 = GetProbabilityFactor();
    double  gapweight = GetWeightsAt( pos );
    double  corweight = gapweight * factor01;   //corrected weight
    double  gapprobab = 2.0 * (( corweight <= 0.5 )? corweight: 1.0 - corweight );
#ifdef __DEBUG__
    if( gapprobab < 0.0 || 1.0 < gapprobab )
        throw myruntime_error(
            mystring( "GapScheme: Gap probability inconsistent." ));
#endif
    return  gapprobab;
}

// -------------------------------------------------------------------------
// GetDeleteOpenProbAt: returns deletion opening probability at the given
//     position
// -------------------------------------------------------------------------
inline
double GapScheme::GetDeleteOpenProbAt( int pos ) const
{
    double  factor01 = GetProbabilityFactor();
    double  opnweight = GetDeleteOpenWeightAt( pos ) * factor01;
    double  delopenpr = 2.0 * (( opnweight <= 0.5 )? opnweight: 1.0 - opnweight );
#ifdef __DEBUG__
    if( delopenpr < 0.0 || 1.0 < delopenpr )
        throw myruntime_error(
            mystring( "GapScheme: Deletion open probability inconsistent." ));
#endif
    return  delopenpr;
}

// -------------------------------------------------------------------------
// GetDeleteExtensionProbAt: returns extension probability of deletion at
//     the given position and time moment
// -------------------------------------------------------------------------
inline
double GapScheme::GetDeleteExtensionProbAt( int pos, int time ) const
{
    double  factor01 = GetProbabilityFactor();
    double  extweight = GetDeleteExtendWeightAt( pos, time ) * factor01;
    double  delprob = 2.0 * (( extweight <= 0.5 )? extweight: 1.0 - extweight );
#ifdef __DEBUG__
    if( delprob < 0.0 || 1.0 < delprob )
        throw myruntime_error(
            mystring( "GapScheme: Deletion extension probability inconsistent." ));
#endif
    return  delprob;
}

// -------------------------------------------------------------------------
// GetPosACcorrectionAt: gets correction value for autocorrelation
//     function at position pos
// -------------------------------------------------------------------------
inline
double GapScheme::GetPosACcorrectionAt( int pos ) const
{
#ifdef __DEBUG__
    if( length <= pos || pos < 0 || positionalaccorr == NULL )
        throw myruntime_error(
            mystring( "GapScheme: Memory access error." ));
#endif
    return positionalaccorr[pos];
}

// /////////////////////////////////////////////////////////////////////////
// Set methods
//
// SetOrgTrProbsAt: set specified target transition probability at the
//     position
//
inline
void GapScheme::SetOrgTrProbsAt( double value, int trans, int pos )
{
    if( orgtrobs_ == NULL || 
        length <= pos || pos < -1 || P_NSTATES <= trans || trans < 0 )
        throw myruntime_error( "GapScheme: Memory access error." );
    orgtrobs_[pos+1][trans] = value;
}

// -------------------------------------------------------------------------
// SetOrgTrProbsAt: set target transition probabilities at the position
//
inline
void GapScheme::SetOrgTrProbsAt( const double( *values )[P_NSTATES], int pos )
{
    if( orgtrobs_ == NULL || values == NULL || length <= pos || pos < -1 ||
        orgtrobs_[pos+1] == NULL || *values == NULL )
        throw myruntime_error( "GapScheme: Memory access error." );
    memcpy( orgtrobs_[pos+1], *values, sizeof( double ) * P_NSTATES );
}


// -------------------------------------------------------------------------
// SetTransProbsAt: set specified recalc. target transition 
//     probability at the position
//
inline
void GapScheme::SetTransProbsAt( double value, int trans, int pos )
{
    if( trnprobs_ == NULL || 
        length <= pos || pos < -1 || P_NSTATES <= trans || trans < 0 )
        throw myruntime_error( "GapScheme: Memory access error." );
    trnprobs_[pos+1][trans] = value;
}

// -------------------------------------------------------------------------
// SetTransProbsAt: set recalc. target transition probabilities at the 
//     position
//
inline
void GapScheme::SetTransProbsAt( const double( *values )[P_NSTATES], int pos )
{
    if( trnprobs_ == NULL || values == NULL || length <= pos || pos < -1 ||
        trnprobs_[pos+1] == NULL || *values == NULL )
        throw myruntime_error( "GapScheme: Memory access error." );
    memcpy( trnprobs_[pos+1], *values, sizeof( double ) * P_NSTATES );
}


// -------------------------------------------------------------------------
// SetLogTransProbsAt: set log of specified target transition 
//     probability at the position
//
inline
void GapScheme::SetLogTransProbsAt( double value, int trans, int pos )
{
    if( logtrnps_ == NULL || 
        length <= pos || pos < -1 || P_NSTATES <= trans || trans < 0 )
        throw myruntime_error( "GapScheme: Memory access error." );
    logtrnps_[pos+1][trans] = value;
}

// -------------------------------------------------------------------------
// SetLogTransProbsAt: set log of target transition 
//     probabilities at the position
//
inline
void GapScheme::SetLogTransProbsAt( const double( *values )[P_NSTATES], int pos )
{
    if( logtrnps_ == NULL || values == NULL || length <= pos || pos < -1 ||
        logtrnps_[pos+1] == NULL || *values == NULL )
        throw myruntime_error( "GapScheme: Memory access error." );
    memcpy( logtrnps_[pos+1], *values, sizeof( double ) * P_NSTATES );
}


// -------------------------------------------------------------------------
// SetNSLogTransProbsAt: set non-scaled log of specified target transition 
//     probability at the position
//
inline
void GapScheme::SetNSLogTransProbsAt( double value, int trans, int pos )
{
    if( nsltrnps_ == NULL || 
        length <= pos || pos < -1 || P_NSTATES <= trans || trans < 0 )
        throw myruntime_error( "GapScheme: Memory access error." );
    nsltrnps_[pos+1][trans] = value;
}

// -------------------------------------------------------------------------
// SetNSLogTransProbsAt: set non-scaled log of target transition 
//     probabilities at the position
//
inline
void GapScheme::SetNSLogTransProbsAt( const double( *values )[P_NSTATES], int pos )
{
    if( nsltrnps_ == NULL || values == NULL || length <= pos || pos < -1 ||
        nsltrnps_[pos+1] == NULL || *values == NULL )
        throw myruntime_error( "GapScheme: Memory access error." );
    memcpy( nsltrnps_[pos+1], *values, sizeof( double ) * P_NSTATES );
}


// -------------------------------------------------------------------------
// SetTransProbsBeg: set beginning target transition probabilities
//
inline
void GapScheme::SetOrgTrProbsBeg( const double( *values )[P_NSTATES])
{
    SetOrgTrProbsAt( values, -1 );
}

// -------------------------------------------------------------------------
// SetPosOpenAt: sets positional gap opening cost at the position
//     (upper bound)
// -------------------------------------------------------------------------
inline
void GapScheme::SetPosOpenAt( double value, int n )
{
#ifdef __DEBUG__
    if( length <= n || n < 0 || posvec == NULL )
        throw myruntime_error(
            mystring( "GapScheme: Memory access error." ));
#endif
    posvec[n] = value;
}

// -------------------------------------------------------------------------
// SetPosExtendAt: sets positional gap extension cost at the position
//     (upper bound)
// -------------------------------------------------------------------------

inline
void GapScheme::SetPosExtendAt( double value, int n )
{
#ifdef __DEBUG__
    if( length <= n || n < 0 || posext == NULL)
        throw myruntime_error(
            mystring( "GapScheme: Memory access error." ));
#endif
    posext[n] = value;
}

// -------------------------------------------------------------------------
// SetOpenAt: sets gap open cost at the position
// -------------------------------------------------------------------------

inline
void GapScheme::SetOpenAt( double value, int n )
{
#ifdef __DEBUG__
    if( length <= n || n < 0 || vector == NULL )
        throw myruntime_error(
            mystring( "GapScheme: Memory access error." ));
#endif
    vector[n] = value;
}

// -------------------------------------------------------------------------
// SetExtendAt: sets cost to extend a gap
// -------------------------------------------------------------------------

inline
void GapScheme::SetExtendAt( double value, int n )
{
#ifdef __DEBUG__
    if( length <= n || n < 0 || extvec == NULL )
        throw myruntime_error(
            mystring( "GapScheme: Memory access error." ));
#endif
    extvec[n] = value;
}

// -------------------------------------------------------------------------
// SetWeightsAt: sets gap weight at the position
// -------------------------------------------------------------------------

inline
void GapScheme::SetWeightsAt( double value, int n )
{
#ifdef __DEBUG__
    if( length <= n || n < 0 || weights == NULL )
        throw myruntime_error(
            mystring( "GapScheme: Memory access error." ));
#endif
    weights[n] = value;
}

// -------------------------------------------------------------------------
// SetDeletesBegAt: sets starting deletion weight at the position
// -------------------------------------------------------------------------

inline
void GapScheme::SetDeletesBegAt( double value, int n )
{
#ifdef __DEBUG__
    if( length <= n || n < 0 || deletes[DBeg] == NULL )
        throw myruntime_error(
            mystring( "GapScheme: Memory access error." ));
#endif
    deletes[DBeg][n] = value;
}

// -------------------------------------------------------------------------
// SetDeletesEndAt: sets final deletion weight at the position
// -------------------------------------------------------------------------

inline
void GapScheme::SetDeletesEndAt( double value, int n )
{
#ifdef __DEBUG__
    if( length <= n || n < 0 || deletes[DEnd] == NULL )
        throw myruntime_error(
            mystring( "GapScheme: Memory access error." ));
#endif
    deletes[DEnd][n] = value;
}

// -------------------------------------------------------------------------
// ComputeDeletesSlopeAt: compute slope of deletion variability at the position
// -------------------------------------------------------------------------

inline
void GapScheme::ComputeDeletesSlopeAt( int n )
{
#ifdef __DEBUG__
    if( length <= n || n < 0 || deletes[DSlope] == NULL )
        throw myruntime_error(
            mystring( "GapScheme: Memory access error." ));
#endif
    double  beg = GetDeletesBegAt( n );
    double  end = GetDeletesEndAt( n );
    int     delta = GetDeletesIntervalAt( n );

    if( delta <= 0 ) {
        deletes[DSlope][n] = 0.0;
        return;
    }

    deletes[DSlope][n] = ( end - beg ) / ( double ) delta;
}

// -------------------------------------------------------------------------
// SetDeletesIntervalAt: sets interval of deletion variability at the
//     position
// -------------------------------------------------------------------------

inline
void GapScheme::SetDeletesIntervalAt( int value, int n )
{
#ifdef __DEBUG__
    if( length <= n || n < 0 || deleteInt == NULL )
        throw myruntime_error(
            mystring( "GapScheme: Memory access error." ));
#endif
    deleteInt[n] = value;
}

// -------------------------------------------------------------------------
// SetDeleteStateAt: sets delete-state information at the position including
//     slope computation
// -------------------------------------------------------------------------

inline
void GapScheme::SetDeleteStateAt( double beg, double end, int delta, int n )
{
#ifdef __DEBUG__
    if( length <= n || n < 0 ||
        deletes[DBeg] == NULL ||
        deletes[DEnd] == NULL ||
        deletes[DSlope] == NULL ||
        deleteInt == NULL
    )
        throw myruntime_error(
            mystring( "GapScheme: Memory access error." ));
#endif
    deletes[DBeg][n] = beg;
    deletes[DEnd][n] = end;
    deleteInt[n] = delta;

    if( delta <= 0 ) {
        deletes[DSlope][n] = 0.0;
        return;
    }

    deletes[DSlope][n] = ( end - beg ) / ( double ) delta;
}

// -------------------------------------------------------------------------
// SetPosACcorrectionAt: sets correction value for autocorrelation
//     function at position pos
// -------------------------------------------------------------------------
inline
void GapScheme::SetPosACcorrectionAt( double value, int pos )
{
#ifdef __DEBUG__
    if( length <= pos || pos < 0 || positionalaccorr == NULL )
        throw myruntime_error(
            mystring( "GapScheme: Memory access error." ));
#endif
    positionalaccorr[pos] = value;
}

// -------------------------------------------------------------------------
// IsCompatible: verifies whether the vector containing gap opening costs is 
//     compatible with another one
// -------------------------------------------------------------------------

inline
bool GapScheme::IsCompatible( const GapScheme& one ) const
{
    return  GetOpenCost()   == one.GetOpenCost() &&
            GetExtendCost() == one.GetExtendCost();
}


#endif//__GapScheme__
