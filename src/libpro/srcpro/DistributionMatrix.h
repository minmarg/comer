/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/



#ifndef __DistributionMatrix__
#define __DistributionMatrix__

#include "compdef.h"
#include "debug.h"
#include "rc.h"
#include "data.h"
#include "datapro.h"
#include "libHDP/HDPbase.h"


class Serializer;
class GapScheme;

// _________________________________________________________________________
// Class DistributionMatrix
//
class DistributionMatrix
{
public:
    DistributionMatrix();
    virtual ~DistributionMatrix();

    int             GetColumns() const { return columns; }
    double          operator()( int m, int a ) const    { return GetValueAt( m, a ); }
    char            operator[]( int m ) const           { return GetResidueAt( m );  }

    double&         operator()( int m, int a );         //modification of values
    char&           operator[]( int m );                //modification of aacids

    char            GetResidueAt( int m ) const;        //returns amino acid character at specified position
    double          GetValueAt( int m, int a ) const;   //returns value at specified position for specified amino acid type
    const double ( *GetVectorAt( int m ) const )[NUMALPH];
    const double ( *GetVector() const )[NUMALPH] { return values; }

    virtual void    Push( const double posvalues[NUMALPH], char );                  //push vector of values
    virtual void    PushAt( const double posvalues[NUMALPH], char, int pos );       //push vector of values at the given position

    virtual void    Serialize( Serializer& ) const;
    virtual void    Deserialize( Serializer& );

    // checks whether this matrix is compatible with another one
    bool        IsCompatible( const DistributionMatrix& one ) const;
    bool        IsCompatible( const GapScheme& goc ) const;

    virtual void    OutputMatrix( const char* = NULL ) const;

    void            Reserve( int amount )   { reallocate( amount ); }

    void            CheckForAllZeros();                     //verify wether extists positions with all values of zero

    virtual void    Clear();                                //erase all information contained in this class

    const char*     GetResidues() const     { return aacids; }

protected:
    double          (*GetVector())[NUMALPH] { return values; }

    virtual void    destroy();                              //deallocate memory and reset values
    virtual void    reallocate( int howmuch );              //memory allocation
    virtual void    init();                                 //initialization method

    void            SetColumns( int col )   { columns = col; }
    void            CheckIntegrity() const;                 //check integrity of the structure

protected:
    double   ( *values )[NUMALPH];      //matrix of values
    char*       aacids;                 //sequence of amino acids for which profile was computed
    int         columns;                //number of columns of matrix
    int         allocated;              //how many positions allocated

};

// _________________________________________________________________________
// Class ExtendedDistributionMatrix
//
class ExtendedDistributionMatrix: public DistributionMatrix
{
public:
    ExtendedDistributionMatrix();
    virtual ~ExtendedDistributionMatrix();
                                    //push vector of values by appending to the end and by inserting at the given position
    virtual void    Push( const double posvalues[NUMALPH], char, double weight, double info, double[PS_NSTATES]);
    virtual void    PushAt( const double posvalues[NUMALPH], char, double weight, double info, double[PS_NSTATES], int pos );
    virtual void    PushAt( const double posvalues[NUMALPH], char, int pos );


    //{{
    double          GetBckPPProbAt( int n ) const;
    void            SetBckPPProbAt( double value, int n );
    const double*   GetPPProbsAt( int n ) const;
    const int*      GetPPPIndsAt( int n ) const;
    void            SetPPProbsAt( int n, const double* probs, const int* ndxs, int size );
    size_t          GetNoPPProbsAt( int n ) const;
    //}}

    //{{
    bool            GetctPsSet() const { return ctpsset_; }
    void            SetctPsSet( bool value ) { ctpsset_ = value; }
    double          GetctBckPPProbAt( int n ) const;
    void            SetctBckPPProbAt( double value, int n );
    const double*   GetctPPProbsAt( int n ) const;
    const int*      GetctPPPIndsAt( int n ) const;
    void            SetctPPProbsAt( int n, const double* probs, const int* ndxs, int size );
    size_t          GetNoctPPProbsAt( int n ) const;
    //}}

    //{{methods for processing context vectors
    bool            GetCtxVecSet() const { return ctxvecset_; }
    void            SetCtxVecSet( bool value ) { ctxvecset_ = value; }
    double          GetLogPProbFactor() const { return lppfact_; }
    void            CalcLogPProbFactor();
    double          GetCtxVecNorm2At( int n ) const;
    double          GetCtxVecLpprobAt( int n ) const;
    const double*   GetCtxVecAt( int n ) const;
    int             GetCtxVecSize() const { return ctxvecsize_; }
    void            SetCtxVecPlusAt( int n, double norm2, double pprob, const double* vec, int vsize );
    //}}

    //{{
    bool            GetSSSSet() const { return sssset_; }
    void            SetSSSSet( bool value ) { sssset_ = value; }
    bool            GetSSSP3() const { return sssp3_; }
    void            SetSSSP3( bool value ) { sssp3_ = value; }
    char            GetSSStateAt( int n ) const;
    double          GetSSStateProbAt( int n ) const;
    double          GetSSStateProbAt( int n, char st ) const;
    void            SetSSStateAt( int n, char st, double prob );
    void            SetSSStateAt( int n, char st, double[SS_NSTATES]);
    //}}

    void            Finalize();
    void            CalcTrgFrequencies();
    void            CalcScores();
    void            MixTrgFrequencies( const HDPbase* );
    void            CalcTFPostPredictives( const HDPbase* );
    void            Calc_ctTFPostPredictives( const HDPbase* );
    void            CalcCtxVector();

    virtual void    Serialize( Serializer& ) const;
    virtual void    Deserialize( Serializer& );
    virtual void    OutputMatrix( const char* = NULL ) const;

    void            SetNoSequences( size_t value )          { nosequences = value; }
    void            SetEffNoSequences( double value )       { effnosequences = value; }

    void            SetRefLambda( double value )            { referenceLambda = value; }
    void            SetRefK( double value )                 { referenceK = value; }
    void            SetLambda( double value )               { lambda = value; }
    void            SetEntropy( double value )              { entropy = value; }
    void            SetK( double value )                    { parameterK = value; }
    void            SetExpectedScore( double value )        { expscore = value; }


    size_t          GetNoSequences() const          { return nosequences; }
    double          GetEffNoSequences() const       { return effnosequences; }

    const double*   GetBackProbs() const { return backprobs_; }
    double          GetBackProbsAt( char res ) const;
    void            SetBackProbsAt( char res, double value );

    const double*   GetPostProbs() const { return postprobs_; }
    double          GetPostProbsAt( char res ) const;
    void            SetPostProbsAt( char res, double value );

    double          GetRefLambda() const            { return referenceLambda; }
    double          GetRefK() const                 { return referenceK; }
    double          GetLambda() const               { return lambda; }
    double          GetEntropy() const              { return entropy; }
    double          GetK() const                    { return parameterK; }
    double          GetExpectedScore() const        { return expscore; }

    const double ( *GetTrgFreqsAt( int m ) const )[NUMAA];
    const double ( *GetTrgFreqs() const )[NUMAA] { return trgfrqs_; }

    const double*   GetInformation() const          { return information; }

    double          GetFrequencyWeightAt( int ) const;      //alpha
    double          GetInformationAt( int ) const;

    double          GetMIDExpNoObservationsBeg( int st ) const { return GetMIDExpNoObservationsAt( -1, st ); }
    double          GetMIDExpNoObservationsAt( int, int st ) const;
    void            SetMIDExpNoObservationsBeg( double vv[PS_NSTATES]) { SetMIDExpNoObservationsAt( -1, vv ); }
    void            SetMIDExpNoObservationsAt( int, double[PS_NSTATES]);

    size_t          GetNameSize() const         { return szname; }
    size_t          GetDescriptionSize() const  { return szdescription; }

    const char*     GetName() const         { return name; }
    const char*     GetDescription() const  { return description; }

    void            SetNameSize( size_t sz )        { szname = sz; }
    void            SetDescriptionSize( size_t sz ) { szdescription = sz; }

    void            SetName( const char* );                 //sets name
    void            SetDescription( const char* );          //sets description

    void            PrintAnnotation( char* sp ) const;      //print short annotation to string stream
    void            PrintAnnotation( FILE* fp ) const;      //print short annotation to file
    void            PrintAnnotation( TPrintFunction, void* vpn ) const;//print short annotation of the pssm

    void            PrintDescriptionFirst( FILE* fp ) const;//format and print to file
    void            PrintDescriptionFirst( TPrintFunction, void* ) const;//format and print name and description

    void            PrintDescription( char* sp ) const;     //format and print to string stream
    void            PrintDescription( FILE* fp ) const;     //format and print to file
    void            PrintDescription( TPrintFunction, void* vpn ) const;//format and print name and description

    size_t          GetMaxAnnotationWidth() const;                  //maximum annotation width
    size_t          GetMinimumRequiredSizeForDescription() const;   //minimum size required to contain description

    void            PrintParameters( FILE* ) const;         //print statistical parameters

    virtual void    Clear();                                //erase all information contained in this class

protected:
    double          (*GetTrgFreqs())[NUMAA] { return trgfrqs_; }
    double*         GetInformation() { return information; }

    void    PrintDescriptionHelper(                         //helper method for formating and printing the description
        TPrintFunction print_func, void* vpn,
        size_t preamble, size_t textwidth, size_t width,
        size_t max_rows, size_t max_length, bool annotation ) const;


    virtual void    destroy();                              //deallocate memory and reset values
    virtual void    reallocate( int howmuch );              //memory allocation
    virtual void    init();                                 //initialization method

    size_t          GetPrivateBufferSize() const    { return MAX_DESCRIPTION_LENGTH; }
    char*           GetPrivateBuffer() const        { return private_buffer; }

    void            FormatBuffer( char*&, const char*,      //auxiliary method to format character buffer
                size_t&, size_t&, const size_t, const size_t, const size_t ) const;

private:
    double   ( *trgfrqs_ )[NUMAA];      //target frequencies
    double*     freqweights;            //frequency weights at each position (known as alpha coefficient)
    double*     information;            //information content vector
    double   ( *expMIDs_ )[PS_NSTATES]; //expected number of MID state observations
    //
    double*     bppprob_;               //background posterior predictive
    double**    ppprobs_;               //posterior predictive probabilities for each cluster
    int**       pppndxs_;               //indices of clusters p.p.probabilities were calculated for
    size_t*     noppps_;                //number of p.p.probability values
    //
    bool        ctpsset_;               //HDP ctx probabilities set
    double*     ctbppprob_;             //background posterior predictive
    double**    ctppprobs_;             //posterior predictive probabilities for each cluster
    int**       ctpppndxs_;             //indices of clusters p.p.probabilities were calculated for
    size_t*     noctppps_;              //number of p.p.probability values
    //
    bool        ctxvecset_;             //whether ctx vector calculation is on
    double      lppfact_;               //log of prior probability factor used in processing context vectors
    double*     ctxvecnorm2_;           //squared norm of vector at each position
    double*     ctxveclpprb_;           //log of prior probabilities at each position
    double**    ctxvecs_;               //context vectors at each position
    int         ctxvecsize_;            //size of context vectors
    //
    bool        sssset_;                //SS states are set
    bool        sssp3_;                 //probabilities for all SS states are set
    char*       sss_;                   //SS state
    double    (*sssprob_)[SS_NSTATES];  //probabilities for SS states
    //
    char*       name;                   //name of the multiple alignment this matrix was transformed from
    char*       description;            //description of the representative sequence of the multiple alignment
    size_t      szname;                 //size of name
    size_t      szdescription;          //size of description
                                        //character buffer for description string storage
    static char private_buffer[MAX_DESCRIPTION_LENGTH];

    size_t      nosequences;            //number of sequences in alignment
    double      effnosequences;         //effective number of sequences

    double      backprobs_[NUMALPH];    //background probabilities
    double      postprobs_[NUMALPH];    //posterior (generalized target) probabilities

    double      referenceLambda;        //reference lambda parameter
    double      referenceK;             //reference parameter K

    double      lambda;                 //computed statistical parameter, lambda
    double      entropy;                //computed entropy given lambda
    double      parameterK;             //computed Karlin's parameter K
    double      expscore;               //expected score per column pair

    friend class RndMSAGenerator;
};



typedef DistributionMatrix          FrequencyMatrix;
typedef ExtendedDistributionMatrix  LogOddsMatrix;

// -------------------------------------------------------------------------
// Functions for output of profile data in the text format
//
static const int    gcpDMscaling = 65536;

void OutputProfile( const char* filename, const FrequencyMatrix&, const LogOddsMatrix&, const GapScheme& );
void TextWriteProfile( FILE*, const FrequencyMatrix&, const LogOddsMatrix&, const GapScheme&, int = gcpDMscaling );
void TextReadProfile( FILE*, FrequencyMatrix&, LogOddsMatrix&, GapScheme& );

// -------------------------------------------------------------------------
// GetValue: used to access score value at the specified position and for
//     the specified amino acid
// -------------------------------------------------------------------------

inline
double DistributionMatrix::GetValueAt( int m, int a ) const
{
#ifdef __DEBUG__
    if( columns <= m || m < 0 )
        throw myruntime_error(
            mystring( "DistributionMatrix: Memory access error." ));

    if( NUMALPH <= a || a < 0 )
        throw myruntime_error(
            mystring( "DistributionMatrix: Memory access error." ));
#endif
    return values[m][a];
}

// -------------------------------------------------------------------------
// GetResidue: returns amino acid at the specified position
// -------------------------------------------------------------------------

inline
char DistributionMatrix::GetResidueAt( int m ) const
{
#ifdef __DEBUG__
    if( columns <= m || m < 0 )
        throw myruntime_error(
            mystring( "DistributionMatrix: Memory access error." ));
#endif
    return aacids[m];
}

// -------------------------------------------------------------------------
// GetVectorAt: Return vector of values at the given position
// -------------------------------------------------------------------------

inline
const double ( *DistributionMatrix::GetVectorAt( int m ) const )[NUMALPH]
{
#ifdef __DEBUG__
    if( !values || columns <= m || m < 0 )
        throw myruntime_error(
            mystring( "DistributionMatrix: Memory access error." ));
#endif
    return values + m;
}
// -------------------------------------------------------------------------
// operator(): used to modify score value at the specified position and for
//     the specified amino acid
// -------------------------------------------------------------------------

inline
double& DistributionMatrix::operator()( int m, int a )
{
#ifdef __DEBUG__
    if( columns <= m || m < 0 )
        throw myruntime_error(
            mystring( "DistributionMatrix: Memory access error." ));

    if( NUMALPH <= a || a < 0 )
        throw myruntime_error(
            mystring( "DistributionMatrix: Memory access error." ));
#endif
    return values[m][a];
}

// -------------------------------------------------------------------------
// operator[]: returns amino acid to be modified
// -------------------------------------------------------------------------

inline
char& DistributionMatrix::operator[]( int m )
{
#ifdef __DEBUG__
    if( columns <= m || m < 0 )
        throw myruntime_error(
            mystring( "DistributionMatrix: Memory access error." ));
#endif
    return aacids[m];
}

////////////////////////////////////////////////////////////////////////////
// ExtendedDistributionMatrix inlines
//
// GetTrgFreqsAt: return vector of target frequency values at given position
//
inline
const double ( *ExtendedDistributionMatrix::GetTrgFreqsAt( int m ) const )[NUMAA]
{
#ifdef __DEBUG__
    if( !trgfrqs_ || columns <= m || m < 0 )
        throw myruntime_error("ExtendedDistributionMatrix: GetTrgFreqsAt: Memory access error." );
#endif
    return trgfrqs_ + m;
}

// -------------------------------------------------------------------------
// GetFrequencyWeightAt: get frequency weight at specified position
//
inline
double ExtendedDistributionMatrix::GetFrequencyWeightAt( int m ) const
{
#ifdef __DEBUG__
    if( columns <= m || m < 0 )
        throw myruntime_error(
            mystring( "ExtendedDistributionMatrix: Memory access error." ));
#endif
    return freqweights[m];
}

// GetInformationAt: returns information content at the position
// -------------------------------------------------------------------------

inline
double ExtendedDistributionMatrix::GetInformationAt( int m ) const
{
#ifdef __DEBUG__
    if( columns <= m || m < 0 )
        throw myruntime_error(
            mystring( "ExtendedDistributionMatrix: Memory access error." ));
#endif
    return information[m];
}

// -------------------------------------------------------------------------
// GetMIDExpNoObservationsAt: return expected number of MID state
//     observations at the position
//
inline
double ExtendedDistributionMatrix::GetMIDExpNoObservationsAt( int m, int st ) const
{
#ifdef __DEBUG__
    if( m < -1 || columns <= m || st < 0 || PS_NSTATES <= st )
        throw myruntime_error( "ExtendedDistributionMatrix: Memory access error." );
#endif
    return expMIDs_[m+1][st];
}

// -------------------------------------------------------------------------
// SetMIDExpNoObservationsAt: set expected number of MID state
//     observations at the position
//
inline
void ExtendedDistributionMatrix::SetMIDExpNoObservationsAt( int m, double mids[PS_NSTATES])
{
#ifdef __DEBUG__
    if( m < -1 || columns <= m )
        throw myruntime_error( "ExtendedDistributionMatrix: Memory access error." );
#endif
    expMIDs_[m+1][PS_M] = mids[PS_M];
    expMIDs_[m+1][PS_I] = mids[PS_I];
    expMIDs_[m+1][PS_D] = mids[PS_D];
}

// -------------------------------------------------------------------------
// GetBackProbsAt: get background probability of residue res
//
inline
double ExtendedDistributionMatrix::GetBackProbsAt( char res ) const
{
    if( res < 0 || NUMALPH <= res )
        throw myruntime_error( "ExtendedDistributionMatrix: Memory access error." );
    return backprobs_[res];
}

// SetBackProbsAt: set background probability for residue res
//
inline
void ExtendedDistributionMatrix::SetBackProbsAt( char res, double value )
{
    if( res < 0 || NUMALPH <= res )
        throw myruntime_error( "ExtendedDistributionMatrix: Memory access error." );
    backprobs_[res] = value;
}

// -------------------------------------------------------------------------
// GetPostProbsAt: get posterior probability of residue res
//
inline
double ExtendedDistributionMatrix::GetPostProbsAt( char res ) const
{
    if( res < 0 || NUMALPH <= res )
        throw myruntime_error( "ExtendedDistributionMatrix: GetPostProbsAt: Memory access error." );
    return postprobs_[res];
}

// SetPostProbsAt: set posterior probability for residue res
//
inline
void ExtendedDistributionMatrix::SetPostProbsAt( char res, double value )
{
    if( res < 0 || NUMALPH <= res )
        throw myruntime_error( "ExtendedDistributionMatrix: SetPostProbsAt: Memory access error." );
    postprobs_[res] = value;
}

// -------------------------------------------------------------------------
// GetBckPPProbAt: get background posterior predictive probability at the 
//  position
//
inline
double ExtendedDistributionMatrix::GetBckPPProbAt( int n ) const
{
#ifdef __DEBUG__
    if( !bppprob_ || n < 0 || columns <= n )
        throw myruntime_error("ExtendedDistributionMatrix: GetBckPPProbAt: Memory access error.");
#endif
    return bppprob_[n];
}

// SetBckPPProbAt: set background posterior predictive probability at the 
//  position
inline
void ExtendedDistributionMatrix::SetBckPPProbAt( double value, int n )
{
#ifdef __DEBUG__
    if( !bppprob_ || n < 0 || columns <= n )
        throw myruntime_error("ExtendedDistributionMatrix: GetBckPPProbAt: Memory access error.");
#endif
    bppprob_[n] = value;
}

// GetPPProbsAt: get posterior predictive probabilities at the position
inline
const double* ExtendedDistributionMatrix::GetPPProbsAt( int n ) const
{
#ifdef __DEBUG__
    if( !ppprobs_ || n < 0 || columns <= n )
        throw myruntime_error("ExtendedDistributionMatrix: GetPPProbsAt: Memory access error.");
#endif
    return ppprobs_[n];
}

// GetPPPIndsAt: get indices of posterior predictive probabilities at the 
//  position
inline
const int* ExtendedDistributionMatrix::GetPPPIndsAt( int n ) const
{
#ifdef __DEBUG__
    if( !pppndxs_ || n < 0 || columns <= n )
        throw myruntime_error("ExtendedDistributionMatrix: GetPPPIndsAt: Memory access error.");
#endif
    return pppndxs_[n];
}

// SetPPProbsAt: set posterior predictive probabilities and their 
//  indices at the position
inline
void ExtendedDistributionMatrix::SetPPProbsAt( int n, const double* probs, const int* ndxs, int size )
{
    int k;
    if( !ppprobs_ || !pppndxs_ || !noppps_ || n < 0 || columns <= n )
        throw myruntime_error("ExtendedDistributionMatrix: SetPPProbsAt: Memory access error.");
    if( size < 0 )
        throw myruntime_error("ExtendedDistributionMatrix: SetPPProbsAt: Invalid size of data.");
    if( ppprobs_[n]) { free( ppprobs_[n]); ppprobs_[n] = NULL; }
    if( pppndxs_[n]) { free( pppndxs_[n]); pppndxs_[n] = NULL; }
    noppps_[n] = (size_t)size;
    if( size < 1 )
        return;
    ppprobs_[n] = ( double* )malloc( size * sizeof( double ));
    pppndxs_[n] = ( int* )malloc( size * sizeof( int ));
    if( ppprobs_[n] == NULL || pppndxs_[n] == NULL )
        throw myruntime_error("ExtendedDistributionMatrix: SetPPProbsAt: Not enough memory.");
    if( probs == NULL || ndxs == NULL )
        throw myruntime_error("ExtendedDistributionMatrix: SetPPProbsAt: Null parameters.");
    for( k = 0; k < size; k ++ ) {
        ppprobs_[n][k] = probs[k];
        pppndxs_[n][k] = ndxs[k];
    }
}

// GetNoPPProbsAt: get number of posterior predictive probabilities at the 
//  position
inline
size_t ExtendedDistributionMatrix::GetNoPPProbsAt( int n ) const
{
#ifdef __DEBUG__
    if( !noppps_ || n < 0 || columns <= n )
        throw myruntime_error("ExtendedDistributionMatrix: GetNoPPProbsAt: Memory access error.");
#endif
    return noppps_[n];
}

// =========================================================================
// GetctBckPPProbAt: get background posterior predictive probability at the 
//  position
//
inline
double ExtendedDistributionMatrix::GetctBckPPProbAt( int n ) const
{
#ifdef __DEBUG__
    if( !ctbppprob_ || n < 0 || columns <= n )
        throw myruntime_error("ExtendedDistributionMatrix: GetctBckPPProbAt: Memory access error.");
#endif
    return ctbppprob_[n];
}

// SetctBckPPProbAt: set background posterior predictive probability at the 
//  position
inline
void ExtendedDistributionMatrix::SetctBckPPProbAt( double value, int n )
{
#ifdef __DEBUG__
    if( !ctbppprob_ || n < 0 || columns <= n )
        throw myruntime_error("ExtendedDistributionMatrix: GetctBckPPProbAt: Memory access error.");
#endif
    ctbppprob_[n] = value;
}

// GetctPPProbsAt: get posterior predictive probabilities at the position
inline
const double* ExtendedDistributionMatrix::GetctPPProbsAt( int n ) const
{
#ifdef __DEBUG__
    if( !ctppprobs_ || n < 0 || columns <= n )
        throw myruntime_error("ExtendedDistributionMatrix: GetctPPProbsAt: Memory access error.");
#endif
    return ctppprobs_[n];
}

// GetctPPPIndsAt: get indices of posterior predictive probabilities at the 
//  position
inline
const int* ExtendedDistributionMatrix::GetctPPPIndsAt( int n ) const
{
#ifdef __DEBUG__
    if( !ctpppndxs_ || n < 0 || columns <= n )
        throw myruntime_error("ExtendedDistributionMatrix: GetctPPPIndsAt: Memory access error.");
#endif
    return ctpppndxs_[n];
}

// SetctPPProbsAt: set posterior predictive probabilities and their 
//  indices at the position
inline
void ExtendedDistributionMatrix::SetctPPProbsAt( int n, const double* probs, const int* ndxs, int size )
{
    int k;
    if( !ctppprobs_ || !ctpppndxs_ || !noctppps_ || n < 0 || columns <= n )
        throw myruntime_error("ExtendedDistributionMatrix: SetctPPProbsAt: Memory access error.");
    if( size < 0 )
        throw myruntime_error("ExtendedDistributionMatrix: SetctPPProbsAt: Invalid size of data.");
    if( ctppprobs_[n]) { free( ctppprobs_[n]); ctppprobs_[n] = NULL; }
    if( ctpppndxs_[n]) { free( ctpppndxs_[n]); ctpppndxs_[n] = NULL; }
    noctppps_[n] = (size_t)size;
    if( size < 1 )
        return;
    ctppprobs_[n] = ( double* )malloc( size * sizeof( double ));
    ctpppndxs_[n] = ( int* )malloc( size * sizeof( int ));
    if( ctppprobs_[n] == NULL || ctpppndxs_[n] == NULL )
        throw myruntime_error("ExtendedDistributionMatrix: SetctPPProbsAt: Not enough memory.");
    if( probs == NULL || ndxs == NULL )
        throw myruntime_error("ExtendedDistributionMatrix: SetctPPProbsAt: Null parameters.");
    for( k = 0; k < size; k ++ ) {
        ctppprobs_[n][k] = probs[k];
        ctpppndxs_[n][k] = ndxs[k];
    }
}

// GetNoctPPProbsAt: get number of posterior predictive probabilities at the 
//  position
inline
size_t ExtendedDistributionMatrix::GetNoctPPProbsAt( int n ) const
{
#ifdef __DEBUG__
    if( !noctppps_ || n < 0 || columns <= n )
        throw myruntime_error("ExtendedDistributionMatrix: GetNoctPPProbsAt: Memory access error.");
#endif
    return noctppps_[n];
}

// =========================================================================
// GetCtxVecNorm2At: get the squared norm of the vector at the position
//
inline
double ExtendedDistributionMatrix::GetCtxVecNorm2At( int n ) const
{
#ifdef __DEBUG__
    if( !ctxvecnorm2_ || n < 0 || columns <= n )
        throw myruntime_error("ExtendedDistributionMatrix: GetCtxVecNorm2At: Memory access error.");
#endif
    return ctxvecnorm2_[n];
}

// GetCtxVecLpprobAt: get the prior probability of the vector at the position
inline
double ExtendedDistributionMatrix::GetCtxVecLpprobAt( int n ) const
{
#ifdef __DEBUG__
    if( !ctxveclpprb_ || n < 0 || columns <= n )
        throw myruntime_error("ExtendedDistributionMatrix: GetCtxVecLpprobAt: Memory access error.");
#endif
    return ctxveclpprb_[n];
}

// GetCtxVecAt: get the vector at the position
inline
const double* ExtendedDistributionMatrix::GetCtxVecAt( int n ) const
{
#ifdef __DEBUG__
    if( !ctxvecs_ || n < 0 || columns <= n )
        throw myruntime_error("ExtendedDistributionMatrix: GetCtxVecAt: Memory access error.");
#endif
    return ctxvecs_[n];
}

// SetCtxVecPlusAt: set vector, its norm and log of prior probability at the position
inline
void ExtendedDistributionMatrix::SetCtxVecPlusAt( int n, 
    double norm2, double pprob, const double* vec, int vsize )
{
    int k;
    if( !ctxvecnorm2_ || !ctxveclpprb_ || !ctxvecs_ || n < 0 || columns <= n )
        throw myruntime_error("ExtendedDistributionMatrix: SetCtxVecPlusAt: Memory access error.");
    if( vsize < 0 )
        throw myruntime_error("ExtendedDistributionMatrix: SetCtxVecPlusAt: Invalid vector size.");
    if( ctxvecs_[n]) { free( ctxvecs_[n]); ctxvecs_[n] = NULL; }
    if( ctxvecsize_ && ctxvecsize_ != vsize )
        throw myruntime_error("ExtendedDistributionMatrix: SetCtxVecPlusAt: Inconsistent vector size.");
    if( !ctxvecsize_ )
        ctxvecsize_ = vsize;
    if( vsize < 1 )
        return;
    if( vec == NULL )
        throw myruntime_error("ExtendedDistributionMatrix: SetCtxVecPlusAt: Null vector given.");
    ctxvecs_[n] = ( double* )malloc( vsize * sizeof( double ));
    if( ctxvecs_[n] == NULL )
        throw myruntime_error("ExtendedDistributionMatrix: SetCtxVecPlusAt: Not enough memory.");
    for( k = 0; k < vsize; k++ )
        ctxvecs_[n][k] = vec[k];
    ctxvecnorm2_[n] = norm2;
    ctxveclpprb_[n] = pprob;
}

// =========================================================================
// GetSSStateAt: get SS state at the position
//
inline
char ExtendedDistributionMatrix::GetSSStateAt( int n ) const
{
    if( !sss_ || n < 0 || columns <= n )
        throw myruntime_error("ExtendedDistributionMatrix: GetSSStateAt: Memory access error.");
    return sss_[n];
}

// GetSSStateProbAt: get SS state probability at the position
//
inline
double ExtendedDistributionMatrix::GetSSStateProbAt( int n ) const
{
    if( !sss_ || !sssprob_ || n < 0 || columns <= n || sss_[n] < 0 || SS_NSTATES <= sss_[n] )
        throw myruntime_error("ExtendedDistributionMatrix: GetSSStateProbAt: Memory access error.");
    return sssprob_[n][sss_[n]];
}
inline
double ExtendedDistributionMatrix::GetSSStateProbAt( int n, char st ) const
{
    if( !sssprob_ || n < 0 || columns <= n || st < 0 || SS_NSTATES <= st )
        throw myruntime_error("ExtendedDistributionMatrix: GetSSStateProbAt: Memory access error.");
    return sssprob_[n][st];
}

// SetSSStateAt: set SS state and probability at the position
inline
void ExtendedDistributionMatrix::SetSSStateAt( int n, char st, double prob )
{
    int k;
    if( !sss_ || !sssprob_ || n < 0 || columns <= n || st < 0 || SS_NSTATES <= st )
        throw myruntime_error("ExtendedDistributionMatrix: SetSSStateAt: Memory access error.");
    sss_[n] = st;
    sssprob_[n][st] = prob;
}
inline
void ExtendedDistributionMatrix::SetSSStateAt( int n, char st, double probs[SS_NSTATES])
{
    int k, s;
    if( !sss_ || !sssprob_ || n < 0 || columns <= n || st < 0 || SS_NSTATES <= st )
        throw myruntime_error("ExtendedDistributionMatrix: SetSSStateAt: Memory access error.");
    sss_[n] = st;
    for( s = 0; s < SS_NSTATES; s++ )
        sssprob_[n][s] = probs[s];
}


#endif//__DistributionMatrix__
