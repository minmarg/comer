/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __DescriptionVector__
#define __DescriptionVector__

#include <stdlib.h>
#include <math.h>


#include "compdef.h"
#include "debug.h"
#include "rc.h"
#include "mystring.h"
#include "myexcept.h"
#include "datapro.h"


#define UNUSED  0
#define USED    1
#define USED_IN_EXTENT  3

//size of epxected number of observations
#define SIZEOFEXPNOBSERVS   ( 400 )

// extern double rint( double x );

// _________________________________________________________________________
// Class PosDescriptionVector
//
class PosDescriptionVector {
public:
    // typedefs ...
    //
    PosDescriptionVector( unsigned reservation );
    virtual ~PosDescriptionVector();

//     unsigned char&      operator[]( int n ) { return residues[n]; };
    unsigned char   operator[]( int n ) const { return residues[n]; };

    size_t          size() const { return length_; }
    size_t          capacity() const { return capacity_; }

    unsigned int    GetCluster() const { return cluster_; }
    void            SetCluster( unsigned int value ) { cluster_ = value; }
    void            IncCluster() { cluster_++; }
    void            ResetCluster() { cluster_ = ( unsigned int ) -1; }


    size_t          GetEffectiveSize() const { return efflength_; }
    void            SetEffectiveSize( size_t value ) { efflength_ = value; }

    void            Reset()     { counter_ = 0; }
    bool            Eos()       { return length_ <= counter_; }
    void            Inc()       { counter_++; }
    void            BackReset() { counter_ = length_ - 1; }
    bool            Geb()       { return 1 <= counter_ + 1; }   //manage situations of possible overflows
    void            Dec()       { --counter_; }

    const char*     GetDescription() const      { return description.c_str(); }
    const unsigned char*    GetResidues() const { return residues; }

    unsigned char   Residue() const         { return ResidueAt( counter_ ); }
    bool            IsUsed () const         { return IsUsedAt ( counter_ ); }
    bool            IsUsedInExtent () const { return IsUsedInExtentAt( counter_ ); }
    double          Weight () const         { return WeightAt ( counter_ ); }

    unsigned char   ResidueAt( size_t n ) const;
    void            SetResidueAt( unsigned char r, size_t );

    bool            IsUsedAt ( size_t n ) const;
    void            SetUnusedAt ( size_t );
    void            SetUsedAt   ( size_t );
    bool            IsUsedInExtentAt ( size_t n ) const;
    void            SetUsedInExtentAt( size_t );
    void            UnsetUsedInExtentAt( size_t );

    bool            IsUsedAndInfAt( size_t p ) const;

    double          WeightAt( size_t n ) const;
    void            SetWeightAt ( double w, size_t );

    double          GetGlbWeight() const { return glbwght_; }
    void            SetGlbWeight( double value ) { glbwght_ = value; }

    void            AppendDescription( const char* desc, size_t pos, size_t len ) { description.append( desc, pos, len ); }

    void            SetResidue( unsigned char r )   { SetResidueAt( r, counter_ ); }
    void            SetUnused ()    { SetUnusedAt( counter_ ); }
    void            SetUsed   ()    { SetUsedAt  ( counter_ ); }


    bool            GetUsed() const { return used; }
    void            SetUsed( bool u ) { used = u; }

    size_t          GetFirstUsed() const { return firstused_; }
    void            SetFirstUsed( size_t value ) { firstused_ = value; }
    size_t          GetLastUsed() const { return lastused_; }
    void            SetLastUsed( size_t value ) { lastused_ = value; }

    virtual void    push( unsigned char r, unsigned char f = USED, double w = 0.0 );
    virtual void    clear();                        //clears all the positions contained in the sequence

protected:
    explicit        PosDescriptionVector( const PosDescriptionVector& );
    explicit        PosDescriptionVector();

    virtual PosDescriptionVector& operator=( const PosDescriptionVector& one );

    virtual void    Realloc( int newcap );
    virtual void    Init();

protected:
    mystring            description;    //string description of vector
    unsigned char*      residues;       //residue codes at the positions
    unsigned char*      flags;          //whether residues is used at the positions
    double*             weights;        //weight of the sequence at the positions
    double              glbwght_;       //global weight
    size_t              length_;        //length of the sequence
    size_t              efflength_;     //effective length of the sequence (length with gaps excluded)
    size_t              capacity_;      //current capacity of the sequence
    bool                used;           //whether or not the sequence is used
    size_t              firstused_;     //index of the first used position
    size_t              lastused_;      //index of the last used position
    //
    size_t              counter_;       //initial counter to iterate over all positions in the sequence
    unsigned int        cluster_;       //cluster number if clustering is in effect
};


// _________________________________________________________________________
// Class ExtendedDescriptionVector
//
class ExtendedDescriptionVector: public PosDescriptionVector {
public:
    // typedefs ...
    //
    enum Direction {    //extent direction type
        xLeft,          //left boundary
        xRight,         //right boundary
        xInterval,      //interval of valid positions in the extent
        xNoSequences,   //number of sequences participating in the extent
        xNoSymbols,     //number of different symbols (sum) occuring in the extent;NOTE:previously
        xCount
    };
    ExtendedDescriptionVector( unsigned reservation );
    ExtendedDescriptionVector( const PosDescriptionVector& sequence );
    virtual ~ExtendedDescriptionVector();

    virtual ExtendedDescriptionVector& operator=( const PosDescriptionVector& sequence );

    //output methods
    void            PrintMatchWeights( FILE* );
    void            PrintTransWeights( FILE* );
    void            PrintTargetTranst( FILE* );
    void            PrintMIDExpNoObservations( FILE* );
    void            PrintPSSMatrix( FILE* );
    void            PrintSuppressedPSSMandWeights( FILE* );
    void            PrintProfile( FILE* );

    //get/set methods...
    size_t*         GetIndicesAt( size_t n ) const;
    size_t          GetIndicesSizeAt( size_t n ) const;
    void            PushIndexAt( size_t value, size_t n );

    const double*   GetBackProbs() const { return backprobs_; }
    double          GetBackProbsAt( unsigned char res ) const;
    void            SetBackProbsAt( unsigned char res, double value );

    const double*   GetPostProbs() const { return postprobs_; }
    double          GetPostProbsAt( unsigned char res ) const;
    void            SetPostProbsAt( unsigned char res, double value );

    int             GetStateAt( size_t n ) const;
    void            SetStateAt( size_t n );



    size_t  GetLeftExtentAt( size_t n ) const  {return GetLeftExtentAt( n, PS_M );}
    void    SetLeftExtentAt( size_t value, size_t n )  {SetLeftExtentAt( n, PS_M, value );}

    size_t  GetRightExtentAt( size_t n ) const  {return GetRightExtentAt( n, PS_M );}
    void    SetRightExtentAt( size_t value, size_t n )  {SetRightExtentAt( n, PS_M, value );}

    size_t  GetExtentIntervalAt( size_t n ) const  {return GetExtentIntervalAt( n, PS_M );}
    void    SetExtentIntervalAt( size_t value, size_t n )  {SetExtentIntervalAt( n, PS_M, value );}

    size_t  GetNoSequencesInExtentAt( size_t n ) const  {return GetNoSequencesInExtentAt( n, PS_M );}
    void    IncNoSequencesInExtentAt( size_t n )  {IncNoSequencesInExtentAt( n, PS_M );}

    double  GetNoSymbolsInExtentAt( size_t n ) const  {return GetNoSymbolsInExtentAt( n, PS_M );}
    void    SetNoSymbolsInExtentAt( double value, size_t n )  {SetNoSymbolsInExtentAt( n, PS_M, value );}


    size_t          GetLeftExtentAt( size_t n, int st ) const;
    void            SetLeftExtentAt( size_t n, int st, size_t value );

    size_t          GetRightExtentAt( size_t n, int st ) const;
    void            SetRightExtentAt( size_t n, int st, size_t value );

    size_t          GetExtentIntervalAt( size_t n, int st ) const;
    void            SetExtentIntervalAt( size_t n, int st, size_t value );

    size_t          GetNoSequencesInExtentAt( size_t n, int st ) const;
    void            IncNoSequencesInExtentAt( size_t n, int st );

    double          GetNoSymbolsInExtentAt( size_t n, int st ) const;
    void            SetNoSymbolsInExtentAt( size_t n, int st, double value );



    size_t          GetLeftMSExtentAt( size_t n ) const;
    void            SetLeftMSExtentAt( size_t value, size_t n );

    size_t          GetRightMSExtentAt( size_t n ) const;
    void            SetRightMSExtentAt( size_t value, size_t n );

    size_t          GetMSExtentIntervalAt( size_t n ) const;
    void            SetMSExtentIntervalAt( size_t value, size_t n );

    size_t          GetNoSequencesInMSExtentAt( size_t n ) const;
    void            IncNoSequencesInMSExtentAt( size_t n );

    size_t          GetNoSymbolsInMSExtentAt( size_t n ) const;
    void            SetNoSymbolsInMSExtentAt( size_t value, size_t n );


    size_t          GetCountAt( size_t n ) const;
    void            IncCountAt( size_t n );

    double          ComputeObsFrequencyWeightAt( size_t n ) const;

    size_t          GetDistributionAt( unsigned char res, size_t n ) const;
    void            IncDistributionAt( unsigned char res, size_t n );

    const double*   GetSqnWeightsAt( size_t p, int st ) const;
    double*         GetSqnWeightsAt( size_t p, int st );
    void            SetSqnWeightsAt( size_t p, int st, size_t pp );
    void            NewSqnWeightsAt( size_t n, int st, size_t size );
    void            FreeSqnWeightsAt( size_t n, int st );
    void            FreeSqnWeightsAt();

    double          GetMatchWeightsAt( unsigned char res, size_t n ) const;
    const double ( *GetMatchWeightsAt( size_t n ) const )[NUMALPH];
    void            SetMatchWeightsAt( double value, unsigned char res, size_t n );
    void            IncMatchWeightsAt( double value, unsigned char res, size_t n );

    double          GetTransWeightsBeg( int trans ) const;
    const double ( *GetTransWeightsBeg() const )[P_NSTATES];
    void            SetTransWeightsBeg( double value, int trans );
    void            IncTransWeightsBeg( double value, int trans );

    double          GetTransWeightsAt( int trans, ssize_t n ) const;
    const double ( *GetTransWeightsAt( ssize_t n ) const )[P_NSTATES];
    void            SetTransWeightsAt( double value, int trans, ssize_t n );
    void            IncTransWeightsAt( double value, int trans, ssize_t n );

    double          GetGapWeightsAt(   size_t n ) const;
    void            SetGapWeightsAt(   double value, size_t n );


    size_t  GetDistinctHistAt( unsigned char res, ssize_t n ) const  {return GetDistinctHistAt(n,PS_M,res);}
    void    SetDistinctHistAt( size_t value, unsigned char res, ssize_t n )  {SetDistinctHistAt(n,PS_M,res,value);}
    void    IncDistinctHistAt( unsigned char res, ssize_t n )  {IncDistinctHistAt(n,PS_M,res);}

    size_t          GetDistinctHistAt( ssize_t n, int st, unsigned char res ) const;
    void            SetDistinctHistAt( ssize_t n, int st, unsigned char res, size_t value );
    void            IncDistinctHistAt( ssize_t n, int st, unsigned char res );

    size_t          GetMSDistinctHistAt( unsigned char res, size_t n ) const;
    void            SetMSDistinctHistAt( size_t value, unsigned char res, size_t n );
    void            IncMSDistinctHistAt( unsigned char res, size_t n );


    double          GetTargetFreqnsAt( unsigned char res, size_t n ) const;
    void            SetTargetFreqnsAt( double value, unsigned char res, size_t n );
    void            NullTargetFreqnsAt(size_t n );

    double          GetTargetTranstBeg( int trans ) const;
    const double ( *GetTargetTranstBeg() const )[P_NSTATES];
    void            SetTargetTranstBeg( double value, int trans );
    void            SetTargetTranstBeg( const double(*)[P_NSTATES]);
    void            NullTargetTranstBeg();

    double          GetTargetTranstAt( int trans, ssize_t n ) const;
    const double ( *GetTargetTranstAt( ssize_t n ) const )[P_NSTATES];
    void            SetTargetTranstAt( double value, int trans, ssize_t n );
    void            SetTargetTranstAt( const double(*)[P_NSTATES], ssize_t n );
    void            NullTargetTranstAt(ssize_t n );

    const double ( *GetPSSMVector() const )[NUMALPH] { return rawPSSMatrix; }
    double          GetPSSMEntryAt( unsigned char res, size_t n ) const;
    void            SetPSSMEntryAt( double value, unsigned char res, size_t n );
    const double ( *GetPSSMVectorAt(  size_t n ) const )[NUMALPH];

    double          GetInformationAt( size_t n ) const;
    void            SetInformationAt( double value, size_t n );

    
    double          GetBckPPProbAt( size_t n ) const;
    void            SetBckPPProbAt( double value, size_t n );

    const double*   GetPPProbsAt( size_t n ) const;
    const int*      GetPPPIndsAt( size_t n ) const;
    void            SetPPProbsAt( size_t n, const double* probs, const int* ndxs, int size );

    size_t          GetNoPPProbsAt( size_t n ) const;


    size_t          GetNoSequencesAt( size_t n ) const;
    void            SetNoSequencesAt( size_t value, size_t n );

    double          GetExpNoObservationsAt( size_t n ) const;
    void            SetExpNoObservationsAt( double value, size_t n );


    double          GetMIDExpNoObservationsAt( int n, int st ) const;
    void            SetMIDExpNoObservationsAt( double value, int n, int st );
    void            IncMIDExpNoObservationsAt( double value, int n, int st );
    void            MulMIDExpNoObservationsAt( double value, int n, int st );


    void            SetSequenceReservation( unsigned amnt ) { reservation_ = amnt; }

    virtual void    clear();//clear all data

    static int      GetSizeOfExpNoDistinctRes() { return SIZEOFEXPNOBSERVS; };
    static void     InitPrexpNoDistinctRes( const double* = NULL );
    static void     DestroyPrexpNoDistinctRes();
    static double   GetExpNoObservations( double avgnodistres );

protected:
    explicit        ExtendedDescriptionVector();
    virtual void    Realloc( int newcap );
    void            ReallocIndices( size_t newcap, size_t n );
    virtual void    Init();

    void            InitRightExtents( size_t from = 0, size_t to = 0 );

private:
    size_t**            indices;                    //indices of the sequences outside the extent
    size_t*             indicesLengths_;            //lengths of indices for each position
    size_t*             indicesCapacities_;         //capacities of indices for each position
    size_t              reservation_;
    //
    int                *states;                     //encodes states for each position
    double              backprobs_[NUMALPH];        //background probabilities for this description vector
    double              postprobs_[NUMALPH];        //posterior probabilities for this description vector
    size_t           ( *extents )[PS_NSTATES][xCount];//left and right boundaries of the extents
    double           ( *nodifrs_ )[PS_NSTATES];     //no. different residues of extents
    size_t           ( *MSextents )[xCount];        //match-state version of extents
    size_t*             counts;                     //number of matched residues (gap inc.) at each position
    size_t           ( *distribution )[NUMALPH];    //residue distribution at each of the positions
    double**            sqnweights_[TIMES2(PS_NSTATES)];//MID sequence weights with allocation mark vectors
    double           ( *matchWeights )[NUMALPH];    //residue match weights
    double           ( *transWeights )[P_NSTATES];  //transition weights
    double             *gapWeights;                 //gap weights
    size_t           ( *distinctHist )[PS_NSTATES][NUMALPH];//histogram of distinct residues for each position
    size_t           ( *MSdistinctHist )[NUMALPH];  //histogram of distinct residues for each match-state position
    double           ( *targetFreqns )[NUMALPH];    //estimated frequencies
    double           ( *targetTranst )[P_NSTATES];  //estimated transition  frequencies
    double           ( *rawPSSMatrix )[NUMALPH];    //not scaled PSSM matrix
    double*             bppprob_;                   //background posterior predictive
    double**            ppprobs_;                   //posterior predictive probabilities for each cluster
    int**               pppndxs_;                   //indices of clusters p.p.probabilities were calculated for
    size_t*             noppps_;                    //number of p.p.probability values
    double*             information;                //information content for each position
    size_t*             noseqs_;                    //number of sequences at each position
    double*             expnoseqs_;                 //expected number of observations (sequences) at each position
    double           ( *expnosmid_ )[PS_NSTATES];   //expected no. observations for each state

    static double*      prexpnores_;                //precomputed expected numbers of distinct residues
};



////////////////////////////////////////////////////////////////////////////
// Class PosDescriptionVector inlines
//
// get residue at the given position
//
inline
unsigned char PosDescriptionVector::ResidueAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !residues || length_ <= n )
        throw myruntime_error("PosDescriptionVector: Memory access error.");
#endif
    return residues[n];
}

// set residue at the position
//
inline
void PosDescriptionVector::SetResidueAt( unsigned char r, size_t n )
{
#ifdef __DEBUG__
    if( residues == NULL || length_ <= n )
        throw myruntime_error("PosDescriptionVector: Memory access error.");
#endif
    residues[n] = r;
}

// -------------------------------------------------------------------------
// get `used' flag at the given position
//
inline
bool PosDescriptionVector::IsUsedAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !flags || length_ <= n )
        throw myruntime_error("PosDescriptionVector: Memory access error.");
#endif
    return GetUsed() && flags[n] != UNUSED;
}

// Set `unused' flag at the position
//
inline
void PosDescriptionVector::SetUnusedAt( size_t n )
{
#ifdef __DEBUG__
    if( flags == NULL || length_ <= n )
        throw myruntime_error("PosDescriptionVector: Memory access error.");
#endif
    flags[n] = UNUSED;
}

// set `used' flag at the position
//
inline
void PosDescriptionVector::SetUsedAt( size_t n )
{
#ifdef __DEBUG__
    if( flags == NULL || length_ <= n )
        throw myruntime_error("PosDescriptionVector: Memory access error.");
#endif
    flags[n] = USED;
}

// check whether the sequence is used and informative at the position
//
inline
bool PosDescriptionVector::IsUsedAndInfAt( size_t p ) const
{
    if( !IsUsedAt( p ))
        return false;
    return GetFirstUsed() <= p && p <= GetLastUsed();
}

// -------------------------------------------------------------------------
// get whether the sequence is used in extent computed for the given position
//
inline
bool PosDescriptionVector::IsUsedInExtentAt ( size_t n ) const
{
#ifdef __DEBUG__
    if( !flags || length_ <= n )
        throw myruntime_error("PosDescriptionVector: Memory access error.");
#endif
    return flags[n] == USED_IN_EXTENT;
}

// set the flag of usage in the extent for the given position
//
inline
void PosDescriptionVector::SetUsedInExtentAt( size_t n )
{
#ifdef __DEBUG__
    if( flags == NULL || length_ <= n )
        throw myruntime_error( "PosDescriptionVector: Memory access error." );
#endif
    flags[n] = USED_IN_EXTENT;
}

// unset the flag of usage in the extent for the given position
//
inline
void PosDescriptionVector::UnsetUsedInExtentAt( size_t n )
{
#ifdef __DEBUG__
    if( flags == NULL || length_ <= n )
        throw myruntime_error( "PosDescriptionVector: Memory access error." );
#endif
    flags[n] &= USED;
}

// -------------------------------------------------------------------------
// get weight at the given position
//
inline
double PosDescriptionVector::WeightAt ( size_t n ) const
{
#ifdef __DEBUG__
    if( !weights || length_ <= n )
        throw myruntime_error("PosDescriptionVector: Memory access error.");
#endif
    return weights[n];
}

// set weight at the position
//
inline
void PosDescriptionVector::SetWeightAt( double w, size_t n )
{
#ifdef __DEBUG__
    if( weights == NULL || length_ <= n )
        throw myruntime_error("PosDescriptionVector: Memory access error.");
#endif
    weights[n] = w;
}





////////////////////////////////////////////////////////////////////////////
// Class ExtendedDescriptionVector inlines
//
// GetBackProbsAt: get background probability for residue res
//
inline
double ExtendedDescriptionVector::GetBackProbsAt( unsigned char res ) const
{
    if( backprobs_ == NULL || NUMALPH <= res )
        throw myruntime_error( "ExtendedDescriptionVector: Memory access error." );
    return backprobs_[res];
}
// SetBackProbsAt: set background probability for residue res
//
inline
void ExtendedDescriptionVector::SetBackProbsAt( unsigned char res, double value )
{
    if( backprobs_ == NULL || NUMALPH <= res )
        throw myruntime_error( "ExtendedDescriptionVector: Memory access error." );
    backprobs_[res] = value;
}


// -------------------------------------------------------------------------
// GetPostProbsAt: get posterior probability for residue res
//
inline
double ExtendedDescriptionVector::GetPostProbsAt( unsigned char res ) const
{
    if( postprobs_ == NULL || NUMALPH <= res )
        throw myruntime_error( "ExtendedDescriptionVector: GetPostProbsAt: Memory access error." );
    return postprobs_[res];
}
// SetPostProbsAt: set posterior probability for residue res
//
inline
void ExtendedDescriptionVector::SetPostProbsAt( unsigned char res, double value )
{
    if( postprobs_ == NULL || NUMALPH <= res )
        throw myruntime_error( "ExtendedDescriptionVector: GetPostProbsAt: Memory access error." );
    postprobs_[res] = value;
}


// -------------------------------------------------------------------------
// Return array of indices at the given position
//
inline
size_t* ExtendedDescriptionVector::GetIndicesAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !indices || length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
        return indices[n];
}
// Returns number of entries in the array of indices at the position
//
inline
size_t ExtendedDescriptionVector::GetIndicesSizeAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !indicesLengths_ || length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    return indicesLengths_[n];
}


// -------------------------------------------------------------------------
// get state at the position
inline
int ExtendedDescriptionVector::GetStateAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !states || length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    return states[n];
}
// set match state at the position
//
inline
void ExtendedDescriptionVector::SetStateAt( size_t n )
{
#ifdef __DEBUG__
    if( states == NULL || length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    states[n] = 1;
}


// -------------------------------------------------------------------------
// get left extent computed at the given position
//
inline
size_t ExtendedDescriptionVector::GetLeftExtentAt( size_t n, int st ) const
{
#ifdef __DEBUG__
    if( !extents || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    return extents[n][st][xLeft];
}
// set left extent value at the position
//
inline
void ExtendedDescriptionVector::SetLeftExtentAt( size_t n, int st, size_t value )
{
#ifdef __DEBUG__
    if( extents == NULL || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    extents[n][st][xLeft] = value;
}


// -------------------------------------------------------------------------
// get right extent computed at the given position
//
inline
size_t ExtendedDescriptionVector::GetRightExtentAt( size_t n, int st ) const
{
#ifdef __DEBUG__
    if( !extents || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
        return extents[n][st][xRight];
}
// set right extent value at the position
//
inline
void ExtendedDescriptionVector::SetRightExtentAt( size_t n, int st, size_t value )
{
#ifdef __DEBUG__
    if( extents == NULL || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    extents[n][st][xRight] = value;
}


// -------------------------------------------------------------------------
// get interval size of extent computed at the given position
//
inline
size_t ExtendedDescriptionVector::GetExtentIntervalAt( size_t n, int st ) const
{
#ifdef __DEBUG__
    if( !extents || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    return extents[n][st][xInterval];
}
// set interval size of extent computed at the position
//
inline
void ExtendedDescriptionVector::SetExtentIntervalAt( size_t n, int st, size_t value )
{
#ifdef __DEBUG__
    if( extents == NULL || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    extents[n][st][xInterval] = value;
}


// -------------------------------------------------------------------------
// get number of sequences in extent at the position
//
inline
size_t ExtendedDescriptionVector::GetNoSequencesInExtentAt( size_t n, int st ) const
{
#ifdef __DEBUG__
    if( !extents || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    return extents[n][st][xNoSequences];
}
// increment number of sequences in extent at the position
//
inline
void ExtendedDescriptionVector::IncNoSequencesInExtentAt( size_t n, int st )
{
#ifdef __DEBUG__
    if( extents == NULL || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    extents[n][st][xNoSequences]++;
}


// -------------------------------------------------------------------------
// get number of symbols present in extent at the position
//
inline
double ExtendedDescriptionVector::GetNoSymbolsInExtentAt( size_t n, int st ) const
{
#ifdef __DEBUG__
    if( !nodifrs_ || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    return nodifrs_[n][st];
}
// set number of symbols in extent at the position
//
inline
void ExtendedDescriptionVector::SetNoSymbolsInExtentAt( size_t n, int st, double value )
{
#ifdef __DEBUG__
    if( nodifrs_ == NULL || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    nodifrs_[n][st] = value;
}



// -------------------------------------------------------------------------
// get left boundary of match-state extent at the given position
//
inline
size_t ExtendedDescriptionVector::GetLeftMSExtentAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !MSextents || length_ <= n )
        throw myruntime_error( "ExtendedDescriptionVector: GetLeftMSExtentAt: Memory access error." );
#endif
    return MSextents[n][xLeft];
}
// set left boundary of match-state extent at the position
//
inline
void ExtendedDescriptionVector::SetLeftMSExtentAt( size_t value, size_t n )
{
#ifdef __DEBUG__
    if( MSextents == NULL || length_ <= n )
        throw myruntime_error( "ExtendedDescriptionVector: SetLeftMSExtentAt: Memory access error." );
#endif
    MSextents[n][xLeft] = value;
}
// -------------------------------------------------------------------------
// get right boundary of match-state extent at the given position
//
inline
size_t ExtendedDescriptionVector::GetRightMSExtentAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !MSextents || length_ <= n )
        throw myruntime_error( "ExtendedDescriptionVector: GetRightMSExtentAt: Memory access error." );
#endif
    return MSextents[n][xRight];
}
// set right boundary of match-state extent at the position
//
inline
void ExtendedDescriptionVector::SetRightMSExtentAt( size_t value, size_t n )
{
#ifdef __DEBUG__
    if( MSextents == NULL || length_ <= n )
        throw myruntime_error( "ExtendedDescriptionVector: SetRightMSExtentAt: Memory access error." );
#endif
    MSextents[n][xRight] = value;
}
// -------------------------------------------------------------------------
// get interval size of match-state extent at the given position
//
inline
size_t ExtendedDescriptionVector::GetMSExtentIntervalAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !MSextents || length_ <= n )
        throw myruntime_error( "ExtendedDescriptionVector: Memory access error." );
#endif
    return MSextents[n][xInterval];
}
// set interval size of match-state extent at the position
//
inline
void ExtendedDescriptionVector::SetMSExtentIntervalAt( size_t value, size_t n )
{
#ifdef __DEBUG__
    if( MSextents == NULL || length_ <= n )
        throw myruntime_error( "ExtendedDescriptionVector: SetMSExtentIntervalAt: Memory access error." );
#endif
    MSextents[n][xInterval] = value;
}
// -------------------------------------------------------------------------
// get number of sequences in match-state extent at the position
//
inline
size_t ExtendedDescriptionVector::GetNoSequencesInMSExtentAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !MSextents || length_ <= n )
        throw myruntime_error( "ExtendedDescriptionVector: GetNoSequencesInMSExtentAt: Memory access error." );
#endif
    return MSextents[n][xNoSequences];
}
// increment number of sequences in match-state extent at the position
//
inline
void ExtendedDescriptionVector::IncNoSequencesInMSExtentAt( size_t n )
{
#ifdef __DEBUG__
    if( MSextents == NULL || length_ <= n )
        throw myruntime_error( "ExtendedDescriptionVector: IncNoSequencesInMSExtentAt: Memory access error." );
#endif
    MSextents[n][xNoSequences]++;
}
// -------------------------------------------------------------------------
// get number of symbols in match-state extent at the position
//
inline
size_t ExtendedDescriptionVector::GetNoSymbolsInMSExtentAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !MSextents || length_ <= n )
        throw myruntime_error( "ExtendedDescriptionVector: GetNoSymbolsInMSExtentAt: Memory access error." );
#endif
    return MSextents[n][xNoSymbols];
}
// set number of symbols in match-state extent at the position
//
inline
void ExtendedDescriptionVector::SetNoSymbolsInMSExtentAt( size_t value, size_t n )
{
#ifdef __DEBUG__
    if( MSextents == NULL || length_ <= n )
        throw myruntime_error( "ExtendedDescriptionVector: SetNoSymbolsInMSExtentAt: Memory access error." );
#endif
    MSextents[n][xNoSymbols] = value;
}




// -------------------------------------------------------------------------
// Compute observed frequency weight known as alpha coefficient;
//     expresses mean value of different residues per column (of extent)
//
inline
double ExtendedDescriptionVector::ComputeObsFrequencyWeightAt( size_t n ) const
{
//     size_t  interval = GetExtentIntervalAt( n );
//     size_t  diffsyms = GetNoSymbolsInExtentAt( n );
    size_t  interval = GetMSExtentIntervalAt( n );
    double  diffsyms = GetNoSymbolsInExtentAt( n, PS_M );
    double  weight = diffsyms;

//     if( !interval || !diffsyms )
//         return 0.0;
//     //NOTE:already processed over all interval
//     if( 1 < interval )
//         weight /= ( double )interval;
    if( weight < 0.0 )
        weight = 0.0;
    return weight;
}




// -------------------------------------------------------------------------
// get number of residues evaluated at the position
//
inline
size_t ExtendedDescriptionVector::GetCountAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !counts || length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    return counts[n];
}
// set number of residues observed at the position
//
inline
void ExtendedDescriptionVector::IncCountAt( size_t n )
{
#ifdef __DEBUG__
    if( counts == NULL || length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    counts[n]++;
}


// -------------------------------------------------------------------------
// get observed frequency of residue `res' at the given position
//
inline
size_t ExtendedDescriptionVector::GetDistributionAt( unsigned char res, size_t n ) const
{
#ifdef __DEBUG__
    if( !distribution || length_ <= n || NUMALPH <= res )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    return distribution[n][res];
}
// set observed frequency value of residue `res' at the given position
//
inline
void ExtendedDescriptionVector::IncDistributionAt( unsigned char res, size_t n )
{
#ifdef __DEBUG__
    if( distribution == NULL || length_ <= n || NUMALPH <= res )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    distribution[n][res]++;
}





// -------------------------------------------------------------------------
// GetSqnWeightsAt: get weights of sequences for position `n' and state `st'
//
inline
const double* ExtendedDescriptionVector::GetSqnWeightsAt( size_t n, int st ) const
{
#ifdef __DEBUG__
    if( sqnweights_ == NULL || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw myruntime_error( "ExtendedDescriptionVector: GetSqnWeightsAt: Memory access error." );
#endif
    return sqnweights_[st][n];
}
// GetSqnWeightsAt: get weights of sequences for position `n' and state `st'
//
inline
double* ExtendedDescriptionVector::GetSqnWeightsAt( size_t n, int st )
{
#ifdef __DEBUG__
    if( sqnweights_ == NULL || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw myruntime_error( "ExtendedDescriptionVector: GetSqnWeightsAt: Memory access error." );
#endif
    return sqnweights_[st][n];
}
// SetSqnWeightsAt: set weights of sequences for position `n' and state `st'
//  equal to those for position `pn'
//
inline
void ExtendedDescriptionVector::SetSqnWeightsAt( size_t n, int st, size_t pn )
{
#ifdef __DEBUG__
    if( sqnweights_ == NULL || length_ <= n || n <= pn || st < 0 || PS_NSTATES <= st )
        throw myruntime_error( "ExtendedDescriptionVector: SetSqnWeightsAt: Memory access error." );
#endif
    sqnweights_[st][n] = sqnweights_[st][pn];
}
// SetSqnWeightsAt: allocate space for weights of sequences for 
//  position `n' and state `st'
//
inline
void ExtendedDescriptionVector::NewSqnWeightsAt( size_t n, int st, size_t size )
{
#ifdef __DEBUG__
    if( sqnweights_ == NULL || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw myruntime_error( "ExtendedDescriptionVector: NewSqnWeightsAt: Memory access error." );
#endif
    if( size < 1 )
        return;
    double* wghts = NULL;
    //check if not occupied
    FreeSqnWeightsAt( n, st );
    wghts = sqnweights_[st][n] = ( double* )malloc( sizeof( double ) * size );
    if( wghts == NULL )
        throw myruntime_error( "ExtendedDescriptionVector: NewSqnWeightsAt: Not enough memory." );
    //initialize memory and mark the allocated block as belonging to this position;
    //the latter is for deallocation
    memset( wghts, 0, sizeof( double ) * size );
    sqnweights_[st+PS_NSTATES][n] = wghts;
}
// FreeSqnWeightsAt: deallocate memory for weights of sequences for 
//  position `n' and state `st'
//
inline
void ExtendedDescriptionVector::FreeSqnWeightsAt( size_t n, int st )
{
#ifdef __DEBUG__
    if( sqnweights_ == NULL || length_ <= n || st < 0 || PS_NSTATES <= st )
        throw myruntime_error( "ExtendedDescriptionVector: NewSqnWeightsAt: Memory access error." );
#endif
    double* marks = sqnweights_[st+PS_NSTATES][n];
    if( marks == NULL )
        return;
    free( marks );
    sqnweights_[st+PS_NSTATES][n] = NULL;
    sqnweights_[st][n] = NULL;
}





// -------------------------------------------------------------------------
// get match weight computed for residue of type `res' at the given position
//
inline
double ExtendedDescriptionVector::GetMatchWeightsAt( unsigned char res, size_t n ) const
{
#ifdef __DEBUG__
    if( !matchWeights || length_ <= n || NUMALPH <= res )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    return matchWeights[n][res];
}
// get match weight vector at the given position
//
inline
const double ( *ExtendedDescriptionVector::GetMatchWeightsAt( size_t n ) const )[NUMALPH]
{
#ifdef __DEBUG__
    if( !matchWeights || length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    return matchWeights + n;
}
// set weight value for residue of type `res' at the given position
//
inline
void ExtendedDescriptionVector::SetMatchWeightsAt( double value, unsigned char res, size_t n )
{
#ifdef __DEBUG__
    if( matchWeights == NULL || length_ <= n || NUMALPH <= res )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    matchWeights[n][res] = value;
}
// increment weight value for residue of type `res' at the given position
//
inline
void ExtendedDescriptionVector::IncMatchWeightsAt( double value, unsigned char res, size_t n )
{
#ifdef __DEBUG__
    if( matchWeights == NULL || length_ <= n || NUMALPH <= res )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    matchWeights[n][res] += value;
}





// -------------------------------------------------------------------------
// get beginning transition weight: 0th pos.
//
inline
double ExtendedDescriptionVector::GetTransWeightsBeg( int trans ) const
{
#ifdef __DEBUG__
    if( !transWeights || !length_ || P_NSTATES <= trans )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    return transWeights[0][trans];
}
// get beginning transition weight vector: 0th pos.
//
inline
const double ( *ExtendedDescriptionVector::GetTransWeightsBeg() const )[P_NSTATES]
{
#ifdef __DEBUG__
    if( !transWeights || !length_ )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    return transWeights;
}
// set beginning transition weight: 0th pos.
//
inline
void ExtendedDescriptionVector::SetTransWeightsBeg( double value, int trans )
{
#ifdef __DEBUG__
    if( transWeights == NULL || !length_ || P_NSTATES <= trans )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    transWeights[0][trans] = value;
}
// increment beginning transition weight: 0th pos.
//
inline
void ExtendedDescriptionVector::IncTransWeightsBeg( double value, int trans )
{
#ifdef __DEBUG__
    if( transWeights == NULL || !length_ || P_NSTATES <= trans )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    transWeights[0][trans] += value;
}


// -------------------------------------------------------------------------
// get transition weight computed at the position; 0th pos. is reserved
//
inline
double ExtendedDescriptionVector::GetTransWeightsAt( int trans, ssize_t n ) const
{
#ifdef __DEBUG__
    if( transWeights == NULL || length_+1 <= n+1 || P_NSTATES <= trans )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    if( n < 0 )
        return GetTransWeightsBeg( trans );
    return transWeights[n+1][trans];
}
// get transition weight vector at the position
//
inline
const double ( *ExtendedDescriptionVector::GetTransWeightsAt( ssize_t n ) const )[P_NSTATES]
{
#ifdef __DEBUG__
    if( transWeights == NULL || length_+1 <= n+1 )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    if( n < 0 )
        return GetTransWeightsBeg();
    return transWeights +( n+1 );
}
// set transition weight at the position; 0th pos. is reserved
//
inline
void ExtendedDescriptionVector::SetTransWeightsAt( double value, int trans, ssize_t n )
{
#ifdef __DEBUG__
    if( transWeights == NULL || length_+1 <= n+1 || P_NSTATES <= trans )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    if( n < 0 )
        SetTransWeightsBeg( value, trans );
    else
        transWeights[n+1][trans] = value;
}
// increment transition weight at the position; 0th pos. is reserved
//
inline
void ExtendedDescriptionVector::IncTransWeightsAt( double value, int trans, ssize_t n )
{
#ifdef __DEBUG__
    if( transWeights == NULL || length_+1 <= n+1 || P_NSTATES <= trans )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    if( n < 0 )
        IncTransWeightsBeg( value, trans );
    else
        transWeights[n+1][trans] += value;
}


// -------------------------------------------------------------------------
// get gap weight computed at the given position
//
inline
double ExtendedDescriptionVector::GetGapWeightsAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !gapWeights || length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    return gapWeights[n];
}
// set gap weight value at the given position
//
inline
void ExtendedDescriptionVector::SetGapWeightsAt( double value, size_t n )
{
#ifdef __DEBUG__
    if( gapWeights == NULL || length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    gapWeights[n] = value;
}



// -------------------------------------------------------------------------
// GetDistinctHistAt: get number of occurences for residue `res'
//
inline
size_t ExtendedDescriptionVector::GetDistinctHistAt( ssize_t n, int st, unsigned char res ) const
{
#ifdef __DEBUG__
    if( !distinctHist || n < -1 || length_+1 <= n+1 || st < 0 || PS_NSTATES <= st || NUMALPH <= res )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    return distinctHist[n+1][st][res];
}
// SetDistinctHistAt: set number of occurrences for residue `res'
//
inline
void ExtendedDescriptionVector::SetDistinctHistAt( ssize_t n, int st, unsigned char res, size_t value )
{
#ifdef __DEBUG__
    if( distinctHist == NULL || n < -1 || length_+1 <= n+1 || st < 0 || PS_NSTATES <= st || NUMALPH <= res )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    distinctHist[n+1][st][res] = value;
}
// IncDistinctHistAt: increase number of occurrences for residue `res'
//
inline
void ExtendedDescriptionVector::IncDistinctHistAt( ssize_t n, int st, unsigned char res )
{
#ifdef __DEBUG__
    if( distinctHist == NULL || n < -1 || length_+1 <= n+1 || st < 0 || PS_NSTATES <= st || NUMALPH <= res )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    distinctHist[n+1][st][res]++;
}


// -------------------------------------------------------------------------
// GetMSDistinctHistAt: get number of occurences for residue `res' at 
//  match-state position
//
inline
size_t ExtendedDescriptionVector::GetMSDistinctHistAt( unsigned char res, size_t n ) const
{
#ifdef __DEBUG__
    if( !MSdistinctHist || length_ <= n || NUMALPH <= res )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    return MSdistinctHist[n][res];
}
// SetMSDistinctHistAt: set number of occurrences for residue `res' at 
//  match-state position
//
inline
void ExtendedDescriptionVector::SetMSDistinctHistAt( size_t value, unsigned char res, size_t n )
{
#ifdef __DEBUG__
    if( MSdistinctHist == NULL || length_ <= n || NUMALPH <= res )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    MSdistinctHist[n][res] = value;
}
// IncMSDistinctHistAt: increase number of occurrences for residue `res' at 
//  match-state position
//
inline
void ExtendedDescriptionVector::IncMSDistinctHistAt( unsigned char res, size_t n )
{
#ifdef __DEBUG__
    if( MSdistinctHist == NULL || length_ <= n || NUMALPH <= res )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    MSdistinctHist[n][res]++;
}


// -------------------------------------------------------------------------
// get estimated probability for residue of type `res' at the given position
//
inline
double ExtendedDescriptionVector::GetTargetFreqnsAt( unsigned char res, size_t n ) const
{
#ifdef __DEBUG__
    if( !targetFreqns || length_ <= n && NUMALPH <= res )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    return targetFreqns[n][res];
}
// set estimated probability for residue of type `res' at the given position
//
inline
void ExtendedDescriptionVector::SetTargetFreqnsAt( double value, unsigned char res, size_t n )
{
#ifdef __DEBUG__
    if( targetFreqns == NULL || length_ <= n || NUMALPH <= res )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    targetFreqns[n][res] = value;
}
// make all target frequencies at the given position equal to zero
//
inline
void ExtendedDescriptionVector::NullTargetFreqnsAt( size_t n )
{
#ifdef __DEBUG__
    if( targetFreqns == NULL || length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    memset( targetFreqns + n, 0, sizeof( double ) * NUMALPH );
}


// -------------------------------------------------------------------------
// get estimated beginning transition probability (at 0th pos.)
//
inline
double ExtendedDescriptionVector::GetTargetTranstBeg( int trans ) const
{
#ifdef __DEBUG__
    if( !targetTranst || !length_ || P_NSTATES <= trans )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    return targetTranst[0][trans];
}
// get estimated beginning transition probabilites: 0th pos.
//
inline
const double ( *ExtendedDescriptionVector::GetTargetTranstBeg() const )[P_NSTATES]
{
#ifdef __DEBUG__
    if( !targetTranst || !length_ )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    return targetTranst;
}
// -------------------------------------------------------------------------
// set estimated beginning transition probability: 0th pos.
//
inline
void ExtendedDescriptionVector::SetTargetTranstBeg( double value, int trans )
{
#ifdef __DEBUG__
    if( targetTranst == NULL || !length_ || P_NSTATES <= trans )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    targetTranst[0][trans] = value;
}
// set estimated beginning transition probabilities: 0th pos.
//
inline
void ExtendedDescriptionVector::SetTargetTranstBeg( const double ( *values )[P_NSTATES])
{
#ifdef __DEBUG__
    if( targetTranst == NULL || values == NULL || *values == NULL || !length_ )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    memcpy( targetTranst[0], *values, sizeof( double ) * P_NSTATES );
}
// make all beginning target transition frequencies equal to zero: 0th pos.
//
inline
void ExtendedDescriptionVector::NullTargetTranstBeg()
{
#ifdef __DEBUG__
    if( targetTranst == NULL || !length_ )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    memset( targetTranst, 0, sizeof( double ) * P_NSTATES );
}


// -------------------------------------------------------------------------
// get estimated transition probability at the given position;
//  0th pos. is reserved
//
inline
double ExtendedDescriptionVector::GetTargetTranstAt( int trans, ssize_t n ) const
{
#ifdef __DEBUG__
    if( targetTranst == NULL || length_+1 <= n+1 || P_NSTATES <= trans )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    if( n < 0 )
        return GetTargetTranstBeg( trans );
    return targetTranst[n+1][trans];
}
// get estimated transition probabilities at the given position;
//  0th pos. is reserved
//
inline
const double ( *ExtendedDescriptionVector::GetTargetTranstAt( ssize_t n ) const )[P_NSTATES]
{
#ifdef __DEBUG__
    if( targetTranst == NULL || length_+1 <= n+1 )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    if( n < 0 )
        return GetTargetTranstBeg();
    return targetTranst + ( n+1 );
}
// -------------------------------------------------------------------------
// set estimated transition probability at the position; 0th pos. is reserved
//
inline
void ExtendedDescriptionVector::SetTargetTranstAt( double value, int trans, ssize_t n )
{
#ifdef __DEBUG__
    if( targetTranst == NULL || length_+1 <= n+1 || P_NSTATES <= trans )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    if( n < 0 )
        SetTargetTranstBeg( value, trans );
    else
        targetTranst[n+1][trans] = value;
}
// set estimated transition probabilities at the position; 0th pos. is reserved
//
inline
void ExtendedDescriptionVector::SetTargetTranstAt( const double ( *values )[P_NSTATES], ssize_t n )
{
#ifdef __DEBUG__
    if( targetTranst == NULL || values == NULL || *values == NULL || length_+1 <= n+1 )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    if( n < 0 )
        SetTargetTranstBeg( values );
    else
        memcpy( targetTranst[n+1], *values, sizeof( double ) * P_NSTATES );
}
// make all target transition frequencies at the position equal to zero;
//  0th pos. is reserved
//
inline
void ExtendedDescriptionVector::NullTargetTranstAt( ssize_t n )
{
#ifdef __DEBUG__
    if( targetTranst == NULL || length_+1 <= n+1 )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    if( n < 0 )
        NullTargetTranstBeg();
    else
        memset( targetTranst +( n+1 ), 0, sizeof( double ) * P_NSTATES );
}


// -------------------------------------------------------------------------
// get PSSM matrix entry given residue type and position
//
inline
double ExtendedDescriptionVector::GetPSSMEntryAt( unsigned char res, size_t n ) const
{
#ifdef __DEBUG__
    if( !rawPSSMatrix || length_ <= n || NUMALPH <= res )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    return rawPSSMatrix[n][res];
}
// get vecor of values from the PSSM matrix given position
//
inline
const double ( *ExtendedDescriptionVector::GetPSSMVectorAt( size_t n ) const )[NUMALPH]
{
#ifdef __DEBUG__
    if( !rawPSSMatrix || length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    return rawPSSMatrix + n;
}
// set PSSM matrix entry given residue type and position
//
inline
void ExtendedDescriptionVector::SetPSSMEntryAt( double value, unsigned char res, size_t n )
{
#ifdef __DEBUG__
    if( rawPSSMatrix == NULL || length_ <= n || NUMALPH <= res )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    rawPSSMatrix[n][res] = value;
}


// -------------------------------------------------------------------------
// get information content calculated for the position
//
inline
double ExtendedDescriptionVector::GetInformationAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !information || length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    return information[n];
}
// set information content at the position
//
inline
void ExtendedDescriptionVector::SetInformationAt( double value, size_t n )
{
#ifdef __DEBUG__
    if( information == NULL || length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    information[n] = value;
}

// -------------------------------------------------------------------------
// GetBckPPProbAt: get background posterior predictive probability at the 
//  position
//
inline
double ExtendedDescriptionVector::GetBckPPProbAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !bppprob_ || length_ <= n )
        throw myruntime_error("ExtendedDescriptionVector: GetBckPPProbAt: Memory access error.");
#endif
    return bppprob_[n];
}

// SetBckPPProbAt: set background posterior predictive probability at the 
//  position
inline
void ExtendedDescriptionVector::SetBckPPProbAt( double value, size_t n )
{
#ifdef __DEBUG__
    if( !bppprob_ || length_ <= n )
        throw myruntime_error("ExtendedDescriptionVector: GetBckPPProbAt: Memory access error.");
#endif
    bppprob_[n] = value;
}

// GetPPProbsAt: get posterior predictive probabilities at the position
inline
const double* ExtendedDescriptionVector::GetPPProbsAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !ppprobs_ || length_ <= n )
        throw myruntime_error("ExtendedDescriptionVector: GetPPProbsAt: Memory access error.");
#endif
    return ppprobs_[n];
}

// GetPPPIndsAt: get indices of posterior predictive probabilities at the 
//  position
inline
const int* ExtendedDescriptionVector::GetPPPIndsAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !pppndxs_ || length_ <= n )
        throw myruntime_error("ExtendedDescriptionVector: GetPPPIndsAt: Memory access error.");
#endif
    return pppndxs_[n];
}

// SetPPProbsAt: set posterior predictive probabilities and their 
//  indices at the position
inline
void ExtendedDescriptionVector::SetPPProbsAt( size_t n, const double* probs, const int* ndxs, int size )
{
    int k;
    if( !ppprobs_ || !pppndxs_ || !noppps_ || length_ <= n )
        throw myruntime_error("ExtendedDescriptionVector: SetPPProbsAt: Memory access error.");
    if( size < 0 )
        throw myruntime_error("ExtendedDescriptionVector: SetPPProbsAt: Invalid size of data.");
    if( ppprobs_[n]) { free( ppprobs_[n]); ppprobs_[n] = NULL; }
    if( pppndxs_[n]) { free( pppndxs_[n]); pppndxs_[n] = NULL; }
    noppps_[n] = (size_t)size;
    if( size < 1 )
        return;
    ppprobs_[n] = ( double* )malloc( size * sizeof( double ));
    pppndxs_[n] = ( int* )malloc( size * sizeof( int ));
    if( ppprobs_[n] == NULL || pppndxs_[n] == NULL )
        throw myruntime_error("ExtendedDescriptionVector: SetPPProbsAt: Not enough memory.");
    if( probs == NULL || ndxs == NULL )
        throw myruntime_error("ExtendedDescriptionVector: SetPPProbsAt: Null parameters.");
    for( k = 0; k < size; k ++ ) {
        ppprobs_[n][k] = probs[k];
        pppndxs_[n][k] = ndxs[k];
    }
}

// GetNoPPProbsAt: get number of posterior predictive probabilities at the 
//  position
inline
size_t ExtendedDescriptionVector::GetNoPPProbsAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !noppps_ || length_ <= n )
        throw myruntime_error("ExtendedDescriptionVector: GetNoPPProbsAt: Memory access error.");
#endif
    return noppps_[n];
}

// -------------------------------------------------------------------------
// get number of sequences at the position
//
inline
size_t ExtendedDescriptionVector::GetNoSequencesAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !noseqs_ || length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    return noseqs_[n];
}
// set number of sequences at the position
//
inline
void ExtendedDescriptionVector::SetNoSequencesAt( size_t value, size_t n )
{
#ifdef __DEBUG__
    if( noseqs_ == NULL || length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    noseqs_[n] = value;
}


// -------------------------------------------------------------------------
// get expected number of observations at the position
//
inline
double ExtendedDescriptionVector::GetExpNoObservationsAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !expnoseqs_ || length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    return expnoseqs_[n];
}
// set expected number of observations at the position
//
inline
void ExtendedDescriptionVector::SetExpNoObservationsAt( double value, size_t n )
{
#ifdef __DEBUG__
    if( expnoseqs_ == NULL || length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    expnoseqs_[n] = value;
}


// -------------------------------------------------------------------------
// Return expected number of observations for state `st' at the position
//
inline
double ExtendedDescriptionVector::GetMIDExpNoObservationsAt( int n, int st ) const
{
#ifdef __DEBUG__
    if( expnosmid_ == NULL || n < -1 || length_+1 <= n+1 || st < 0 || PS_NSTATES <= st )
        throw myruntime_error( mystring( "GetMIDExpNoObservationsAt: Memory access error." ));
#endif
    //data of beginning position is at 0th place
    return expnosmid_[n+1][st];
}
// set expected number of observations for state `st' at the position
//
inline
void ExtendedDescriptionVector::SetMIDExpNoObservationsAt( double value, int n, int st )
{
#ifdef __DEBUG__
    if( expnosmid_ == NULL || n < -1 || length_+1 <= n+1 || st < 0 || PS_NSTATES <= st )
        throw myruntime_error( mystring( "SetMIDExpNoObservationsAt: Memory access error." ));
#endif
    //save beginning position at 0th place
    expnosmid_[n+1][st] = value;
}
// add expected number of observations for state `st' at the position
//
inline
void ExtendedDescriptionVector::IncMIDExpNoObservationsAt( double value, int n, int st )
{
#ifdef __DEBUG__
    if( expnosmid_ == NULL || n < -1 || length_+1 <= n+1 || st < 0 || PS_NSTATES <= st )
        throw myruntime_error( mystring( "IncMIDExpNoObservationsAt: Memory access error." ));
#endif
    //save beginning position at 0th place
    expnosmid_[n+1][st] += value;
}
// multiply expected number of observations for state `st' at position `n' by `value'
//
inline
void ExtendedDescriptionVector::MulMIDExpNoObservationsAt( double value, int n, int st )
{
#ifdef __DEBUG__
    if( expnosmid_ == NULL || n < -1 || length_+1 <= n+1 || st < 0 || PS_NSTATES <= st )
        throw myruntime_error( mystring( "MulMIDExpNoObservationsAt: Memory access error." ));
#endif
    //save beginning position at 0th place
    expnosmid_[n+1][st] *= value;
}


#endif//__DescriptionVector__
