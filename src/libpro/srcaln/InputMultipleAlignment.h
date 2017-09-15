/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __InputMultipleAlignment__
#define __InputMultipleAlignment__

#include "compdef.h"
#include "debug.h"
#include "rc.h"
#include "data.h"
#include "mystring.h"
#include "myexcept.h"
#include "libpro/srcpro/DescriptionVector.h"
#include "libpro/srcpro/DistributionMatrix.h"
#include "libHDP/HDPbase.h"

//{{functional interface
void MedianValue( size_t* histdata, size_t szdat, double* median, size_t scale = 1 );
//}}

class AlignmentSimulation;
class RndMSAGenerator;

// _________________________________________________________________________
// CLASS InputMultipleAlignment
//
class InputMultipleAlignment
{
public:
    // typedefs ...
    //
    InputMultipleAlignment( unsigned = ALLOCSEQ );
    ~InputMultipleAlignment();

    void                        PlainPreprocess();
    void                        SelectSequences();
    void                        ConstructProfile();

    size_t                      size() const { return length_; };
    size_t                      capacity() const { return capacity_; }

    double                      GetIdentityLevel() const { return identity_level; }
    void                        SetIdentityLevel( double value ) { identity_level = value; };

    const char*                 GetName() const     { return name; }
    const char*                 GetTitle() const    { return titletext; }

    void                        SetName( const char* );             //sets name
    void                        SetTitle( const char* );            //sets title information
                                                                    //appends title information to the existing one
    void                        AppendTitle( const char*, size_t = 0, size_t = UINT_MAX );

    void                        ReadAlignment( const char* );
    void                        ReadFASTA( const char* filename );
    void                        ReadSTOCKHOLM1( const char* filename );
    void                        PutAlignment( const char* = NULL );

    bool                        GetKeepTitles() const   { return keeptitles_; }
    void                        SetKeepTitles( bool value ) { keeptitles_ = value; }

    void                        SetIgnoreGapsInQuery( bool value ) { ignoregapsinquery_ = value; }
    bool                        GetIgnoreGapsInQuery() const   { return ignoregapsinquery_; }

    void                        SetComputeDELETEstates( bool value );
    bool                        GetComputeDELETEstates() const { return deletestateson_; }


    int                         SetTarFrMix( int value ) { tfrmix_ = value; }
    bool                        GetTarFrMixHDPCtx() const { return tfrmix_ == tfrmixHDPCtx; }

    int                         SetScoAdjment( int value ) { scoadj_ = value; }
    bool                        GetScoAdjmentHDPCtx() const {
                                    return scoadj_ == scoadjHDPCtx || scoadj_ == scoadjHDPsco; 
                                }
    void                        SetHDPbase( const HDPbase* value ) { HDPbase_ = value; }
    const HDPbase*              GetHDPbase() const { return HDPbase_; }
    void                        SetHDPctBase( const HDPbase* value ) { HDPctbase_ = value; }
    const HDPbase*              GetHDPctBase() const { return HDPctbase_; }


    bool                        GetUsingSEGFilter() const           { return usingsegfilt_; }
    void                        SetUsingSEGFilter( bool value )     { usingsegfilt_ = value; }

    size_t                      GetSEGWindow() const                { return segfiltwinlenval_; }
    void                        SetSEGWindow( size_t value )        { segfiltwinlenval_ = value; SetUsingSEGFilter( true ); }

    double                      GetSEGLowEntropy() const            { return segfiltlowentval_; }
    void                        SetSEGLowEntropy( double value )    { segfiltlowentval_ = value; SetUsingSEGFilter( true ); }

    double                      GetSEGHighEntropy() const           { return segfilthighentval_; }
    void                        SetSEGHighEntropy( double value )   { segfilthighentval_ = value;SetUsingSEGFilter( true ); }


    bool                        GetUsingSeqSEGFilter() const        { return usingseqseg_; }
    void                        SetUsingSeqSEGFilter( bool value )  { usingseqseg_ = value; }

    size_t                      GetSeqSEGWindow() const             { return seqsegwinlenval_; }
    void                        SetSeqSEGWindow( size_t value )     { seqsegwinlenval_ = value; SetUsingSeqSEGFilter( true ); }

    double                      GetSeqSEGLowEntropy() const         { return seqseglowentval_; }
    void                        SetSeqSEGLowEntropy( double value ) { seqseglowentval_ = value; SetUsingSeqSEGFilter( true ); }

    double                      GetSeqSEGHighEntropy() const        { return seqseghighentval_; }
    void                        SetSeqSEGHighEntropy( double value ){ seqseghighentval_ = value;SetUsingSeqSEGFilter( true ); }


    size_t                      GetExtentMinWindow() const                  { return extminwindow_; }
    void                        SetExtentMinWindow( size_t value )          { extminwindow_ = value; }

    double                      GetExtentMinSeqPercentage() const           { return extminseqperc_; }
    void                        SetExtentMinSeqPercentage( double value )   { extminseqperc_ = value; }

    double                      GetPseudoCountWeight() const            { return pseudocntweight_; }
    void                        SetPseudoCountWeight( double value )    { pseudocntweight_ = value; }

#if 0
    std::ifstream&              operator<<( std::ifstream& );
    std::ofstream&              operator>>( std::ofstream& );
#endif

    void                        ExportFrequenciesTo( FrequencyMatrix& ) const;
    void                        ExportPSSMTo( LogOddsMatrix&, bool = true ) const;
    void                        ExportGapWeightsTo( GapScheme& ) const;

    void                        OutputWeightedFrequencies( const char* = NULL );
    void                        OutputPSSM( const char* = NULL );
    void                        OutputSuppressedProfile( const char* = NULL );
    void                        OutputProfile( const char* = NULL );

    double                      GetEffNoSequences() const { return effnoseqs_; }

    PosDescriptionVector*       SequenceAt( size_t n ) const;

protected:
    void                        InitQueryDescription();
    void                        PreprocessAlignment();
    void                        PurgeAtSequenceIdentity();
    void                        PurgeAtSequenceIdentityObs();

    void                        RefineWithSEG();
    void                        FilterSequencesWithSEG();

    void                        SetPosNoSequences();

    void                        SetStates();
    void                        SetBackgroundProbabilities();
    void                        RecalcBackgroundProbabilities();
    void                        CalcPosteriorProbabilities();

    void                        ComputeExtents();

    void                        ComputeSequenceWeights();
    void                        ComputeMstateSequenceWeights();
    void                        ComputeIstateSequenceWeights();
    void                        ComputeDstateSequenceWeights();

    void                        CalculateExpNoResiduesAt( 
                                      size_t p, bool usdgwght, size_t left, size_t right, 
                                      const double* wghts, double* express, size_t scale = 10 );
    void                        ComputeGlobSequenceWeights( double* );
    void                        ComputeMIDstateSequenceWeights();
    void                        ComputeMIDstateSequenceWeightsNoExtents();
    void                        ComputeGWMIDstateSequenceWeights( double );

    void                        ComputeTransitionFrequencies( bool usegwgts, bool expmids = false );

    void                        AdjustWeights();
    void                        AdjustWeightsAt( size_t, double (*)[NUMALPH] = NULL );

    double                      ExpectedNoObservationsAt( size_t pos, int st ) const;
    void                        DeriveExpectedNoObservations();
    void                        CalculateEffNoSequences();

    void                        ComputeTargetTransFrequencies();

    void                        ComputeTargetFrequenciesMDL();
    void                        ComputeTargetFrequenciesMDLVar();
    void                        ComputeTargetFrequencies();
    void                        MixTargetFrequenciesHDPCtx();
    void                        CalcTFPosteriorPredictives();
    void                        ComputePSSM();

    void            Reset() { counter_ = 0; }
    bool            Eof() { return length_ <= counter_; }
    void            Inc() { counter_++; }

    bool                        NotNull() const { return SequenceAt( counter_ ) != NULL; }
//     bool            NotNull( size_t n ) const { return SequenceAt( n ) != NULL; }
    PosDescriptionVector*       Sequence() const { return SequenceAt( counter_ ); }

    void                        DestroyAt( size_t n );

    void                        SetSequenceAt( size_t, PosDescriptionVector* );
    void                        push( PosDescriptionVector* );
    void                        clear();                            //clears all sequences in the alignment matrix
    void                        Realloc( int newcap );

    size_t                      GetNoSequences() const              { return noseqs_; }
    void                        SetNoSequences( size_t val )        { noseqs_ = val;  }

    void                        SetEffNoSequences( double val )     { effnoseqs_ = val;  }

    PosDescriptionVector*       NewPositionVector( size_t = ALLOCPOS ) const;
    void                        DeletePositionVector( PosDescriptionVector* );

    ExtendedDescriptionVector*  NewExtendedDescriptionVector( size_t = ALLOCPOS ) const;
    ExtendedDescriptionVector*  NewExtendedDescriptionVector( const PosDescriptionVector& ) const;
    void                        DeleteExtendedDescriptionVector( ExtendedDescriptionVector* );

    void TranslateSTOConsSequence( const mystring&, PosDescriptionVector* svs );
    void TranslateSTOStates( const mystring&, PosDescriptionVector* svs );

private:
    PosDescriptionVector**      alignmentMatrix;        //input alignmemt information
    ExtendedDescriptionVector*  queryDescription;       //used to contain results of computations for query profile

    size_t                  length_;                //number of sequences contained in the alignment matrix
    size_t                  capacity_;              //current capacity of the alignment matrix
    double                  identity_level;         //sequence identity level; one of both sequences with mutual
                                                    //  similarity above this value will be ignored
    size_t                  noseqs_;                //number of sequences
    double                  effnoseqs_;             //effective number of sequences
    //
    size_t                  counter_;               //initial counter to iterate over all sequences
    //
    char*                   name;                   //name of the multiple alignment
    char*                   titletext;              //text description of the first sequence in the multiple alignment

    bool                    keeptitles_;            //whether to keep title descriptions of fasta sequences read
    bool                    ignoregapsinquery_;     //whether to ignore gaps in query sequence while reading multiple alignment
    bool                    deletestateson_;        //whether to compute delete states

    int                     tfrmix_;                //mixing of target frequencies
    int                     scoadj_;                //score adjustment
    const HDPbase*          HDPbase_;               //HDP base structure with data
    const HDPbase*          HDPctbase_;             //HDP base structure for environmental context data

    bool                    usingsegfilt_;          //flag of using filter
    size_t                  segfiltwinlenval_;      //window length
    double                  segfiltlowentval_;      //low entropy value
    double                  segfilthighentval_;     //high entropy value

    bool                    usingseqseg_;           //whether using SEG for sequences in multiple alignment
    size_t                  seqsegwinlenval_;       //window length of SEG for sequences
    double                  seqseglowentval_;       //low entropy value of SEG for sequences
    double                  seqseghighentval_;      //high entropy value of SEG for sequences

    size_t                  extminwindow_;          //minimum required window length in positions of extent
    double                  extminseqperc_;         //minimum required sequence percentage an extent must cover

    double                  pseudocntweight_;       //weight for pseudo count frequencies

    //FRIENDS...
    friend class AlignmentSimulation;
    friend class RndMSAGenerator;
};



////////////////////////////////////////////////////////////////////////////
// Class InputMultipleAlignment inlines
//
// -------------------------------------------------------------------------
// creates new position vector

inline
PosDescriptionVector* InputMultipleAlignment::NewPositionVector( size_t reserve ) const
{
    return new PosDescriptionVector( reserve );
}

// deallocates position vector

inline
void InputMultipleAlignment::DeletePositionVector( PosDescriptionVector* v )
{
    if( v )
        delete v;
}

// new default extended description vector

inline
ExtendedDescriptionVector* InputMultipleAlignment::NewExtendedDescriptionVector( size_t reserve ) const
{
    return new ExtendedDescriptionVector( reserve );
}

// new extended description vector

inline
ExtendedDescriptionVector* InputMultipleAlignment::NewExtendedDescriptionVector( const PosDescriptionVector& s ) const
{
    return new ExtendedDescriptionVector( s );
}

// deallocates extended description vector

inline
void InputMultipleAlignment::DeleteExtendedDescriptionVector( ExtendedDescriptionVector* v )
{
    if( v )
        delete v;
}

// -------------------------------------------------------------------------
// SetComputeDeleteStates: sets flag for computing of delete states
//
inline
void InputMultipleAlignment::SetComputeDELETEstates( bool value )
{
    if( value )
        SetIgnoreGapsInQuery( !value );
    deletestateson_ = value;
}

// -------------------------------------------------------------------------
// SequenceAt: picks up sequence by position

inline
PosDescriptionVector* InputMultipleAlignment::SequenceAt( size_t n ) const
{
#ifdef __DEBUG__
    if( n < length_ )
#endif
        return alignmentMatrix[n];

    throw myruntime_error( mystring( "InputMultipleAlignment: Memory access error." ));
}
// -------------------------------------------------------------------------
// SetSsequenceAt: set sequence at position n
//  NOTE: no memory management takes place!
//
inline
void InputMultipleAlignment::SetSequenceAt( size_t n, PosDescriptionVector* seq )
{
    if( length_ <= n )
        throw myruntime_error("InputMultipleAlignment: SetSsequenceAt: Memory access error.");
    alignmentMatrix[n] = seq;
}



#endif//__InputMultipleAlignment__
