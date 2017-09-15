/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __ProFreqWriter__
#define __ProFreqWriter__

#include <stdio.h>

#include "ext/rng.h"
#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "rc.h"
#include "mystring.h"
#include "myexcept.h"

#include "BinarySearchStructure.h"
#include "CtxtFrequencies.h"
#include "CtxtCoefficients.h"
#include "libpro/srcpro/DistributionMatrix.h"

// _________________________________________________________________________
// Class ProFreqWriter
//

class ProFreqWriter
{
    enum {
        RndNumHeader = 0xfffffff3,       //header of the formatted random numbers
        TargetHeader = 0xfffffff1        //header for the target function values
    };

public:
    //typedefs...
    typedef void ( ProFreqWriter::*PMETHOD )( const char*, int, void* );
    //
    ProFreqWriter( const char* name );
    ProFreqWriter( const char* output, const char* directory );
    ProFreqWriter( const char* output, char* arguments[], int no_args );
    ~ProFreqWriter();

    virtual unsigned long Run();

    void                SetMinEffThickness( int value ) { mineffthickn_ = value; }
    int                 GetMinEffThickness() const      { return mineffthickn_; }

    void                SetNoSamples( int value )       { no_samples_ = value; }
    int                 GetNoSamples() const            { return no_samples_; }

    void                SetContextSize( int value )     { context_size_ = value; }
    int                 GetContextSize() const          { return context_size_; }

    void                SetContextStep( int value )     { ctxt_step_ = value; }
    int                 GetContextStep() const          { return ctxt_step_; }

    void                SetMix( bool value )            { mix_ = value; }
    bool                GetMix() const                  { return mix_; }

    void                SetCentralWeight( double value ){ cwght_ = value; }
    double              GetCentralWeight() const        { return cwght_; }

    void                SetDishPerGroup( bool value )   { dishpergroup_ = value; }
    bool                GetDishPerGroup() const         { return dishpergroup_; }

    void                SetSeed( unsigned long value )  { seed_ = value; }
    unsigned long       GetSeed() const                 { return seed_; }

    const char*         GetOutputName() const           { return outname_; }

protected:
    explicit ProFreqWriter();

    void                ProcessInput( PMETHOD, void* );                     //iterate over all files
    void                CountFiles( const char*, int ser, void* args );     //count files
    void                CountPositions( const char* file, int, void* args );//count positions
    void                ProcessFile( const char* filename, int ser, void* );//process one file

    void                MixSample( CtxtFrequencies** p_sample, bool tonormal );
    void                MixSample3( CtxtFrequencies** p_sample, bool tonormal );

    void                Sample( const char* filename, FrequencyMatrix&, LogOddsMatrix&, 
                                int pos, CtxtFrequencies** p_sample );
    void                SampleHlp( const char* filename, FrequencyMatrix&, LogOddsMatrix&, 
                                int pos, CtxtFrequencies** p_sample, bool tonormal = true );
    bool                SampleErrCorrection( CtxtFrequencies* sample, int n );
    const char*         SampleToNormal( CtxtFrequencies* sample, int n );

    void                Output();
    void                PrintFrequenciesDishPerGroup( FILE* );
    void                PrintFrequenciesDish1( FILE* );

    const char*         GetDataDir() const          { return data_dir_; }
    const char* const*  GetFiles() const            { return files_; }
    const char*         GetFileAt( int ) const;
    int                 GetNoPositions() const      { return no_posits_; }
    void                SetNoPositions( int value ) { no_posits_ = value; }
    int*                GetNoPositionsAddr()        { return &no_posits_; }
    const int           GetNoArguments() const      { return no_arguments_; }

    unsigned long*      GetSeedAddr()               { return &seed_; }
    MTRng&              GetRng()                    { return rng_; }

    int                 GetSampleCounter() const    { return samplecounter_; }
    void                ResetSampleCounter()        { samplecounter_ = 0; }
    void                IncSampleCounter()          { samplecounter_++; }

    void                DestroyRndNumbers();
    void                NewRndNumbers();
    int                 GetNoAuxFields() const      { return 3; };
    int                 GetNoPreambleFields() const { return 2; };
    double              GetNextRndNumber();
    double*             GetRndNumbers()             { return rndnumbers_; }
    const double*       GetRndNumbers() const       { return rndnumbers_; }
    size_t              GetSizeOfRndNumbers() const;
    void                MakeRndNumbers();
    void                GenerateRndNumbers();
    void                VerifyRndNumbers();

    const CtxtCoefficients* GetMixCoeffs() const    { return coeffs_; }
    CtxtCoefficients*   GetMixCoeffs()              { return coeffs_; }
    void                DestroyCoeffs();
    void                NewCoeffs( size_t length, double cweight );

    const SimpleVector* GetSamples() const          { return samples_; }
    SimpleVector*       GetSamples()                { return samples_; }
    void                DestroySamples();
    void                NewSamples( size_t );

    int                 GetActNoSamples() const     { return actual_nosmpls_; }
    void                IncActNoSamples()           { actual_nosmpls_++; }
    void                ResetActNoSamples()         { actual_nosmpls_ = 0; }

    int                 GetNoProcessedFiles() const { return no_processedfiles_; }
    void                IncNoProcessedFiles()       { no_processedfiles_++; }
    void                ResetNoProcessedFiles()     { no_processedfiles_ = 0; }

    static size_t       GetDefGroupSize()           { return s_defgrsize; }
    static size_t       GetDefNoGroups()            { return s_defnogrps; }

private:
    const char*         outname_;           //output name
    const char*         data_dir_;          //directory that contains observed frequency files to be processed
    const char* const*  files_;             //names of files to be processed
    int                 no_posits_;         //number of positions
    const int           no_arguments_;      //number of arguments
    int                 no_processedfiles_; //number of processed files
    int                 mineffthickn_;      //minimum effective thickness
    int                 no_samples_;        //number of samples to produce
    int                 actual_nosmpls_;    //actual number of samples generated
    int                 context_size_;      //sample context size
    int                 ctxt_step_;         //step between adjacent context positions
    bool                mix_;               //whether to mix context positions
    double              cwght_;             //weight of central context position (when mixing)
    bool                dishpergroup_;      //arrange dishes one per group
    unsigned long       seed_;              //seed for random number generator
    CtxtCoefficients*   coeffs_;            //coefficients for mixing positions
    MTRng               rng_;               //random number generator
    int                 samplecounter_;     //sample counter for accessing of random numbers
    double*             rndnumbers_;        //random numbers
    SimpleVector*       samples_;           //all samples recorded to array
    static size_t       s_defgrsize;        //default group size
    static size_t       s_defnogrps;        //default number of groups
};


// /////////////////////////////////////////////////////////////////////////
// INLINES
//
inline
const char* ProFreqWriter::GetFileAt( int n ) const
{
#ifdef __DEBUG__
    if( no_arguments_ <= n )
        throw myruntime_error( mystring( "ProFreqWriter: Memory access error." ));
#endif
    return files_[ n ];
}

// -------------------------------------------------------------------------
// GetNextRndNumber: get next random number
//
inline
double ProFreqWriter::GetNextRndNumber()
{
    double  value;
#ifdef __DEBUG__
    if( samplecounter_ < 0 || rndnumbers_ == NULL )
        throw myruntime_error( mystring( "ProFreqWriter: Memory access error." ));
#endif
    if( GetNoSamples() <= samplecounter_ )
        return -1.0;

    value = rndnumbers_[ GetNoPreambleFields() + samplecounter_ ];
    IncSampleCounter();
    return value;
}

// -------------------------------------------------------------------------
// GetSizeOfRndNumbers: get size of random numbers in bytes
//
inline
size_t ProFreqWriter::GetSizeOfRndNumbers() const
{
    if( GetNoSamples() <= 0 )
        return 0;
    return ( GetNoSamples() + GetNoAuxFields()) * sizeof( double );
}

// -------------------------------------------------------------------------
// DestroyCoeffs: destroy coefficients
//
inline
void ProFreqWriter::DestroyCoeffs()
{
    if( coeffs_ )
        delete coeffs_;
    coeffs_ = NULL;
}

// -------------------------------------------------------------------------
// NewCoeffs: create coefficients
//
inline
void ProFreqWriter::NewCoeffs( size_t length, double cweight )
{
    DestroyCoeffs();
    coeffs_ = new CtxtCoefficients( length, cweight );
    if( coeffs_ == NULL )
        throw myruntime_error("ProFreqWriter: NewCoeffs: Not enough memory.");
}


#endif//__ProFreqWriter__
