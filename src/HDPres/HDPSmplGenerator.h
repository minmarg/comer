/***************************************************************************
 *   Copyright (C) 2012 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __HDPSmplGenerator__
#define __HDPSmplGenerator__

#include <stdio.h>
#include <stdlib.h>

#include "ext/rng.h"
#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "rc.h"
#include "mystring.h"
#include "myexcept.h"

#include "BinarySearchStructure.h"

static const char*  _sch_ext_grp_ = ".grp";
static const char*  _sch_ext_crt_ = ".crt";
static const char*  _sch_ext_mix_ = ".mix";
static const char*  _sch_ext_one_ = ".one";
static const char*  _sch_ext_two_ = ".two";
static const char*  _sch_ext_th3_ = ".th3";

// _________________________________________________________________________
// Class HDPSmplGenerator
//

class HDPSmplGenerator
{
public:
    HDPSmplGenerator( const char* patoutput );
    ~HDPSmplGenerator();

    void                Run();

    const char*         GetOutputPattern() const        { return patoutname_; }

    int                 GetDimensions() const           { return dim_; }
    void                SetDimensions( int value )      { dim_ = value; }

    int                 GetNoVars() const           { return novariates_; }
    void                SetNoVars( int value )      { novariates_ = value; }

    int                 GetContextLength() const        { return ctxtlen_; }
    void                SetContextLength( int value )   { ctxtlen_ = value; }

    int                 GetNoClusters() const           { return noclusters_; }
    void                SetNoClusters( int value )      { noclusters_ = value; }

    int                 GetNoSamplesC() const           { return nocsamples_; }
    void                SetNoSamplesC( int value )      { nocsamples_ = value; }

    bool                GetLGSTNSampling() const        { return lgstnrmsampling_; }
    void                SetLGSTNSampling( bool value )  { lgstnrmsampling_ = value; }

    bool                GetDoOverlap() const            { return overlap_; }
    void                SetDoOverlap( bool value )      { overlap_ = value; }

protected:
    explicit HDPSmplGenerator();

    void                GenerateRec( int, int, int* = NULL );
    void                GenSamplesForCluster( int* varndxs, int vsize );
    void                GenSamplesForClusterMeanVar( int* varndxs, int vsize );
    void                GenSamplesForClusterLGSTNRM( int* varndxs, int vsize );

    void                Output();
    void                PrintCorrect( FILE* );
    void                PrintMixed( FILE* );
    void                PrintOne( FILE* );
    void                PrintAny( FILE*, int noclsts );

    int                 GetSampleCounter() const    { return smplcounter_; }
    void                ResetSampleCounter()        { smplcounter_ = 0; }
    void                IncSampleCounter()          { smplcounter_++; }

    const SimpleVector* GetSamples() const          { return samples_; }
    SimpleVector*       GetSamples()                { return samples_; }
    void                DestroySamples();
    void                NewSamples( int );

    const int*          GetOLBuffer() const         { return olbuff_; }
    int*                GetOLBuffer()               { return olbuff_; }
    void                DestroyOLBuffer()           { if( olbuff_ ) free( olbuff_ ); olbuff_ = NULL; }
    void                NewOLBuffer( int size );

private:
    const char*         patoutname_;        //pattern of output files
    int                 dim_;               //dimensions of vectors
    int                 novariates_;        //number of variates determining one cluster
    int                 ctxtlen_;           //sample context length
    int                 noclusters_;        //number of clusters to produce
    int                 nocsamples_;        //number of samples to output per cluster
    bool                lgstnrmsampling_;   //sample from logistic-normal distribution
    bool                overlap_;           //overlap of support variates in generating data
    SimpleVector*       samples_;           //all generated samples
    int                 smplcounter_;       //sample counter
    int*                olbuff_;            //buffer for overlap
};


// /////////////////////////////////////////////////////////////////////////
// INLINES
//
// NewOLBuffer: create new overlap buffer
//
inline
void HDPSmplGenerator::NewOLBuffer( int size )
{
    if( size <= 0 )
        throw myruntime_error("HDPSmplGenerator: NewOLBuffer: Invalid size.");
    DestroyOLBuffer();
    olbuff_ = ( int* )malloc( size * sizeof( int ));
    if( olbuff_ == NULL )
        throw myruntime_error("HDPSmplGenerator: NewOLBuffer: Not enough memory.");
    memset( olbuff_, 0, size * sizeof( int ));
}

// -------------------------------------------------------------------------


#endif//__HDPSmplGenerator__
