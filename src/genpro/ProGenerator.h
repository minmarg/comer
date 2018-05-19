/***************************************************************************
 *   Copyright (C) 2014 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __ProGenerator__
#define __ProGenerator__

#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "rc.h"
#include "mystring.h"
#include "myexcept.h"
#include "ext/rng.h"
#include "libpro/srcpro/Serializer.h"


// _________________________________________________________________________
// Class ProGenerator
//

class ProGenerator {
    enum {
        pgnFREQ,
        pgnPSSM,
        pgnGAPS,
        pgnCnt
    };
public:
    //typedefs...
    typedef void ( ProGenerator::*PROCMET )( const char*, void* );
    //
    ProGenerator( const char* outdir, const char* namepat, const char* directory );
    ProGenerator( const char* outdir, const char* namepat, char* arguments[], int no_args );
    ~ProGenerator();

    void        Generate(); //generate random profiles

    size_t      GetNoProfiles() const { return no_profs_; }
    Uint64      GetTotalLength() const { return totlen_; }

    int         GetLengthToGenerate() const { return gen_lenpro_; }
    void        SetLengthToGenerate( int value ) { gen_lenpro_ = value; }

    int         GetNumberToGenerate() const { return gen_nopros_; }
    void        SetNumberToGenerate( int value ) { gen_nopros_ = value; }

    int         GetFragLength() const { return gen_fraglen_; }
    void        SetFragLength( int value ) { gen_fraglen_ = value; }

    bool        GetOperateWithSSE() const { return gen_opsse_; }
    void        SetOperateWithSSE( bool value ) { gen_opsse_ = value; }

protected:
    explicit ProGenerator();

    void    Init();                                     //necessary initialization

    void    ProcessInput( PROCMET, void* );             //process profiles given explicitly or found in directory
    void    ReadInProfile( const char* filename, void* );//read in profile
    void    GenerateProfile( const char* fname, const char* name );
    void    GenerateProfileFrag( const char* name, 
                FrequencyMatrix*, LogOddsMatrix*, GapScheme* );
    void    GenerateProfileSSE( const char* name, 
                FrequencyMatrix*, LogOddsMatrix*, GapScheme* );
    void    GenerateProfileSSEtest( const char* name, 
                FrequencyMatrix*, LogOddsMatrix*, GapScheme* );
    void    AddFrag( const char* name, FrequencyMatrix*, LogOddsMatrix*, GapScheme* );
    void    AddSSE( const char* name, char sse, FrequencyMatrix*, LogOddsMatrix*, GapScheme* );
    void    AddSSEtest( const char* name, char sse, FrequencyMatrix*, LogOddsMatrix*, GapScheme*,
                        size_t* sslen );
    void    AddColumn( const char* name, int pp, 
                const FrequencyMatrix*, const LogOddsMatrix*, const GapScheme*, 
                FrequencyMatrix*, LogOddsMatrix*, GapScheme* );
    void    RecalcBgProbs( FrequencyMatrix*, LogOddsMatrix*, GapScheme* );
    void    RecalcPostProbs( FrequencyMatrix*, LogOddsMatrix*, GapScheme* );

    const char* GetOutputDir() const { return outdir_; }
    const char* GetNamePattern() const { return nampat_; }
    const char* GetInputDir() const { return inpdir_; }
    const char* const* GetProfileNames() const { return pronames_; }
    const char* GetProfileName( int ) const;
    const int   GetNoProNames() const { return no_names_; }

    void    SetNoProfiles( size_t value ) { no_profs_ = value; }
    void    SetTotalLength( Uint64 value ) { totlen_ = value; }
    void    IncTotalLength( size_t length ) { totlen_ += length; }

    void        PushProfile( FrequencyMatrix*, LogOddsMatrix*, GapScheme* );
    const FrequencyMatrix*  GetProFREQAt( int ) const;
    const LogOddsMatrix*    GetProPSSMAt( int ) const;
    const GapScheme*        GetProGAPSAt( int ) const;

private:
    void        Realloc( int newsize );
    void        DestroyProfiles();

private:
    const char*         outdir_;            //output directory name
    const char*         nampat_;            //filename pattern for new generated profiles
    const char*         inpdir_;            //input directory of profiles
    const char* const*  pronames_;          //names of profiles used to generate new profiles
    const int           no_names_;          //number of profiles provided with arguments

    size_t              no_profs_;          //number of profiles read
    Uint64              totlen_;            //total length of profiles read

    void*             (*profiles_)[pgnCnt]; //profiles
    int                 allocated_;         //allocated cells for profiles_

    int                 gen_lenpro_;        //length of profiles to be generated
    int                 gen_nopros_;        //number of profiles to generate
    int                 gen_fraglen_;       //length of indivisible profile fragment
    bool                gen_opsse_;         //operate with SS elements

    MTRng               rng_p_;//RNG for profile indices
    MTRng               rng_pp_;//RNG for profile positions
    MTRng               rng_sse_p_;//RNG for profile indices in generating SSE
    MTRng               rng_sse_pp_;//RNG for profile positions in generating SSE

    Serializer          serializer;         //object to serialize data to file
};


////////////////////////////////////////////////////////////////////////////
// INLINES
//
inline
const char* ProGenerator::GetProfileName( int n ) const
{
#ifdef __DEBUG__
    if( no_names_ <= n )
        throw myruntime_error("ProGenerator: Memory access error.");
#endif
    return pronames_[n];
}

// -------------------------------------------------------------------------
// PushProfile: save profile in memory
//
inline
void ProGenerator::PushProfile( FrequencyMatrix* freq, LogOddsMatrix* pssm, GapScheme* gaps )
{
    if( allocated_ <= no_profs_ )
        Realloc( 1 + TIMES2( allocated_ ));
    if( profiles_ == NULL )
        throw myruntime_error("ProGenerator: PushProfile: Memory access error.");
    profiles_[no_profs_][pgnFREQ] = freq;
    profiles_[no_profs_][pgnPSSM] = pssm;
    profiles_[no_profs_][pgnGAPS] = gaps;
    no_profs_++;
}

// -------------------------------------------------------------------------
// GetPSSMAt: get profile terms at the position
//
inline
const FrequencyMatrix* ProGenerator::GetProFREQAt( int n ) const
{
    if( profiles_ == NULL || allocated_ <= n || no_profs_ <= n || !profiles_[n])
        throw myruntime_error("ProGenerator: GetProFREQAt: Memory access error.");
    return (const FrequencyMatrix*)profiles_[n][pgnFREQ];
}

inline
const LogOddsMatrix* ProGenerator::GetProPSSMAt( int n ) const
{
    if( profiles_ == NULL || allocated_ <= n || no_profs_ <= n || !profiles_[n])
        throw myruntime_error("ProGenerator: GetProPSSMAt: Memory access error.");
    return (const LogOddsMatrix*)profiles_[n][pgnPSSM];
}

inline
const GapScheme* ProGenerator::GetProGAPSAt( int n ) const
{
    if( profiles_ == NULL || allocated_ <= n || no_profs_ <= n || !profiles_[n])
        throw myruntime_error("ProGenerator: GetProGAPSAt: Memory access error.");
    return (const GapScheme*)profiles_[n][pgnGAPS];
}

// -------------------------------------------------------------------------
//

#endif//__ProGenerator__
