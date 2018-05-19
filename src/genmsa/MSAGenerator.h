/***************************************************************************
 *   Copyright (C) 2014 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __MSAGenerator__
#define __MSAGenerator__

#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "rc.h"
#include "mystring.h"
#include "myexcept.h"
#include "ext/rng.h"
#include "ext/ivector.h"


// _________________________________________________________________________
// Class MSAGenerator
//

class MSAGenerator {
public:
    //typedefs...
    typedef void ( MSAGenerator::*PROCMET )( const char*, void* );
    //
    MSAGenerator( const char* outdir, const char* namepat, const char* directory );
    MSAGenerator( const char* outdir, const char* namepat, char* arguments[], int no_args );
    ~MSAGenerator();

    void        Generate(); //generate random MSAs

    size_t      GetNoMSAs() const { return no_msas_; }
    Uint64      GetTotalLength() const { return totlen_; }

    int         GetLengthToGenerate() const { return gen_lenmsa_; }
    void        SetLengthToGenerate( int value ) { gen_lenmsa_ = value; }

    int         GetNumberToGenerate() const { return gen_nomsas_; }
    void        SetNumberToGenerate( int value ) { gen_nomsas_ = value; }

    int         GetFragLength() const { return gen_fraglen_; }
    void        SetFragLength( int value ) { gen_fraglen_ = value; }

protected:
    explicit MSAGenerator();

    void    Init();//initialization

    void    ProcessInput( PROCMET, void* );//process alignments given
    void    ReadInMSA( const char* filename, void* );//read in alignment
    void    GenerateMSA( const char* fname, const char* name );
    void    GenerateMSAFrag( const char* name, mystring**, int width );
    void    AddFrag( const char* name, mystring**, int width, int frag );
    void    AddColumn( const char* name, int pp, const mystring**, const int, mystring**, int );

    const char* GetOutputDir() const { return outdir_; }
    const char* GetNamePattern() const { return nampat_; }
    const char* GetInputDir() const { return inpdir_; }
    const char* const* GetMSANames() const { return msanames_; }
    const char* GetMSAName( int ) const;
    const int   GetNoMSANames() const { return no_names_; }

    void    SetNoMSAs( size_t value ) { no_msas_ = value; }
    void    SetTotalLength( Uint64 value ) { totlen_ = value; }
    void    IncTotalLength( size_t length ) { totlen_ += length; }

    void    MakeMState( const mystring**, const int, Ivector* ) const;

    void    PushMSA( mystring**, int, Ivector* );
    const mystring**    GetMSAAt( int ) const;
    mystring**          GetMSAAt( int );
    int                 GetMSAWidthAt( int ) const;
    const Ivector*      GetMSAStatesAt( int ) const;

    int     GetMaxWidth() const { return maxwdt_; }
    void    SetMaxWidth( int value ) { if( maxwdt_ < value ) maxwdt_ = value; }

private:
    void        Realloc( int newsize );
    void        DestroyMSAs();

private:
    const char*         outdir_;            //output directory name
    const char*         nampat_;            //filename pattern for new generated files
    const char*         inpdir_;            //input directory of alignments
    const char* const*  msanames_;          //names of files used to generate new alignments
    const int           no_names_;          //number of files provided with arguments

    size_t              no_msas_;           //number of files read
    Uint64              totlen_;            //total length of alignments read

    mystring***         msas_;              //alignments
    int*                msawdts_;           //alignment widths (#sequences within each alignment)
    int                 maxwdt_;            //maximum width of alignments
    Ivector**           mstates_;           //sequence of indices of m states
    int                 allocated_;         //allocated cells for msas_

    int                 gen_lenmsa_;        //length of alignments to be generated
    int                 gen_nomsas_;        //number of alignments to generate
    int                 gen_fraglen_;       //length of indivisible alignment fragment

    MTRng               rng_p_;//RNG for file indices
    MTRng               rng_pp_;//RNG for alignment positions
};


////////////////////////////////////////////////////////////////////////////
// INLINES
//
inline
const char* MSAGenerator::GetMSAName( int n ) const
{
#ifdef __DEBUG__
    if( n < 0 || no_names_ <= n )
        throw myruntime_error("MSAGenerator: GetMSAName: Memory access error.");
#endif
    return msanames_[n];
}

// -------------------------------------------------------------------------
// PushProfile: save alignment in memory
//
inline
void MSAGenerator::PushMSA( mystring** msa, int width, Ivector* mstate )
{
    if( allocated_ <= no_msas_ )
        Realloc( 1 + TIMES2( allocated_ ));
    if( msas_ == NULL || msawdts_ == NULL || mstates_ == NULL )
        throw myruntime_error("MSAGenerator: PushMSA: Memory access error.");
    msas_[no_msas_] = msa;
    msawdts_[no_msas_] = width;
    mstates_[no_msas_] = mstate;
    no_msas_++;
}

// -------------------------------------------------------------------------
// GetMSAAt: get alignment at the position
//
inline
const mystring** MSAGenerator::GetMSAAt( int n ) const
{
    if( msas_ == NULL || n < 0 || allocated_ <= n || no_msas_ <= n || !msas_[n])
        throw myruntime_error("MSAGenerator: GetMSAAt: Memory access error.");
    return (const mystring**)msas_[n];
}

inline
mystring** MSAGenerator::GetMSAAt( int n )
{
    if( msas_ == NULL || n < 0 || allocated_ <= n || no_msas_ <= n || !msas_[n])
        throw myruntime_error("MSAGenerator: GetMSAAt: Memory access error.");
    return msas_[n];
}

inline
int MSAGenerator::GetMSAWidthAt( int n ) const
{
    if( msawdts_ == NULL || n < 0 || allocated_ <= n || no_msas_ <= n)
        throw myruntime_error("MSAGenerator: GetMSAWidthAt: Memory access error.");
    return msawdts_[n];
}

inline
const Ivector* MSAGenerator::GetMSAStatesAt( int n ) const
{
    if( mstates_ == NULL || n < 0 || allocated_ <= n || no_msas_ <= n )
        throw myruntime_error("MSAGenerator: GetMSAStatesAt: Memory access error.");
    return mstates_[n];
}

// -------------------------------------------------------------------------
//

#endif//__MSAGenerator__
