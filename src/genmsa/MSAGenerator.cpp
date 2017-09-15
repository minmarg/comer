/***************************************************************************
 *   Copyright (C) 2014 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "ext/psl.h"
#include "MSAGenerator.h"


const char* slg_EXTN = ".fa";

// -------------------------------------------------------------------------
// CLASS MSAGenerator
//
// Constructors
//
MSAGenerator::MSAGenerator( const char* outdir, const char* namepat, const char* directory )
:
    outdir_( outdir ),
    nampat_( namepat ),
    inpdir_( directory ),
    msanames_( NULL ),
    no_names_( 0 ),

    no_msas_( 0 ),
    totlen_( 0 ),
    maxwdt_( 0 ),
    msawdts_( NULL ),
    msas_( NULL ),
    mstates_( NULL ),
    allocated_( 0 ),
    gen_lenmsa_( 0 ),
    gen_nomsas_( 0 ),
    gen_fraglen_( 0 )
{
    Init();
}

// Constructor
//
MSAGenerator::MSAGenerator( const char* outdir, const char* namepat, char* arguments[], int no_args )
:
    outdir_( outdir ),
    nampat_( namepat ),
    inpdir_( NULL ),
    msanames_( arguments ),
    no_names_( no_args ),

    no_msas_( 0 ),
    totlen_( 0 ),
    maxwdt_( 0 ),
    msawdts_( NULL ),
    msas_( NULL ),
    mstates_( NULL ),
    allocated_( 0 ),
    gen_lenmsa_( 0 ),
    gen_nomsas_( 0 ),
    gen_fraglen_( 0 )
{
    Init();
}

// Default constructor
//
MSAGenerator::MSAGenerator()
:
    outdir_( NULL ),
    nampat_( NULL ),
    inpdir_( NULL ),
    msanames_( NULL ),
    no_names_( 0 ),
    no_msas_( 0 ),
    totlen_( 0 ),
    maxwdt_( 0 ),
    msawdts_( NULL ),
    msas_( NULL ),
    mstates_( NULL ),
    allocated_( 0 ),
    gen_lenmsa_( 0 ),
    gen_nomsas_( 0 ),
    gen_fraglen_( 0 )
{
    throw( myruntime_error("MSAGenerator: Default initialization is not allowed."));
}

// Destructor
//
MSAGenerator::~MSAGenerator()
{
    DestroyMSAs();
}

// -------------------------------------------------------------------------
// Init: initialization
//
void MSAGenerator::Init()
{
    time_t  tm;
    time( &tm );
    rng_p_.Set((unsigned long)(size_t)( &rng_p_ ) +(unsigned long)tm );
    rng_pp_.Set((unsigned long)(size_t)( &rng_pp_ ) +(unsigned long)tm );
    Realloc( 1000 );
}

// -------------------------------------------------------------------------
// Realloc: reallocate memory to contain alignments
//
void MSAGenerator::Realloc( int newsize )
{
    Ivector**   tmp_mstates;
    mystring*** tmp_msas;
    int*        tmp_msawdts;

    if( newsize <= allocated_ )
        return;

    if( allocated_ < 1 ) {
        tmp_mstates = ( Ivector**)malloc( newsize * sizeof(void*));
        tmp_msas = ( mystring***)malloc( newsize * sizeof(void*));
        tmp_msawdts = ( int*)malloc( newsize * sizeof(int));
    } else {
        tmp_mstates = ( Ivector**)realloc( mstates_, newsize * sizeof(void*));
        tmp_msas = ( mystring***)realloc( msas_, newsize * sizeof(void*));
        tmp_msawdts = ( int*)realloc( msawdts_, newsize * sizeof(int));
    }

    if( tmp_mstates == NULL || tmp_msas == NULL || tmp_msawdts == NULL )
        throw myruntime_error("MSAGenerator: Realloc: Not enough memory.");

    mstates_ = tmp_mstates;
    msas_ = tmp_msas;
    msawdts_ = tmp_msawdts;

    //fill uninitialized memory with zeros
    memset( mstates_ + allocated_, 0, (newsize-allocated_)*sizeof(void*));
    memset( msas_ + allocated_, 0, (newsize-allocated_)*sizeof(void*));
    memset( msawdts_ + allocated_, 0, (newsize-allocated_)*sizeof(int));

    allocated_ = newsize;
}

// -------------------------------------------------------------------------
// DestroyMSAs: deallocate memory for MSAs
//
void MSAGenerator::DestroyMSAs()
{
    const Ivector* sv;
    const mystring** msa;
    int n, i, w;
    for( n = 0; n < no_msas_ && n < allocated_; n++ ) {
        sv = GetMSAStatesAt(n);
        msa = (const mystring**)GetMSAAt(n);
        w = GetMSAWidthAt(n);
        if( sv ) {
            delete sv;
            mstates_[n] = NULL;
        }
        if( msa ) {
            for( i = 0; i < w; i++ )
                if( msa[i]) {
                    delete msa[i];
                    msa[i] = NULL;
                }
            free( msa );
            msas_[n] = NULL;
        }
    }
    if( mstates_ ) { free( mstates_ ); mstates_ = NULL; }
    if( msas_ ) { free( msas_ ); msas_ = NULL; }
    if( msawdts_ ) { free( msawdts_ ); msawdts_ = NULL; }
    allocated_ = 0;
    no_msas_ = 0;
}

// -------------------------------------------------------------------------
// Generate: generate random MSAs
//
void MSAGenerator::Generate()
{
    mystring    outfname, name;
    char        namebuf[BUF_MAX];
    int n;

    message("Reading...");

    //read in all alignments
    ProcessInput( &MSAGenerator::ReadInMSA, NULL );

    try {
        if( 0 < GetNumberToGenerate()) {
            sprintf( namebuf, "Generating %d MSA(s)...", GetNumberToGenerate());
            message( namebuf );
        }

        for( n = 0; n < GetNumberToGenerate(); n++ ) {
            sprintf( namebuf, "_%05d", n );
            outfname = mystring( GetOutputDir()) + DIRSEP + GetNamePattern() + namebuf + slg_EXTN;
            name = mystring( GetNamePattern()) + namebuf;
            GenerateMSA( outfname.c_str(), name.c_str());
        }

    } catch( myexception const& ex )
    {
        throw myruntime_error( ex.what(), ex.eclass());
    }

    message( "Finished." );
}

// -------------------------------------------------------------------------
// ProcessInput: process alignments given or found in directory
//
void MSAGenerator::ProcessInput( PROCMET proc_method, void* args )
{
    if( proc_method == NULL )
        throw myruntime_error("MSAGenerator: ProcessInput: Null processing method.");

    if( GetInputDir()) {
        struct stat     info;
        dirent*         entry = NULL;
        DIR*            direct = opendir( GetInputDir());
        mystring        dirname = GetInputDir() + mystring( DIRSEP );
        mystring        fullname;

        if( !direct )
            throw myruntime_error("MSAGenerator: ProcessInput: Unable to open input directory.");

        while(( entry = readdir( direct )))
        {
            fullname = dirname + entry->d_name;

            if( stat( fullname.c_str(), &info ) == -1 ) {
                error( mystring( mystring("Cannot stat '") + mystring(entry->d_name) + "'; skipping.").c_str());
                continue;
            }
            if(( S_IFMT & info.st_mode ) != S_IFREG )
//             if( !S_ISREG( info.st_mode ))
                continue;

            ( this->*proc_method )( fullname.c_str(), args );
        }
        closedir( direct );

    } else
        if( GetMSANames()) {
            bool option = false;

            for( int n = 0; n < GetNoMSANames(); n++ ) {
                if( option ) {
                    option = false;
                    continue;
                }
                if( GetMSAName(n) && GetMSAName(n)[0] == '-' && GetMSAName(n)[2] == 0 ) {
                    if( GetMSAName(n)[1] != 'v' )
                        option = true;//assume an option
                    continue;
                }
                ( this->*proc_method )( GetMSAName(n), args );
            }
        }
}

// -------------------------------------------------------------------------
// MakeMState: make m-state vector given the alignment
//
void MSAGenerator::MakeMState( const mystring** msa, const int width, Ivector* msvec ) const
{
    int p;
    if( msa == NULL || width < 1 || msvec == NULL )
        throw myruntime_error("MSAGenerator: MakeMState: Null arguments.");
    if( msa[0] == NULL )
        throw myruntime_error("MSAGenerator: MakeMState: Null leading sequence.");
    for( p = 0; p < msa[0]->length(); p++ ) {
        if((*msa[0])[p] != '-' )
            msvec->Push(p);
    }
}

// -------------------------------------------------------------------------
// ReadInMSA: read alignment and save it in memory
//
void MSAGenerator::ReadInMSA( const char* filename, void* )
{
    FILE* fp = NULL;
    int emsg;
    const char* p;
    mystring    line;
    mystring**  msa = NULL;
    Ivector*    svec = NULL;
    int length = 1000;//default sequence length
    int width = 500;
    int nn = -1;//number of sequences within msa
    int i = 0;

    if(( fp = fopen( filename, "r")) == NULL ) {
        error( mystring( mystring("Failed to open '") + filename + mystring("'. Skipped.")).c_str());
        return;
    }

    try {
        svec = new Ivector;
        msa = (mystring**)malloc( width*sizeof(void*));
        if( svec == NULL || msa == NULL )
            throw myruntime_error("MSAGenerator: ReadInMSA: Not enough memory.");
        memset(msa,0,width*sizeof(void*));

        //READ
        while(!feof(fp)) {
            if(( emsg = skip_comments( fp, line )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( !line.length())
                continue;

            p = line.c_str();
            for( ; *p && (*p==' ' || *p=='\t' || *p=='\n' || *p=='\r'); p++ );

            if( *p=='>') {
                if( 0 <= nn && msa[nn] == NULL )
                    throw myruntime_error("MSAGenerator: ReadInMSA: Invalid file format.");
                if( 0 < nn && msa[nn]->length() != msa[nn-1]->length())
                    throw myruntime_error("MSAGenerator: ReadInMSA: Inconsistent sequence lengths.");
                nn++;
                if( width <= nn ) {
                    msa = (mystring**)realloc( msa, (width+width)*sizeof(void*));
                    if( msa == NULL )
                        throw myruntime_error("MSAGenerator: ReadInMSA: Not enough memory.");
                    memset(msa+width,0,width*sizeof(void*));
                    width += width;
                }
                continue;
            }

            if( nn < 0 )
                throw myruntime_error("MSAGenerator: ReadInMSA: No sequence description.");

            if( msa[nn] == NULL ) {
                msa[nn] = new mystring;
                if( msa[nn] == NULL )
                    throw myruntime_error("MSAGenerator: ReadInMSA: Not enough memory.");
                if( 0 < nn )
                    msa[nn]->reserve( msa[nn-1]->length()+10);
                else
                    msa[nn]->reserve( length );
            }

            for( ; *p && *p!='\n' && *p!='\r'; p++ )
                if( *p!=' ' && *p!='\t' )
                    msa[nn]->append(*p);
        }

        MakeMState( (const mystring**)msa, nn+1, svec );
        IncTotalLength( svec->GetSize());

        if( 0 <= nn ) {
            PushMSA( msa, nn+1, svec );
            SetMaxWidth( nn+1 );
        }
    } catch( myexception const& ex )
    {
        if( svec ) {
            delete svec;
            svec = NULL;
        }
        if( msa ) {
            for( i = 0; i < width; i++ )
                if( msa[i]) {
                    delete msa[i];
                    msa[i] = NULL;
                }
            free( msa );
            msa = NULL;
        }
        error( mystring( ex.what() + mystring(" Skipping '") + filename + mystring("'.")).c_str());
    }

    if( fp )
        fclose( fp );
}

// -------------------------------------------------------------------------
// GenerateMSA: generate one alignment and write it to file
//
void MSAGenerator::GenerateMSA( const char* outfilename, const char* name )
{
    FILE* fp = NULL;
    mystring**  msa;
    int width = GetMaxWidth();
    int length = GetLengthToGenerate();
    int i;

    if( GetNoMSAs() < 1 )
        return;
    if( width < 1 )
        return;

    try {
        msa = (mystring**)malloc( width*sizeof(void*));
        if( msa == NULL )
            throw myruntime_error("MSAGenerator: GenerateMSA: Not enough memory.");
        memset(msa,0,width*sizeof(void*));

        for( i = 0; i < width; i++ ) {
            msa[i] = new mystring;
            if( msa[i] == NULL )
                throw myruntime_error("MSAGenerator: GenerateMSA: Not enough memory.");
            msa[i]->reserve( length+10 );
        }

        GenerateMSAFrag( name, msa, width );

        if(( fp = fopen( outfilename, "w")) == NULL )
            throw myruntime_error("MSAGenerator: GenerateMSA: Failed to open file for writing.");

        for( i = 0; i < width; i++ ) {
            if(i)
                fprintf( fp, ">SGEN_%06d\n", i+1 );
            else
                fprintf( fp, ">%s Generated multiple sequence alignment\n", name );
            fprintf( fp, "%s\n", msa[i]->c_str());
        }

    } catch( myexception const& ex )
    {
        error( ex.what());
    }
    if( fp )
        fclose( fp );
    if( msa ) {
        for( i = 0; i < width; i++ )
            if( msa[i]) {
                delete msa[i];
                msa[i] = NULL;
            }
        free( msa );
        msa = NULL;
    }
}

// -------------------------------------------------------------------------
// GenerateMSAFrag: generate an MSA by shuffling and laying down 
//  alignment fragments
//
void MSAGenerator::GenerateMSAFrag( const char* name, mystring** msa, int width )
{
    if( msa == NULL )
        throw myruntime_error("MSAGenerator: GenerateMSAFrag: Null target alignment.");

    if( GetNoMSAs() < 1 )
        return;

    const int   lenerr = 5;
    const int   lcnFRAG = GetFragLength();//fragment length
    int         length = GetLengthToGenerate();
    int lc = -lcnFRAG, lp = 0, cnt = 0;

    while(( lc += lcnFRAG ) < length - lenerr ) {
        if( lc == lp )
            cnt++;
        else
            cnt = 0;
        if( 10 <= cnt )
            throw myruntime_error("MSAGenerator: GenerateMSAFrag: Live lock in generating an MSA.");
        lp = lc;
        AddFrag( name, msa, width );
    }
}

// -------------------------------------------------------------------------
// AddFrag: add alignment fragment to the MSA being generated
//
void MSAGenerator::AddFrag( const char* name, mystring** msa, int width )
{
    if( msa == NULL )
        throw myruntime_error("MSAGenerator: AddFrag: Null arguments.");

    const int   lcnFRAG = GetFragLength();//fragment length
    const mystring** srcmsa = NULL;//source MSA
    const Ivector*   srcsv = NULL;//source state vector 
    int         srcwdt = 0;//width of source MSA

    if( GetNoMSAs() < 1 )
        throw myruntime_error("MSAGenerator: AddFrag: No alignments.");
    if( lcnFRAG < 1 )
        throw myruntime_error("MSAGenerator: AddFrag: Invalid fragment length.");

    const int   length = GetLengthToGenerate();
    size_t  nomsas = GetNoMSAs();
    double  urv = rng_p_.GetDouble();
    size_t  ind = ( size_t )rint( urv *(double)(nomsas-1));

    if( nomsas <= ind )
        throw myruntime_error("MSAGenerator: AddFrag: Invalid file index generated.");

    srcmsa = (const mystring**)GetMSAAt(ind);
    srcsv = GetMSAStatesAt(ind);
    srcwdt = GetMSAWidthAt(ind);

    if( srcmsa == NULL || srcsv == NULL )
        throw myruntime_error("MSAGenerator: AddFrag: Null source MSA.");
    if( srcwdt < 1 )
        throw myruntime_error("MSAGenerator: AddFrag: Invalid width of source MSA.");
    if( width < srcwdt )
        throw myruntime_error("MSAGenerator: AddFrag: Too small target width.");

    int     lenmst = srcsv->GetSize();
    int     left = SLC_MIN( 40, lenmst * 2 / 10 );
    int     rght = lenmst - 1 - left;
    double  urv2 = rng_pp_.GetDouble();
    int     pp = ( int )rint( urv2 *(double)(lenmst-1));
    mystring gaps;
    int     n, ii, iii;

    if( rght < pp + lcnFRAG )
        pp = rght - lcnFRAG;
    if( pp < left )
        pp = left;

    if( lenmst <= pp + lcnFRAG )
        pp = lenmst - lcnFRAG;

    if( pp < 0 )
        return;//too short alignment

    ii = srcsv->GetValueAt(pp);
    iii = srcsv->GetValueAt(pp+lcnFRAG-1);

    if( iii < ii || srcmsa[0]->length() <= iii )
        throw myruntime_error("MSAGenerator: AddFrag: Invalid positions from state vector.");

    for( n = 0; n < iii-ii+1; n++ )
        gaps.append('-');

    //copy alignment fragment
    for( n = 0; n < srcwdt; n++ )
        msa[n]->append( srcmsa[n]->substr(ii,iii-ii+1));
    for( ; n < width; n++ )
        msa[n]->append( gaps );
}
