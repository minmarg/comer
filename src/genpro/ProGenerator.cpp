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
#include <time.h>
#include "libpro/srcpro/GapScheme.h"
#include "libpro/srcpro/DistributionMatrix.h"
#include "ProGenerator.h"


const char* slg_EXTN = ".pro";

// -------------------------------------------------------------------------
// CLASS ProGenerator
//
// Constructors
//
ProGenerator::ProGenerator( const char* outdir, const char* namepat, const char* directory )
:
    outdir_( outdir ),
    nampat_( namepat ),
    inpdir_( directory ),
    pronames_( NULL ),
    no_names_( 0 ),

    no_profs_( 0 ),
    totlen_( 0 ),
    profiles_( NULL ),
    allocated_( 0 ),
    gen_lenpro_( 0 ),
    gen_nopros_( 0 ),
    gen_fraglen_( 0 ),
    gen_opsse_( false )
{
    Init();
}

// Constructor
//
ProGenerator::ProGenerator( const char* outdir, const char* namepat, char* arguments[], int no_args )
:
    outdir_( outdir ),
    nampat_( namepat ),
    inpdir_( NULL ),
    pronames_( arguments ),
    no_names_( no_args ),

    no_profs_( 0 ),
    totlen_( 0 ),
    profiles_( NULL ),
    allocated_( 0 ),
    gen_lenpro_( 0 ),
    gen_nopros_( 0 ),
    gen_fraglen_( 0 ),
    gen_opsse_( false )
{
    Init();
}

// Default constructor
//
ProGenerator::ProGenerator()
:
    outdir_( NULL ),
    nampat_( NULL ),
    inpdir_( NULL ),
    pronames_( NULL ),
    no_names_( 0 ),
    no_profs_( 0 ),
    totlen_( 0 ),
    profiles_( NULL ),
    allocated_( 0 ),
    gen_lenpro_( 0 ),
    gen_nopros_( 0 ),
    gen_fraglen_( 0 ),
    gen_opsse_( false )
{
    throw( myruntime_error("ProGenerator: Default initialization is not allowed."));
}

// Destructor
//
ProGenerator::~ProGenerator()
{
    DestroyProfiles();
}

// -------------------------------------------------------------------------
// Init: initialization
//
void ProGenerator::Init()
{
    time_t  tm;
    time( &tm );
    rng_p_.Set((unsigned long)(size_t)( &rng_p_ ) +(unsigned long)tm );
    rng_pp_.Set((unsigned long)(size_t)( &rng_pp_ ) +(unsigned long)tm );
    rng_sse_p_.Set((unsigned long)(size_t)( &rng_sse_p_ ) +(unsigned long)tm );
    rng_sse_pp_.Set((unsigned long)(size_t)( &rng_sse_pp_ ) +(unsigned long)tm );
    Realloc( 1000 );
}

// -------------------------------------------------------------------------
// Realloc: reallocate memory to contain profiles
//
void ProGenerator::Realloc( int newsize )
{
    void* (*tmp_profiles )[pgnCnt];

    if( newsize <= allocated_ )
        return;

    if( allocated_ < 1 ) {
        tmp_profiles = ( void*(*)[pgnCnt] )malloc( newsize * pgnCnt * sizeof(void*));

    } else {
        tmp_profiles = ( void*(*)[pgnCnt] )realloc( profiles_, newsize * pgnCnt * sizeof(void*));
    }

    if( tmp_profiles == NULL )
        throw myruntime_error("ProGenerator: Realloc: Not enough memory.");

    profiles_ = tmp_profiles;

    // fill uninitialized memory with zeros
    memset( profiles_ + allocated_, 0, sizeof(void*) * pgnCnt * ( newsize - allocated_ ));

    allocated_ = newsize;
}

// -------------------------------------------------------------------------
// DestroyProfiles: deallocate memory for profiles
//
void ProGenerator::DestroyProfiles()
{
    int n;
    for( n = 0; n < no_profs_ && n < allocated_; n++ ) {
        if( profiles_[n][pgnFREQ]) {
            delete (FrequencyMatrix*)profiles_[n][pgnFREQ];
            profiles_[n][pgnFREQ] = NULL;
        }
        if( profiles_[n][pgnPSSM]) {
            delete (LogOddsMatrix*)profiles_[n][pgnPSSM];
            profiles_[n][pgnPSSM] = NULL;
        }
        if( profiles_[n][pgnGAPS]) {
            delete (GapScheme*)profiles_[n][pgnGAPS];
            profiles_[n][pgnGAPS] = NULL;
        }
    }
    free( profiles_ );
    profiles_ = NULL;
    allocated_ = 0;
    no_profs_ = 0;
}

// -------------------------------------------------------------------------
// Generate: generate random profiles
//
void ProGenerator::Generate()
{
    mystring    outfname, name;
    char        namebuf[BUF_MAX];
    int n;

    message("Reading...");

    //read in all profiles
    ProcessInput( &ProGenerator::ReadInProfile, NULL );

    try {
        if( 0 < GetNumberToGenerate()) {
            sprintf( namebuf, "Generating %d profile(s)...", GetNumberToGenerate());
            message( namebuf );
        }

        for( n = 0; n < GetNumberToGenerate(); n++ ) {
            sprintf( namebuf, "_%05d", n );
            outfname = mystring( GetOutputDir()) + DIRSEP + GetNamePattern() + namebuf + slg_EXTN;
            name = mystring( GetNamePattern()) + namebuf;
            GenerateProfile( outfname.c_str(), name.c_str());
        }

    } catch( myexception const& ex )
    {
        throw myruntime_error( ex.what(), ex.eclass());
    }

    message( "Finished." );
}

// -------------------------------------------------------------------------
// ProcessInput: process profiles given explicitly or found in directory
//
void ProGenerator::ProcessInput( PROCMET proc_method, void* args )
{
    if( proc_method == NULL )
        throw myruntime_error("ProGenerator: ProcessInput: Null processing method.");

    if( GetInputDir()) {
        struct stat     info;
        dirent*         entry = NULL;
        DIR*            direct = opendir( GetInputDir());
        mystring        dirname = GetInputDir() + mystring( DIRSEP );
        mystring        fullname;

        if( !direct )
            throw myruntime_error("ProGenerator: ProcessInput: Unable to open input directory.");

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
        if( GetProfileNames()) {
            bool    option = false;

            for( int n = 0; n < GetNoProNames(); n++ ) {
                if( option ) {
                    option = false;
                    continue;
                }
                if( GetProfileName(n) && GetProfileName(n)[0] == '-' && GetProfileName(n)[2] == 0 ) {
                    if( GetProfileName(n)[1] != 'v' )
                        option = true;//assume an option
                    continue;
                }
                ( this->*proc_method )( GetProfileName(n), args );
            }
        }
}

// -------------------------------------------------------------------------
// ReadInProfile: read profile and save it in memory
//
void ProGenerator::ReadInProfile( const char* filename, void* )
{
    FrequencyMatrix*    freq = new FrequencyMatrix;
    LogOddsMatrix*      pssm = new LogOddsMatrix;
    GapScheme*          gaps = new GapScheme;
    double              prob = 0.0;

    if( freq == NULL || pssm == NULL || gaps == NULL )
        throw myruntime_error("ProGenerator: ReadInProfile: Not enough memory.");

    try {
        serializer.ReadProfile( filename, *freq, *pssm, *gaps );
        IncTotalLength( pssm->GetColumns());

//         freq.CheckForAllZeros();//to ensure valid scores

        if( GetOperateWithSSE() && !pssm->GetSSSSet())
            throw myruntime_error("ProGenerator: ReadInProfile: No SS predictions.");

        PushProfile( freq, pssm, gaps );

    } catch( myexception const& ex )
    {
        if( freq ) { delete freq; freq = NULL; } 
        if( pssm ) { delete pssm; pssm = NULL; }
        if( gaps ) { delete gaps; gaps = NULL; }
        error( mystring( ex.what() + mystring(" Skipping '") + filename + mystring("'.")).c_str());
    }
}

// -------------------------------------------------------------------------
// GenerateProfile: generate one profile and write it to file
//
void ProGenerator::GenerateProfile( const char* outfilename, const char* name )
{
    FrequencyMatrix genfreq;
    LogOddsMatrix   genpssm;
    GapScheme       gengaps;
    int length = GetLengthToGenerate();

    if( GetNoProfiles() < 1 )
        return;

    genfreq.Reserve( length+10 );
    genpssm.Reserve( length+10 );
    gengaps.Reserve( length+10 );

    if( GetOperateWithSSE())
        GenerateProfileSSE( name, &genfreq, &genpssm, &gengaps );
    else
        GenerateProfileFrag( name, &genfreq, &genpssm, &gengaps );

    RecalcBgProbs( &genfreq, &genpssm, &gengaps );

    serializer.WriteProfile( outfilename, genfreq, genpssm, gengaps );
}

// -------------------------------------------------------------------------
// GenerateProfileFrag: generate a profile by shuffling and laying down 
//  profile fragments
//
void ProGenerator::GenerateProfileFrag( const char* name, 
          FrequencyMatrix* genfreq, LogOddsMatrix* genpssm, GapScheme* gengaps  )
{
    if( genfreq == NULL || genpssm == NULL || gengaps == NULL )
        throw myruntime_error("ProGenerator: GenerateProfileFrag: Null target profile.");

    if( GetNoProfiles() < 1 )
        return;

    const int   lenerr = 5;
    int         length = GetLengthToGenerate();
    int lc = 0, lp = 0, cnt = 0;

    while(( lc = genpssm->GetColumns()) < length - lenerr ) {
        if( lc == lp )
            cnt++;
        else
            cnt = 0;
        if( 10 <= cnt )
            throw myruntime_error("ProGenerator: GenerateProfileFrag: Live lock in generating a profile.");
        lp = lc;
        AddFrag( name, genfreq, genpssm, gengaps );
    }
}

// -------------------------------------------------------------------------
// AddFrag: add profile fragment to profile being generated
//
void ProGenerator::AddFrag( const char* name, 
          FrequencyMatrix* genfreq, LogOddsMatrix* genpssm, GapScheme* gengaps )
{
    if( genfreq == NULL || genpssm == NULL || gengaps == NULL )
        throw myruntime_error("ProGenerator: AddFrag: Null arguments.");

    const int   lcnFRAG = GetFragLength();//fragment length
    const FrequencyMatrix*  freq = NULL;
    const LogOddsMatrix*    pssm = NULL;//profile providing SS element
    const GapScheme*        gaps = NULL;

    if( GetNoProfiles() < 1 )
        throw myruntime_error("ProGenerator: AddFrag: No profiles.");
    if( lcnFRAG < 1 )
        throw myruntime_error("ProGenerator: AddFrag: Invalid fragment length.");

    const int   length = GetLengthToGenerate();
    size_t  nopros = GetNoProfiles();
    double  urv = rng_sse_p_.GetDouble();
    size_t  ind = ( size_t )rint( urv *(double)(nopros-1));

    if( nopros <= ind )
        throw myruntime_error("ProGenerator: AddFrag: Invalid profile index generated.");

    freq = GetProFREQAt(ind);
    pssm = GetProPSSMAt(ind);
    gaps = GetProGAPSAt(ind);

    if( freq == NULL || pssm == NULL || gaps == NULL )
        throw myruntime_error("ProGenerator: AddFrag: Null profile.");

    int     lenpssm = pssm->GetColumns();
    int     left = SLC_MIN( 40, lenpssm * 2 / 10 );
    int     rght = lenpssm - 1 - left;
    double  urv2 = rng_sse_pp_.GetDouble();
    int     pp = ( int )rint( urv2 *(double)(lenpssm-1));
    int     n, ii;

    if( rght < pp + lcnFRAG )
        pp = rght - lcnFRAG;
    if( pp < left )
        pp = left;

    if( lenpssm <= pp + lcnFRAG )
        pp = lenpssm - lcnFRAG;

    if( pp < 0 )
        return;//too short profile

    //copy profile fragment
    for( ii = 0; pp < lenpssm && ii < lcnFRAG; pp++, ii++ )
        AddColumn( name, pp, freq, pssm, gaps, genfreq, genpssm, gengaps );
}

// -------------------------------------------------------------------------
// GenerateProfileSSE: generate a profile operating with SS elements
//
void ProGenerator::GenerateProfileSSE( const char* name, 
          FrequencyMatrix* genfreq, LogOddsMatrix* genpssm, GapScheme* gengaps  )
{
    if( genfreq == NULL || genpssm == NULL || gengaps == NULL )
        throw myruntime_error("ProGenerator: GenerateProfileSSE: Null target profile.");

    const LogOddsMatrix* pssm = NULL;//profiles providing SS chain

    if( GetNoProfiles() < 1 )
        return;

    const int   lenerr = 5;
    int         length = GetLengthToGenerate();
    size_t  nopros = GetNoProfiles();
    double  urv = rng_p_.GetDouble();
    size_t  ind = ( size_t )rint( urv *(double)(nopros-1));
    char    sse;

    if( nopros <= ind )
        throw myruntime_error("ProGenerator: GenerateProfileSSE: Invalid profile index generated.");

    pssm = GetProPSSMAt(ind);

    if( pssm == NULL )
        throw myruntime_error("ProGenerator: GenerateProfileSSE: Null starting profile.");

    int     lenpssm = pssm->GetColumns();
    double  urv2 = rng_pp_.GetDouble();
    int     pp = ( int )rint( urv2 *(double)(lenpssm-1));
    int     lc = 0, lp = 0, cnt = 0;

    if( pp < 0 || lenpssm <= pp )
        throw myruntime_error("ProGenerator: GenerateProfileSSE: Invalid profile position generated.");

    while(( lc = genpssm->GetColumns()) < length - lenerr ) {
        if( lc == lp )
            cnt++;
        else
            cnt = 0;
        if( 10 <= cnt )
            throw myruntime_error("ProGenerator: GenerateProfileSSE: Live lock in generating a profile.");
        lp = lc;
        sse = pssm->GetSSStateAt(pp);
        AddSSE( name, sse, genfreq, genpssm, gengaps );
        for( pp++; pp < pssm->GetColumns() && sse == pssm->GetSSStateAt(pp); pp++ );
        if( pssm->GetColumns() <= pp ) {
            ind++;
            if( nopros <= ind )
                ind = 0;
            pssm = GetProPSSMAt(ind);
            if( pssm == NULL )
                throw myruntime_error("ProGenerator: GenerateProfileSSE: Null profile.");
            pp = 0;
        }
    }
}

// -------------------------------------------------------------------------
// AddSSE: add SS element to profile being generated
//
void ProGenerator::AddSSE( const char* name, char sse, 
          FrequencyMatrix* genfreq, LogOddsMatrix* genpssm, GapScheme* gengaps )
{
    if( genfreq == NULL || genpssm == NULL || gengaps == NULL )
        throw myruntime_error("ProGenerator: AddSSE: Null arguments.");
    if( sse < 0 || SS_NSTATES <= sse )
        throw myruntime_error("ProGenerator: AddSSE: Invalid SS state.");

    const int   lcnCLEN = 10;//14[obs];//length of a long coil: 75% of real coils up to this value
    const int   lcnCFRAG = 6;//9[obs];//length of a coil fragment: real median length
    const FrequencyMatrix*  freq = NULL;
    const LogOddsMatrix*    pssm = NULL;//profile providing SS element
    const GapScheme*        gaps = NULL;

    if( GetNoProfiles() < 1 )
        throw myruntime_error("ProGenerator: AddSSE: No profiles.");

    const int   length = GetLengthToGenerate();
    size_t  nopros = GetNoProfiles();
    double  urv = rng_sse_p_.GetDouble();
    size_t  ind = ( size_t )rint( urv *(double)(nopros-1));

    if( nopros <= ind )
        throw myruntime_error("ProGenerator: AddSSE: Invalid profile index generated.");

    freq = GetProFREQAt(ind);
    pssm = GetProPSSMAt(ind);
    gaps = GetProGAPSAt(ind);

    if( freq == NULL || pssm == NULL || gaps == NULL )
        throw myruntime_error("ProGenerator: AddSSE: Null profile.");

    int     lenpssm = pssm->GetColumns();
    int     left = SLC_MIN( 40, lenpssm * 2 / 10 );
    int     rght = lenpssm - 1 - left;
    double  urv2 = rng_sse_pp_.GetDouble();
    int     pp = ( int )rint( urv2 *(double)(lenpssm-1));
    int     pe, clen = 0;
    double  expMIDs[PS_NSTATES];
    double  sssprob[SS_NSTATES];
    char    sssp;
    int     n, ii;

    if( rght < pp )
        pp = rght;
    if( pp < left )
        pp = left;

    if( pp < 0 || lenpssm <= pp )
        throw myruntime_error("ProGenerator: AddSSE: Invalid profile position generated.");

    if( !pssm->GetSSSSet())
        throw myruntime_error("ProGenerator: AddSSE: No SS prediction in a profile.");

    for(; pp <= rght && (sssp=pssm->GetSSStateAt(pp)) != sse; pp++ );
    if( rght < pp ) pp = left;
    for(; pp <= rght && (sssp=pssm->GetSSStateAt(pp)) != sse; pp++ );
    if( rght < pp )
        return;//cycled

    //get to the beginning of sse
    for(; 0 <= pp && (sssp=pssm->GetSSStateAt(pp)) == sse; pp-- );
    sssp = pssm->GetSSStateAt(++pp);
    if( sssp != sse )
        throw myruntime_error("ProGenerator: AddSSE: Unexpected SS state.");

    //if sse is a long coil, sample a fragment of it
    if( sse == SS_C ) {
        //learn the end position of sse
        for( pe = pp; pe < lenpssm && pssm->GetSSStateAt(pe) == sse; pe++ );
        if( pssm->GetSSStateAt(--pe) != sse )
            throw myruntime_error("ProGenerator: AddSSE: Unexpected SS end state.");
        clen = pe - pp + 1;
        if( lcnCLEN < clen ) {
            urv2 = rng_sse_pp_.GetDouble();
            pp += ( int )rint( urv2 *(double)(clen-1));
            if( lenpssm <= pp || ( sssp = pssm->GetSSStateAt(pp)) != sse )
                throw myruntime_error("ProGenerator: AddSSE: Positioning SSE failed.");
        }
    }

    //start filling SS element
    for( ii = 0; pp < lenpssm && pssm->GetSSStateAt(pp) == sse &&
       (( sse == SS_C )? genpssm->GetColumns() < length: 1 ); pp++, ii++ )
    {
        //sample just a fragment of a long coil
        if( sse == SS_C && lcnCLEN < clen && lcnCFRAG < ii )
            break;
        //
        //rewrite the whole SS element
        AddColumn( name, pp, freq, pssm, gaps, genfreq, genpssm, gengaps );
    }
}

// -------------------------------------------------------------------------
// AddColumn: add profile column to profile being generated
//
void ProGenerator::AddColumn( const char* name, int pp, 
                const FrequencyMatrix* freq, const LogOddsMatrix* pssm, const GapScheme* gaps, 
                FrequencyMatrix* genfreq, LogOddsMatrix* genpssm, GapScheme* gengaps )
{
    if( freq == NULL || pssm == NULL || gaps == NULL )
        throw myruntime_error("ProGenerator: AddColumn: Null source profile.");
    if( genfreq == NULL || genpssm == NULL || gengaps == NULL )
        throw myruntime_error("ProGenerator: AddColumn: Null target profile.");

    int     lenpssm = pssm->GetColumns();
    double  expMIDs[PS_NSTATES];
    double  sssprob[SS_NSTATES];
    int     n;

    if( pp < 0 || lenpssm <= pp )
        throw myruntime_error("ProGenerator: AddColumn: Invalid profile position.");

    //fill starting position
    if( gengaps->GetColumns() < 1 )
        gengaps->SetOrgTrProbsBeg( gaps->GetOrgTrProbsAt(pp-1));//-1 is ok
    if( genpssm->GetColumns() < 1 ) {
        expMIDs[PS_M] = pssm->GetMIDExpNoObservationsAt( pp-1, PS_M );
        expMIDs[PS_I] = pssm->GetMIDExpNoObservationsAt( pp-1, PS_I );
        expMIDs[PS_D] = pssm->GetMIDExpNoObservationsAt( pp-1, PS_D );
        genpssm->SetMIDExpNoObservationsBeg( expMIDs );
        //NOTE: bg probs
        for( n = 0; n < NUMALPH; n++ )
            genpssm->SetBackProbsAt( n, pssm->GetBackProbsAt(n));
        //
        genpssm->SetNoSequences( pssm->GetNoSequences());
        genpssm->SetEffNoSequences( pssm->GetEffNoSequences());
        genpssm->SetRefLambda( pssm->GetRefLambda());
        genpssm->SetRefK( pssm->GetRefK());
        genpssm->SetLambda( pssm->GetLambda());
        genpssm->SetEntropy( pssm->GetEntropy());
        genpssm->SetK( pssm->GetK());
        genpssm->SetExpectedScore( pssm->GetExpectedScore());
        genpssm->SetName((mystring(name)+slg_EXTN).c_str());
        genpssm->SetDescription((mystring(name)+" Generated profile").c_str());
        //
        genpssm->SetctPsSet( pssm->GetctPsSet());
        genpssm->SetCtxVecSet( pssm->GetCtxVecSet());
        genpssm->SetSSSSet( pssm->GetSSSSet());
        genpssm->SetSSSP3( pssm->GetSSSP3());
    }

    //copy profile column
    if( freq->GetVectorAt(pp) == NULL )
        throw myruntime_error("ProGenerator: AddColumn: Memory access error.");
    genfreq->Push( *freq->GetVectorAt(pp), freq->GetResidueAt(pp));
    ////
    n = genpssm->GetColumns();
    expMIDs[PS_M] = pssm->GetMIDExpNoObservationsAt( pp, PS_M );
    expMIDs[PS_I] = pssm->GetMIDExpNoObservationsAt( pp, PS_I );
    expMIDs[PS_D] = pssm->GetMIDExpNoObservationsAt( pp, PS_D );
    if( pssm->GetVectorAt(pp) == NULL )
        throw myruntime_error("ProGenerator: AddColumn: Memory access error: Null scores.");
    genpssm->Push( *pssm->GetVectorAt(pp), pssm->GetResidueAt(pp), 
                    pssm->GetFrequencyWeightAt(pp), pssm->GetInformationAt(pp), 
                    expMIDs );
    genpssm->SetBckPPProbAt( pssm->GetBckPPProbAt(pp), n );
    genpssm->SetPPProbsAt( n, pssm->GetPPProbsAt(pp), pssm->GetPPPIndsAt(pp), 
                          (int)pssm->GetNoPPProbsAt(pp));
    //
    if( genpssm->GetctPsSet() != pssm->GetctPsSet())
        throw myruntime_error("ProGenerator: AddColumn: Incompatible profiles (ct).");
    //genpssm->SetctPsSet( pssm->GetctPsSet());
    if( pssm->GetctPsSet()) {
        genpssm->SetctBckPPProbAt( pssm->GetctBckPPProbAt(pp), n );
        genpssm->SetctPPProbsAt( n, pssm->GetctPPProbsAt(pp), pssm->GetctPPPIndsAt(pp), 
                                (int)pssm->GetNoctPPProbsAt(pp));
    }
    //
    if( genpssm->GetCtxVecSet() != pssm->GetCtxVecSet())
        throw myruntime_error("ProGenerator: AddColumn: Incompatible profiles (CtxVec).");
    //genpssm->SetCtxVecSet( pssm->GetCtxVecSet());
    if( pssm->GetCtxVecSet()) {
        genpssm->SetCtxVecPlusAt( n, pssm->GetCtxVecNorm2At(pp), pssm->GetCtxVecLpprobAt(pp), 
                                  pssm->GetCtxVecAt(pp), pssm->GetCtxVecSize());
    }
    //
    if( genpssm->GetSSSSet() != pssm->GetSSSSet() || genpssm->GetSSSP3() != genpssm->GetSSSP3())
        throw myruntime_error("ProGenerator: AddColumn: Incompatible profiles (SS).");
    //genpssm->SetSSSSet( pssm->GetSSSSet());
    //genpssm->SetSSSP3( pssm->GetSSSP3());
    if( pssm->GetSSSSet()) {
        sssprob[SS_C] = pssm->GetSSStateProbAt( pp, SS_C );
        sssprob[SS_E] = pssm->GetSSStateProbAt( pp, SS_E );
        sssprob[SS_H] = pssm->GetSSStateProbAt( pp, SS_H );
        genpssm->SetSSStateAt( n, pssm->GetSSStateAt(pp), sssprob );
    }
    ////
    gengaps->Push( gaps->GetOrgTrProbsAt(pp), gaps->GetWeightsAt(pp), 
                    0.0/*delbeg*/, 0.0/*delend*/, 0/*delta*/, gaps->AAcid(pp));
}

// -------------------------------------------------------------------------
// RecalcBgProbs: recalculate profile background probabilities
// NOTE: consistent with 
//  InputMultipleAlignment::RecalcBackgroundProbabilities
//
void ProGenerator::RecalcBgProbs( 
          FrequencyMatrix* genfreq, LogOddsMatrix* genpssm, GapScheme* gengaps )
{
    if( genfreq == NULL || genpssm == NULL || gengaps == NULL )
        throw myruntime_error("ProGenerator: RecalcBgProbs: Null target profile.");

    const size_t    noress = NUMALPH;
    const size_t    noeffress = NUMAA;
    const double    effnos = genpssm->GetEffNoSequences();
    const double    maxtgtobss = 20.0;
    const double    tadj = 1.;//5.
    const double    bckpseudocounts = SLC_MAX(1., tadj * effnos );//10.0;
    const double    errtol = 1.0e-3;
    char            errbuf[KBYTE];

    double          bpp[noeffress], mixp[noeffress];
    double*         tgprobs = NULL;
    double          denm, consv;
    double          bckp, tgtp;
    double          expM, tgtpc;
    size_t          p, r;

    //Get target frequencies calculated
    genpssm->Finalize();

    for( r = 0; r < noeffress; r++ ) {
        bckp = LOSCORES.PROBABility( r );
        bpp[r] = bckpseudocounts * bckp;
        mixp[r] = 0.0;
    }

    consv = 0.0;
    for( p = 0; p < genpssm->GetColumns(); p++ ) 
    {
        expM = genpssm->GetMIDExpNoObservationsAt( p, PS_M );
        tgtpc = SLC_MIN( maxtgtobss, expM );
        denm = 1.0 /( bckpseudocounts + tgtpc );

        if( genpssm->GetTrgFreqsAt(p) == NULL )
            throw myruntime_error( "ProGenerator: RecalcBgProbs: "
                "Null profile target probabilities.");
        tgprobs = (double*)*genpssm->GetTrgFreqsAt(p);
        if( tgprobs == NULL )
            throw myruntime_error("ProGenerator: RecalcBgProbs: Null target probabilities.");

        for( r = 0; r < noeffress; r++ ) {
            tgtp = ( bpp[r] + tgtpc * tgprobs[r]) * denm;
            mixp[r] += tgtp;
            consv += tgtp;
        }
    }

    denm = 1.0;
    if( consv )
        denm /= consv;

    consv = 0.0;
    for( r = 0; r < noeffress; r++ ) {
        mixp[r] *= denm;
        consv += mixp[r];
        genpssm->SetBackProbsAt( r, mixp[r]);
    }
    for( ; r < noress; r++ )
        genpssm->SetBackProbsAt( r, 0.0 );
    if( consv < 1.0 - errtol || consv > 1.0 + errtol ) {
        sprintf( errbuf, "ProGenerator: RecalcBgProbs: Probabilities not conserved: %g", consv );
        throw myruntime_error( errbuf );
    }
}
