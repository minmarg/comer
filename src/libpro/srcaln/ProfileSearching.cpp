/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <stdlib.h>
#include <string.h>

#include "mystring.h"
#include "myexcept.h"
#include "libseg/SEGSequence.h"
#include "libpro/srcpro/Serializer.h"
#include "libpro/srcaln/InputMultipleAlignment.h"
#include "ProfileSearching.h"


// -------------------------------------------------------------------------
// CLASS HitInformation
//
// Constructor
//
HitInformation::HitInformation(
        double sc,
        double eval,
        double refeval,
        char* annot,
        char* fullinfo,
        bool destroy )
:
    score( sc ),
    evalue( eval ),
    ref_evalue( refeval ),
    annotation( annot ),
    fullalignment( fullinfo ),
    autodestroy( destroy )
{
}

// Default constructor
//
HitInformation::HitInformation()
:   score( 0.0 ),
    evalue( -1.0 ),
    ref_evalue( -1.0 ),
    annotation( NULL ),
    fullalignment( NULL )
{
    throw( myruntime_error(
            mystring( "Default initialization of the HitInformation objects is prohibited." )));
}

// Destructor
//
HitInformation::~HitInformation()
{
    if( !autodestroy )
        return;

    if( annotation )
        free( annotation );
    if( fullalignment )
        free( fullalignment );
}

////////////////////////////////////////////////////////////////////////////
// CLASS ProfileSearching
//
// Constructor
//
ProfileSearching::ProfileSearching(
        const char* paramconfigfile,
        const char* input,
        const char* database,
        const char* output,
        double eval_threshold,
        int no_hits,
        int no_alns,
        double ident_threshold,
        double infrm_threshold,
        double mask_percnt,
        int gapopen,
        int gapextend,
        bool showpars,
        bool fixedcosts,
        AbstractScoreMatrix::TType  met,
        AbstractScoreMatrix::TBehaviour behaviour,
        AbstractScoreMatrix::TScaling precstrategy,
        TMask maskapp
    )
:
    paramConfigFile( paramconfigfile ),
    input_name( input ),
    database_name( database ),
    output_name( output ),
    scoreSystem( NULL ),

    max_evalue( eval_threshold ),
    max_no_hits( no_hits ),
    max_no_alns( no_alns ),
    seqidentity( ident_threshold ),
    information_content( infrm_threshold ),
    deletestates_( true ),

    thickness_number( 0 ),
    thickness_percnt( 0.0 ),

    maskscale_percnt( mask_percnt ),
    gapopencost( gapopen ),
    gapextncost( gapextend ),
    show_pars( showpars ),
    fixed_costs( fixedcosts ),

    hitListing( NULL ),
    capacity( 0 ),
    size( 0 ),

    profile_db( database ),

    no_sequences( 0 ),
    db_size( 0 ),

    method( met ),
    statbehaviour( behaviour ),
    scaling( precstrategy ),
    masking( maskapp ),

    tfrmix_( tfrmixNo ),
    scoadj_( scoadjNo ),
    HDPbase_( NULL ),
    HDPctbase_( NULL ),

    ssemodel_( 0 ),

    autogapcosts( false ),
    autocorrwinsize( -1 ),

    usingseg( false ),
    segwinlen( 0 ),
    seglowentropy( -1.0 ),
    seghighentropy( -1.0 ),
    segdistance( 0.0 ),

    usingseqseg( false ),
    seqsegwinlen( DEFAULT_SEGSEQ_WIN_LENGTH ),
    seqseglowent( DEFAULT_SEGSEQ_LOW_ENTROPY ),
    seqseghighent( DEFAULT_SEGSEQ_HIGH_ENTROPY ),

    deletioncoeff( 1.0 ),

    gapprobfactevalue( 1.0 ),
    gapprobfactweight( 0.0 ),
    gapprobfactshift( 0.0 ),

    acorrnumerator1st( 1.0 ),
    acorrnumerator2nd( 1.0 ),
    acorrlogscale( 1.0 ),
    acorrdenomscale( 1.0 ),

    ainfoupper2nd( 0.0 ),
    ainfonumerator2nd( 1.0 ),
    ainfoscale2nd( 1.0 ),
    ainfonumeratoralt( 1.0 ),
    ainfoscalealt( 1.0 ),

    hsplength_( 0 ),
    hspscore_( 0 ),
    hspdistance_( 0 ),
    hspnohsps_( 0 ),

    evalue_alnlength( -1.0 ),
    positcorrects( false ),
    autoaccorrection( true )
{
    Realloc( ALLOCHITS );

    if( GetGapOpenCost() == 0 ) {   // 0 means ungapped allignments
        SetGapOpenCost( Configuration::NOGAPVAL );
        SetGapExtnCost( Configuration::NOGAPVAL );
    }

    ProfileAlignment::SetInformationThreshold( infrm_threshold );
}

// Default constructor
//
ProfileSearching::ProfileSearching()
:
    scoreSystem( NULL ),
    information_content( 0.0 ),
    deletestates_( true ),
    thickness_number( 0 ),
    thickness_percnt( 0.0 ),
    maskscale_percnt( 0.0 ),
    profile_db(( const char* )NULL ),

    tfrmix_( tfrmixNo ),
    scoadj_( scoadjNo ),
    HDPbase_( NULL ),
    HDPctbase_( NULL ),

    ssemodel_( 0 ),

    autogapcosts( false ),
    autocorrwinsize( -1 ),

    usingseg( false ),
    segwinlen( 0 ),
    seglowentropy( -1.0 ),
    seghighentropy( -1.0 ),
    segdistance( 0.0 ),

    usingseqseg( false ),
    seqsegwinlen( DEFAULT_SEGSEQ_WIN_LENGTH ),
    seqseglowent( DEFAULT_SEGSEQ_LOW_ENTROPY ),
    seqseghighent( DEFAULT_SEGSEQ_HIGH_ENTROPY ),

    evalue_alnlength( -1.0 ),
    positcorrects( false ),
    autoaccorrection( true )
{
    throw( myruntime_error(
            mystring( "Default initialization of the ProfileSearching objects is prohibited." )));
}

// Destructor
//
ProfileSearching::~ProfileSearching()
{
    if( hitListing ) {
        for( size_t n = 0; n < size; n++ ) {
            delete hitListing[n];
        }
        free( hitListing );
        hitListing = NULL;
        capacity = 0;
        size = 0;
    }
    DestroyScoreSystem();
}

// -------------------------------------------------------------------------
// Realloc: reallocates memory fo hit list
//
void ProfileSearching::Realloc( size_t newcap )
{
    if( capacity == 0 ) {
        hitListing = ( HitInformation** )malloc( sizeof( void* ) * newcap );
    } else {
        hitListing = ( HitInformation** )realloc( hitListing, sizeof( void* ) * newcap );
    }

    if( !hitListing )
        throw myruntime_error( mystring( "ProfileSearching: Not enough memory." ));

    HitInformation**  listing = hitListing;

    if( capacity != 0 ) {
        listing = hitListing + capacity;
    }

    memset( listing, 0, sizeof( void* ) * ( newcap - capacity ));

    capacity = newcap;
}

// -------------------------------------------------------------------------
// Push: pushe hit information into the hit list
//
void ProfileSearching::Push( HitInformation* hit )
{
    if( capacity  <  size + 1 ) {
        Realloc( capacity * 2 );
    }

    hitListing[size] = hit;

    size++;
}

// -------------------------------------------------------------------------
// Run: start searching of the query profile against the database
//
void ProfileSearching::Run()
{
    FrequencyMatrix dbfreq;
    LogOddsMatrix   dbpssm;
    GapScheme       dbgaps;

    GetConfiguration( ProcomUngapped ).SetFilename( GetParamConfigFile());
    GetConfiguration( ProcomGapped ).SetFilename( GetParamConfigFile());

    //read parameters
    GetConfiguration( ProcomUngapped ).ReadUngapped();  //read ungapped configuration
    GetConfiguration( ProcomGapped ).SetAutoGapOpenCost( GetAutoGapCosts());
    GetConfiguration( ProcomGapped ).SetGapOpenCost( GetGapOpenCost());
    GetConfiguration( ProcomGapped ).SetGapExtendCost( GetGapExtnCost());
    GetConfiguration( ProcomGapped ).Read();            //read reference parameters

    char    strbuf[BUF_MAX];

    sprintf( strbuf, "Lambda (u/g), %8.4f/%8.4f",   GetConfiguration( ProcomUngapped ).GetLambda(),
                                                    GetConfiguration( ProcomGapped ).GetLambda());   message( strbuf, false );
    sprintf( strbuf, "K      (u/g), %8.4f/%8.4f",   GetConfiguration( ProcomUngapped ).GetK(),
                                                    GetConfiguration( ProcomGapped ).GetK());        message( strbuf, false );
    sprintf( strbuf, "H      (u/g), %8.4f/%8.4f",   GetConfiguration( ProcomUngapped ).GetH(),
                                                    GetConfiguration( ProcomGapped ).GetH());        message( strbuf, false );
    sprintf( strbuf, "Alpha  (u/g), %8.4f/%8.4f",   GetConfiguration( ProcomUngapped ).GetAlpha(),
                                                    GetConfiguration( ProcomGapped ).GetAlpha());    message( strbuf, false );
    sprintf( strbuf, "Beta   (u/g), %8.4f/%8.4f",   GetConfiguration( ProcomUngapped ).GetBeta(),
                                                    GetConfiguration( ProcomGapped ).GetBeta());     message( strbuf );

    if( scoreSystem )
        throw myruntime_error( mystring( "ProfileSearching: Score system unexpectedly initialized." ));

    query_gaps.SetGapProbabFactorEvalue( GetGapProbabFactorEvalue());
    query_gaps.SetGapProbabFactorWeight( GetGapProbabFactorWeight());
    query_gaps.SetGapProbabFactorShift( GetGapProbabFactorShift());

    query_gaps.SetAutocorrectionNumerator1st( GetAutocorrectionNumerator1st());
    query_gaps.SetAutocorrectionNumerator2nd( GetAutocorrectionNumerator2nd());
    query_gaps.SetAutocorrectionLogScale( GetAutocorrectionLogScale());
    query_gaps.SetAutocorrectionDenomScale( GetAutocorrectionDenomScale());
    query_gaps.SetConfiguration( GetConfiguration());


    dbgaps.SetGapProbabFactorEvalue( GetGapProbabFactorEvalue());
    dbgaps.SetGapProbabFactorWeight( GetGapProbabFactorWeight());
    dbgaps.SetGapProbabFactorShift( GetGapProbabFactorShift());

    dbgaps.SetAutocorrectionNumerator1st( GetAutocorrectionNumerator1st());
    dbgaps.SetAutocorrectionNumerator2nd( GetAutocorrectionNumerator2nd());
    dbgaps.SetAutocorrectionLogScale( GetAutocorrectionLogScale());
    dbgaps.SetAutocorrectionDenomScale( GetAutocorrectionDenomScale());
    dbgaps.SetConfiguration( GetConfiguration());

    ProcessQuery();

    try{
        //try read database header next
        profile_db.Open();

    } catch( myexception const& ex )
    {
        try{
            //try profile or multiple alignment in fasta
            FillMatrices( dbfreq, dbpssm, dbgaps, GetDatabase());
        } catch( myexception const& ex1 ) {
            throw myruntime_error(( mystring( ex.what()) + "\n" )+ ex1.what(), ex.eclass());
        }
        

        message( "Searching..." );

        SetNoSequences( 1 );
        SetDbSize( dbpssm.GetColumns());

        AbstractScoreMatrix::ComputeLengthAdjustment(
                        GetConfiguration( ProcomGapped ),
                        query_pssm.GetColumns(),
                        GetDbSize(),
                        GetNoSequences()
        );

        ComputationLogicWithProfiles( dbfreq, dbpssm, dbgaps );

        PostComputationLogic();
        message( "Finished." );
        return;
    }

#ifdef  UNIVERSALSCORES
    message( "Reading frequencies..." );
    //read the distinct frequency vectors first
    profile_db.ReadInFrequencies();

    message( "Computing scores..." );
    CreateScoreSystem();
#endif

    message( "Searching..." );

    SetNoSequences( profile_db.GetNoSequences());
    SetDbSize( profile_db.GetDbSize());

    AbstractScoreMatrix::ComputeLengthAdjustment(
                    GetConfiguration( ProcomGapped ),
                    query_pssm.GetColumns(),
                    GetDbSize(),
                    GetNoSequences()
    );

    //while not having reached the end of database
    while( profile_db.Next( dbfreq, dbpssm, dbgaps, GetGapOpenCost(), GetGapExtnCost(), GetFixedCosts())) {
        ComputationLogicWithProfiles( dbfreq, dbpssm, dbgaps );
    }

    PostComputationLogic();

    profile_db.Close();
    message( "Finished." );
}

// -------------------------------------------------------------------------
// ComputationLogicWithProfiles: performs computations with profiles; this
//     includes profile alignment, statistical significance calculations and
//     gathering of information for output
// -------------------------------------------------------------------------

void ProfileSearching::ComputationLogicWithProfiles( 
    FrequencyMatrix& freq, LogOddsMatrix& pssm, GapScheme& gaps )
{
//  gaps.OutputGapScheme();
//  freq.OutputMatrix();
//  pssm.OutputMatrix();

    double      adjscore = 0.0;
    double      expscore = 0.0;
    double      explength = 0.0;
    mystring    errstr;
    int         errcls;

    try{
        if( !PreprocessSubject( freq, pssm, gaps, true/*first call*/ ))
            return;
    } catch( myexception const& ex )
    {
        if( 1/*any*/ || ex.eclass() == SCALING ) {
            warning( ex.what());
            return;
        }
        throw myruntime_error( ex.what(), ex.eclass());
    }

    ProfileAlignment*   proaln = NULL;

    if( GetAlnAlgorithm() == ProfileAlignment::Optimizational )
        proaln = new ProfileAlignment(
            query_freq, query_pssm, query_gaps,
            freq, pssm, gaps, GetScoreSystem(), UngappedAlignments()
        );
    else if( GetAlnAlgorithm() == ProfileAlignment::Hybrid )
        proaln = new HybridAlignment(
            query_freq, query_pssm, query_gaps,
            freq, pssm, gaps, GetScoreSystem(), UngappedAlignments()
        );

    if( proaln == NULL )
        throw myruntime_error("ProfileSearching: Not enough memory.");

    try{
        proaln->SetModelSS18( GetSSEModel());
        proaln->SetProProbs( profile_db.GetProProbs());
        //run alignment algorithm...
        proaln->Run( profile_db.GetNoSequences());
        proaln->MARealign();

// explength = GetScoreSystem()->GetMeanLengthGivenExpectation(proaln->GetExpectation(),&expscore);
// int difflen = proaln->GetNoMatched() - (int)explength;
// if( 0 < difflen ) {
//     adjscore = proaln->GetScore() - (double)difflen * GetScoreSystem()->GetEntropy() * 0.1;
//     if( adjscore < 0.0 )
//         adjscore = 1.0;
//     proaln->AdjustScore( adjscore );
// }

        //{{2nd pass if needed
        if( GetAutoACcorrection() && GetScoreSystem()) {
            gaps.AdjustContextByEval(
//                     proaln.GetRawExpectation(),
                    proaln->GetExpectPerAlignment(),
                    GetScoreSystem()->GetEntropy(),
                    GetScoreSystem()->GetSbjctInfContent());
            query_gaps.AdjustContextByEval(
//                     proaln.GetRawExpectation(),
                    proaln->GetExpectPerAlignment(),
                    GetScoreSystem()->GetEntropy(),
                    GetScoreSystem()->GetQueryInfContent());
//             GetScoreSystem()->SetInfoThresholdByEval( proaln.GetRawExpectation());
            GetScoreSystem()->SetInfoThresholdByEval( proaln->GetExpectPerAlignment());
            PreprocessSubject( freq, pssm, gaps, false );
            proaln->Run( profile_db.GetNoSequences());
        }
        //}}

        //if alignment length is shorter than the expected mean length given
        //computed statistical parameters
        if(0)
        if( 0.0 < GetExpectForAlnLength() &&
            proaln->GetExpectation() < GetExpectForAlnLength())
            if( proaln->GetAlnLength() <
                GetScoreSystem()->GetMeanLengthGivenExpectation( GetExpectForAlnLength(), &expscore ))
                    if( expscore < proaln->GetScore())
                        if( expscore < proaln->GetScore() - expscore )
                            proaln->AdjustScore( proaln->GetScore() - expscore );
                        else
                            proaln->AdjustScore( expscore );


        //if expectation is above threshold used
        if( proaln->GetExpectation() <= GetEvalueThreshold() &&
            proaln->GetReferenceExpectation() <= GetEvalueThreshold() &&
            0 < proaln->GetScore())
        {
            size_t  title_size  = pssm.GetMinimumRequiredSizeForDescription();
            size_t  aln_size = proaln->GetMaxReqSizeForAlignment();
            size_t  total_size = title_size + aln_size;//size of full alignment information (desc.inc.)

            char*   annotatn = ( char* )malloc( pssm.GetMaxAnnotationWidth());
            char*   fullinfo = ( char* )malloc( total_size );

            if( !annotatn || !fullinfo )
                throw myruntime_error( mystring( "Not enough memory." ));

            pssm.PrintAnnotation( annotatn );
            proaln->Print( fullinfo, ToShowPars());

            HitInformation* hit =
                new HitInformation( proaln->GetScore(), proaln->GetExpectation(),
                            proaln->GetReferenceExpectation(), annotatn, fullinfo );

            Push( hit );
    //      proaln->OutputScoringMatrix(); //
        }
    } catch( myexception const& ex )
    {
        errstr = ex.what();
        errcls = ex.eclass();
    }

    if( proaln ) {
        delete proaln;
        proaln = NULL;
    }

    if( !errstr.empty())
        throw myruntime_error( errstr.c_str(), errcls );
}

// -------------------------------------------------------------------------
// PostComputationLogic: performs final output operations
// -------------------------------------------------------------------------

void ProfileSearching::PostComputationLogic()
{
    FILE*   fp = stdout;

    if( GetOutput() && strlen( GetOutput()))
        fp = fopen( GetOutput(), "w" );

    if( fp == NULL )
        throw myruntime_error( mystring( "Failed to open file for writing." ));

    QSortHits();

    PrintSearchingHeader( fp );
    PrintHits( fp );
    PrintSearchingFooter( fp );

    if( fp != stdout )
        fclose( fp );
}

// -------------------------------------------------------------------------
// PrintSearchingHeader: Prints general information related to the query and
//     database
// -------------------------------------------------------------------------

void ProfileSearching::PrintSearchingHeader( FILE* fp )
{
    if( fp == NULL )
        return;

    const size_t    padding = OUTPUTINDENT;
    char            blanks[padding+1];
    size_t          n = 0;

    for( n = 0; n < padding; blanks[n++] = 32 );
    blanks[n] = 0;

    progname_and_version( fp );

    if( GetFixedCosts()) {
        if( GetAutoGapCosts())
            fprintf( fp, " -positional gap costs unaffected by probabilities-\n" );
        else
            fprintf( fp, " -positionally invariable gap costs-\n" );
    }

    fprintf( fp, "\n" );
    fprintf( fp, " Query (%d positions):\n", query_pssm.GetColumns());

    query_pssm.PrintDescriptionFirst( fp );

    fprintf( fp, "\n\n" );
    fprintf( fp, " Database:\n" );
    fprintf( fp, "%s\n", my_basename( profile_db.GetDbName()));

    
    fprintf( fp, "%s%d profiles\n%s%llu total positions\n\n\n",
            blanks, GetNoSequences(),
            blanks, GetDbSize());
}

// -------------------------------------------------------------------------
// PrintHits: Prints all hits were found while searching query against the
//     database
// -------------------------------------------------------------------------

void ProfileSearching::PrintHits( FILE* fp )
{
    if( fp == NULL )
        return;

    size_t      width = OUTPUTINDENT + OUTPUTWIDTH;
    const char* title = " Profiles found below the e-value threshold:";
    const char* nohitstitle = " No profiles found below the e-value threshold";
    size_t      padding = width - strlen( title );

    if( !GetHitlistSize()) {
        fprintf( fp, "%s of %g.\n\n\n", nohitstitle, GetEvalueThreshold());
        return;
    }

    fprintf( fp, "%s", title );
    for( size_t n = 0; n < padding; n++ ) fprintf( fp, " " );
    fprintf( fp, " %7s %7s\n\n", "Score", "E-value" );

    for( size_t h = 0; h < GetHitlistSize() && ( int )h < GetNoHitsThreshold(); h++ )
    {
        const HitInformation*   hit = GetHitAt( h );
        //annotation will never be null
        size_t                  alength = strlen( hit->GetAnnotation());

        fprintf( fp, "%s", hit->GetAnnotation());
        for( size_t n = 0; n < width - alength; n++ ) fprintf( fp, " " );

        if( hit->GetEvalue() < 0.0 )
//             fprintf( fp, " %7d %7s\n", ( int )hit->GetScore(), "n/a" );
            fprintf( fp, " %7d %7.1g\n", ( int )hit->GetScore(), hit->GetRefEvalue());
        else
            fprintf( fp, " %7d %7.1g\n", ( int )hit->GetScore(), hit->GetEvalue());
    }
    fprintf( fp, "\n\n" );


    for( size_t h = 0; h < GetHitlistSize() && ( int )h < GetNoAlnsThreshold(); h++ )
    {
        const HitInformation*   hit = GetHitAt( h );
        fprintf( fp, "%s\n", hit->GetFullAlignment());
    }
}

// -------------------------------------------------------------------------
// PrintSearchingFooter: Prints some general statistics related to profile
//     searching and computation of statistical significance
// -------------------------------------------------------------------------

void ProfileSearching::PrintSearchingFooter( FILE* fp )
{
    if( fp == NULL )
        return;

    size_t  length = ScoringMatrix::GetDeltaLength();       //length adjustment
    Uint64  sspace = ScoringMatrix::GetSearchSpace();       //effective search space

    int     eff_query_length = query_pssm.GetColumns() - length;
    Int64   eff_db_length = GetDbSize() - GetNoSequences() * length;

    fprintf( fp, "\n\n" );
    PrintMethodName( fp );
    if( UngappedAlignments())
        fprintf( fp, "No gap costs used: Ungapped alignments\n\n" );
    else if(0) {
        if( GetAutoGapCosts()) {
            fprintf( fp, "Gap open cost, Computed (window size, %d)\n", GetAutocorrWinsize());
            fprintf( fp, "Initial gap extension cost, %d\n", abs( GetGapExtnCost()));
        } else {
            fprintf( fp, "Gap open cost, %d\n", abs( GetGapOpenCost()));
            fprintf( fp, "Gap extension cost, %d\n", abs( GetGapExtnCost()));
        }
        if( !GetFixedCosts())
            fprintf( fp, "Deletion probability weight, %.2f\n", GetDeletionCoefficient());
        fprintf( fp, "\n" );
    }
    PrintParameterTable( fp );
    fprintf( fp, "\n" );
    fprintf( fp, "Length of query, %d\n", query_pssm.GetColumns());
    fprintf( fp, "Length of database, %llu\n", GetDbSize());
    fprintf( fp, "Effective length of query, %d\n", eff_query_length );
    fprintf( fp, "Effective length of database, %lld\n", eff_db_length );
    fprintf( fp, "Effective search space, %llu\n\n", sspace );
}

// -------------------------------------------------------------------------
// ProcessQuery: Reads query profile from file and processes it
// -------------------------------------------------------------------------

void ProfileSearching::ProcessQuery()
{
    FillMatrices( query_freq, query_pssm, query_gaps, GetInput());
    if( GetTarFrMixHDPCtx())
        query_pssm.MixTrgFrequencies( GetHDPbase());

    if( GetUsingSeg()) {
        //SEG logic
        SEGProfile  segpro(
                query_freq,
                query_pssm,
                GetSegWinLength(),
                GetSegLowEntropy(),
                GetSegHighEntropy()
        );
        segpro.SetDistance( GetSegDistance());
        segpro.Run();
        segpro.MaskSeggedPositions( query_freq, query_pssm, query_gaps );
    }
}

// -------------------------------------------------------------------------
// FillMatrices: fills matrices describing profile; this includes frequency
//     matrix, pssm, position-specific gap costs
// -------------------------------------------------------------------------
#include "libpro/srcsco/ProfileMatrix.h"

void ProfileSearching::FillMatrices(
        FrequencyMatrix&    freq,
        LogOddsMatrix&      pssm,
        GapScheme&          gaps,
        const char*         filename )
{
    Serializer  serializer;

    try {
        //read profile
        serializer.ReadProfile( filename, freq, pssm, gaps );

    } catch( myexception const& ex )
    {
        //make profile
        InputMultipleAlignment* inmaln = new InputMultipleAlignment;

        if( !inmaln )
            throw( myruntime_error( "Not enough memory." ));

        if( GetHCSeg()) {
            inmaln->SetSEGWindow( GetHCWinLength());
            inmaln->SetSEGLowEntropy( GetHCLowEntropy());
            inmaln->SetSEGHighEntropy( GetHCHighEntropy());
        }
        else
            inmaln->SetUsingSEGFilter( false );

        if( GetUsingSeqSeg()) {
            inmaln->SetSeqSEGWindow( GetSeqSegWinLength());
            inmaln->SetSeqSEGLowEntropy( GetSeqSegLowEntropy());
            inmaln->SetSeqSEGHighEntropy( GetSeqSegHighEntropy());
        }
        else
            inmaln->SetUsingSeqSEGFilter( false );

        inmaln->SetIdentityLevel( GetIdentityThreshold());
        inmaln->SetComputeDELETEstates( GetComputeDELETEstates());

        inmaln->SetExtentMinWindow( GetExtentMinWindow());
        inmaln->SetExtentMinSeqPercentage( GetExtentMinSeqPercentage());
        inmaln->SetPseudoCountWeight( GetPseudoCountWeight());

        inmaln->ReadAlignment( filename );
        inmaln->ConstructProfile();

        inmaln->ExportFrequenciesTo( freq );
        inmaln->ExportPSSMTo( pssm );
        inmaln->ExportGapWeightsTo( gaps );

        delete inmaln;

        // SEG logic
        if( GetUsingSeg()) {
            SEGProfile  segpro(
                    freq,
                    pssm,
                    GetSegWinLength(),
                    GetSegLowEntropy(),
                    GetSegHighEntropy()
            );
            segpro.SetDistance( GetSegDistance());
            segpro.Run();
            segpro.MaskSeggedPositions( freq, pssm, gaps );
        }
        // --
//         message( "Data read." );
    }

    freq.CheckForAllZeros();    //necessary to ensure valid scores

//     //set gap open and extension costs as given by the program arguments
//     gaps.Prepare( GetGapOpenCost(), GetGapExtnCost(), GetFixedCosts());
}

// -------------------------------------------------------------------------
// qsort: this quick sort implementation is from 'Numerical recipes in C'
// Avoiding recursion makes the implementation much faster.
// QS_MIN is the size of subarrays sorted by straight insertion.
// QNSTACK is the required auxiliary storage whose size should be at least
// 2 log2 N
// -------------------------------------------------------------------------

#define SWAP( a, b ) temp = ( a ); ( a ) = ( b ); ( b ) = temp;
#define QS_MIN 7
#define QNSTACK 100

void ProfileSearching::QSortHits()
{
    size_t  n = GetHitlistSize();

    if( !n )
        return;

    int     lpos = 0;                           //left position
    int     rpos = ( int )n - 1;                //right position
    int     i, j, k;
    int     s = 0;                              //stack iterator
    int     stack[QNSTACK+1];                   //position stack
    HitInformation* hit, *temp;

    for(;;) {
        if( 1|| rpos - lpos < QS_MIN ) {                    //insertion sort when subarray small enough.
            for( j = lpos + 1; j <= rpos; j++ ) {
                hit = GetHitAt( j );
                for( i = j - 1; i >= lpos; i-- ) {
                    if( !( *hit < *GetHitAt( i )))
                        break;
                    GetHitAt( i+1 ) = GetHitAt( i );
                }
                GetHitAt( i+1 ) = hit;
            }
            if( s == 0 ) break;
            rpos = stack[s--];                          //pop stack and begin a new round of partitioning
            lpos = stack[s--];

        } else {
            k = ( lpos + rpos ) >> 1;       //choose median of left, center, and right elements as partitioning element hit
            SWAP( GetHitAt( k ), GetHitAt( lpos+1 ));   //also rearrange so that hit[lpos] <= hit[lpos+1] <= hit[rpos]

            if( *GetHitAt( rpos ) < *GetHitAt( lpos )) {
                SWAP( GetHitAt( lpos ), GetHitAt( rpos ));
            }
            if( *GetHitAt( rpos ) < *GetHitAt( lpos+1 )) {
                SWAP( GetHitAt( lpos+1 ), GetHitAt( rpos ));
            }
            if( *GetHitAt( lpos+1 ) < *GetHitAt( lpos )) {
                SWAP( GetHitAt( lpos ), GetHitAt( lpos+1 ));
            }

            i = lpos + 1;                           //initialize pointers for partitioning
            j = rpos;
            hit = GetHitAt( lpos+1 );               //partitioning element

            for(;;) {
                do { i++; } while ( *GetHitAt( i ) < *hit );//scan up to find element > hit
                do { j--; } while ( *hit < *GetHitAt( j )); //scan down to find element < hit
                if( j < i ) break;                          //pointers crossed; partitioning complete
                SWAP( GetHitAt( i ), GetHitAt( j ));
            }

            GetHitAt( lpos+1 ) = GetHitAt( j );             //insert partitioning element
            GetHitAt( j ) = hit;
            s += 2;

            if( s > QNSTACK ) {
                error( "Sort of e-values failed: Pre-allocated stack size is too small." );
                break;
            }

            //push pointers to larger subarray on stack, process smaller subarray immediately
            if( rpos - i + 1 >= j - lpos ) {
                stack[s] = rpos;
                stack[s-1] = i;
                rpos = j - 1;
            } else {
                stack[s] = j - 1;
                stack[s-1] = lpos;
                lpos = i;
            }
        }
    }
}
