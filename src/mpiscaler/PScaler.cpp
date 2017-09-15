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
#include "libpro/srcpro/Serializer.h"
#include "PScaler.h"


////////////////////////////////////////////////////////////////////////////
// CLASS PScaler
//
// Constructor
//
PScaler::PScaler(
        const char* paramconfigfile,
        const char* database,
        const char* output,
        bool no_scaling,
        double infrm_threshold,
        double lambda,
        AbstractScoreMatrix::TScaling a_scaling,
        TMask c_masking,
        AbstractScoreMatrix::TType type )
:
    paramConfigFile( paramconfigfile ),
    database_name( database ),
    output_name( output ),
    scoreSystem( NULL ),

    not_to_scale( no_scaling ),
    infcontent( infrm_threshold ),
    lambda_spec( lambda ),
    profile_db( database, FrequencyStore::TFreqBSEnosStructure ),

    no_sequences( 0 ),
    db_size( 0 ),

    scaling( a_scaling ),
    masking( c_masking ),
    type_( type )
{
    SetFMinMsgSize( &AbstractUniversalScoreMatrix::GetMinSizeOfProbCalculator );
    SetFIsMsgValid( &AbstractUniversalScoreMatrix::IsFormattedDataValid );
}

// Default constructor
//
PScaler::PScaler()
:
    scoreSystem( NULL ),
    infcontent( 0.0 ),
    profile_db(( const char* )NULL )
{
    throw( myruntime_error(
            mystring( "PScaler: Default initialization is not allowed." )));
}

// Destructor
//
PScaler::~PScaler()
{
    DestroyScoreSystem();
}

// -------------------------------------------------------------------------
// Run: Begins searching of the query profile against the database
// -------------------------------------------------------------------------

long PScaler::Run()
{
    PrivateMPIInit();

    bool        needsending = false;            //no need for sending of message for a while
    static char s_errmsg[MSG_MAX];              //string buffer for error message
    int         s_szemsg = 0;                   //size of message
    int         szlimit = MSG_MAX;
    mystring    throwstr;

    sprintf( s_errmsg, "PScaler(%d): ", mMPIMyRank());
    szlimit -= strlen( s_errmsg );

    try {
        Configuration&  config_ungapped = GetConfiguration( ProcomUngapped );

        config_ungapped.SetFilename( GetParamConfigFile());
        //set parameters
        SetUngappedParams( config_ungapped );
        if( UsingLambdaSpecific()) {
            config_ungapped.SetLambda( GetLambdaSpecific());
        }

        //not using gapped configuration...
        GetConfiguration( ProcomGapped ) = GetConfiguration( ProcomUngapped );


        if( scoreSystem )
            throw myruntime_error( mystring( "PScaler: Score system has unexpectedly been initialized." ));

        //read database header
        profile_db.Open();

#ifdef  UNIVERSALSCORES
        MasterMessage( "Reading frequencies..." );
        //read the distinct frequency vectors
        profile_db.ReadInFrequencies();

        SetNoSequences( profile_db.GetNoSequences());
        SetDbSize( profile_db.GetDbSize());

        MasterMessage( "Scaling of scores..." );

        CreateScoreSystem();
        //from this point on the master should be notified on error occured in slave nodes
        needsending = true;
        ScaleScoreSystem();
#endif
        //no need here of sending notification message to the master if
        //error has occured in the slave nodes
        needsending = false;

        PostComputationLogic();

        profile_db.Close();

        MasterMessage( "Writing configuration..." );

        if( mMPIIAmMaster()) {
            config_ungapped.SetLambda(      scoreSystem->GetLambda());
            config_ungapped.SetK(           scoreSystem->GetK());
            config_ungapped.SetH(           scoreSystem->GetEntropy());
            config_ungapped.SetScaleFactor( scoreSystem->GetMultiplier());
            config_ungapped.WriteUngapped();
        }

        MasterMessage( "Finished." );

    } catch( myexception const& ex )
    {
        if( mMPIIAmMaster()) {
            throwstr = ex.what();
            //inform the slave processes if needed about the error occurred
            if( needsending )
                Terminate();
        } else {
            //inform the master if needed about the error occurred and exit silently
            if( needsending ) {
                s_szemsg = strlen( ex.what());
                if( szlimit <= s_szemsg )
                    s_szemsg = szlimit - 1;

                memcpy( s_errmsg + strlen( s_errmsg ), ex.what(), s_szemsg );
                s_errmsg[s_szemsg] = 0;
                SendMPIMessage( s_errmsg, s_szemsg, false/*throw_on_error*/ );
            }
        }
    }

    PrivateMPITerminate();

    if( !throwstr.empty())
        throw myruntime_error( throwstr );

    return 0;
}

// -------------------------------------------------------------------------
// PostComputationLogic: performs final output operations
// -------------------------------------------------------------------------

void PScaler::PostComputationLogic()
{
    if( ! mMPIIAmMaster())
        //nothing to do if not the master
        return;

    FILE*   fp = stdout;

    if( GetOutput() && strlen( GetOutput()))
        fp = fopen( GetOutput(), "w" );

    if( fp == NULL )
        throw myruntime_error( mystring( "PScaler: Failed to open file for writing." ));

    PrintHeader( fp );
    PrintFooter( fp );

    if( fp != stdout )
        fclose( fp );
}

// -------------------------------------------------------------------------
// PrintSearchingHeader: Prints general information related to the query and
//     database
// -------------------------------------------------------------------------

void PScaler::PrintHeader( FILE* fp )
{
    if( fp == NULL )
        return;

    const size_t    padding = OUTPUTINDENT;
    char            blanks[padding+1];
    size_t          n = 0;

    for( n = 0; n < padding; blanks[n++] = 32 );
    blanks[n] = 0;

    progname_and_version( fp );

    fprintf( fp, "\n" );
    fprintf( fp, " Database:\n" );
    fprintf( fp, "%s\n", my_basename( profile_db.GetDbName()));

    fprintf( fp, "%s%d profiles\n%s%llu total positions\n\n\n",
            blanks, GetNoSequences(),
            blanks, GetDbSize());
}

// -------------------------------------------------------------------------
// PrintSearchingFooter: Prints some general statistics related to profile
//     searching and computation of statistical significance
// -------------------------------------------------------------------------

void PScaler::PrintFooter( FILE* fp )
{
    if( fp == NULL )
        return;

    fprintf( fp, "\n\n" );
    PrintMethodName( fp );
    PrintParameterTable( fp );
    fprintf( fp, "\n" );
}

