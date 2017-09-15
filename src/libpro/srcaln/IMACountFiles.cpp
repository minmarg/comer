/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdlib.h>
#include <string.h>


#include "IMACountFiles.h"


// -------------------------------------------------------------------------
// CLASS IMACountFiles
//
// Constructor
//
IMACountFiles::IMACountFiles( const char* name )
:
    outname_( name ),
    data_dir_( NULL ),
    files_( NULL ),
    no_files_( 0 ),
    no_processedfiles_( 0 )
{
}

// Constructor
//
IMACountFiles::IMACountFiles( const char* output, const char* directory )
:
    outname_( output ),
    data_dir_( directory ),
    files_( NULL ),
    no_files_( 0 ),
    no_processedfiles_( 0 )
{
}

// Constructor
//
IMACountFiles::IMACountFiles( const char* output, char* arguments[], int no_args )
:
    outname_( output ),
    data_dir_( NULL ),
    files_( arguments ),
    no_files_( no_args ),
    no_processedfiles_( 0 )
{
}

// Default constructor
//
IMACountFiles::IMACountFiles()
:
    outname_( NULL ),
    data_dir_( NULL ),
    files_( NULL ),
    no_files_( 0 ),
    no_processedfiles_( 0 )
{
    throw( myruntime_error( mystring( "IMACountFiles: Default construction is not allowed." )));
}

// Destructor
//
IMACountFiles::~IMACountFiles()
{
}

// -------------------------------------------------------------------------
// Make: creates the database by processing profiles
// -------------------------------------------------------------------------

void IMACountFiles::Make()
{
    char    strbuf[BUF_MAX];

    ResetNoProcessedFiles();

    //process all files
    ProcessInput( &IMACountFiles::ProcessFile, NULL );

    GetOvrCounts().ComputeScores();
    try {
        GetOvrCounts().ComputeParameters();
    } catch( myexception const& ex ) {
        error( ex.what());
    }

    OutputScores();

    sprintf( strbuf, "Finished: Processed %d file(s).", GetNoProcessedFiles());
    message( strbuf );
}

// -------------------------------------------------------------------------
// Make: creates the database by processing profiles
// -------------------------------------------------------------------------

void IMACountFiles::OutputScores()
{
    FILE*   fp = stdout;

    if( GetOutputName() && strlen( GetOutputName()))
        fp = fopen( GetOutputName(), "w" );

    if( fp == NULL ) {
        error( mystring( mystring( "IMAClusters: Failed to open file '" ) + GetOutputName() + mystring( "'." )).c_str());
        fp = stdout;
    }

    fprintf( fp, "## Processed %d file(s)", GetNoProcessedFiles());
    if( GetDataDir())
        fprintf( fp, " from directory: %s", GetDataDir()? my_basename( GetDataDir()): "" );
    fprintf( fp, "\n" );

    GetOvrCounts().WriteScores( fp );

    if( fp != stdout )
        fclose( fp );
}

// -------------------------------------------------------------------------
// ProcessInput: process profiles given explicitly or found in directory
// -------------------------------------------------------------------------

void IMACountFiles::ProcessInput( PMETHOD processing_method, void* args )
{
    if( GetDataDir()) {
        struct stat     info;
        dirent*         entry = NULL;
        DIR*            direct = opendir( GetDataDir());
        mystring        dirname = GetDataDir() + mystring( DIRSEP );

        if( !direct )
            throw myruntime_error( mystring( "IMACountFiles: Unable to process directory specified." ));

        while(( entry = readdir( direct ))) {

            mystring    fullname = dirname + entry->d_name;

            if( stat( fullname.c_str(), &info ) == -1 ) {
                error( mystring( mystring( "Cannot stat '" ) + mystring( entry->d_name ) + "'; Skipping." ).c_str());
                continue;
            }
            if(( S_IFMT & info.st_mode ) != S_IFREG )
//             if( !S_ISREG( info.st_mode ))
                continue;

            ( this->*processing_method )( fullname.c_str()/*entry->d_name*/, args );
        }
        closedir( direct );

    } else
        if( GetFiles()) {
            bool    option = false;

            for( int n = 0; n < GetNoFiles(); n++ ) {
                if( option ) {
                    option = false;
                    continue;
                }
                if( GetFileAt( n )[0] == '-' && GetFileAt( n )[2] == 0 ) {
                    option = true;  //assume it is an option
                    continue;
                }
                ( this->*processing_method )( GetFileAt( n ), args );
            }
        }
}

// -------------------------------------------------------------------------
// ProcessFile: extracts needed information from file and processes it
// -------------------------------------------------------------------------

void IMACountFiles::ProcessFile( const char* filename, void* )
{
    try {
        IMACounts   observation;
        observation.Read( filename );

        GetOvrCounts().SumUp( observation );

        IncNoProcessedFiles();

    } catch( myexception const& ex )
    {
        error( mystring( ex.what() + mystring( " Skipping '" ) + filename + mystring( "'." )).c_str());
    }
}
