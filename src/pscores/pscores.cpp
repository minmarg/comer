/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <errno.h>
#include <getopt.h>
#include <stdlib.h>

#include "rc.h"
#include "libpro/srcaln/IMACountFiles.h"
#include "pscores.h"



int main( int argc, char *argv[] )
{
    //string values of the parameters
    mystring        output;
    mystring        directory;
    mystring        scale;
    bool            suppress = true;    //suppress output
    int             c;

    SetGlobalProgName( argv[0], version );

    if( argc <= 1 ) {
        fprintf( stdout, "%s", usage( argv[0], instructions, version, verdate ).c_str());
        return EXIT_SUCCESS;
    }

    while( 1 ) {
#ifdef USE_GETOPT_LONG
        int option_index = 0;
        static struct option long_options[] = {
            {"o",       required_argument, 0, 'o'},
            {"d",       required_argument, 0, 'd'},
            {"n",       required_argument, 0, 'n'},
            {"v",       no_argument,       0, 'v'},
            { 0, 0, 0, 0 }
        };
        if(( c = getopt_long_only( argc, argv, "hvo:d:n:", long_options, &option_index )) == -1 )
            break;
#else
        if(( c = getopt( argc, argv, "hvo:d:n:" )) == -1 )
            break;
#endif

        switch( c ) {
            case ':':   fprintf( stdout, "No argument value.\n%s",
                                                usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_FAILURE;
            case '?':   fprintf( stdout, "%s",  usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_FAILURE;
            case 'h':   fprintf( stdout, "%s",  usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_SUCCESS;

            case 'o':   output       = optarg;      break;
            case 'd':   directory    = optarg;      break;
            case 'n':   scale        = optarg;      break;

            case 'v':   suppress     = false;       break;
            default:    break;
        }
    }

    char*   p;
    double  tmpval;
    int     intscale;

    if( !scale.empty()) {
        intscale = strtol( scale.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Scale factor specified is invalid." );
            return EXIT_FAILURE;
        }
        if( intscale < 1 || intscale > 100 ) {
            error( "Scale factor is not within interval [1-100]." );
            return EXIT_FAILURE;
        }
    }

    SetQuiet( suppress );
    message( NULL );


    try {
        IMACountFiles*  files;

        if( !directory.empty())
            files = new IMACountFiles( output.c_str(), directory.c_str());
        else
            files = new IMACountFiles( output.c_str(), argv + 1, argc - 1 );

        if( !files ) {
            error( "Not enough memory." );
            return EXIT_FAILURE;
        }

        if( !scale.empty()) {
            IMACounts::SetScoreScale(( double )intscale );
            IMACounts::SetAutoScoreScale( false );
        }

        files->Make();
        delete files;

    } catch( myexception const& ex )
    {
        error( ex.what());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
