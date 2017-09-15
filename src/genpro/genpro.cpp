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
#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/datapro.h"
#include "ProGenerator.h"
#include "genpro.h"



int main( int argc, char *argv[] )
{
    //string values of the parameters
    mystring        outdir;
    mystring        outfnamepat;
    mystring        directory;
    mystring        lenpro;
    mystring        nopros;
    mystring        fraglen;//length of indivisible profile fragment
    bool            opsse = false;//operate with SS elements
    bool            suppress = true;//suppress output
    int             c;

    SetArguments( &argc, &argv );
    SetGlobalProgName( argv[0], version );

    if( argc <= 1 ) {
        fprintf( stdout, "%s", usage( argv[0], makeinst, version, verdate ).c_str());
        return EXIT_SUCCESS;
    }

    while( 1 ) {
#ifdef USE_GETOPT_LONG
        int option_index = 0;
        static struct option long_options[] = {
            {"o",       required_argument, 0, 'o'},
            {"n",       required_argument, 0, 'n'},
            {"d",       required_argument, 0, 'd'},
            {"l",       required_argument, 0, 'l'},
            {"N",       required_argument, 0, 'N'},
            {"f",       required_argument, 0, 'f'},
            {"s",       no_argument,       0, 's'},
            {"v",       no_argument,       0, 'v'},
            {"h",       no_argument,       0, 'h'},
            { 0, 0, 0, 0 }
        };
        if(( c = getopt_long_only( argc, argv, "hvo:n:d:l:N:f:s", long_options, &option_index )) == -1 )
            break;
#else
        if(( c = getopt( argc, argv, "hvo:n:d:l:N:f:s" )) == -1 )
            break;
#endif

        switch( c ) {
            case ':':   fprintf( stdout, "No parameter value.\n%s",
                                                usage( argv[0], makeinst, version, verdate ).c_str());  return EXIT_FAILURE;
            case '?':   fprintf( stdout, "%s",  usage( argv[0], makeinst, version, verdate ).c_str());  return EXIT_FAILURE;
            case 'h':   fprintf( stdout, "%s",  usage( argv[0], makeinst, version, verdate ).c_str());  return EXIT_SUCCESS;

            case 'o':   outdir       = optarg;      break;
            case 'n':   outfnamepat  = optarg;      break;
            case 'd':   directory    = optarg;      break;

            case 'l':   lenpro       = optarg;      break;
            case 'N':   nopros       = optarg;      break;
            case 'f':   fraglen      = optarg;      break;
            case 's':   opsse        = true;        break;

            case 'v':   suppress     = false;       break;
            default:    break;
        }
    }

    char*   p;
    double  tmpval;
    int     intvalue;
    char    strbuf[BUF_MAX];

    mystring    valSUBMAT = defSUBMAT;
    int         intlenpro = 0;
    int         intnopros = DEFAULT_NO_PROFILES;
    int         intfraglen = DEFAULT_FRAGLEN;

    // -- sub. matrix --
    try {
        SetLOSCORES( valSUBMAT, DEFAULT_PSCORES_FILENAME );
    } catch( myexception const& ex ) {
        error( ex.what());
        return EXIT_FAILURE;
    }

    if( outdir.empty()) {
        error( "Output directory name is not specified." );
        return EXIT_FAILURE;
    }

    // -- options --
    if( lenpro.empty()) {
        error( "Profile length should be provided." );
        return EXIT_FAILURE;
    }
    intvalue = strtol( lenpro.c_str(), &p, 10 );
    if( errno || *p ) {
        error( "Profile length is invalid." );
        return EXIT_FAILURE;
    }
    if( intvalue < 1 ) {
        error( "Invalid profile length." );
        return EXIT_FAILURE;
    }
    intlenpro = intvalue;


    if( !nopros.empty()) {
        intvalue = strtol( nopros.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Number of profiles is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue < 1 ) {
            error( "Invalid number of profiles." );
            return EXIT_FAILURE;
        }
        intnopros = intvalue;
    }


    if( !fraglen.empty()) {
        intvalue = strtol( fraglen.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Length of indivisible profile fragment is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue < 1 ) {
            error( "Invalid length of indivisible profile fragment." );
            return EXIT_FAILURE;
        }
        intfraglen = intvalue;
    }


    // --

    SetQuiet( suppress );
    message( NULL );

    sprintf( strbuf, "GEN: %d profile(s) of length %d", intnopros, intlenpro );
    message( strbuf );

    try {
        ProGenerator*   progen = NULL;

        if( !directory.empty())
            progen = new ProGenerator( outdir.c_str(), outfnamepat.c_str(), directory.c_str());
        else
            progen = new ProGenerator( outdir.c_str(), outfnamepat.c_str(), argv+1, argc-1 );

        if( progen == NULL ) {
            error( "Not enough memory." );
            return EXIT_FAILURE;
        }

        progen->SetLengthToGenerate( intlenpro );
        progen->SetNumberToGenerate( intnopros );
        progen->SetFragLength( intfraglen );
        progen->SetOperateWithSSE( opsse );
        progen->Generate();
        delete progen;

    } catch( myexception const& ex )
    {
        error( ex.what());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
