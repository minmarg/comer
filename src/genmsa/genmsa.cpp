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
#include "MSAGenerator.h"
#include "genmsa.h"



int main( int argc, char *argv[] )
{
    //string values of the parameters
    mystring        outdir;
    mystring        outfnamepat;
    mystring        directory;
    mystring        lenmsa;
    mystring        nomsas;
    mystring        fraglen;//length of indivisible profile fragment
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
            {"v",       no_argument,       0, 'v'},
            {"h",       no_argument,       0, 'h'},
            { 0, 0, 0, 0 }
        };
        if(( c = getopt_long_only( argc, argv, "hvo:n:d:l:N:f:", long_options, &option_index )) == -1 )
            break;
#else
        if(( c = getopt( argc, argv, "hvo:n:d:l:N:f:" )) == -1 )
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

            case 'l':   lenmsa       = optarg;      break;
            case 'N':   nomsas       = optarg;      break;
            case 'f':   fraglen      = optarg;      break;

            case 'v':   suppress     = false;       break;
            default:    break;
        }
    }

    char*   p;
    double  tmpval;
    int     intvalue;
    char    strbuf[BUF_MAX];

    mystring    valSUBMAT = defSUBMAT;
    int         intlenmsa = 0;
    int         intnomsas = DEFAULT_NO_MSAS;
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
    if( lenmsa.empty()) {
        error( "Alignment length should be provided." );
        return EXIT_FAILURE;
    }
    intvalue = strtol( lenmsa.c_str(), &p, 10 );
    if( errno || *p ) {
        error( "Alignment length is invalid." );
        return EXIT_FAILURE;
    }
    if( intvalue < 1 ) {
        error( "Invalid alignment length." );
        return EXIT_FAILURE;
    }
    intlenmsa = intvalue;


    if( !nomsas.empty()) {
        intvalue = strtol( nomsas.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Number of alignments is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue < 1 ) {
            error( "Invalid number of alignments." );
            return EXIT_FAILURE;
        }
        intnomsas = intvalue;
    }


    if( !fraglen.empty()) {
        intvalue = strtol( fraglen.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Length of indivisible alignment fragment is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue < 1 ) {
            error( "Invalid length of indivisible alignment fragment." );
            return EXIT_FAILURE;
        }
        intfraglen = intvalue;
    }


    // --

    SetQuiet( suppress );
    message( NULL );

    sprintf( strbuf, "GEN: %d MSA(s) of length %d", intnomsas, intlenmsa );
    message( strbuf );

    try {
        MSAGenerator* msagen = NULL;

        if( !directory.empty())
            msagen = new MSAGenerator( outdir.c_str(), outfnamepat.c_str(), directory.c_str());
        else
            msagen = new MSAGenerator( outdir.c_str(), outfnamepat.c_str(), argv+1, argc-1 );

        if( msagen == NULL ) {
            error( "Not enough memory." );
            return EXIT_FAILURE;
        }

        msagen->SetLengthToGenerate( intlenmsa );
        msagen->SetNumberToGenerate( intnomsas );
        msagen->SetFragLength( intfraglen );
        msagen->Generate();
        delete msagen;

    } catch( myexception const& ex )
    {
        error( ex.what());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
