/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
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
#include "fsample.h"
#include "ProFreqWriter.h"


int main( int argc, char *argv[] )
{
    //string values of the parameters
    mystring        directory;          //directory to read profiles from
    mystring        outfile;            //output file of samples
    mystring        mineffthickn;       //minimum effective thickness
    mystring        nosamples;          //number of samples
    mystring        ctxtlength;         //context length
    mystring        step;               //step between adjacent context positions
    mystring        cwght;              //weight of central context position
    mystring        seed;               //seed for random number generator
    bool            mix = false;        //mix context positions
    bool            dishpergroup = false;//arrange dishes one per group
    bool            suppress = true;    //suppress warnings
    int             c;

    SetArguments( &argc, &argv );
    SetGlobalProgName( argv[0], version );

    if( argc <= 1 ) {
        fprintf( stdout, "%s", usage( argv[0], instructions, version, verdate ).c_str());
        return EXIT_SUCCESS;
    }

    while( 1 ) {
#ifdef USE_GETOPT_LONG
        int option_index = 0;
        static struct option long_options[] = {
            {"d",       required_argument, 0, 'd'},
            {"o",       required_argument, 0, 'o'},
            {"t",       required_argument, 0, 't'},
            {"n",       required_argument, 0, 'n'},
            {"c",       required_argument, 0, 'c'},
            {"p",       required_argument, 0, 'p'},
            {"m",       no_argument,       0, 'm'},
            {"w",       required_argument, 0, 'w'},
            {"a",       no_argument,       0, 'a'},
            {"s",       required_argument, 0, 's'},
            {"v",       no_argument,       0, 'v'},
            {"h",       no_argument,       0, 'h'},
            { 0, 0, 0, 0 }
        };
        if(( c = getopt_long_only( argc, argv, "hvd:o:t:n:c:p:mw:as:", 
                  long_options, &option_index )) == -1 )
            break;
#else
        if(( c = getopt( argc, argv, "hvd:o:t:n:c:p:mw:as:" )) == -1 )
            break;
#endif

        switch( c ) {
            case ':':   fprintf( stdout, "No parameter value.\n%s",
                                                usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_FAILURE;
            case '?':   fprintf( stdout, "%s",  usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_FAILURE;
            case 'h':   fprintf( stdout, "%s",  usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_SUCCESS;

            case 'd':   directory       = optarg;   break;
            case 'o':   outfile         = optarg;   break;
            case 't':   mineffthickn    = optarg;   break;
            case 'n':   nosamples       = optarg;   break;
            case 'c':   ctxtlength      = optarg;   break;
            case 'p':   step            = optarg;   break;
            case 'm':   mix             = true;     break;
            case 'w':   cwght           = optarg;   break;
            case 'a':   dishpergroup    = true;     break;
            case 's':   seed            = optarg;   break;
            case 'v':   suppress        = false;    break;
            default:    break;
        }
    }

    char*   p;
    int     intvalue;
    int     intmineffthickn = DEFAULT_MIN_EFFECTIVE_THICKNESS;
    int     intnosamples = DEFAULT_NUMBER_OF_SAMPLES;
    int     intctxtlen = DEFAULT_CTXT_LENGTH;
    int     intstep = DEFAULT_CTXT_STEP;
    double  dblcwght = DEFAULT_CTXT_CWGHT;
    int     seed_value = 0;//seed for random number generator
    double  tmpval;


    // -- sub. matrix --
    //TODO: save sub.matrix name in profile
    try {
        SetLOSCORES( "Gonnet", DEFAULT_PSCORES_FILENAME );
    } catch( myexception const& ex ) {
        error( ex.what());
        return EXIT_FAILURE;
    }

    if( !mineffthickn.empty()) {
        intvalue = strtol( mineffthickn.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Minimum effective number of sequences is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue <= 0 ) {
            error( "Minimum effective number of sequences is non-positive." );
            return EXIT_FAILURE;
        }
        intmineffthickn = intvalue;
    }
    if( !nosamples.empty()) {
        intvalue = strtol( nosamples.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Number of samples specified is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue < 1 ) {
            error( "Wrong number of samples specified." );
            return EXIT_FAILURE;
        }
        intnosamples = intvalue;
    }
    if( !ctxtlength.empty()) {
        intvalue = strtol( ctxtlength.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Context length specified is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue < 1 || 100 < intvalue ) {
            error( "Wrong context length specified." );
            return EXIT_FAILURE;
        }
        intctxtlen = intvalue;
    }
    if( !step.empty()) {
        intvalue = strtol( step.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Context step specified is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue < 1 || 10 < intvalue ) {
            error( "Wrong context step specified." );
            return EXIT_FAILURE;
        }
        intstep = intvalue;
    }
    if( !cwght.empty()) {
        tmpval = strtod( cwght.c_str(), &p );
        if( errno || *p ) {
            error( "Weight of central context position is invalid." );
            return EXIT_FAILURE;
        }
        if( tmpval <= 0.0 || 1.0 < tmpval ) {
            error( "Wrong weight of central context position." );
            return EXIT_FAILURE;
        }
        dblcwght = tmpval;
    }

    if( !seed.empty()) {
        intvalue = strtol( seed.c_str(), &p, 10 );
        if( errno || *p || intvalue < 0 ) {
            error( "Seed specified is invalid." );
            return EXIT_FAILURE;
        }
        seed_value = intvalue;
    }

    SetQuiet( suppress );


    ProFreqWriter*  frqwriter = NULL;
    int             retcode = EXIT_SUCCESS;

    try {
        //Sampling frequency vectors given profiles
        if( !directory.empty())
            frqwriter = new ProFreqWriter( outfile.c_str(), directory.c_str());
        else
            frqwriter = new ProFreqWriter( outfile.c_str(), argv + 1, argc - 1 );

        if( !frqwriter ) {
            error( "Not enough memory." );
            return EXIT_FAILURE;
        }

        frqwriter->SetMinEffThickness( intmineffthickn );
        frqwriter->SetNoSamples( intnosamples );
        frqwriter->SetContextSize( intctxtlen );
        frqwriter->SetContextStep( intstep );
        frqwriter->SetMix( mix );
        frqwriter->SetCentralWeight( dblcwght );
        frqwriter->SetDishPerGroup( dishpergroup );
        frqwriter->SetSeed( seed_value );
        frqwriter->Run();

    } catch( myexception const& ex )
    {
        error( ex.what());
        retcode = EXIT_FAILURE;
    }

    if( frqwriter ) {
        delete frqwriter;
        frqwriter = NULL; 
    }

    return retcode;
}

//////////////////////////////////////////////////////////////////////
