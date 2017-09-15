/***************************************************************************
 *   Copyright (C) 2012 by Mindaugas Margelevicius                         *
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
#include "HDPSmplGenerator.h"
#include "hdpsim.h"


int main( int argc, char *argv[] )
{
    //string values of the parameters
    mystring        patoutfile;         //pattern of output files
    mystring        dimensions;         //dimensions of vectors
    mystring        novariates;         //number of variates determining one cluster
    mystring        ctxtlength;         //context length
    mystring        noclusters;         //number of clusters
    mystring        nocsamples;         //number of samples per cluster
    bool            lgstnrmsampling = false;//sample from logistic-normal distribution
    bool            overlap = true;     //overlap of support variates
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
            {"o",       required_argument, 0, 'o'},
            {"d",       required_argument, 0, 'd'},
            {"a",       required_argument, 0, 'a'},
            {"c",       required_argument, 0, 'c'},
            {"k",       required_argument, 0, 'k'},
            {"n",       required_argument, 0, 'n'},
            {"p",       no_argument,       0, 'p'},
            {"l",       no_argument,       0, 'l'},
            {"v",       no_argument,       0, 'v'},
            {"h",       no_argument,       0, 'h'},
            { 0, 0, 0, 0 }
        };
        if(( c = getopt_long_only( argc, argv, "hvo:d:a:c:k:n:pl", 
                  long_options, &option_index )) == -1 )
            break;
#else
        if(( c = getopt( argc, argv, "hvo:d:a:c:k:n:pl" )) == -1 )
            break;
#endif

        switch( c ) {
            case ':':   fprintf( stdout, "No parameter value.\n%s",
                                                usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_FAILURE;
            case '?':   fprintf( stdout, "%s",  usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_FAILURE;
            case 'h':   fprintf( stdout, "%s",  usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_SUCCESS;

            case 'o':   patoutfile      = optarg;   break;
            case 'd':   dimensions      = optarg;   break;
            case 'a':   novariates      = optarg;   break;
            case 'c':   ctxtlength      = optarg;   break;
            case 'k':   noclusters      = optarg;   break;
            case 'n':   nocsamples      = optarg;   break;
            case 'p':   lgstnrmsampling = true;     break;
            case 'l':   overlap         = false;    break;
            case 'v':   suppress        = false;    break;
            default:    break;
        }
    }

    char*   p;
    int     intvalue;
    int     intdimensions = DEFAULT_DIM;
    int     intnovariates = DEFAULT_NUMBER_OF_VARS;
    int     intctxtlen = DEFAULT_CTXT_LENGTH;
    int     intnoclusters = DEFAULT_NUMBER_OF_CLUSTERS;
    int     intnocsamples = DEFAULT_NUMBER_OF_SAMPLES;
    double  tmpval;


    if( !dimensions.empty()) {
        intvalue = strtol( dimensions.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Invalid value of dimensions." );
            return EXIT_FAILURE;
        }
        if( intvalue < 1 || 1000 < intvalue ) {
            error( "Dimensions value is either non-positive or too large." );
            return EXIT_FAILURE;
        }
        intdimensions = intvalue;
    }
    if( !novariates.empty()) {
        intvalue = strtol( novariates.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Invalid number of variates." );
            return EXIT_FAILURE;
        }
        intnovariates = intvalue;
    }
    if( intnovariates < 1 || intdimensions < intnovariates ) {
        error( "Number of variates is either non-positive or greater than dimensions." );
        return EXIT_FAILURE;
    }
    if( !ctxtlength.empty()) {
        intvalue = strtol( ctxtlength.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Invalid context length." );
            return EXIT_FAILURE;
        }
        if( intvalue < 1 || 100 < intvalue ) {
            error( "Invalid value of context length." );
            return EXIT_FAILURE;
        }
        intctxtlen = intvalue;
    }
    if( !noclusters.empty()) {
        intvalue = strtol( noclusters.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Invalid number of clusters." );
            return EXIT_FAILURE;
        }
        if( intvalue < 1 ) {
            error( "Non-positive number of clusters." );
            return EXIT_FAILURE;
        }
        intnoclusters = intvalue;
    }
    if( !nocsamples.empty()) {
        intvalue = strtol( nocsamples.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Invalid number of samples per cluster." );
            return EXIT_FAILURE;
        }
        if( intvalue < 1 ) {
            error( "Non-positive number of samples per cluster." );
            return EXIT_FAILURE;
        }
        intnocsamples = intvalue;
    }
    if( intnoclusters < intnocsamples ) {
        error( "Number of samples should be upper-bounded by number of clusters." );
        return EXIT_FAILURE;
    }

    SetQuiet( suppress );


    HDPSmplGenerator*   genor = NULL;
    int                 retcode = EXIT_SUCCESS;

    try {
        genor = new HDPSmplGenerator( patoutfile.c_str());

        if( !genor ) {
            error( "Not enough memory." );
            return EXIT_FAILURE;
        }

        genor->SetDimensions( intdimensions );
        genor->SetNoVars( intnovariates );
        genor->SetContextLength( intctxtlen );
        genor->SetNoClusters( intnoclusters );
        genor->SetNoSamplesC( intnocsamples );
        genor->SetLGSTNSampling( lgstnrmsampling );
        genor->SetDoOverlap( overlap );
        genor->Run();

    } catch( myexception const& ex )
    {
        error( ex.what());
        retcode = EXIT_FAILURE;
    }

    if( genor ) {
        delete genor;
        genor = NULL; 
    }

    return retcode;
}

//////////////////////////////////////////////////////////////////////
