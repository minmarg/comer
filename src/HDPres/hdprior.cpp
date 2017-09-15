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
#include "HDPCalcPrior.h"
#include "hdprior.h"


int main( int argc, char *argv[] )
{
    //string values of the parameters
    mystring        unscafac;           //scale factor in uninformative configuration
    mystring        outfile;            //name of output file
    mystring        inpfile;            //input file of samples
    mystring        adjdegf;            //adjustment to the degrees of freedom
    bool            fileprior = false;  //calculate prior parameters from file
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
            {"i",       required_argument, 0, 'i'},
            {"o",       required_argument, 0, 'o'},
            {"n",       required_argument, 0, 'n'},
            {"s",       required_argument, 0, 's'},
            {"p",       no_argument,       0, 'p'},
            {"v",       no_argument,       0, 'v'},
            {"h",       no_argument,       0, 'h'},
            { 0, 0, 0, 0 }
        };
        if(( c = getopt_long_only( argc, argv, "hvi:o:n:s:p", 
                  long_options, &option_index )) == -1 )
            break;
#else
        if(( c = getopt( argc, argv, "hvi:o:n:s:p" )) == -1 )
            break;
#endif

        switch( c ) {
            case ':':   fprintf( stdout, "No parameter value.\n%s",
                                                usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_FAILURE;
            case '?':   fprintf( stdout, "%s",  usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_FAILURE;
            case 'h':   fprintf( stdout, "%s",  usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_SUCCESS;

            case 'i':   inpfile         = optarg;   break;
            case 'o':   outfile         = optarg;   break;
            case 'n':   adjdegf         = optarg;   break;
            case 's':   unscafac        = optarg;   break;
            case 'p':   fileprior       = true;     break;
            case 'v':   suppress        = false;    break;
            default:    break;
        }
    }

    char*   p;
    int     intvalue;
    double  dblunscafac = DEFAULT_UNINF_SCALEFAC;
    double  dbladjdegf = DEFAULT_DEGF_ADJUSTMENT;
    double  tmpval;

    if( !unscafac.empty()) {
        tmpval = strtod( unscafac.c_str(), &p );
        if( errno || *p ) {
            error( "Scale factor specified is invalid." );
            return EXIT_FAILURE;
        }
        if( tmpval <= 0.0 ) {
            error( "Wrong value of scale factor." );
            return EXIT_FAILURE;
        }
        dblunscafac = tmpval;
    }

    if( outfile.empty()) {
        error( "No output file given." );
        return EXIT_FAILURE;
    }

    if( !adjdegf.empty()) {
        tmpval = strtod( adjdegf.c_str(), &p );
        if( errno || *p ) {
            error( "Adjustment to deg. of freedom specified is invalid." );
            return EXIT_FAILURE;
        }
        if( tmpval < -1.0 ) {
            error( "Wrong value of adjustment to deg. of freedom." );
            return EXIT_FAILURE;
        }
        dbladjdegf = tmpval;
    }


    SetQuiet( suppress );


    HDPCalcPrior*   calcor = NULL;
    int             retcode = EXIT_SUCCESS;

    try {
        //Sampling frequency vectors given profiles
        calcor = new HDPCalcPrior();

        if( !calcor ) {
            error( "Not enough memory." );
            return EXIT_FAILURE;
        }

        if( !inpfile.empty())
            calcor->SetFilename( inpfile.c_str());
        calcor->SetOutputFile( outfile.c_str());
        calcor->SetDegFAdjustment( dbladjdegf );
        calcor->SetUninfScaleFactor( dblunscafac );
        calcor->SetUninfPrior( !fileprior );
        calcor->Run();

    } catch( myexception const& ex )
    {
        error( ex.what());
        retcode = EXIT_FAILURE;
    }

    if( calcor ) { delete calcor; calcor = NULL; }

    return retcode;
}

//////////////////////////////////////////////////////////////////////
