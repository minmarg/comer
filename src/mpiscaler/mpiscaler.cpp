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

#include "needconfig.h"

#include "rc.h"
#include "libpro/srcpro/MOptions.h"

#include "mpiscaler.h"
#include "PScaler.h"


// -------------------------------------------------------------------------

int main( int argc, char *argv[] )
{
    //names of input file and database
    mystring        database;           //name of database
    mystring        output;             //name of output file
    mystring        information;        //information content threshold
    mystring        lambda;             //value of parameter lambda
    mystring        scomet;             //scoring method
    bool            no_scaling = true;  //do not perform scaling of matrix
    int             c;

    SetArguments( &argc, &argv );
    SetGlobalProgName( argv[0], version );

    while( 1 ) {
#ifdef USE_GETOPT_LONG
        int option_index = 0;
        static struct option long_options[] = {
            {"d",       required_argument, 0, 'd'},
            {"o",       required_argument, 0, 'o'},

            {"I",       required_argument, 0, 'I'},
            {"l",       required_argument, 0, 'l'},
            {"s",       required_argument, 0, 's'},
            { 0, 0, 0, 0 }
        };
        if(( c = getopt_long_only( argc, argv, "hd:I:l:s:o:", long_options, &option_index )) == -1 )
            break;
#else
        if(( c = getopt( argc, argv, "hd:I:l:s:o:" )) == -1 )
            break;
#endif
        switch( c ) {
            case ':':   error( "You must specify a value for the arguments." );
                        fprintf( stdout, "%s", usage( argv[0], instructions, version, verdate ).c_str());   return EXIT_FAILURE;
            case '?':   fprintf( stdout, "%s", usage( argv[0], instructions, version, verdate ).c_str());   return EXIT_FAILURE;
            case 'h':   fprintf( stdout, "%s", usage( argv[0], instructions, version, verdate ).c_str());   return EXIT_SUCCESS;

            case 'd':   database    = optarg;       break;
            case 'o':   output      = optarg;       break;
            case 'I':   information = optarg;       break;
            case 'l':   lambda      = optarg; 
                        no_scaling  = false;        break;
            case 's':   scomet      = optarg;       break;
            default:    break;
        }
    }

    char*   p;
    double  tmpval;
    mystring valSUBMAT = defSUBMAT;
    double  infrm_value = DEFAULT_INFORMATION_CONTENT;  //information content threshold
    double  lambdaval   = -1.0;                         //indication of absence of the value specification

    AbstractScoreMatrix::TScaling
            scaling     = DEFAULT_PRECISION;            //precision strategy

    TMask   masking      = DEFAULT_MASKING;             //how to mask positions if needed
    AbstractScoreMatrix::TType
            scoring     = DEFAULT_SCORING_METHOD;       //scoring method
    bool    usingmasking = true;                        //whether masking is in use


    // -- sub. matrix --
    try {
        SetLOSCORES( valSUBMAT, DEFAULT_PSCORES_FILENAME );
    } catch( myexception const& ex ) {
        error( ex.what());
        return EXIT_FAILURE;
    }

    if( !information.empty()) {
        tmpval = strtod( information.c_str(), &p );
        if( errno || *p ) {
            error( "Information threshold specified is invalid." );
            return EXIT_FAILURE;
        }
        if( tmpval < 0.0 ) {
            error( "Information threshold specified is negative." );
            return EXIT_FAILURE;
        }
        infrm_value = tmpval;
    }


    if( !lambda.empty()) {
        tmpval = strtod( lambda.c_str(), &p );
        if( errno || *p ) {
            error( "Parameter lambda specified is invalid." );
            return EXIT_FAILURE;
        }
        if( tmpval <= 0.0 || 1.0 <= tmpval ) {
            error( "Parameter lambda specified is not within interval [0-1]." );
            return EXIT_FAILURE;
        }
        lambdaval = tmpval;
    }

    if( !scomet.empty()) {
        if( scomet == "pro" )
            scoring = AbstractScoreMatrix::ParallelUniversal;
        else if( scomet == "LSO" ) {
            scoring = AbstractScoreMatrix::ParallelLSOUniversal;
            scaling = AbstractScoreMatrix::AutoScalling;
        }
        else if( scomet == "HDP" ) {
            scoring = AbstractScoreMatrix::ParallelHDPUniversal;
            scaling = AbstractScoreMatrix::AutoScalling;
        }
        else {
            error( "Invalid scoring method." );
            return EXIT_FAILURE;
        }
    }

    if( database.empty()) {
        error( "Database is not specified." );
        return EXIT_FAILURE;
    }


    int         ret = EXIT_SUCCESS;
    PScaler*    dispatcher = NULL;


    try {
        dispatcher = new PScaler(
                        GetFullParamFilename(),
                        database.c_str(),
                        output.c_str(),
                        no_scaling,
                        infrm_value,
                        lambdaval,
                        scaling,
                        masking,
                        scoring
        );

        if( !dispatcher ) {
            error( "Not enough memory." );
            return EXIT_FAILURE;
        }

        dispatcher->Run();

    } catch( myexception const& ex )
    {
        error( ex.what());
        ret = EXIT_FAILURE;
    }

    if( dispatcher )
        delete dispatcher;

    return ret;
}
