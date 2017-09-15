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
#include "libpro/srcpro/Database.h"
#include "makedb.h"



int main( int argc, char *argv[] )
{
    //string values of the parameters
    mystring        output;
    mystring        distribution;
    mystring        directory;
    mystring        segwindow;
    mystring        seglowent;
    mystring        seghighent;
    mystring        segdistance;
    bool            suppress = true;    //suppress output
    bool            usingseg = false;   //whether using seg
    int             c;

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
            {"t",       required_argument, 0, 't'},
            {"d",       required_argument, 0, 'd'},
            {"v",       no_argument,       0, 'v'},

            {"U",       no_argument,       0, 'U'},
            {"w",       required_argument, 0, 'w'},
            {"f",       required_argument, 0, 'f'},
            {"F",       required_argument, 0, 'F'},
            {"D",       required_argument, 0, 'D'},
            { 0, 0, 0, 0 }
        };
        if(( c = getopt_long_only( argc, argv, "hvo:t:d:Uw:f:F:D:", long_options, &option_index )) == -1 )
            break;
#else
        if(( c = getopt( argc, argv, "hvo:t:d:Uw:f:F:D:" )) == -1 )
            break;
#endif

        switch( c ) {
            case ':':   fprintf( stdout, "No parameter value.\n%s",
                                                usage( argv[0], makeinst, version, verdate ).c_str());  return EXIT_FAILURE;
            case '?':   fprintf( stdout, "%s",  usage( argv[0], makeinst, version, verdate ).c_str());  return EXIT_FAILURE;
            case 'h':   fprintf( stdout, "%s",  usage( argv[0], makeinst, version, verdate ).c_str());  return EXIT_SUCCESS;

            case 'o':   output       = optarg;      break;
            case 't':   distribution = optarg;      break;
            case 'd':   directory    = optarg;      break;

            case 'U':   usingseg    = true;                     break;
            case 'w':   segwindow   = optarg; usingseg = true;  break;
            case 'f':   seglowent   = optarg; usingseg = true;  break;
            case 'F':   seghighent  = optarg; usingseg = true;  break;
            case 'D':   segdistance = optarg; usingseg = true;  break;

            case 'v':   suppress     = false;       break;
            default:    break;
        }
    }

    char*   p;
    double  tmpval;
    int     intvalue;
    char    strbuf[BUF_MAX];

    TFVectorProbabilities   distribtype = DEFAULT_DISTRIBUTION_TYPE;
    mystring                description;
    mystring                valSUBMAT = defSUBMAT;

    size_t                  segwinlenval = DEFAULT_SEGPRO_WIN_LENGTH;
    double                  seglowentval = DEFAULT_SEGPRO_LOW_ENTROPY;
    double                  seghighentval = DEFAULT_SEGPRO_HIGH_ENTROPY;
    double                  segdistanceval = DEFAULT_SEGPRO_EQUALITY_DISTANCE;


    // -- sub. matrix --
    try {
        SetLOSCORES( valSUBMAT, DEFAULT_PSCORES_FILENAME );
    } catch( myexception const& ex ) {
        error( ex.what());
        return EXIT_FAILURE;
    }

    if( output.empty()) {
        error( "Database name as output is not specified." );
        return EXIT_FAILURE;
    }

    if( ! distribution.empty()) {
        if( distribution == "simple" ) {
            distribtype = DISCRETE;
        }
        else
        if( distribution == "multin" ) {
            distribtype = MULTINOMIAL;
        }
        else
        if( distribution == "profile" ) {
            distribtype = PROVECTOR;
        }
        else {
            error( "Unknown distribution type." );
            return EXIT_FAILURE;
        }
    }

    switch( distribtype ) {
        case DISCRETE:
            description = "Simple discrete distribution.";
            break;

        case MULTINOMIAL:
            description = "Multinomial distribution.";
            break;

        case PROVECTOR:
            description = "Profile vector distribution.";
            break;

        default:
            error( "Unknown distribution type." );
            return EXIT_FAILURE;
    }

    // SEG options --
    if( !segwindow.empty()) {
        intvalue = strtol( segwindow.c_str(), &p, 10 );

        if( errno || *p ) {
            error( "Window length is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue < 2 ) {
            error( "Window length must be >1." );
            return EXIT_FAILURE;
        }
        segwinlenval = intvalue;
    }

    if( !seglowent.empty()) {
        tmpval = strtod( seglowent.c_str(), &p );
        if( errno || *p ) {
            error( "Low entropy threshold is invalid." );
            return EXIT_FAILURE;
        }
        if( tmpval < 0.0 ) {
            error( "Low entropy threshold is negative." );
            return EXIT_FAILURE;
        }
        seglowentval = tmpval;
    }

    if( !seghighent.empty()) {
        tmpval = strtod( seghighent.c_str(), &p );
        if( errno || *p ) {
            error( "High entropy threshold is invalid." );
            return EXIT_FAILURE;
        }
        if( tmpval < 0.0 ) {
            error( "High entropy threshold is negative." );
            return EXIT_FAILURE;
        }
        seghighentval = tmpval;
    }

    if( !segdistance.empty()) {
        tmpval = strtod( segdistance.c_str(), &p );
        if( errno || *p ) {
            error( "Distance of equivalence is invalid." );
            return EXIT_FAILURE;
        }
        if( tmpval < 0.0 ) {
            error( "Distance of equivalence is negative." );
            return EXIT_FAILURE;
        }
        segdistanceval = tmpval;
    }
    // --

    SetQuiet( suppress );
    message( description.c_str());
    message( NULL );

    if( usingseg ) {
        sprintf( strbuf, "SEG:  %4d, %.2f, %.2f; Distance %.2f",
                segwinlenval, seglowentval, seghighentval, segdistanceval );
        message( strbuf );
    }


    try {
        Database*   database;

        if( !directory.empty())
            database = new Database( output.c_str(), directory.c_str());
        else
            database = new Database( output.c_str(), argv + 1, argc - 1 );

        if( !database ) {
            error( "Not enough memory." );
            return EXIT_FAILURE;
        }

        if( usingseg )
            database->SetSegParameters(
                segwinlenval,
                seglowentval,
                seghighentval,
                segdistanceval
            );

        database->SetDistributionType( distribtype );
        database->Make();
        delete database;

    } catch( myexception const& ex )
    {
        error( ex.what());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
