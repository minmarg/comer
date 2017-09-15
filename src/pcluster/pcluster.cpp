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
#include "libpro/srcaln/IMAClusters.h"
#include "pcluster.h"


int main( int argc, char *argv[] )
{
    //string values of the parameters
    mystring        input;
    mystring        output;
    mystring        clusthresh;         //clusterization threshold
    mystring        seqidlevel;         //sequence identity level
    mystring        posgapthresh;       //gap percentage per position
    mystring        segwindow;
    mystring        seglowent;
    mystring        seghighent;
    bool            ignoreDs = false;       //ignore positions with gaps in query
    bool            suppress = true;    //suppress warnings
    bool            usingseg = false;   //whether using seg
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
            {"i",       required_argument, 0, 'i'},
            {"o",       required_argument, 0, 'o'},

            {"c",       required_argument, 0, 'c'},
            {"t",       required_argument, 0, 't'},
            {"g",       no_argument,       0, 'g'},
            {"p",       required_argument, 0, 'p'},

            {"h",       no_argument,       0, 'h'},
            {"v",       no_argument,       0, 'v'},

            {"U",       no_argument,       0, 'U'},
            {"w",       required_argument, 0, 'w'},
            {"f",       required_argument, 0, 'f'},
            {"F",       required_argument, 0, 'F'},
            { 0, 0, 0, 0 }
        };
        if(( c = getopt_long_only( argc, argv, "hvi:o:c:t:gp:Uw:f:F:", long_options, &option_index )) == -1 )
            break;
#else
        if(( c = getopt( argc, argv, "hvi:o:c:t:gp:Uw:f:F:" )) == -1 )
            break;
#endif

        switch( c ) {
            case ':':   fprintf( stdout, "No parameter value.\n%s",
                                                usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_FAILURE;
            case '?':   fprintf( stdout, "%s",  usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_FAILURE;
            case 'h':   fprintf( stdout, "%s",  usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_SUCCESS;

            case 'i':   input       = optarg;       break;
            case 'o':   output      = optarg;       break;

            case 'c':   clusthresh  = optarg;       break;
            case 't':   seqidlevel  = optarg;       break;
            case 'g':   ignoreDs    = true;         break;
            case 'p':   posgapthresh= optarg;       break;

            case 'U':   usingseg    = true;                     break;
            case 'w':   segwindow   = optarg; usingseg = true;  break;
            case 'f':   seglowent   = optarg; usingseg = true;  break;
            case 'F':   seghighent  = optarg; usingseg = true;  break;

            case 'v':   suppress   = false;         break;
            default:    break;
        }
    }

    char*   p;
    double  tmpval;
    int     intvalue;
    char    strbuf[BUF_MAX];

    double  sidlevel = CLUSTERING_THRESHOLD;
    double  seqidlevel_val = DEFAULT_SEQIDENTITY_PERCENTAGE;
    double  posgap_ignore = DEFAULT_POSGAPIGNORE_PERCENTAGE;

    size_t  segwinlenval = MAC_SEGSEQ_WIN_LENGTH;
    double  seglowentval = MAC_SEGSEQ_LOW_ENTROPY;
    double  seghighentval = MAC_SEGSEQ_HIGH_ENTROPY;


    if( input.empty()) {
        error( "Input multiple alignment is not specified." );
        return EXIT_FAILURE;
    }

    if( !clusthresh.empty()) {
        intvalue = strtol( clusthresh.c_str(), &p, 10 );

        if( errno || *p ) {
            error( "Clustering threshold specified is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue < 1 || intvalue > 100 ) {
            error( "Clustering threshold is not within interval [1-100]." );
            return EXIT_FAILURE;
        }
        sidlevel = ( double )intvalue / 100.0;
    }

    if( !seqidlevel.empty()) {
        intvalue = strtol( seqidlevel.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Sequence identity level specified is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue < 0 || intvalue > 100 ) {
            error( "Sequence identity level specified is not within interval [0-100]." );
            return EXIT_FAILURE;
        }
        seqidlevel_val = ( double )intvalue / 100.0;
    }

    if( !posgapthresh.empty()) {
        intvalue = strtol( posgapthresh.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Gap percentage specified is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue < 0 || intvalue > 100 ) {
            error( "Gap percentage specified is not within interval [0-100]." );
            return EXIT_FAILURE;
        }
        posgap_ignore = ( double )intvalue / 100.0;
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
    // --

    SetQuiet( suppress );
    message( NULL );

    if( usingseg ) {
        sprintf( strbuf, "SEG Filter:  %4d, %.2f, %.2f", segwinlenval, seglowentval, seghighentval );
        message( strbuf );
    }

    int             ret = EXIT_SUCCESS;
    IMAClusters*    clusters = NULL;

    try {
        clusters = new IMAClusters;

        if( !clusters ) {
            error( "Not enough memory." );
            return EXIT_FAILURE;
        }

        if( usingseg ) {
            clusters->SetSEGWindow( segwinlenval );
            clusters->SetSEGLowEntropy( seglowentval );
            clusters->SetSEGHighEntropy( seghighentval );
        }
        else
            clusters->SetUsingSEGFilter( false );

        clusters->SetIgnoreGapsInQuery( ignoreDs );
        clusters->SetIdentityLevel( seqidlevel_val );
        clusters->SetClusteringThreshold( sidlevel );
        clusters->SetIgnoreGapPercentage( posgap_ignore );
        clusters->ReadAlignment( input.c_str());
        clusters->Make();
        clusters->OutputCounts( output.c_str());

    } catch( myexception const& ex )
    {
        error( ex.what());
        ret = EXIT_FAILURE;
    }

    if( clusters )
        delete clusters;

    return ret;
}
