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
#include "needconfig.h"
#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/datapro.h"
#include "RndMSAGenerator.h"
#include "rndmsa.h"



int main( int argc, char *argv[] )
{
    //string values of the parameters
    mystring        outdir;
    mystring        outfnamepat;
    mystring        lenmsa;
    mystring        nomsas;
    mystring        trprobs;
    mystring        rsprobs;
    mystring        profile;
    mystring        noise;
    mystring        msaeffnos;
    mystring        msawidth;
    bool            suppress = true;    //suppress output
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
            {"l",       required_argument, 0, 'l'},
            {"N",       required_argument, 0, 'N'},
            {"t",       required_argument, 0, 't'},
            {"p",       required_argument, 0, 'p'},
            {"m",       required_argument, 0, 'm'},
            {"r",       required_argument, 0, 'r'},
            {"e",       required_argument, 0, 'e'},
            {"b",       required_argument, 0, 'b'},
            {"v",       no_argument,       0, 'v'},
            {"h",       no_argument,       0, 'h'},
            { 0, 0, 0, 0 }
        };
        if(( c = getopt_long_only( argc, argv, "hvo:n:l:N:t:p:m:r:e:b:", long_options, &option_index )) == -1 )
            break;
#else
        if(( c = getopt( argc, argv, "hvo:n:l:N:t:p:m:r:e:b:" )) == -1 )
            break;
#endif

        switch( c ) {
            case ':':   fprintf( stdout, "No parameter value.\n%s",
                                                usage( argv[0], makeinst, version, verdate ).c_str());  return EXIT_FAILURE;
            case '?':   fprintf( stdout, "%s",  usage( argv[0], makeinst, version, verdate ).c_str());  return EXIT_FAILURE;
            case 'h':   fprintf( stdout, "%s",  usage( argv[0], makeinst, version, verdate ).c_str());  return EXIT_SUCCESS;

            case 'o':   outdir       = optarg;      break;
            case 'n':   outfnamepat  = optarg;      break;

            case 'l':   lenmsa       = optarg;      break;
            case 'N':   nomsas       = optarg;      break;

            case 't':   trprobs      = optarg;      break;
            case 'p':   rsprobs      = optarg;      break;
            case 'm':   profile      = optarg;      break;
            case 'r':   noise        = optarg;      break;

            case 'e':   msaeffnos    = optarg;      break;
            case 'b':   msawidth     = optarg;      break;

            case 'v':   suppress     = false;       break;
            default:    break;
        }
    }

    char*   p;
    double  tmpval;
    int     intvalue, i;
    char    strbuf[BUF_MAX];

    mystring    valSUBMAT = defSUBMAT;
    int         intlenmsa = DEFAULT_MSA_LEN;
    int         intnomsas = DEFAULT_NO_MSAs;
    int         intmsawidth = DEFAULT_MSA_WDT;
    int         intmsaeffnos = 0;

    double      mytrprobs[P_NSTATES];
    double      myrsprobs[NUMAA];
    double      valnoise = DEFAULT_MSA_NOISE;
    const char* pp, *bp;
    int         emsg;
    size_t      length, rbts;

    // -- options --
    MOptions        options;
    mystring        insoptfile = GetFullOptionsFilename();
    mystring        altinsoptfile =
                        mystring( my_dirname( argv[0])) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetOptionsFilename();

    if( !file_exists( insoptfile.c_str()))
        insoptfile = altinsoptfile;

    try {
        options.SetFilename( insoptfile.c_str());
        options.Read();
    } catch( myexception const& ex ) {
        error( ex.what());
        return EXIT_FAILURE;
    }

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
    // probabilities
    for( i = 0; i < P_NSTATES; i++ )
        mytrprobs[i] = (*TRANSPROBS.GetPriors())[i];
    for( i = 0; i < NUMAA; i++ )
        myrsprobs[i] = LOSCORES.PROBABility(i);

    if( !trprobs.empty()) {
        pp = bp = trprobs.c_str();
        length = trprobs.length();
        for( i = 0; i < P_NSTATES; i++, pp+=rbts+1 ) {//1 for separator
            if( length <= size_t( pp - bp )) {
                error("Short of transition probabilities.");
                return EXIT_FAILURE;
            }
            if(( emsg = read_double( pp, length - size_t( pp - bp ), &tmpval, &rbts )) != 0 ) {
                error( TranslateReadError( emsg ));
                return EXIT_FAILURE;
            }
            if( tmpval < 0.0 || 1.0 < tmpval ) {
                error("Transition probabilities not within a range of [0-1].");
                return EXIT_FAILURE;
            }
            mytrprobs[i] = tmpval;
        }
    }
    if( !rsprobs.empty()) {
        pp = bp = rsprobs.c_str();
        length = rsprobs.length();
        for( i = 0; i < NUMAA; i++, pp+=rbts+1 ) {//1 for separator
            if( length <= size_t( pp - bp )) {
                error("Short of residue probabilities.");
                return EXIT_FAILURE;
            }
            if(( emsg = read_double( pp, length - size_t( pp - bp ), &tmpval, &rbts )) != 0 ) {
                error( TranslateReadError( emsg ));
                return EXIT_FAILURE;
            }
            if( tmpval < 0.0 || 1.0 < tmpval ) {
                error("Residue probabilities not within a range of [0-1].");
                return EXIT_FAILURE;
            }
            myrsprobs[i] = tmpval;
        }
    }

    //length
    if( !lenmsa.empty()) {
        intvalue = strtol( lenmsa.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Length is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue < 1 ) {
            error( "Invalid length." );
            return EXIT_FAILURE;
        }
        intlenmsa = intvalue;
    }

    //sample size
    if( !nomsas.empty()) {
        intvalue = strtol( nomsas.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Number of MSAs is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue < 1 ) {
            error( "Invalid number of MSAs." );
            return EXIT_FAILURE;
        }
        intnomsas = intvalue;
    }

    //alignment width
    if( !msawidth.empty()) {
        intvalue = strtol( msawidth.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Number of sequences within MSA is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue < 1 ) {
            error( "Invalid number of sequences within MSA." );
            return EXIT_FAILURE;
        }
        intmsawidth = intvalue;
    }

    //eff. no. sequences of alignment
    if( !msaeffnos.empty()) {
        intvalue = strtol( msaeffnos.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Eff. number of sequences of MSA is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue < 1 || 20 < intvalue ) {
            error( "Invalid eff. number of sequences of MSA." );
            return EXIT_FAILURE;
        }
        intmsaeffnos = intvalue;
    }

    //noise level
    if( !noise.empty()) {
        tmpval = strtod( noise.c_str(), &p );
        if( errno || *p ) {
            error( "Noise level is invalid." );
            return EXIT_FAILURE;
        }
        if( tmpval < 0.0 || 1.0 < tmpval ) {
            error( "Invalid noise level." );
            return EXIT_FAILURE;
        }
        valnoise = tmpval;
    }

    // --

    SetQuiet( suppress );
    message( NULL );

    if( intmsaeffnos )
        sprintf( strbuf, "GEN: %d MSA(s); length %d; effnos %d", intnomsas, intlenmsa, intmsaeffnos );
    else
        sprintf( strbuf, "GEN: %d MSA(s); length %d; width %d", intnomsas, intlenmsa, intmsawidth );
    message( strbuf );

    try {
        RndMSAGenerator* msagen = NULL;

        msagen = new RndMSAGenerator( outdir.c_str(), outfnamepat.c_str());
        msagen = new RndMSAGenerator( outdir.c_str(), outfnamepat.c_str());

        if( msagen == NULL ) {
            error( "Not enough memory." );
            return EXIT_FAILURE;
        }

        msagen->SetOptions( &options );
        msagen->SetGenLength( intlenmsa );
        msagen->SetGenNumber( intnomsas );
        msagen->SetGenWidth( intmsawidth );
        msagen->SetGenEffnos( intmsaeffnos );
        msagen->SetGenTrProbs( mytrprobs );
        msagen->SetGenReProbs( myrsprobs );
        if( !profile.empty())
            msagen->SetGenProfile( profile );
        msagen->SetGenNoise( valnoise );

        msagen->Generate();
        delete msagen;

    } catch( myexception const& ex )
    {
        error( ex.what());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
