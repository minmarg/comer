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
#include "HDPsampler.h"
#include "hdpclust.h"


int main( int argc, char *argv[] )
{
    //string values of the parameters
    mystring        outfile;            //name pattern of output files
    mystring        restartoutfile;     //file of samples to restart sampling with
    mystring        parclust;           //processing behaviour of parallel clustering
    mystring        noscans;            //number of scans (iterations)
    mystring        norgsscans;         //no. interm. restricted Gibbs sampling scans
    mystring        nomhups;            //no. of M-H updates within a single GS scan
    mystring        smprops;            //split/merge proposals only
    mystring        tau;                //concentration parameter tau
    mystring        gamma;              //concentration parameter gamma
    mystring        kappa0;             //prior value of scale factor
    mystring        nu0;                //prior value of degrees of freedom
    mystring        adjdegf;            //adjustment to the degrees of freedom
    mystring        scalefac;           //scale factor for uninformative prior
    bool            jnmodsmpl = false;  //modified Jain-Neal MCMC sampling algorithm
    bool            clst4each = false;  //a separate cluster for each sample
    bool            notblsmpl = false;  //no table sampling
    bool            unfdish = false;    //uniform dish in M-H updates
    bool            unfvect = false;    //uniform dish vectors in M-H updates
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
            {"o",       required_argument, 0, 'o'},
            {"r",       required_argument, 0, 'r'},
            {"e",       no_argument,       0, 'e'},
            {"a",       required_argument, 0, 'a'},
            {"t",       no_argument,       0, 't'},
            {"k",       required_argument, 0, 'k'},
            {"m",       required_argument, 0, 'm'},
            {"u",       required_argument, 0, 'u'},
            {"l",       required_argument, 0, 'l'},
            {"j",       no_argument,       0, 'j'},
            {"d",       no_argument,       0, 'd'},
            {"f",       no_argument,       0, 'f'},
            {"y",       required_argument, 0, 'y'},
            {"g",       required_argument, 0, 'g'},
            {"b",       required_argument, 0, 'b'},
            {"c",       required_argument, 0, 'c'},
            {"n",       required_argument, 0, 'n'},
            {"s",       required_argument, 0, 's'},
            {"p",       no_argument,       0, 'p'},
            {"v",       no_argument,       0, 'v'},
            {"h",       no_argument,       0, 'h'},
            { 0, 0, 0, 0 }
        };
        if(( c = getopt_long_only( argc, argv, "hvo:r:ea:tk:m:u:l:jdfy:g:b:c:n:s:p", 
                  long_options, &option_index )) == -1 )
            break;
#else
        if(( c = getopt( argc, argv, "hvo:r:ea:tk:m:u:l:jdfy:g:b:c:n:s:p" )) == -1 )
            break;
#endif

        switch( c ) {
            case ':':   fprintf( stdout, "No parameter value.\n%s",
                                                usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_FAILURE;
            case '?':   fprintf( stdout, "%s",  usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_FAILURE;
            case 'h':   fprintf( stdout, "%s",  usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_SUCCESS;

            case 'o':   outfile         = optarg;   break;
            case 'r':   restartoutfile  = optarg;   break;
            case 'e':   clst4each       = true;     break;
            case 'a':   parclust        = optarg;   break;
            case 't':   notblsmpl       = true;     break;
            case 'k':   noscans         = optarg;   break;
            case 'm':   norgsscans      = optarg;   break;
            case 'u':   nomhups         = optarg;   break;
            case 'l':   smprops         = optarg;   break;
            case 'j':   jnmodsmpl       = true;     break;
            case 'd':   unfdish         = true;     break;
            case 'f':   unfvect         = true;     break;
            case 'y':   tau             = optarg;   break;
            case 'g':   gamma           = optarg;   break;
            case 'b':   kappa0          = optarg;   break;
            case 'c':   nu0             = optarg;   break;
            case 'n':   adjdegf         = optarg;   break;
            case 's':   scalefac        = optarg;   break;
            case 'p':   fileprior       = true;     break;
            case 'v':   suppress        = false;    break;
            default:    break;
        }
    }

    char*   p;
    int     intvalue;
    int     intparclust = DEFAULT_PARALLEL_CLUST;
    int     intnoscans = DEFAULT_NUMBER_OF_SCANS;
    int     intnorgsscans = DEFAULT_NUMBER_OF_IRGS_SCANS;
    int     intnomhups = DEFAULT_NUMBER_OF_MH_UPDATES;
    int     intsmprops = 0;
    double  dbltau = DEFAULT_TAU;
    double  dblgamma = DEFAULT_GAMMA;
    double  dblkappa0 = DEFAULT_KAPPA0;
    double  dblnu0 = DEFAULT_NU0;
    double  dbladjdegf = DEFAULT_DEGF_ADJUSTMENT;
    double  dblscalefac = DEFAULT_UNINF_SCALEFAC;
    double  tmpval;

    if( restartoutfile.empty()) {
        error( "No input file given." );
        return EXIT_FAILURE;
    }

    if( outfile.empty()) {
        error( "No output file given." );
        return EXIT_FAILURE;
    }

    if( !parclust.empty()) {
        intvalue = strtol( parclust.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Invalid value of parallel processing option." );
            return EXIT_FAILURE;
        }
        if( intvalue != 1 && intvalue != 2 && intvalue != 3 ) {
            error( "Invalid parallel processing option." );
            return EXIT_FAILURE;
        }
        intparclust = intvalue;
    }

    if( !noscans.empty()) {
        intvalue = strtol( noscans.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Number of scans specified is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue < 0 ) {
            error( "Wrong number of scans specified." );
            return EXIT_FAILURE;
        }
        intnoscans = intvalue;
    }

    if( !norgsscans.empty()) {
        intvalue = strtol( norgsscans.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Number of restricted scans specified is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue < -1 ) {
            error( "Wrong number of restricted scans specified." );
            return EXIT_FAILURE;
        }
        intnorgsscans = intvalue;
    }

    if( !nomhups.empty()) {
        intvalue = strtol( nomhups.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Number of MH updates specified is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue < 0 ) {
            error( "Wrong number of MH updates specified." );
            return EXIT_FAILURE;
        }
        intnomhups = intvalue;
    }

    if( !smprops.empty()) {
        intvalue = strtol( smprops.c_str(), &p, 10 );
        if( errno || *p ) {
            error( "Invalid choice of split/merge proposals." );
            return EXIT_FAILURE;
        }
        if( intvalue != 1 && intvalue != 2 ) {
            error( "Wrong choice of split/merge proposals." );
            return EXIT_FAILURE;
        }
        intsmprops = intvalue;
    }

    if( !tau.empty()) {
        tmpval = strtod( tau.c_str(), &p );
        if( errno || *p ) {
            error( "Concentration parameter tau is invalid." );
            return EXIT_FAILURE;
        }
        if( tmpval < -1.0 ) {
            error( "Wrong value of Concentration parameter tau." );
            return EXIT_FAILURE;
        }
        dbltau = tmpval;
    }
    if( !gamma.empty()) {
        tmpval = strtod( gamma.c_str(), &p );
        if( errno || *p ) {
            error( "Concentration parameter gamma is invalid." );
            return EXIT_FAILURE;
        }
        if( tmpval < -1.0 ) {
            error( "Wrong value of Concentration parameter gamma." );
            return EXIT_FAILURE;
        }
        dblgamma = tmpval;
    }

    if( !kappa0.empty()) {
        tmpval = strtod( kappa0.c_str(), &p );
        if( errno || *p ) {
            error( "Scale factor kappa_0 is invalid." );
            return EXIT_FAILURE;
        }
        if( tmpval < -1.0 ) {
            error( "Wrong value of scale factor kappa_0." );
            return EXIT_FAILURE;
        }
        dblkappa0 = tmpval;
    }
    if( !nu0.empty()) {
        tmpval = strtod( nu0.c_str(), &p );
        if( errno || *p ) {
            error( "Value of degrees of freedom is invalid." );
            return EXIT_FAILURE;
        }
        if( tmpval < -1.0 ) {
            error( "Wrong value of degrees of freedom." );
            return EXIT_FAILURE;
        }
        dblnu0 = tmpval;
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

    if( !scalefac.empty()) {
        tmpval = strtod( scalefac.c_str(), &p );
        if( errno || *p ) {
            error( "Scale factor specified is invalid." );
            return EXIT_FAILURE;
        }
        if( tmpval <= 0.0 ) {
            error( "Non-positive values of scale factor are not allowed." );
            return EXIT_FAILURE;
        }
        dblscalefac = tmpval;
    }


    SetQuiet( suppress );


    HDPsampler*     sampler = NULL;
    int             retcode = EXIT_SUCCESS;

    try {
        //Sampling frequency vectors given profiles
        sampler = new HDPsampler();

        if( !sampler ) {
            error( "Not enough memory." );
            return EXIT_FAILURE;
        }

        sampler->SetFilename( restartoutfile.c_str());
        sampler->SetOutputFile( outfile.c_str());
        sampler->SetNoIterations( intnoscans );
        sampler->SetNoRstrGSScans( intnorgsscans );
        sampler->SetNoMHUpdates( intnomhups );
        sampler->SetSMProposalsOnly( intsmprops );
        sampler->SetCluster4Each( clst4each );
        sampler->SetParallelProc( intparclust );
        sampler->SetModJNMCMC( jnmodsmpl );
        sampler->SetMHSampleDishUnf( unfdish );
        sampler->SetMHSampleVectUnf( unfvect );
        sampler->SetNoTableSampling( notblsmpl );
        sampler->SetTauPreset( dbltau );
        sampler->SetGammaPreset( dblgamma );
        sampler->SetKappa0Preset( dblkappa0 );
        sampler->SetNu0Preset( dblnu0 );
        sampler->SetDegFAdjustment( dbladjdegf );
        sampler->SetS0ScaleFac( dblscalefac );
        sampler->SetUninfPrior( !fileprior );
        sampler->Run();

    } catch( myexception const& ex )
    {
        error( ex.what());
        retcode = EXIT_FAILURE;
    }

    if( sampler ) { delete sampler; sampler = NULL; }

    return retcode;
}

//////////////////////////////////////////////////////////////////////
