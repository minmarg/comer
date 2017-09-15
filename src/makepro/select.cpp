/***************************************************************************
 *   Copyright (C) 2013 by Mindaugas Margelevicius                         *
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
#include "libpro/srcpro/Serializer.h"
#include "libpro/srcpro/GapScheme.h"
#include "libpro/srcpro/DistributionMatrix.h"
#include "libpro/srcpro/SEGProfile.h"
#include "libpro/srcaln/InputMultipleAlignment.h"
#include "libHDP/HDPbase.h"
#include "select.h"



int main( int argc, char *argv[] )
{
    int             c;
    char            strbuf[BUF_MAX];
    //string values of parameters
    mystring        input;
    mystring        output;
    mystring        optfile;
    bool            noselect = false;//do not apply select
    bool            suppress = true;//suppress warnings

    SetArguments( &argc, &argv );
    SetGlobalProgName( argv[0], version );

    if( argc <= 1 ) {
        fprintf( stdout, "%s", usage( argv[0], selinst, version, verdate ).c_str());
        return EXIT_SUCCESS;
    }

    while( 1 ) {
#ifdef USE_GETOPT_LONG
        int option_index = 0;
        static struct option long_options[] = {
            {"i",       required_argument, 0, 'i'},
            {"o",       required_argument, 0, 'o'},
            {"p",       required_argument, 0, 'p'},
            {"t",       no_argument,       0, 't'},
            {"v",       no_argument,       0, 'v'},
            { 0, 0, 0, 0 }
        };
        if(( c = getopt_long_only( argc, argv, "hvi:o:p:t", long_options, &option_index )) == -1 )
            break;
#else
        if(( c = getopt( argc, argv, "hvi:o:p:t" )) == -1 )
            break;
#endif

        switch( c ) {
            case ':':   fprintf( stdout, "Values should follow options.\n%s",
                                                usage( argv[0], selinst, version, verdate ).c_str());  return EXIT_FAILURE;
            case '?':   fprintf( stdout, "%s",  usage( argv[0], selinst, version, verdate ).c_str());  return EXIT_FAILURE;
            case 'h':   fprintf( stdout, "%s",  usage( argv[0], selinst, version, verdate ).c_str());  return EXIT_SUCCESS;

            case 'i':   input       = optarg;       break;
            case 'o':   output      = optarg;       break;

            case 'p':   optfile     = optarg;       break;
            case 't':   noselect    = true;         break;
            case 'v':   suppress    = false;        break;
            default:    break;
        }
    }

    SetQuiet( suppress );


    if( input.empty()) {
        error( "Input multiple alignment is not specified." );
        return EXIT_FAILURE;
    }
    if( output.empty()) {
        error( "Output file is not specified." );
        return EXIT_FAILURE;
    }


    mystring        insparamfile = GetFullParamFilename();
    mystring        altinsparamfile =
                        mystring( my_dirname( argv[0] )) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetParamFilename();

    if( !file_exists( insparamfile.c_str()))
        insparamfile = altinsparamfile;

    mystring        inhdpctxfile = GetFullHDPCtxFilename();
    mystring        altinhdpctxfile =
                        mystring( my_dirname( argv[0] )) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetHDPCtxFilename();

    if( !file_exists( inhdpctxfile.c_str()))
        inhdpctxfile = altinhdpctxfile;

    mystring        inhdpfile = inhdpctxfile;

    mystring        insoptfile = GetFullOptionsFilename();
    mystring        altinsoptfile =
                        mystring( my_dirname( argv[0] )) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetOptionsFilename();

    if( !file_exists( insoptfile.c_str()))
        insoptfile = altinsoptfile;

    if( optfile.empty())
        optfile = insoptfile;

    try {
        OPTIONS.Read( optfile.c_str());
    } catch( myexception const& ex ) {
        error( ex.what());
        return EXIT_FAILURE;
    }


    double          valIDENTITY = OPTIONS.GetIDENTITY();
    int             valDELSTATE = OPTIONS.GetDELSTATE();

    double          valPCFWEIGHT = OPTIONS.GetPCFWEIGHT();
    double          valMINALNFRN = OPTIONS.GetMINALNFRN();
    int             valMINALNPOS = OPTIONS.GetMINALNPOS();
    mystring        valSUBMAT = OPTIONS.GetSUBMAT();
    mystring        valTFRMIX = OPTIONS.GetTFRMIX();
    int             inttfrmix = tfrmixNo;
    double          valMIXWGT = OPTIONS.GetMIXWGT();
    mystring        valSCOADJ = OPTIONS.GetSCOADJ();
    int             intscoadj = scoadjNo;
    int             valSUPCLT = OPTIONS.GetSUPCLT();
    double          valADJWGT = OPTIONS.GetADJWGT();

    int             valHCFILTER = OPTIONS.GetHCFILTER();
    int             valHCWINDOW = OPTIONS.GetHCWINDOW();
    double          valHCLOWENT = OPTIONS.GetHCLOWENT();
    double          valHCHIGHENT = OPTIONS.GetHCHIGHENT();

    int             valINVLCFILTER = OPTIONS.GetINVLCFILTER();
    int             valLCFILTEREACH = OPTIONS.GetLCFILTEREACH();
    int             valLCWINDOW = OPTIONS.GetLCWINDOW();
    double          valLCLOWENT = OPTIONS.GetLCLOWENT();
    double          valLCHIGHENT = OPTIONS.GetLCHIGHENT();
    double          valDISTANCE = OPTIONS.GetDISTANCE();


    //{{ -- sub. matrix --
    try {
        SetLOSCORES( valSUBMAT, DEFAULT_PSCORES_FILENAME );
    } catch( myexception const& ex ) {
        error( ex.what());
        return EXIT_FAILURE;
    }
    //}}


    //{{ -- Print --
    if( noselect ) {
        sprintf( strbuf, "IDENTITY = %.2f, DELSTATE = %s", valIDENTITY, valDELSTATE? "yes": "no" );
        message( strbuf, false );
        sprintf( strbuf, "HCFILTER = %s", valHCFILTER? "yes": "no" );
        message( strbuf, false );
        if( valHCFILTER ){
            sprintf( strbuf, "HCWINDOW = %d, HCLOWENT = %.2f, HCHIGHENT = %.2f",
                        valHCWINDOW, valHCLOWENT, valHCHIGHENT );
            message( strbuf, false );
        }
        sprintf( strbuf, "INVLCFILTER = %s, LCFILTEREACH = %s",
                    valINVLCFILTER? "yes": "no", valLCFILTEREACH? "yes": "no" );
        message( strbuf, false );
        if( valINVLCFILTER || valLCFILTEREACH ){
            sprintf( strbuf, "LCWINDOW = %d, LCLOWENT = %.2f, LCHIGHENT = %.2f",
                        valLCWINDOW, valLCLOWENT, valLCHIGHENT );
            if( valINVLCFILTER )
                sprintf( strbuf + strlen( strbuf ), ", DISTANCE = %.2f", valDISTANCE );
            message( strbuf, false );
        }
        message( NULL );
    }
    //}}

    int                     ret = EXIT_SUCCESS;
    Serializer              serializer;
    InputMultipleAlignment* inmaln = NULL;

    try {
        inmaln = new InputMultipleAlignment();

        if( !inmaln ) {
            error( "Not enough memory." );
            return EXIT_FAILURE;
        }

        if( valHCFILTER ) {
            inmaln->SetSEGWindow( valHCWINDOW );
            inmaln->SetSEGLowEntropy( valHCLOWENT );
            inmaln->SetSEGHighEntropy( valHCHIGHENT );
        }
        else
            inmaln->SetUsingSEGFilter( false );

        if( valLCFILTEREACH ) {
            inmaln->SetSeqSEGWindow( valLCWINDOW );
            inmaln->SetSeqSEGLowEntropy( valLCLOWENT );
            inmaln->SetSeqSEGHighEntropy( valLCHIGHENT );
        }
        else
            inmaln->SetUsingSeqSEGFilter( false );

        inmaln->SetKeepTitles( true );
        inmaln->SetIdentityLevel( valIDENTITY );
        inmaln->SetComputeDELETEstates( valDELSTATE );
        inmaln->SetExtentMinWindow( valMINALNPOS );
        inmaln->SetExtentMinSeqPercentage( valMINALNFRN );
        inmaln->SetPseudoCountWeight( valPCFWEIGHT );


        inmaln->ReadAlignment( input.c_str());
        if( noselect ) {
            inmaln->PlainPreprocess();
        } else {
            message( "Selecting sequences...");
            inmaln->SelectSequences();
        }
        inmaln->PutAlignment( output.c_str());


        message( "Done.");

    } catch( myexception const& ex )
    {
        error( ex.what());
        ret = EXIT_FAILURE;
    }

    if( inmaln )
        delete inmaln;

    return ret;
}
