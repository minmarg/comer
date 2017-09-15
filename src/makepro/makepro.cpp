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
#include "libpro/srcpro/Serializer.h"
#include "libpro/srcpro/GapScheme.h"
#include "libpro/srcpro/DistributionMatrix.h"
#include "libpro/srcpro/SEGProfile.h"
#include "libpro/srcaln/InputMultipleAlignment.h"
#include "libHDP/HDPbase.h"
#include "makepro.h"



int main( int argc, char *argv[] )
{
    int             c;
    char            strbuf[BUF_MAX];
    //string values of parameters
    mystring        input;
    mystring        output;
    mystring        ascii;
    mystring        optfile;
    bool            suppress = true;//suppress warnings

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
            {"i",       required_argument, 0, 'i'},
            {"o",       required_argument, 0, 'o'},
            {"f",       required_argument, 0, 'f'},
            {"p",       required_argument, 0, 'p'},
            {"v",       no_argument,       0, 'v'},
            { 0, 0, 0, 0 }
        };
        if(( c = getopt_long_only( argc, argv, "hvi:o:f:p:", long_options, &option_index )) == -1 )
            break;
#else
        if(( c = getopt( argc, argv, "hvi:o:f:p:" )) == -1 )
            break;
#endif

        switch( c ) {
            case ':':   fprintf( stdout, "Values should follow options.\n%s",
                                                usage( argv[0], makeinst, version, verdate ).c_str());  return EXIT_FAILURE;
            case '?':   fprintf( stdout, "%s",  usage( argv[0], makeinst, version, verdate ).c_str());  return EXIT_FAILURE;
            case 'h':   fprintf( stdout, "%s",  usage( argv[0], makeinst, version, verdate ).c_str());  return EXIT_SUCCESS;

            case 'i':   input       = optarg;       break;
            case 'o':   output      = optarg;       break;

            case 'f':   ascii       = optarg;       break;
            case 'p':   optfile     = optarg;       break;
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
    //// hdp1
    mystring        inhdp1file = GetFullHDP1Filename();
    mystring        altinhdp1file =
                        mystring( my_dirname( argv[0] )) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetHDP1Filename();

    if( !file_exists( inhdp1file.c_str()))
        inhdp1file = altinhdp1file;
    //// hdpctx
    mystring        inhdpctxfile = GetFullHDPCtxFilename();
    mystring        altinhdpctxfile =
                        mystring( my_dirname( argv[0] )) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetHDPCtxFilename();

    if( !file_exists( inhdpctxfile.c_str()))
        inhdpctxfile = altinhdpctxfile;
    ////
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
    double          valcADJWGT = OPTIONS.GetcADJWGT();

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

    //{{ -- mix. parameters --
    try {
        mystring lcval1 = valTFRMIX;
        mystring lcval2 = valSCOADJ;
        lcval1.lower();
        lcval2.lower();
        if( lcval1 == "hdpctx" || 
            lcval2 == "hdpctx" || lcval2 == "hdpsco") {
            SetHDPBASE( inhdp1file.c_str());
            if( lcval2 == "hdpsco")
                intscoadj = scoadjHDPsco;
            else if( lcval2 == "hdpctx")
                intscoadj = scoadjHDPCtx;
            else if( lcval1 == "hdpctx")
                inttfrmix = tfrmixHDPCtx;
            HDPBASE.SetMixWeight( valMIXWGT );
            HDPBASE.SetNoSupClusters( valSUPCLT );
            HDPBASE.SetAdjWeight( valADJWGT );
            //
            if( GetHDPctsUsed()) {
                SetHDPctBASE( inhdpctxfile.c_str());
                HDPctBASE.SetAdjWeight( valcADJWGT );
                HDPctBASE.SetNoSupClusters( valSUPCLT );
            }
        }
        else {
            if( lcval1 != "no") {
                error("Unknown choice value for Mixing of profile target frequencies.");
                return EXIT_FAILURE;
            }
            if( lcval2 != "no") {
                error("Unknown choice value for Score adjustment.");
                return EXIT_FAILURE;
            }
        }
    } catch( myexception const& ex ) {
        error( ex.what());
        return EXIT_FAILURE;
    }
    //}}

    //{{ -- Print --
    sprintf( strbuf, "IDENTITY = %.2f", valIDENTITY );
    message( strbuf, false );

    sprintf( strbuf, "PCFWEIGHT = %.1f, MINALNFRN = %.2f, MINALNPOS = %d", valPCFWEIGHT, valMINALNFRN, valMINALNPOS );
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

        //target frequencies will be mixed, if required, just before searching
//         inmaln->SetTarFrMix( inttfrmix );
        inmaln->SetScoAdjment( intscoadj );
        inmaln->SetHDPbase( &HDPBASE );
        if( GetHDPctsUsed())
            inmaln->SetHDPctBase( &HDPctBASE );

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


        inmaln->SetIdentityLevel( valIDENTITY );
        inmaln->SetComputeDELETEstates( valDELSTATE );
        inmaln->SetExtentMinWindow( valMINALNPOS );
        inmaln->SetExtentMinSeqPercentage( valMINALNFRN );
        inmaln->SetPseudoCountWeight( valPCFWEIGHT );


        inmaln->ReadAlignment( input.c_str());
        message( "Constructing profile...");
        inmaln->ConstructProfile();

        FrequencyMatrix frequencies;
        LogOddsMatrix   pssm;
        GapScheme       gaps;

        inmaln->ExportFrequenciesTo( frequencies );
        inmaln->ExportPSSMTo( pssm );
        inmaln->ExportGapWeightsTo( gaps );

        // SEG logic
        if( valINVLCFILTER ) {
            SEGProfile  segpro(
                    frequencies,
                    pssm,
                    valLCWINDOW,
                    valLCLOWENT,
                    valLCHIGHENT
            );
            segpro.SetDistance( valDISTANCE );
            segpro.Run();
            segpro.MaskSeggedPositions( frequencies, pssm, gaps );
        }
        // --

        serializer.WriteProfile( output.c_str(), frequencies, pssm, gaps );


        if( !ascii.empty())
            //inmaln->OutputProfile( ascii.c_str());
            OutputProfile( ascii.c_str(), frequencies, pssm, gaps );


//         serializer.SerializeProfile( frequencies, pssm, gaps, output.c_str());//obsolete
//         serializer.DeserializeProfile( frequencies, pssm, gaps, output.c_str());//obsolete

//         gaps.ComputeCosts();
//         gaps.OutputGapScheme();
//         frequencies.OutputMatrix();
//         pssm.OutputMatrix();

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
