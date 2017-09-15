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
#include "libpro/srcsco/ProfileMatrix.h"
#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/FrequencyStore.h"
#include "libpro/srcpro/SSSScores.h"
#include "libpro/srcpro/CVS2Scores.h"
#include "libpro/srcaln/ProfileSearching.h"
#include "libHDP/HDPbase.h"
#include "comer.h"

// Functional-interface-----------------------------------------------------

bool GetScoringScheme( const mystring&, AbstractScoreMatrix::TType* );
bool GetStatBehaviour( bool, AbstractScoreMatrix::TBehaviour* );
bool GetMasking( bool, bool, TMask* );

// =========================================================================

int main( int argc, char *argv[] )
{
    int             c;
    char            strbuf[BUF_MAX];
    //names of input file and database
    mystring        input;
    mystring        database;
    mystring        output;
    mystring        optfile;
    bool            suppress = true;    //suppress warnings

    SetGlobalProgName( argv[0], version );

    while( 1 ) {
#ifdef USE_GETOPT_LONG
        int option_index = 0;
        static struct option long_options[] = {
            {"i",       required_argument, 0, 'i'},
            {"d",       required_argument, 0, 'd'},
            {"o",       required_argument, 0, 'o'},
            {"p",       required_argument, 0, 'p'},
            { 0, 0, 0, 0 }
        };
        if(( c = getopt_long_only(
                    argc, argv,
                    "hi:d:o:p:v",
                    long_options,
                    &option_index )) == -1 )
            break;
#else
        if(( c = getopt( argc, argv, "hi:d:o:p:v" )) == -1 )
            break;
#endif
        switch( c ) {
            case ':':   error( "Values should follow options." );
                        fprintf( stdout, "%s", usage( argv[0], instructions, version, verdate ).c_str());   return EXIT_FAILURE;
            case '?':   fprintf( stdout, "%s", usage( argv[0], instructions, version, verdate ).c_str());   return EXIT_FAILURE;
            case 'h':   fprintf( stdout, "%s", usage( argv[0], instructions, version, verdate ).c_str());   return EXIT_SUCCESS;

            case 'i':   input       = optarg;       break;
            case 'd':   database    = optarg;       break;
            case 'o':   output      = optarg;       break;
            case 'p':   optfile     = optarg;       break;

            case 'v':   suppress    = false;        break;
            default:    break;
        }
    }

    SetQuiet( suppress );

    if( input.empty() && !database.empty()) {
        error( "Input multiple alignment file is not specified." );
        return EXIT_FAILURE;
    }
    if( !input.empty() && database.empty()) {
        error( "Database is not specified." );
        return EXIT_FAILURE;
    }
    if( input.empty() && database.empty()) {
        fprintf( stdout, "%s", usage( argv[0], instructions, version, verdate ).c_str());
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

    mystring        inhdpscofile = GetFullHDPscoFilename();
    mystring        altinhdpscofile =
                        mystring( my_dirname( argv[0] )) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetHDPscoFilename();

    if( !file_exists( inhdpscofile.c_str()))
        inhdpscofile = altinhdpscofile;
    //// hdpctx
    mystring        inhdpctxfile = GetFullHDPCtxFilename();
    mystring        altinhdpctxfile =
                        mystring( my_dirname( argv[0] )) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetHDPCtxFilename();

    if( !file_exists( inhdpctxfile.c_str()))
        inhdpctxfile = altinhdpctxfile;

    mystring        inhdpctsfile = GetFullHDPctsFilename();
    mystring        altinhdpctsfile =
                        mystring( my_dirname( argv[0] )) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetHDPctsFilename();

    if( !file_exists( inhdpctsfile.c_str()))
        inhdpctsfile = altinhdpctsfile;
    //// cvs2s
    mystring        incvs2sfile = GetFullCVS2SFilename();
    mystring        altincvs2sfile =
                        mystring( my_dirname( argv[0] )) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetCVS2SFilename();

    if( !file_exists( incvs2sfile.c_str()))
        incvs2sfile = altincvs2sfile;


    //// ssss
    mystring        inssssfile = GetFullSSSSFilename();
    mystring        altinssssfile =
                        mystring( my_dirname( argv[0] )) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetSSSSFilename();

    if( !file_exists( inssssfile.c_str()))
        inssssfile = altinssssfile;


    mystring        inhdpssssfile = GetFullHDPSSSSFilename();
    mystring        altinhdpssssfile =
                        mystring( my_dirname( argv[0] )) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetHDPSSSSFilename();

    if( !file_exists( inhdpssssfile.c_str()))
        inhdpssssfile = altinhdpssssfile;


    mystring        inihdpssssfile = GetFulliHDPSSSSFilename();
    mystring        altinihdpssssfile =
                        mystring( my_dirname( argv[0] )) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetiHDPSSSSFilename();

    if( !file_exists( inihdpssssfile.c_str()))
        inihdpssssfile = altinihdpssssfile;


    //// OPTIONS
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

    double          valEVAL = OPTIONS.GetEVAL();
    int             valNOHITS = OPTIONS.GetNOHITS();
    int             valNOALNS = OPTIONS.GetNOALNS();
    int             valSHOW = OPTIONS.GetSHOW();

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

    double          valCVSWGT = OPTIONS.GetCVSWGT();

    double          valSSSWGT = OPTIONS.GetSSSWGT();
    double          valSSSHDP = OPTIONS.GetSSSHDP();

    double          valINFCON = OPTIONS.GetINFCON();
    int             valMASKAFTER = OPTIONS.GetMASKAFTER();
    double          valSCALEDOWN = OPTIONS.GetSCALEDOWN();

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

    int             fixedOPENCOST = DEFAULTGAPOPENCOST;
    int             intOPENCOST = 4;
    bool            boolAutoOpenCost = true;
    int             valEXTCOST = OPTIONS.GetEXTCOST();
    double          valDELPROBWEIGHT = OPTIONS.GetDELPROBWEIGHT();
    mystring        valSCHEME = OPTIONS.GetSCHEME();
    double          valMINPP = OPTIONS.GetMINPP();
    int             valCOMPSTATS = OPTIONS.GetCOMPSTATS();
    int             valUSEGCPROBS = OPTIONS.GetUSEGCPROBS();

    double          valGPROBEVAL = OPTIONS.GetGPROBEVAL();
    double          valGPFARGWEIGHT = OPTIONS.GetGPFARGWEIGHT();
    double          valGPFARGSHIFT = OPTIONS.GetGPFARGSHIFT();

    double          valAC1NUMER = OPTIONS.GetAC1NUMER();
    double          valAC2UBNUMER = OPTIONS.GetAC2UBNUMER();
    double          valAC2LOGSCALE = OPTIONS.GetAC2LOGSCALE();
    double          valAC2DENOMSCALE = OPTIONS.GetAC2DENOMSCALE();
    int             valANPOSCOR = OPTIONS.GetANPOSCOR();
    int             valPROHIBITCOR = OPTIONS.GetPROHIBITCOR();

    double          valINFCON2UB = OPTIONS.GetINFCON2UB();
    double          valINFCON2NUMER = OPTIONS.GetINFCON2NUMER();
    double          valINFCON2LOGSCALE = OPTIONS.GetINFCON2LOGSCALE();
    double          valINFCONALTNUMER = OPTIONS.GetINFCONALTNUMER();
    double          valINFCONALTLOGSCALE = OPTIONS.GetINFCONALTLOGSCALE();

    int             valHSPLEN = OPTIONS.GetHSPLEN();
    int             valHSPMINSCORE = OPTIONS.GetHSPMINSCORE();
    int             valHSPMAXDIST = OPTIONS.GetHSPMAXDIST();
    int             valNOHSPS = OPTIONS.GetNOHSPS();

    OPTIONS.GetOPENCOST( &intOPENCOST, &boolAutoOpenCost );
    if( !boolAutoOpenCost )
        fixedOPENCOST = intOPENCOST;


    AbstractScoreMatrix::TType      method      = DEFAULT_SCORING_SCHEME;
    AbstractScoreMatrix::TBehaviour behaviour   = DEFAULT_STATISTICAL_BEHAVIOUR;
    AbstractScoreMatrix::TScaling   precision   = DEFAULT_PRECISION;
    ProfileAlignment::TAlgorithm    alnalgo     = DEFAULT_ALNALGORITHM;

    TMask   masking      = DEFAULT_MASKING;     //how to mask positions if needed
    bool    usingmasking = true;                //whether masking is in use

//     {
//         alnalgo = ProfileAlignment::Hybrid;
//         precision = AbstractScoreMatrix::NoScaling;
//     }

    // -- sub. matrix --
    try {
        SetLOSCORES( valSUBMAT, DEFAULT_PSCORES_FILENAME );
    } catch( myexception const& ex ) {
        error( ex.what());
        return EXIT_FAILURE;
    }

    //{{ SS state score parameters
    try {/////////////////////////
        if( valSSSWGT ) {
            SSSSCORES.ReadScores( inssssfile.c_str());
            SSSSCORES.SetSSSWeight( valSSSWGT );
            if( valSSSHDP ) {
#if defined( USEiHDPSSSSCORES )
                iHDPSSSSCORES.ReadScores( inihdpssssfile.c_str());
#else
                HDPSSSSCORES.ReadScores( inhdpssssfile.c_str());
#endif
                iHDPSSSSCORES.SetiHDPSSSWeight( valSSSHDP );
                HDPSSSSCORES.SetHDPSSSWeight( valSSSHDP );
            }
        }
    } catch( myexception const& ex ) {
        error( ex.what());
        return EXIT_FAILURE;
    }
    //}}

    //{{ log-odds score-to-score parameters
    try {/////////////////////////
        if( valCVSWGT ) {
            CVS2SCORES.ReadScores( incvs2sfile.c_str());
            CVS2SCORES.SetCVSWeight( valCVSWGT );
        }
    } catch( myexception const& ex ) {
        error( ex.what());
        return EXIT_FAILURE;
    }
    //}}

    //{{ -- mix. parameters --
    try {
        if( 1||!valSSSHDP ) {
            mystring lcval1 = valTFRMIX;
            mystring lcval2 = valSCOADJ;
            lcval1.lower();
            lcval2.lower();
            if( lcval1 == "hdpctx" || 
                lcval2 == "hdpctx" || lcval2 == "hdpsco") {
                SetHDPBASE( inhdp1file.c_str());
                if( GetHDPctsUsed())
                    SetHDPctBASE( inhdpctxfile.c_str());
                if( lcval2 == "hdpsco") {
                    HDPSCORES.ReadScores( inhdpscofile.c_str());
                    HDPBASE.SetScores( &HDPSCORES );
                    if( GetHDPctsUsed()) {
                        HDPctSCORES.ReadScores( inhdpctsfile.c_str());
                        HDPctBASE.SetScores( &HDPctSCORES );
                    }
                    intscoadj = scoadjHDPsco;
                }
                else if( lcval2 == "hdpctx")
                    intscoadj = scoadjHDPCtx;
                else if( lcval1 == "hdpctx")
                    inttfrmix = tfrmixHDPCtx;
                HDPBASE.SetMixWeight( valMIXWGT );
                HDPBASE.SetNoSupClusters( valSUPCLT );
                HDPBASE.SetAdjWeight( valADJWGT );
                //
                if( GetHDPctsUsed()) {
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
        }
    } catch( myexception const& ex ) {
        error( ex.what());
        return EXIT_FAILURE;
    }
    //}}

    // -- method --

    if( !GetScoringScheme( valSCHEME, &method )) {
        error( "Unknown scoring scheme." );
        return EXIT_FAILURE;
    }
    if( method == AbstractScoreMatrix::Universal &&
        precision == AbstractScoreMatrix::FPScaling )
    {
        error( "Global scoring scheme is not compatible with floating-point precision." );
        return EXIT_FAILURE;
    }

    // -- statistical behaviour --

    if( !GetStatBehaviour( valCOMPSTATS, &behaviour )) {
        error( "Unknown statistics." );
        return EXIT_FAILURE;
    }

    // -- masking approach --

    usingmasking = ( 0.0 < valINFCON );

    if( !GetMasking( valMASKAFTER, usingmasking, &masking )) {
        error( "Unknown masking approach." );
        return EXIT_FAILURE;
    }

    // Print --

    message( NULL );

    sprintf( strbuf, "EVAL = %.1f, NOHITS = %d, NOALNS = %d", valEVAL, valNOHITS, valNOALNS );
    message( strbuf, false );

    sprintf( strbuf, "IDENTITY = %.2f", valIDENTITY );
    message( strbuf, false );

    sprintf( strbuf, "PCFWEIGHT = %.1f, MINALNFRN = %.2f, MINALNPOS = %d", valPCFWEIGHT, valMINALNFRN, valMINALNPOS );
    message( strbuf, false );


    if( 0.0 < valINFCON ) {
        sprintf( strbuf, "INFCON = %.2f, MASKAFTER = %s, SCALEDOWN = %.2f",
                    valINFCON, valMASKAFTER? "yes": "no", valSCALEDOWN );
        message( strbuf, false );
    }

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

//     sprintf( strbuf, "OPENCOST = %s%d, EXTCOST = %d, DELPROBWEIGHT = %.2f",
//                     boolAutoOpenCost? "A": "", boolAutoOpenCost? intOPENCOST: -intOPENCOST,
//                     -valEXTCOST, valDELPROBWEIGHT );
//     message( strbuf, false );
//     sprintf( strbuf, "SCHEME = %s, COMPSTATS = %s, USEGCPROBS = %s",
//                     valSCHEME.c_str(), valCOMPSTATS? "yes": "no", valUSEGCPROBS? "yes": "no" );
    sprintf( strbuf, "SCHEME = %s, MINPP = %.2f", valSCHEME.c_str(), valMINPP );
    message( strbuf, false );

//     sprintf( strbuf, "GPROBEVAL = %.1g, GPFARGWEIGHT = %.2f, GPFARGSHIFT = %.2f",
//                     valGPROBEVAL, valGPFARGWEIGHT, valGPFARGSHIFT );
//     message( strbuf, false );

    if( !valPROHIBITCOR ) {
        sprintf( strbuf, "AC1NUMER = %.1f, AC2UBNUMER = %.1f, AC2LOGSCALE = %.1f, AC2DENOMSCALE = %.2f",
                        valAC1NUMER, valAC2UBNUMER, valAC2LOGSCALE, valAC2DENOMSCALE );
        message( strbuf, false );
        sprintf( strbuf, "ANPOSCOR = %s, PROHIBITCOR = %s",
                        valANPOSCOR? "yes": "no", valPROHIBITCOR? "yes": "no" );
        message( strbuf, false );

        if( 0.0 < valINFCON2UB ) {
            sprintf( strbuf, "INFCON2UB =  %.2f, INFCON2NUMER =  %.1f, INFCON2LOGSCALE = %.1f",
                            valINFCON2UB, valINFCON2NUMER, valINFCON2LOGSCALE );
            message( strbuf, false );
            sprintf( strbuf, "INFCONALTNUMER =  %.1f, INFCONALTLOGSCALE = %.1f",
                            valINFCONALTNUMER, valINFCONALTLOGSCALE );
            message( strbuf, false );
        }
    }

    sprintf( strbuf, "HSPLEN = %d, HSPMINSCORE = %d, HSPMAXDIST = %d, NOHSPS = %d",
                    valHSPLEN, valHSPMINSCORE, valHSPMAXDIST, valNOHSPS );
    message( strbuf );

    // --


    if( method == AbstractScoreMatrix::Universal && valANPOSCOR ) {
        warning( "Analitically computed corrections are ignored in Global scoring scheme." );
    }

    int                 ret = EXIT_SUCCESS;
    ProfileSearching*   searching = NULL;

    try {
        searching = new ProfileSearching(
                    insparamfile.c_str(),
                    input.c_str(),
                    database.c_str(),
                    output.c_str(),
                    valEVAL,
                    valNOHITS,
                    valNOALNS,
                    valIDENTITY,
                    valINFCON,
                    valSCALEDOWN,
                    fixedOPENCOST,
                    valEXTCOST,
                    valSHOW,
                    !valUSEGCPROBS,
                    method,
                    behaviour,
                    precision,
                    masking
        );

        if( !searching ) {
            error( "Not enough memory." );
            return EXIT_FAILURE;
        }

        searching->SetTarFrMix( inttfrmix );
        searching->SetScoAdjment( intscoadj );
        searching->SetHDPbase( &HDPBASE );
        if( GetHDPctsUsed())
            searching->SetHDPctbase( &HDPctBASE );

        searching->SetAutoGapCosts( boolAutoOpenCost, intOPENCOST );
        searching->SetComputeDELETEstates( valDELSTATE );
        searching->SetDeletionCoefficient( valDELPROBWEIGHT );

        searching->SetExtentMinWindow( valMINALNPOS );
        searching->SetExtentMinSeqPercentage( valMINALNFRN );
        searching->SetPseudoCountWeight( valPCFWEIGHT );

        searching->SetGapProbabFactorEvalue( valGPROBEVAL );
        searching->SetGapProbabFactorWeight( valGPFARGWEIGHT );
        searching->SetGapProbabFactorShift( valGPFARGSHIFT );

        searching->SetAutocorrectionPositional( valANPOSCOR );
        searching->SetAutoACcorrection( !valPROHIBITCOR );
        searching->SetAutocorrectionNumerator1st( valAC1NUMER );
        searching->SetAutocorrectionNumerator2nd( valAC2UBNUMER );
        searching->SetAutocorrectionLogScale( valAC2LOGSCALE );
        searching->SetAutocorrectionDenomScale( valAC2DENOMSCALE );

        searching->SetInfoCorrectionUpperBound2nd( valINFCON2UB );
        searching->SetInfoCorrectionNumerator2nd( valINFCON2NUMER );
        searching->SetInfoCorrectionScale2nd( valINFCON2LOGSCALE );
        searching->SetInfoCorrectionNumeratorAlt( valINFCONALTNUMER );
        searching->SetInfoCorrectionScaleAlt( valINFCONALTLOGSCALE );

        searching->SetAlnAlgorithm( alnalgo );

        searching->SetHSPLength( valHSPLEN );
        searching->SetHSPScore( valHSPMINSCORE );
        searching->SetHSPDistance( valHSPMAXDIST );
        searching->SetHSPNoHSPs( valNOHSPS );

        if( valHCFILTER )
            searching->SetHCParameters(
                valHCWINDOW,
                valHCLOWENT,
                valHCHIGHENT
            );
        if( valINVLCFILTER )
            searching->SetSegParameters(
                valLCWINDOW,
                valLCLOWENT,
                valLCHIGHENT,
                valDISTANCE
            );
        if( valLCFILTEREACH )
            searching->SetSeqSegParameters(
                valLCWINDOW,
                valLCLOWENT,
                valLCHIGHENT
            );

        searching->Run();

    } catch( myexception const& ex )
    {
        error( ex.what());
        ret = EXIT_FAILURE;
    }

    if( searching )
        delete searching;

    return ret;
}

// -------------------------------------------------------------------------
// GetScoringScheme: determines code of schoring scheme given its name
// -------------------------------------------------------------------------

bool GetScoringScheme(
    const mystring&             scorscheme,
    AbstractScoreMatrix::TType* method )
{
    mystring    description;

    if( method == NULL )
        return false;

    if( ! scorscheme.empty())
    {
        if( scorscheme == "profile" ) {
            *method = AbstractScoreMatrix::ProfileSpecific;
        }
        else
        if( scorscheme == "psLSO" ) {
            *method = AbstractScoreMatrix::ProfSpecLSO;
        }
        else
        if( scorscheme == "adjusted" ) {
            *method = AbstractScoreMatrix::AdjustedProfileSpecific;
        }
        else
        if( scorscheme == "HDPpos" ) {
            *method = AbstractScoreMatrix::HDPProfileSpecific;
        }
        else
        if( scorscheme == "HDPctx" ) {
            *method = AbstractScoreMatrix::HDP0CtxProfileSpecific;
        }
        else
        if( scorscheme == "global" )
        {
            *method = AbstractScoreMatrix::Universal;
        }
        else {
            return false;
        }
    }

    switch( *method ) {
        case AbstractScoreMatrix::ProfileSpecific:
            description = "Scoring scheme: Profile.";
            break;

        case AbstractScoreMatrix::ProfSpecLSO:
            description = "Scoring scheme: PS-LSO.";
            break;

        case AbstractScoreMatrix::AdjustedProfileSpecific:
            description = "Scoring scheme: Adjusted profile.";
            break;

        case AbstractScoreMatrix::HDPProfileSpecific:
            description = "Scoring scheme: HDP Profile.";
            break;

        case AbstractScoreMatrix::HDP0CtxProfileSpecific:
            description = "Scoring scheme: HDP0-Ctxt Profile.";
            break;

        case AbstractScoreMatrix::Universal:
            description = "Scoring scheme: Global.";
            break;

        default:
            return false;
    }

//     message( description.c_str(), false );
    return true;
}

// -------------------------------------------------------------------------
// GetStatBehaviour: determines code of statistics to use
// -------------------------------------------------------------------------

bool GetStatBehaviour(
    bool                                compstats,
    AbstractScoreMatrix::TBehaviour*    behaviour )
{
    mystring    description;

    if( behaviour == NULL )
        return false;

    if( !compstats )
        *behaviour = AbstractScoreMatrix::StatisticsGiven;

    switch( *behaviour ) {
        case AbstractScoreMatrix::StatisticsGiven:
            description = "Statistics: Non-composition-based.";
            break;

        case AbstractScoreMatrix::ComputeStatistics:
            description = "Statistics: Composition-based.";
            break;

        default:
            return false;
    }

//     message( description.c_str(), false );
    return true;
}

// -------------------------------------------------------------------------
// GetStatBehaviour: determines masking approach
// -------------------------------------------------------------------------

bool GetMasking(
    bool        maskafter,
    bool        usingmasking,
    TMask*      masking )
{
    mystring    description;

    if( masking == NULL )
        return false;

    if( maskafter )
        *masking = MaskToConsider;

    mystring    locdesc;

    switch( *masking ) {
        case Unmasked:
            locdesc = "No masking.";
            break;

        case MaskToIgnore:
            locdesc = "Instant masking.";
            break;

        case MaskToConsider:
            locdesc = "Late masking.";
            break;

        default:
            return false;
    }

    if( usingmasking ) {
        description = locdesc;
//         message( description.c_str(), false);
    }
    return true;
}

