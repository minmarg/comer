/***************************************************************************
 *   Copyright (C) 2010 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __MOptions__
#define __MOptions__

#include "debug.h"
#include "mystring.h"
#include "myexcept.h"


//global OPTIONS object
class MOptions;
extern MOptions OPTIONS;

// _________________________________________________________________________
// Class MOptions
//
class MOptions
{
public:
    MOptions();
    MOptions( const char* fullname );
    ~MOptions();

    void            Read(); //read options from file
    void            Read( const char* filename );

    double          GetEVAL() const { return valEVAL_; }
    int             GetNOHITS() const { return valNOHITS_; }
    int             GetNOALNS() const { return valNOALNS_; }
    int             GetSHOW() const { return valSHOW_; }

    double          GetIDENTITY() const { return ( double )valIDENTITY_ / 100.0; }
    int             GetDELSTATE() const { return valDELSTATE_; }

    double          GetPCFWEIGHT() const { return valPCFWEIGHT_; }
    double          GetMINALNFRN() const { return ( double )valMINALNFRN_ / 100.0; }
    int             GetMINALNPOS() const { return valMINALNPOS_; }
    mystring        GetSUBMAT() const { return valSUBMAT_; }

    mystring        GetTFRMIX() const { return valTFRMIX_; }
    double          GetMIXWGT() const { return valMIXWGT_; }
    mystring        GetSCOADJ() const { return valSCOADJ_; }
    int             GetSUPCLT() const { return valSUPCLT_; }
    double          GetADJWGT() const { return valADJWGT_; }
    double          GetcADJWGT() const { return valcADJWGT_; }

    double          GetCVSWGT() const { return valCVSWGT_; }

    double          GetSSSWGT() const { return valSSSWGT_; }
    double          GetSSSHDP() const { return valSSSHDP_; }

    int             GetSSEMODEL() const { return valSSEMODEL_; }

    double          GetINFCON() const { return valINFCON_; }
    int             GetMASKAFTER() const { return valMASKAFTER_; }
    double          GetSCALEDOWN() const { return 1.0 - ( double )valSCALEDOWN_ / 100.0; }

    int             GetHCFILTER() const { return valHCFILTER_; }
    int             GetHCWINDOW() const { return valHCWINDOW_; }
    double          GetHCLOWENT() const { return valHCLOWENT_; }
    double          GetHCHIGHENT() const { return valHCHIGHENT_; }

    int             GetINVLCFILTER() const { return valINVLCFILTER_; }
    int             GetLCFILTEREACH() const { return valLCFILTEREACH_; }
    int             GetLCWINDOW() const { return valLCWINDOW_; }
    double          GetLCLOWENT() const { return valLCLOWENT_; }
    double          GetLCHIGHENT() const { return valLCHIGHENT_; }
    double          GetDISTANCE() const { return valDISTANCE_; }

    void            GetOPENCOST( int* cost, bool* cauto ) {
                        TranslateOPENCOST();
                        if( cost ) *cost = intOPENCOST_;
                        if( cauto ) *cauto = boolAutoOpenCost_;
                    }
    int             GetEXTCOST() const { return valEXTCOST_; }
    double          GetDELPROBWEIGHT() const { return valDELPROBWEIGHT_; }
    mystring        GetSCHEME() const { return valSCHEME_; }
    double          GetMINPP() const { return valMINPP_; }
    int             GetCOMPSTATS() const { return valCOMPSTATS_; }
    int             GetUSEGCPROBS() const { return valUSEGCPROBS_; }

    double          GetGPROBEVAL() const { return valGPROBEVAL_; }
    double          GetGPFARGWEIGHT() const { return valGPFARGWEIGHT_; }
    double          GetGPFARGSHIFT() const { return valGPFARGSHIFT_; }

    double          GetAC1NUMER() const { return valAC1NUMER_; }
    double          GetAC2UBNUMER() const { return valAC2UBNUMER_; }
    double          GetAC2LOGSCALE() const { return valAC2LOGSCALE_; }
    double          GetAC2DENOMSCALE() const { return valAC2DENOMSCALE_; }
    int             GetANPOSCOR() const { return valANPOSCOR_; }
    int             GetPROHIBITCOR() const { return valPROHIBITCOR_; }

    double          GetINFCON2UB() const { return valINFCON2UB_; }
    double          GetINFCON2NUMER() const { return valINFCON2NUMER_; }
    double          GetINFCON2LOGSCALE() const { return valINFCON2LOGSCALE_; }
    double          GetINFCONALTNUMER() const { return valINFCONALTNUMER_; }
    double          GetINFCONALTLOGSCALE() const { return valINFCONALTLOGSCALE_; }

    int             GetHSPLEN() const { return valHSPLEN_; }
    int             GetHSPMINSCORE() const { return valHSPMINSCORE_; }
    int             GetHSPMAXDIST() const { return valHSPMAXDIST_; }
    int             GetNOHSPS() const { return valNOHSPS_; }

    const char*     GetFilename() const                 { return filename_; }
    void            SetFilename( const char* name )     { filename_ = name; }

protected:
    void            Init();
    void            TranslateOPENCOST();

    void            ReadEVAL();
    void            ReadNOHITS();
    void            ReadNOALNS();
    void            ReadSHOW();

    void            ReadIDENTITY();
    void            ReadDELSTATE();

    void            ReadPCFWEIGHT();
    void            ReadMINALNFRN();
    void            ReadMINALNPOS();
    void            ReadSUBMAT();

    void            ReadTFRMIX();
    void            ReadMIXWGT();
    void            ReadSCOADJ();
    void            ReadSUPCLT();
    void            ReadADJWGT();
    void            ReadcADJWGT();

    void            ReadCVSWGT();

    void            ReadSSSWGT();
    void            ReadSSSHDP();

    void            ReadSSEMODEL();

    void            ReadINFCON();
    void            ReadMASKAFTER();
    void            ReadSCALEDOWN();

    void            ReadHCFILTER();
    void            ReadHCWINDOW();
    void            ReadHCLOWENT();
    void            ReadHCHIGHENT();

    void            ReadINVLCFILTER();
    void            ReadLCFILTEREACH();
    void            ReadLCWINDOW();
    void            ReadLCLOWENT();
    void            ReadLCHIGHENT();
    void            ReadDISTANCE();

    void            ReadOPENCOST();
    void            ReadEXTCOST();
    void            ReadDELPROBWEIGHT();
    void            ReadSCHEME();
    void            ReadMINPP();
    void            ReadCOMPSTATS();
    void            ReadUSEGCPROBS();

    void            ReadGPROBEVAL();
    void            ReadGPFARGWEIGHT();
    void            ReadGPFARGSHIFT();

    void            ReadAC1NUMER();
    void            ReadAC2UBNUMER();
    void            ReadAC2LOGSCALE();
    void            ReadAC2DENOMSCALE();
    void            ReadANPOSCOR();
    void            ReadPROHIBITCOR();

    void            ReadINFCON2UB();
    void            ReadINFCON2NUMER();
    void            ReadINFCON2LOGSCALE();
    void            ReadINFCONALTNUMER();
    void            ReadINFCONALTLOGSCALE();

    void            ReadHSPLEN();
    void            ReadHSPMINSCORE();
    void            ReadHSPMAXDIST();
    void            ReadNOHSPS();

private:
    const char* filename_;

    double      valEVAL_;
    int         valNOHITS_;
    int         valNOALNS_;
    int         valSHOW_;

    int         valIDENTITY_;
    int         valDELSTATE_;

    double      valPCFWEIGHT_;
    int         valMINALNFRN_;
    int         valMINALNPOS_;
    mystring    valSUBMAT_;

    mystring    valTFRMIX_;
    double      valMIXWGT_;
    mystring    valSCOADJ_;
    int         valSUPCLT_;
    double      valADJWGT_;
    double      valcADJWGT_;

    double      valCVSWGT_;

    double      valSSSWGT_;
    double      valSSSHDP_;

    int         valSSEMODEL_;

    double      valINFCON_;
    int         valMASKAFTER_;
    int         valSCALEDOWN_;

    int         valHCFILTER_;
    int         valHCWINDOW_;
    double      valHCLOWENT_;
    double      valHCHIGHENT_;

    int         valINVLCFILTER_;
    int         valLCFILTEREACH_;
    int         valLCWINDOW_;
    double      valLCLOWENT_;
    double      valLCHIGHENT_;
    double      valDISTANCE_;

    mystring    valOPENCOST_;
    int         intOPENCOST_;
    bool        boolAutoOpenCost_;
    int         valEXTCOST_;
    double      valDELPROBWEIGHT_;
    mystring    valSCHEME_;
    double      valMINPP_;
    int         valCOMPSTATS_;
    int         valUSEGCPROBS_;

    double      valGPROBEVAL_;
    double      valGPFARGWEIGHT_;
    double      valGPFARGSHIFT_;

    double      valAC1NUMER_;
    double      valAC2UBNUMER_;
    double      valAC2LOGSCALE_;
    double      valAC2DENOMSCALE_;
    int         valANPOSCOR_;
    int         valPROHIBITCOR_;

    double      valINFCON2UB_;
    double      valINFCON2NUMER_;
    double      valINFCON2LOGSCALE_;
    double      valINFCONALTNUMER_;
    double      valINFCONALTLOGSCALE_;

    int         valHSPLEN_;
    int         valHSPMINSCORE_;
    int         valHSPMAXDIST_;
    int         valNOHSPS_;

};

// INLINES
// -------------------------------------------------------------------------

// DEFINES
// -------------------------------------------------------------------------

#define defEVAL ( 10.0 )
#define defNOHITS ( 700 )
#define defNOALNS ( 700 )
#define defSHOW ( 1 )

#define defIDENTITY ( 90 )
#define defDELSTATE ( 1 )

#define defPCFWEIGHT ( 1.5 )
#define defMINALNFRN ( 5 )
#define defMINALNPOS ( 10 )
#define defSUBMAT ( "Gonnet" )

#define defTFRMIX ( "no" )
#define defMIXWGT ( 0.1 )
#define defSCOADJ ( "hdpsco" )
#define defSUPCLT ( 5 )
#define defADJWGT ( 0.33 )
#define defcADJWGT ( 0.33 )

#define defCVSWGT ( 0.15 )

#define defSSSWGT ( 0.12 )
#define defSSSHDP ( 0.01 )


#define defSSEMODEL ( 1 )


#define defINFCON ( 0 )
#define defMASKAFTER ( 0 )
#define defSCALEDOWN ( 50 )

#define defHCFILTER ( 0 )
#define defHCWINDOW ( 12 )
#define defHCLOWENT ( 3.3 )
#define defHCHIGHENT ( 3.4 )

#define defINVLCFILTER ( 0 )
#define defLCFILTEREACH ( 1 )
#define defLCWINDOW ( 8 )
#define defLCLOWENT ( 1.6 )
#define defLCHIGHENT ( 1.6 )
#define defDISTANCE ( 12.96 )

#define defOPENCOST ( "A4" )
#define defintOPENCOST ( 4 )
#define defboolAutoOpenCost ( true )
#define defEXTCOST ( 1 )
#define defDELPROBWEIGHT ( 0.6 )
#define defSCHEME ( "psLSO" )
#define defMINPP ( 0.3 )
#define defCOMPSTATS ( 0/*0*/ )
#define defUSEGCPROBS ( 1 )

#define defGPROBEVAL ( 1.0e-5 )
#define defGPFARGWEIGHT ( 0.4 )
#define defGPFARGSHIFT ( 0.0 )

#define defAC1NUMER ( 4.0 )
#define defAC2UBNUMER ( 5.0 )
#define defAC2LOGSCALE ( 14.0 )
#define defAC2DENOMSCALE ( 0.12 )
#define defANPOSCOR ( 0 )
#define defPROHIBITCOR ( 1 )

#define defINFCON2UB ( 0 )
#define defINFCON2NUMER ( 4.4 )
#define defINFCON2LOGSCALE ( 4.0 )
#define defINFCONALTNUMER ( 1.0 )
#define defINFCONALTLOGSCALE ( 3.0 )

#define defHSPLEN ( 3 )
#define defHSPMINSCORE ( 0 )
#define defHSPMAXDIST ( 60 )
#define defNOHSPS ( 3 )

#define MIN_AUTOCORR_WINDOW_SIZE (  1 )
#define MAX_AUTOCORR_WINDOW_SIZE ( 50 )


// CONSTANTS
// -------------------------------------------------------------------------

static const char*  SECOPTIONS = "OPTIONS";

static const char*  EVAL = "EVAL";
static const char*  NOHITS = "NOHITS";
static const char*  NOALNS = "NOALNS";
static const char*  SHOW = "SHOW";

static const char*  IDENTITY = "IDENTITY";
static const char*  DELSTATE = "DELSTATE";

static const char*  PCFWEIGHT = "PCFWEIGHT";
static const char*  MINALNFRN = "MINALNFRN";
static const char*  MINALNPOS = "MINALNPOS";
static const char*  SUBMAT = "SUBMAT";

static const char*  TFRMIX = "TFRMIX";
static const char*  MIXWGT = "MIXWGT";
static const char*  SCOADJ = "SCOADJ";
static const char*  SUPCLT = "SUPCLT";
static const char*  ADJWGT = "ADJWGT";
static const char*  cADJWGT = "cADJWGT";

static const char*  CVSWGT = "CVSWGT";

static const char*  SSSWGT = "SSSWGT";
static const char*  SSSHDP = "SSSHDP";

static const char*  SSEMODEL = "SSEMODEL";

static const char*  INFCON = "INFCON";
static const char*  MASKAFTER = "MASKAFTER";
static const char*  SCALEDOWN = "SCALEDOWN";

static const char*  HCFILTER = "HCFILTER";
static const char*  HCWINDOW = "HCWINDOW";
static const char*  HCLOWENT = "HCLOWENT";
static const char*  HCHIGHENT = "HCHIGHENT";

static const char*  INVLCFILTER = "INVLCFILTER";
static const char*  LCFILTEREACH = "LCFILTEREACH";
static const char*  LCWINDOW = "LCWINDOW";
static const char*  LCLOWENT = "LCLOWENT";
static const char*  LCHIGHENT = "LCHIGHENT";
static const char*  DISTANCE = "DISTANCE";

static const char*  OPENCOST = "OPENCOST";
static const char*  EXTCOST = "EXTCOST";
static const char*  DELPROBWEIGHT = "DELPROBWEIGHT";
static const char*  SCHEME = "SCHEME";
static const char*  MINPP = "MINPP";
static const char*  COMPSTATS = "COMPSTATS";
static const char*  USEGCPROBS = "USEGCPROBS";

static const char*  GPROBEVAL = "GPROBEVAL";
static const char*  GPFARGWEIGHT = "GPFARGWEIGHT";
static const char*  GPFARGSHIFT = "GPFARGSHIFT";

static const char*  AC1NUMER = "AC1NUMER";
static const char*  AC2UBNUMER = "AC2UBNUMER";
static const char*  AC2LOGSCALE = "AC2LOGSCALE";
static const char*  AC2DENOMSCALE = "AC2DENOMSCALE";
static const char*  ANPOSCOR = "ANPOSCOR";
static const char*  PROHIBITCOR = "PROHIBITCOR";

static const char*  INFCON2UB = "INFCON2UB";
static const char*  INFCON2NUMER = "INFCON2NUMER";
static const char*  INFCON2LOGSCALE = "INFCON2LOGSCALE";
static const char*  INFCONALTNUMER = "INFCONALTNUMER";
static const char*  INFCONALTLOGSCALE = "INFCONALTLOGSCALE";

static const char*  HSPLEN = "HSPLEN";
static const char*  HSPMINSCORE = "HSPMINSCORE";
static const char*  HSPMAXDIST = "HSPMAXDIST";
static const char*  NOHSPS = "NOHSPS";


#endif//__MOptions__
