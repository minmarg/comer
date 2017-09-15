/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __HDPsampler__
#define __HDPsampler__

#include <stdio.h>
#include <stdlib.h>

#include "rc.h"
#include "myexcept.h"
#include "ext/psl.h"
#include "ext/rng.h"
#include "lmpi/MessageDispatcher.h"
#include "libHDP/HDPbase.h"

extern  MTRng   RNGi, RNGc;
extern  MTRng   RNGsp, RNGsm, RNGsa;
extern  MTRng   RNGjl, RNGja;

extern  const double cdp_LOG_DBL_MIN;
extern  const double cdp_LOG_DBL_MAX;

static const char*  _sch_ext_last_ = ".last";
static const char*  _sch_ext_best_ = ".best";
static const char*  _sch_ext_grp_ = ".grp";
static const char*  _sch_ext_dsh_ = ".dsh";
static const char*  _sch_ext_par_ = ".par";
static const char*  _sch_ext_sum_ = ".sum";

// -------------------------------------------------------------------------
// class HDPsampler: sampler for Hierarchical Dirichlet Process
//
class HDPsampler: public HDPbase, public MessageDispatcher
{
public:
    enum THeaders {
        THeadDat = 0xfffffff7, //data header
        THeadMsg = 0xfffffff9  //message header
    };
    enum TCommands {
        TComTerminate,      //command terminate

        TComInitMCMC,       //initiliaze mixing by MCMC procedure
        TComMCMCaccpd,      //MH update accepted
        TComMCMCcomplete,   //complete listing of vectors to be migrated

        TComInitVec,        //initiliaze of sampling by multivariate probabilities
        TComValsVec,        //values of sampling by multivariate probabilities
        TComMigrVec,        //migrate vectors

        TComInitMtx,        //initiliaze of sampling by matrix probabilities
        TComValsMtx,        //values of sampling by matrix probabilities
        TComMigrTbl,        //migrate tables

        TComInitTest,       //test initiation
        TComValsTest,       //test values from slaves

        TComValsKappa0,     //value of hyperparameter kappa0
        TComValsNu0,        //value of hyperparameter nu0
        TComValsConc,       //values of concentration parameters
        TComNoComs
    };

public:
    HDPsampler();
    ~HDPsampler();

    const char* GetFilename() const { return filename_; }
    void        SetFilename( const char* value ) { filename_ = value; }

    const char* GetOutputFile() const { return outputfile_; }
    void        SetOutputFile( const char* value ) { outputfile_ = value; }

    int         GetParallelProc() const { return parproc_; }
    void        SetParallelProc( int value ) { parproc_ = value; }

    int         GetNoIterations() const { return noiters_; }
    void        SetNoIterations( int value ) { noiters_ = value; }

    int         GetNoRstrGSScans() const { return noMHRGSscans_; }
    void        SetNoRstrGSScans( int value ) { noMHRGSscans_ = value; }

    int         GetNoMHUpdates() const { return noMHups_; }
    void        SetNoMHUpdates( int value ) { noMHups_ = value; }

    bool        GetSplitProposalsOnly() const { return smprops_ == 1; }
    bool        GetMergeProposalsOnly() const { return smprops_ == 2; }
    void        SetSMProposalsOnly( int value ) { smprops_ = value; }

    bool        GetModJNMCMC() const { return jnmodsmpl_; }
    void        SetModJNMCMC( bool value ) { jnmodsmpl_ = value; }

    bool        GetMHSampleDishUnf() const { return unfdishMHups_; }
    void        SetMHSampleDishUnf( bool value ) { unfdishMHups_ = value; }

    bool        GetMHSampleVectUnf() const { return unfvectMHups_; }
    void        SetMHSampleVectUnf( bool value ) { unfvectMHups_ = value; }

    bool        GetNoTableSampling() const { return notblsmpl_; }
    void        SetNoTableSampling( bool value ) { notblsmpl_ = value; }


    double      GetTauPreset() const { return taupreset_; }
    void        SetTauPreset( double value ) { taupreset_ = value; }

    double      GetGammaPreset() const { return gammapreset_; }
    void        SetGammaPreset( double value ) { gammapreset_ = value; }

    double      GetKappa0Preset() const { return kappa0preset_; }
    void        SetKappa0Preset( double value ) { kappa0preset_ = value; }

    double      GetNu0Preset() const { return nu0preset_; }
    void        SetNu0Preset( double value ) { nu0preset_ = value; }


    virtual long    Run();

    void        Read();
    void        SaveBestAndLast( double*, int, bool last = false );
    void        Save( mystring& fname, double* = NULL, int* = NULL );

protected:
    bool        GetParsRead() const { return parsread_; }
    void        SetParsRead( bool value ) { parsread_ = value; }

    bool        MoveVector( Restaurant*, 
                  int n, int& d, Table*& from, const int v, const Pslvector* vec, 
                  int& nn, int& dd, Table*& to, bool adjust = true );
    void        AdjustDishParams( int k, const Pslvector* vec, bool add );
    void        AdjustDishParams( int k, const Table* tbl, bool add );

    bool        MoveTable( Restaurant* rest, int n, int d, Table* from, int& dd, bool adjust = true );


    //{{MPI-RELATED SECTION OF METHODS
    bool        MSSampling();
    bool        MasterProcess();
    bool        SlaveProcess();

    bool        MasterSampleHyperparameters();
    bool        MasterDistMCMCindices( int* noups );
    bool        MasterDistMCMCindicesObs( int* noups );
    bool        MasterGenVecIndex( int dd, int vv, int* vvv );
    bool        MasterGenVecIndexS( int dd, int vv, int* vvv );
    bool        MasterGenVecIndexD( int dd, int vv, int* vvv );
    bool        MasterCmdMCMCaccpd();
    bool        ProcessCmdMCMCcomplete( int* = NULL );
    bool        SplitMergeDish( double* values, int nn );

    bool        MasterProcessRestVectors( int rb, int re );
    void        MasterProcessCmdValsVec( int, int, double*, int*, int* rests, int rsz );
    bool        MasterProcessRestTables( int n, int nn );
    void        MasterProcessCmdValsMtx( int rb, int re, double* values, int* pp );

    bool        MasterCmdInitTest();
    bool        MasterCmdValsTest();

    bool        ProcessCmdNewKappa0( double newk0 );
    bool        MasterprocessSampleKappa0();
    bool        MasterprocessBcastKappa0();

    bool        ProcessCmdNewNu0( double newnu0 );
    bool        MasterprocessSampleNu0();
    bool        MasterprocessBcastNu0();

    bool        MasterprocessSampleTau();
    bool        MasterprocessSampleTauObs();
    bool        MasterprocessSampleGamma();
    bool        MasterprocessSampleGammaObs();
    bool        MasterprocessBcastConcParams();

    friend int ss_hdphyperkappa0( double x, double* lfx, void* params );
    friend int ss_hdphypernu0( double x, double* lfx, void* params );
    friend int ss_hdpconctau( double x, double* lfx, void* params );
    friend int ss_hdpconcgamma( double x, double* lfx, void* params );
    friend int ars_hdpconctau( double x, double* hx, double* hdx, void* params );
    friend int ars_hdpconctau_helper( double x, double* hx, double* hdx, void* params );
    friend int ars_hdpconcgamma( double x, double* hx, double* hdx, void* params );
    friend int ars_hdpconcgamma_helper( double x, double* hx, double* hdx, void* params );

    bool        SlaveCmdInitMCMC();

    bool        SlaveProcessSample( 
                      double*, int, int* newt, int* newk, bool logs = true );
    bool        SlaveProcessSampleFromProbs( 
                      double*, int, int* newt, int* newk, bool = true, int* = NULL, bool logs = true );

    bool        SlaveProcessCmdSampleVec( int ttt = 0, int vvv = 0 );
    bool        SlaveProcessCmdSampleVecObs( int ttt = 0, int vvv = 0 );
    bool        SlaveProcessCmdMigrVec();
    bool        SlaveProcessCmdSampleMtx( int ttt = 0 );
    bool        SlaveProcessCmdMigrTbl();
    bool        SlaveProcessCmdValsKappa0();
    bool        SlaveProcessCmdValsNu0();
    bool        SlaveProcessCmdValsConc();

    bool        SlaveCmdInitTest();

        //{{Metropolis-Hastings update methods
    bool        MHupdate( int d1, int v1, int d2, int v2 );
    bool        MHupdateS( int d1, int v1, int d2, int v2 );
    bool        MHInitLaunchState( int d1, int v1, int d2, int v2, int* pd1, int* pd2, int* pv1, int* pv2 );
    bool        MHInitLaunchStateThroughSampling( int d1, int v1, int d2, int v2, int* pd1, int* pd2, int* pv1, int* pv2 );
    bool        MHIntmRstrGSScans( int dn1, int dn2, int vn1, int vn2 );
    bool        MHIntmRstrGSSingle( int dn1, int dn2, int vn1, int vn2, 
                                    int* varray = NULL, double* = NULL );
    bool        MHTransProbRatio( int d1, int v1, int d2, int v2, 
                                  int dn1, int dn2, int vn1, int vn2, double* );
    bool        MHTransitionProb( int d1, int v1, int d2, int v2,
                                  int dn1, int dn2, int vn1, int vn2, double* );
    bool        MHCalcNoTablesOfDishes( int d1, int v1, int d2, int v2,
                                    int dn1, int dn2, int vn1, int vn2,
                                    int* md1, int* md2, int* mm, int* nots );
    bool        MHPriorProbRatio( int d1, int v1, int d2, int v2,
                                  int dn1, int dn2, int vn1, int vn2, double* );
    bool        MHLikelihoodRatio( int d1, int v1, int d2, int v2,
                                   int dn1, int dn2, int vn1, int vn2, double* );
    bool        MHLogLikelihood( int d1, int d2, double* lvlhood );
    bool        MHLogLikelihoodObs( int d1, int d2, double* lvlhood );
    bool        MHSaveSplitMergeInfo( int d1, int v1, int d2, int v2,
                                      int dn1, int dn2, int vn1, int vn2 );
        //}}

    void        TellSlavesToTerminate( bool tothrow );

    size_t      FormatMessage( char*, int header, int cmd, int novals, double* values );
    size_t      FormatTermMessage();

    size_t          GetMaxSizeOfDataMessage();
    size_t          GetMaxSizeOfResMBuffer();
    static size_t   GetMinSizeOfMessage();
    static bool     AreDataValid( const char* );

    void        InitResMBuffer();
    void        DestroyResMBuffer();
    char*       GetResMBuffer() { return resmbuffer_; }
    const char* GetResMBuffer() const { return resmbuffer_; }
    int         GetResMBufferHead() const { return GetBufferHead( GetResMBuffer()); }
    void        SetResMBufferHead( int value ) { SetBufferHead( GetResMBuffer(), value ); }
    int         GetResMBufferCmd() const { return GetBufferCmd( GetResMBuffer()); }
    void        SetResMBufferCmd( int value ) { SetBufferCmd( GetResMBuffer(), value ); }
    int         GetResMBufferNoVals() const { return GetBufferNoVals( GetResMBuffer()); }
    void        SetResMBufferNoVals( int value ) { SetBufferNoVals( GetResMBuffer(), value ); }
    double      GetResMBufferValueAt( int n ) const { return GetBufferValueAt( GetResMBuffer(), n ); }
    void        SetResMBufferValueAt( int n, double value ) { SetBufferValueAt( GetResMBuffer(), n, value ); }
    const double* GetResMBufferValues() const { return GetBufferValues( GetResMBuffer()); }
    double*     GetResMBufferValues() { return GetBufferValues( GetResMBuffer()); }
    bool        CheckResMBufferCRC() const { return CheckBufferCRC( GetResMBuffer()); }
    int         SetResMBufferCRC() { return SetBufferCRC( GetResMBuffer()); }

    void        InitMBuffer();
    void        DestroyMBuffer();
    char*       GetMBuffer() { return mbuffer_; }
    const char* GetMBuffer() const { return mbuffer_; }
    int         GetMBufferHead() const { return GetBufferHead( GetMBuffer()); }
    void        SetMBufferHead( int value ) { SetBufferHead( GetMBuffer(), value ); }
    int         GetMBufferCmd() const { return GetBufferCmd( GetMBuffer()); }
    void        SetMBufferCmd( int value ) { SetBufferCmd( GetMBuffer(), value ); }
    int         GetMBufferNoVals() const { return GetBufferNoVals( GetMBuffer()); }
    void        SetMBufferNoVals( int value ) { SetBufferNoVals( GetMBuffer(), value ); }
    double      GetMBufferValueAt( int n ) const { return GetBufferValueAt( GetMBuffer(), n ); }
    void        SetMBufferValueAt( int n, double value ) { SetBufferValueAt( GetMBuffer(), n, value ); }
    const double* GetMBufferValues() const { return GetBufferValues( GetMBuffer()); }
    double*     GetMBufferValues() { return GetBufferValues( GetMBuffer()); }
    bool        CheckMBufferCRC() const { return CheckBufferCRC( GetMBuffer()); }
    int         SetMBufferCRC() { return SetBufferCRC( GetMBuffer()); }

    int         GetBufferHead( const char* buffer ) const;
    void        SetBufferHead( char* buffer, int value );
    int         GetBufferCmd( const char* buffer ) const;
    void        SetBufferCmd( char* buffer, int value );
    int         GetBufferNoVals( const char* buffer ) const;
    void        SetBufferNoVals( char* buffer, int value );
    double      GetBufferValueAt( const char* buffer, int n ) const;
    void        SetBufferValueAt( char* buffer, int n, double value );
    const double* GetBufferValues( const char* buffer ) const;
    double*     GetBufferValues( char* buffer );
    bool        CheckBufferCRC( const char* buffer ) const;
    int         SetBufferCRC( char* buffer );
    //}}

    int*        GetLocIntBuf() { return locintbuf_; }
    int         GetSizeLocIntBuf() const { return szlocintbuf_; }
    void        ResetLocIntBuf( int value = 0 );
    void        InitLocIntBuf( int size );
    void        DestroyLocIntBuf() { if( locintbuf_ ) free( locintbuf_ ); locintbuf_ = NULL; }

private:
    void        hdmessage( int );

private:
    const char* filename_;//filename pattern of grouped frequencies and parameters
    const char* outputfile_;//name pattern of output files
    bool        parsread_;//parameters are read

    //{{sampling attributes
    int         parproc_;//parallel processing options
    int         noiters_;//number of Gibbs scans
    int         noMHRGSscans_;//no. interm. restr. Gibbs smpl. scans within a single M-H update
    int         noMHups_;//number of M-H updates within a single complete Gibbs scan
    int         smprops_;//split/merge proposals only in M-H updates 
    bool        jnmodsmpl_;//modified Jain-Neal MCMC sampling algorithm
    bool        unfdishMHups_;//sample dishes uniformly in M-H updates
    bool        unfvectMHups_;//sample dish vectors uniformly in M-H updates
    bool        notblsmpl_;//no whole-table sampling

    double      taupreset_;//tau
    double      gammapreset_;//gamma
    double      kappa0preset_;//kappa0
    double      nu0preset_;//nu0
    //}}

    //{{MPI-messaging
    char*       mbuffer_;   //messaging buffer
    char*       resmbuffer_;//reserved messaging buffer
    //}}

    //{{buffers for local calculations
    int*        locintbuf_;//local buffer of integers
    int         szlocintbuf_;//size of buffer
    //}}
};


// =========================================================================
// GetMaxSizeOfDataMessage: return maximum size of valid message
//
inline
size_t HDPsampler::GetMaxSizeOfDataMessage()
{
    size_t  size = 20 * sizeof( double );//values
    if( 0 < mMPIRingSize())
        size *= mMPIRingSize();
    size += GetMinSizeOfMessage();

    size_t  auxflds = 10 * sizeof( double );//size of auxiliary fields
    size_t  auxsize = TIMES5( sizeof( double ));//table, dish, factor, prob, org.prob
    if( GetBasin())
        auxsize *= GetBasin()->GetActualSize();
    auxsize = GetMinSizeOfMessage() +
            TIMES2( auxsize ) + auxflds;//2*size, values for dishes and tables

    return SLC_MAX( size, auxsize );
}

// -------------------------------------------------------------------------
// GetMaxSizeOfResMBuffer: return maximum size of reserved buffer
//
inline
size_t HDPsampler::GetMaxSizeOfResMBuffer()
{
    return GetMaxSizeOfDataMessage();
}

// -------------------------------------------------------------------------
// GetMinSizeOfMessage: return minimum size of valid message
//
inline
size_t HDPsampler::GetMinSizeOfMessage()
{
    return  ( sizeof( int )) +      //header
            ( sizeof( int )) +      //command
            ( sizeof( int )) +      //number of elements
            ( sizeof( int ));       //CRC
}

// -------------------------------------------------------------------------
// AreDataValid: verify wthether the data message is valid
//
inline
bool HDPsampler::AreDataValid( const char* strstream )
{
    if( !strstream )
        return false;
    if( *( int* )strstream != THeadDat )
        return false;
    return true;
}

// -------------------------------------------------------------------------
// FormatTermMessage: format terminating message
//
inline
size_t HDPsampler::FormatTermMessage()
{
    return FormatMessage( GetMBuffer(), THeadMsg, TComTerminate, 0, NULL );
}

// -------------------------------------------------------------------------
// InitLocIntBuf: initialize local buffer of integers
//
inline
void HDPsampler::InitLocIntBuf( int size )
{
    DestroyLocIntBuf();
    locintbuf_ = ( int* )malloc( size * sizeof( int ));
    if( locintbuf_ == NULL )
        throw myruntime_error("HDPsampler: InitLocIntBuf: Not enough memory.");
    memset( locintbuf_, 0, size * sizeof( int ));
    szlocintbuf_ = size;
}

// -------------------------------------------------------------------------
// ResetLocIntBuf: reset local buffer of integers
//
inline
void HDPsampler::ResetLocIntBuf( int value )
{
    int n;
    if( locintbuf_ && 0 < szlocintbuf_ )
        for( n = 0; n < szlocintbuf_; n++ )
            locintbuf_[n] = value;
}

// -------------------------------------------------------------------------
// InitMBuffer: initialization of messaging buffer
//
inline
void HDPsampler::InitMBuffer()
{
    DestroyMBuffer();
    mbuffer_ = ( char* )malloc( GetMaxSizeOfDataMessage());
    if( mbuffer_ == NULL )
        throw myruntime_error( "HDPsampler: Not enough memory." );
}

// -------------------------------------------------------------------------
// DestroyMBuffer: memory deallocation
//
inline
void HDPsampler::DestroyMBuffer()
{
    if( mbuffer_ ) {
        free( mbuffer_ );
        mbuffer_ = NULL;
    }
}

// -------------------------------------------------------------------------
// InitResMBuffer: initialization of messaging buffer
//
inline
void HDPsampler::InitResMBuffer()
{
    DestroyResMBuffer();
    resmbuffer_ = ( char* )malloc( GetMaxSizeOfResMBuffer());
    if( resmbuffer_ == NULL )
        throw myruntime_error( "HDPsampler: Not enough memory." );
}

// -------------------------------------------------------------------------
// DestroyResMBuffer: memory deallocation
//
inline
void HDPsampler::DestroyResMBuffer()
{
    if( resmbuffer_ ) {
        free( resmbuffer_ );
        resmbuffer_ = NULL;
    }
}

// =========================================================================
// GetBufferHead: get the header present in the buffer
//
inline
int HDPsampler::GetBufferHead( const char* buffer ) const
{
    const char* p = buffer;
    int     value;
    if( p == NULL )
        throw myruntime_error("HDPsampler: Memory access error.");
    value = *( int* )p;
    return value;
}

// -------------------------------------------------------------------------
// SetBufferHead: set the header in the buffer
//
inline
void HDPsampler::SetBufferHead( char* buffer, int value )
{
    char* p = buffer;
    if( p == NULL )
        throw myruntime_error("HDPsampler: Memory access error.");
    *( int* )p = value;
}

// -------------------------------------------------------------------------
// GetBufferCmd: get the command present in the buffer
//
inline
int HDPsampler::GetBufferCmd( const char* buffer ) const
{
    const char* p = buffer;
    int     value;
    if( p == NULL )
        throw myruntime_error("HDPsampler: Memory access error.");
    p += sizeof( int );//pass header
    value = *( int* )p;
    return value;
}

// -------------------------------------------------------------------------
// SetBufferCmd: set the command in the buffer
//
inline
void HDPsampler::SetBufferCmd( char* buffer, int value )
{
    char* p = buffer;
    if( p == NULL )
        throw myruntime_error("HDPsampler: Memory access error.");
    p += sizeof( int );//pass header
    *( int* )p = value;
}

// -------------------------------------------------------------------------
// GetBufferNoVals: get the number of values present in the buffer
//
inline
int HDPsampler::GetBufferNoVals( const char* buffer ) const
{
    const char* p = buffer;
    int     value;
    if( p == NULL )
        throw myruntime_error("HDPsampler: Memory access error.");
    p += sizeof( int );//pass header
    p += sizeof( int );//pass cmd
    value = *( int* )p;
    return value;
}

// -------------------------------------------------------------------------
// SetBufferNoVals: set the number of values in the buffer
//
inline
void HDPsampler::SetBufferNoVals( char* buffer, int value )
{
    char* p = buffer;
    if( p == NULL )
        throw myruntime_error("HDPsampler: Memory access error.");
    p += sizeof( int );//pass header
    p += sizeof( int );//pass cmd
    *( int* )p = value;
}

// -------------------------------------------------------------------------
// GetBufferValueAt: get the value at nth position in the buffer
//
inline
double HDPsampler::GetBufferValueAt( const char* buffer, int n ) const
{
    const char* p = buffer;
    int     novals;
    double  value;
    if( p == NULL || n < 0 )
        throw myruntime_error("HDPsampler: Memory access error.");
    p += sizeof( int );//pass header
    p += sizeof( int );//pass cmd
    novals = *( int* )p;
    if( n < 0 || novals <= n )//|| 
//       ( GetMaxSizeOfDataMessage() - GetMinSizeOfMessage() <= n * sizeof( double )))
        throw myruntime_error("HDPsampler: Memory access error.");
    p += sizeof( int );//pass number of values
    p += sizeof( double ) * n;
    value = *( double* )p;
    return value;
}

// -------------------------------------------------------------------------
// SetBufferValueAt: set the value at nth position in the buffer
//
inline
void HDPsampler::SetBufferValueAt( char* buffer, int n, double value )
{
    char*   p = buffer;
    int     novals;
    if( p == NULL || n < 0 )
        throw myruntime_error("HDPsampler: Memory access error.");
    p += sizeof( int );//pass header
    p += sizeof( int );//pass cmd
    novals = *( int* )p;
    if( n < 0 || novals <= n )//|| 
//       ( GetMaxSizeOfDataMessage() - GetMinSizeOfMessage() <= n * sizeof( double )))
        throw myruntime_error("HDPsampler: Memory access error.");
    p += sizeof( int );//pass number of values
    p += sizeof( double ) * n;
    *( double* )p = value;
}

// -------------------------------------------------------------------------
// GetBufferValues: get the values in the buffer
//
inline
const double* HDPsampler::GetBufferValues( const char* buffer ) const
{
    const char* p = buffer;
    if( p == NULL )
        throw myruntime_error("HDPsampler: Memory access error.");
    p += sizeof( int );//header
    p += sizeof( int );//cmd
    p += sizeof( int );//number of values
    return ( const double* )p;
}

// -------------------------------------------------------------------------
// GetBufferValues: get the values in the buffer
//
inline
double* HDPsampler::GetBufferValues( char* buffer )
{
    char* p = buffer;
    if( p == NULL )
        throw myruntime_error("HDPsampler: Memory access error.");
    p += sizeof( int );//header
    p += sizeof( int );//cmd
    p += sizeof( int );//number of values
    return ( double* )p;
}

// -------------------------------------------------------------------------
// CheckBufferCRC: check calculated and the one present in the buffer CRCs
//
inline
bool HDPsampler::CheckBufferCRC( const char* buffer ) const
{
    const char* p = buffer;
    int     novals;
    int     crc = 0;
    if( p == NULL )
        throw myruntime_error("HDPsampler: Memory access error.");
    crc += *( int* )p; p += sizeof( int );//header
    crc += *( int* )p; p += sizeof( int );//cmd
    crc += novals = *( int* )p; p += sizeof( int );
    if( novals < 0 )
        throw myruntime_error("HDPsampler: Memory access error.");
    p += sizeof( double ) * novals;//pass all values
    if( crc != *( int* )p )
        return false;
    return true;
}

// -------------------------------------------------------------------------
// SetBufferCRC: set calculated CRC in the buffer; return amount of bytes
//  written to buffer
//
inline
int HDPsampler::SetBufferCRC( char* buffer )
{
    char*   p = buffer;
    int     novals;
    int     crc = 0;
    if( p == NULL )
        throw myruntime_error("HDPsampler: Memory access error.");
    crc += *( int* )p; p += sizeof( int );//header
    crc += *( int* )p; p += sizeof( int );//cmd
    crc += novals = *( int* )p; p += sizeof( int );
    if( novals < 0 )
        throw myruntime_error("HDPsampler: Memory access error.");
    p += sizeof( double ) * novals;//pass all values
    *( int* )p = crc; p += sizeof( int );
    return ( int )( p - buffer );
}

// -------------------------------------------------------------------------
// hdmessage: output symbol as message
//
inline
void HDPsampler::hdmessage( int emexp )
{
    static int  maxcntperline = 60;
    static int  countofmarks = 0;
    char        symb = '.';

    switch( emexp ) {
      case 0: symb = '*'; break;
      case 1: symb = '.'; break;
      case 2: symb = '+'; break;
      case 3: symb = '-'; break;
      case 4: symb = '*'; break;
      case 5: symb = '='; break;
    }

    if( emexp < 0 ) {
        message( NULL );
        countofmarks = 0;
        return;
    }

    if( maxcntperline < countofmarks ) {
        message( NULL );
        countofmarks = 0;
    }

    if( !countofmarks )
        messagec( 32, true );

    messagec( symb );
    countofmarks++;
}


#endif//__HDPsampler__
