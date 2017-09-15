/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "libpro/srcpro/Serializer.h"
#include "libpro/srcpro/GapScheme.h"
#include "ProFreqWriter.h"

// -------------------------------------------------------------------------
// static variables
//
size_t ProFreqWriter::s_defgrsize = 200;//default group size
size_t ProFreqWriter::s_defnogrps = 20000;//default number of groups

// -------------------------------------------------------------------------
// CLASS ProFreqWriter
//
// Constructor
//
ProFreqWriter::ProFreqWriter( const char* name )
:
    rng_(),
    outname_( name ),
    data_dir_( NULL ),
    files_( NULL ),
    no_posits_( 0 ),
    no_arguments_( 0 ),
    no_processedfiles_( 0 ),
    mineffthickn_( 0 ),
    no_samples_( 0 ),
    actual_nosmpls_( 0 ),
    context_size_( 0 ),
    ctxt_step_( 0 ),
    mix_( false ),
    cwght_( 0.0 ),
    dishpergroup_( false ),
    seed_( 0 ),
    coeffs_( NULL ),
    samplecounter_( 0 ),
    rndnumbers_( NULL ),
    samples_( NULL )
{
}

// Constructor
//
ProFreqWriter::ProFreqWriter( const char* output, const char* directory )
:
    rng_(),
    outname_( output ),
    data_dir_( directory ),
    files_( NULL ),
    no_posits_( 0 ),
    no_arguments_( 0 ),
    no_processedfiles_( 0 ),
    mineffthickn_( 0 ),
    no_samples_( 0 ),
    actual_nosmpls_( 0 ),
    context_size_( 0 ),
    ctxt_step_( 0 ),
    mix_( false ),
    cwght_( 0.0 ),
    dishpergroup_( false ),
    seed_( 0 ),
    coeffs_( NULL ),
    samplecounter_( 0 ),
    rndnumbers_( NULL ),
    samples_( NULL )
{
}

// Constructor
//
ProFreqWriter::ProFreqWriter( const char* output, char* arguments[], int no_args )
:
    rng_(),
    outname_( output ),
    data_dir_( NULL ),
    files_( arguments ),
    no_posits_( 0 ),
    no_arguments_( no_args ),
    no_processedfiles_( 0 ),
    mineffthickn_( 0 ),
    no_samples_( 0 ),
    actual_nosmpls_( 0 ),
    context_size_( 0 ),
    ctxt_step_( 0 ),
    mix_( false ),
    cwght_( 0.0 ),
    dishpergroup_( false ),
    seed_( 0 ),
    coeffs_( NULL ),
    samplecounter_( 0 ),
    rndnumbers_( NULL ),
    samples_( NULL )
{
}

// Default constructor
//
ProFreqWriter::ProFreqWriter()
:
    rng_(),
    outname_( NULL ),
    data_dir_( NULL ),
    files_( NULL ),
    no_posits_( 0 ),
    no_arguments_( 0 ),
    no_processedfiles_( 0 ),
    mineffthickn_( 0 ),
    no_samples_( 0 ),
    actual_nosmpls_( 0 ),
    context_size_( 0 ),
    ctxt_step_( 0 ),
    mix_( false ),
    cwght_( 0.0 ),
    dishpergroup_( false ),
    seed_( 0 ),
    coeffs_( NULL ),
    samplecounter_( 0 ),
    rndnumbers_( NULL ),
    samples_( NULL )
{
    throw( myruntime_error( mystring( "ProFreqWriter: Default construction is not allowed." )));
}

// Destructor
//
ProFreqWriter::~ProFreqWriter()
{
    DestroyRndNumbers();
    DestroySamples();
    DestroyCoeffs();
}

// -------------------------------------------------------------------------
// DestroyRandNumbers: destroy array of random numbers
//
void ProFreqWriter::DestroyRndNumbers()
{
    if( rndnumbers_ ) {
        free( rndnumbers_ );
        rndnumbers_ = NULL;
    }
}

// -------------------------------------------------------------------------
// DestroySamples: destroy all samples in vector
//
void ProFreqWriter::DestroySamples()
{
    CtxtFrequencies*    frq;
    SimpleVector*       sv;
    size_t          n, m;
    if( samples_ ) {
        for( n = 0; n < samples_->GetSize(); n++ ) {
            sv = ( SimpleVector* )samples_->GetValueAt( n );
            if( sv == NULL )
                continue;
            for( m = 0; m < sv->GetSize(); m++ ) {
                frq = ( CtxtFrequencies* )sv->GetValueAt( m );
                if( frq == NULL )
                    continue;
                delete frq;
            }
            delete sv;
            sv = NULL;
        }
        delete samples_;
        samples_ = NULL;
    }
}

// -------------------------------------------------------------------------
// NewSamples: allocate required memory for new samples
//
void ProFreqWriter::NewRndNumbers()
{
    int no_auxfields = GetNoAuxFields();   //header, no. elems, crc
    if( 0 < GetNoSamples()) {
        DestroyRndNumbers();
        rndnumbers_ = ( double* )malloc( sizeof( double ) * ( GetNoSamples() + no_auxfields ));
        if( !rndnumbers_ )
            throw( myruntime_error( mystring( "ProFreqWriter: Not enough memory." )));
    }
}

// -------------------------------------------------------------------------
// GenerateRndNumbers: generate random numbers
//
void ProFreqWriter::GenerateRndNumbers()
{
    double  p;
    int     no_auxfields = GetNoAuxFields();   //header, no. elems, crc
    size_t  n = 0, i, crc = 0;

    if( rndnumbers_ && 0 < GetNoSamples()) {
        GetRng().Set( GetSeed());
        rndnumbers_[n] = ( double )RndNumHeader;  crc += ( size_t )rndnumbers_[n++];
        rndnumbers_[n] = ( double )GetNoSamples();crc += ( size_t )rndnumbers_[n++];
        for( ; n < GetNoSamples() + no_auxfields - 1; n++ ) {
            p = GetRng().GetDouble();//p in [0,1)
            rndnumbers_[n] = p;
            for( i = 0; i < sizeof( double ); i++ ) crc += (( char* )&p )[i];
        }
        rndnumbers_[n] = ( double )crc;
    }
}

// -------------------------------------------------------------------------
// VerifyRndNumbers: verify random numbers
//
void ProFreqWriter::VerifyRndNumbers()
{
    double  p;
    int     no_auxfields = GetNoAuxFields();   //header, no. elems, crc
    size_t  n = 0, i, crc = 0;

    if( rndnumbers_ && 0 < GetNoSamples()) {
        if( rndnumbers_[n] != RndNumHeader )
            throw( myruntime_error( "ProFreqWriter: VerifyRndNumbers: Wrong header." ));

        crc += ( size_t )rndnumbers_[n++];

        if( rndnumbers_[n] != GetNoSamples())
            throw( myruntime_error( "ProFreqWriter: VerifyRndNumbers: Wrong number of cells." ));

        crc += ( size_t )rndnumbers_[n++];

        for( ; n < GetNoSamples() + no_auxfields - 1; n++ ) {
            p = rndnumbers_[n];
            if( p < 0.0 || 1.0 <= p )
                throw( myruntime_error( "ProFreqWriter: VerifyRndNumbers: Wrong random number." ));

            for( i = 0; i < sizeof( double ); i++ ) crc += (( char* )&p )[i];
        }

        if( rndnumbers_[n] != crc )
            throw( myruntime_error( "ProFreqWriter: VerifyRndNumbers: Wrong CRC." ));
    }
}

// -------------------------------------------------------------------------
// MakeRndNumbers: make random numbers
//
void ProFreqWriter::MakeRndNumbers()
{
    DestroyRndNumbers();
    if( GetNoSamples() <= 0 )
        return;

    NewRndNumbers();
    GenerateRndNumbers();

//     if( mMPIIAmMaster()) {
//         //master
//         MasterMessage( "Generating random numbers..." );
//         GenerateRndNumbers();
//     }
//     BcastMPIMessage(( char* )GetRndNumbers(), GetSizeOfRndNumbers(), true/*throw on error*/ );
// 
//     if( !mMPIIAmMaster()) {
//         //slaves
//         VerifyRndNumbers();
//     }
}

// -------------------------------------------------------------------------
// NewSamples: allocate required memory for new samples
//
void ProFreqWriter::NewSamples( size_t size )
{
    samples_ = new SimpleVector( size );
    if( !samples_ )
        throw( myruntime_error( mystring( "ProFreqWriter: Not enough memory." )));
}

// -------------------------------------------------------------------------
// Run: read profiles and sample frequency vectors
// -------------------------------------------------------------------------

unsigned long ProFreqWriter::Run()
{
    char        strbuf[BUF_MAX];
    double      nosperpos = 0.0;    //number of samples per position
    mystring    throwstr;

    try {
        if( GetMix() && 1 < GetContextSize()) {
            if( GetCentralWeight() <= 0.0 || 1.0 <= GetCentralWeight())
                throw( myruntime_error("ProFreqWriter: Invalid central weight."));
            NewCoeffs( GetContextSize(), GetCentralWeight());
            GetMixCoeffs()->FindCoefficients();
        }

        ResetSampleCounter();
        message( "Counting positions...", false );
        SetNoPositions( 0 );
//         ProcessInput( &ProFreqWriter::CountFiles, GetNoFilesAddr());//obsolete
        ProcessInput( &ProFreqWriter::CountPositions, GetNoPositionsAddr());
        sprintf( strbuf, "  total positions, %d", GetNoPositions());
        message( strbuf );

        MakeRndNumbers();

        if( 0 < GetNoPositions() && 0 < GetNoSamples())
            nosperpos = (( double )GetNoSamples()) / GetNoPositions();

        if( nosperpos <= 0.0 )
            throw( myruntime_error( "ProFreqWriter: Unable to sample." ));
    
        DestroySamples();
        if( 0 < GetNoSamples())
            NewSamples( GetDefNoGroups());

        message( "Sampling..." );
        ResetActNoSamples();
        ResetNoProcessedFiles();
        ProcessInput( &ProFreqWriter::ProcessFile, &nosperpos );

        sprintf( strbuf, "Processed %d file(s); Generated %d sample(s).",
            GetNoProcessedFiles(), GetActNoSamples());
        message( strbuf );

        message( "Writing samples..." );
        Output();

    } catch( myexception const& ex )
    {
        throwstr = ex.what();
    }

    if( !throwstr.empty())
        throw myruntime_error( throwstr );

    return 0;
}

// -------------------------------------------------------------------------
// OutputParameters: print optimized parameters to file
// -------------------------------------------------------------------------

void ProFreqWriter::Output()
{
    FILE*   fp = NULL;
    const double* wghts = NULL;
    int n;

    if( !GetOutputName() || !strlen( GetOutputName()))
        throw myruntime_error( "ProFreqWriter: Unable to write: No filename." );

    fp = fopen( GetOutputName(), "w" );

    if( fp == NULL ) {
        throw myruntime_error( 
            mystring( "ProFreqWriter: Failed to open file '" ) + 
            GetOutputName() + mystring( "'." ));
    }

    fprintf( fp, "## Processed %d file(s)\n", GetNoProcessedFiles());
    if( GetDataDir())
        fprintf( fp, "## Directory: %s\n", /*my_basename( */GetDataDir());
    fprintf( fp, "## Minimum effective profile thickness, %d\n", GetMinEffThickness());
    fprintf( fp, "## Context length = %d\n", GetContextSize());
    fprintf( fp, "## Context step = %d\n", GetContextStep());
    fprintf( fp, "## Mixing context positions = %s\n", GetMix()? "true": "false" );
    if( GetMix() && 1 < GetContextSize()) {
        fprintf( fp, "## Weight of central context position = %g\n", GetCentralWeight());
        fprintf( fp, "## Coefficients:\n##  ");
        wghts = GetMixCoeffs()? GetMixCoeffs()->GetCoefficients(): NULL;
        if( wghts )
            for( n = 0; n < GetContextSize(); n++ )
                fprintf( fp, " %.2f", wghts[n]);
        fprintf( fp, "\n");
    }

    if( GetDishPerGroup())
        PrintFrequenciesDishPerGroup( fp );
    else
        PrintFrequenciesDish1( fp );
    fclose( fp );
}

// -------------------------------------------------------------------------
// PrintFrequenciesDishPerGroup: print frequencies to file with arrangement of
//  having one dish per group
//
void ProFreqWriter::PrintFrequenciesDishPerGroup( FILE* fp )
{
    SimpleVector*       group;
    CtxtFrequencies*    frq;
    size_t  nogrps, nosmps;
    size_t  n, m;

    if( fp == NULL )
        return;

    if( !GetSamples())
        return;

    const char*     patstrnog = "Number of groups =";
    const char*     patstrnds = "Number of dishes =";
    const char*     patstrtns = "Total number of samples =";
    const char*     patstrssz = "Context length =";
    const char*     patstrgrp = "Group";
    const char*     patstrnot = "Number of tables =";
    const char*     patstrtbl = "Table";
    const char*     patstrtid = "Id =";
    const char*     patstrdid = "Dish id =";
    const char*     patstrnos = "Number of samples =";
    const char*     patstrendofsmp = "*";
    const char*     patstrendoftbl = ";;";
    const char*     patstrendofgrp = "////";
    int tt, dd;
    int ctxtsz = 0;

    if( GetSamples()->GetSize()) {
        group = ( SimpleVector* )GetSamples()->GetValueAt(0);
        if( group ) {
            frq = ( CtxtFrequencies* )group->GetValueAt(0);
            if( frq )
                ctxtsz = frq->GetLength();
        }
    }

    fprintf( fp, "%s %d; %s %d; %s %d; %s %d\n", 
              patstrnog, GetSamples()->GetSize(), 
              patstrnds, GetSamples()->GetSize(),
              patstrtns, GetActNoSamples(), patstrssz, ctxtsz );

    nogrps = 0;
    for( dd = n = 0; n < GetSamples()->GetSize(); n++ ) {
        group = ( SimpleVector* )GetSamples()->GetValueAt( n );
        if( !group )
            continue;
        nosmps = group->GetSize();
        fprintf( fp, "%s %d; %s %d\n", patstrgrp, nogrps++, patstrnot, nosmps );
        for( tt = m = 0; m < nosmps; m++ ) {
            frq = ( CtxtFrequencies* )group->GetValueAt( m );
            if( !frq )
                continue;
            fprintf( fp, "%s %d (%s %d; %s %d); %s %d\n", 
                     patstrtbl, tt++, patstrtid, m, patstrdid, nogrps-1,
                     patstrnos, 1 );
            frq->WriteDoubles1( fp, n == 0 );
            fprintf( fp, "%s\n", patstrendofsmp );
            fprintf( fp, "%s\n", patstrendoftbl );//end of table
        }
        fprintf( fp, "%s\n\n", patstrendofgrp );
    }
}

// -------------------------------------------------------------------------
// PrintFrequenciesDish1: print frequencies to file with arrangement of
//  having single dish and one table per group
//
void ProFreqWriter::PrintFrequenciesDish1( FILE* fp )
{
    SimpleVector*       group;
    CtxtFrequencies*    frq;
    size_t  nogrps, nosmps;
    size_t  n, m;

    if( fp == NULL )
        return;

    if( !GetSamples())
        return;

    const char*     patstrnog = "Number of groups =";
    const char*     patstrnds = "Number of dishes =";
    const char*     patstrtns = "Total number of samples =";
    const char*     patstrssz = "Context length =";
    const char*     patstrgrp = "Group";
    const char*     patstrnot = "Number of tables =";
    const char*     patstrtbl = "Table";
    const char*     patstrtid = "Id =";
    const char*     patstrdid = "Dish id =";
    const char*     patstrnos = "Number of samples =";
    const char*     patstrendofsmp = "*";
    const char*     patstrendoftbl = ";;";
    const char*     patstrendofgrp = "////";
    int tt;
    int ctxtsz = 0;

    if( GetSamples()->GetSize()) {
        group = ( SimpleVector* )GetSamples()->GetValueAt(0);
        if( group ) {
            frq = ( CtxtFrequencies* )group->GetValueAt(0);
            if( frq )
                ctxtsz = frq->GetLength();
        }
    }

    fprintf( fp, "%s %d; %s %d; %s %d; %s %d\n", 
              patstrnog, GetSamples()->GetSize(), 
              patstrnds, 1,
              patstrtns, GetActNoSamples(), patstrssz, ctxtsz );

    nogrps = 0;
    for( n = 0; n < GetSamples()->GetSize(); n++ ) {
        group = ( SimpleVector* )GetSamples()->GetValueAt( n );
        if( !group )
            continue;
        nosmps = group->GetSize();
        fprintf( fp, "%s %d; %s %d\n", patstrgrp, nogrps++, patstrnot, 1 );
        fprintf( fp, "%s %d (%s %d; %s %d); %s %d\n", 
                  patstrtbl, 0, patstrtid, 0, patstrdid, 0,
                  patstrnos, nosmps );
        for( m = 0; m < nosmps; m++ ) {
            frq = ( CtxtFrequencies* )group->GetValueAt( m );
            if( !frq )
                continue;
            frq->WriteDoubles1( fp, n == 0 );
            fprintf( fp, "%s\n", patstrendofsmp );
        }
        fprintf( fp, "%s\n", patstrendoftbl );//end of table
        fprintf( fp, "%s\n\n", patstrendofgrp );
    }
}

// -------------------------------------------------------------------------
// ProcessInput: process profiles given explicitly or found in directory
// -------------------------------------------------------------------------

void ProFreqWriter::ProcessInput( PMETHOD processing_method, void* args )
{
    int     n = 0, m = 0;
    struct stat info;

    if( GetDataDir()) {
        dirent*         entry = NULL;
        DIR*            direct = opendir( GetDataDir());
        mystring        dirname = GetDataDir() + mystring( DIRSEP );

        if( !direct )
            throw myruntime_error( "ProFreqWriter: Unable to process directory specified." );

        while(( entry = readdir( direct ))) {

            mystring    fullname = dirname + entry->d_name;

            if( stat( fullname.c_str(), &info ) == -1 ) {
                error( mystring( mystring( "Cannot stat '" ) + 
                    mystring( entry->d_name ) + "'; Skipping." ).c_str());
                continue;
            }
            if(( S_IFMT & info.st_mode ) != S_IFREG )
//             if( !S_ISREG( info.st_mode ))
                continue;

            ( this->*processing_method )( fullname.c_str()/*entry->d_name*/, m++, args );
        }
        closedir( direct );

    } else
        if( GetFiles()) {
            bool    option = false;

            for( n = 0; n < GetNoArguments(); n++ ) {
                if( option ) {
                    option = false;
                    continue;
                }
                if( GetFileAt( n )[0] == '-' && GetFileAt( n )[2] == 0 ) {
                    option = true;  //assume it is an option
                    continue;
                }
                if( stat( GetFileAt( n ), &info ) == -1 ) {
                    error( mystring( mystring( "Cannot stat '" ) + 
                        mystring( GetFileAt( n )) + "'; Skipping." ).c_str());
                    continue;
                }
                if(( S_IFMT & info.st_mode ) != S_IFREG )
                    continue;

                ( this->*processing_method )( GetFileAt( n ), m++, args );
            }
        }
}

// -------------------------------------------------------------------------
// CountFiles: count number of profiles
//
void ProFreqWriter::CountFiles( const char*, int, void* args )
{
    int*    counter = ( int* )args;

    if( counter )
        (*counter)++;
}

// -------------------------------------------------------------------------
// CountPositions: count total number of positions of all profiles
//
void ProFreqWriter::CountPositions( const char* filename, int, void* args )
{
    int*                counter = NULL;
    FrequencyMatrix     freq;
    LogOddsMatrix       pssm;
    GapScheme           gaps;
    Serializer          serializer;

    if( !args )
        throw myruntime_error( "ProFreqWriter: CountPositions: Argument expected." );

    counter = ( int* )args;

    try {
//         serializer.DeserializeProfile( freq, pssm, gaps, filename );//obsolete
        serializer.ReadProfile( filename, freq, pssm, gaps );

//         freq.CheckForAllZeros();    //necessary to ensure valid scores

        if( pssm.GetEffNoSequences() < GetMinEffThickness())
            return;

        if( counter )
            (*counter) += pssm.GetColumns();

    } catch( myexception const& ex )
    {
        error( mystring( ex.what() + 
            mystring( " Skipping '" ) + filename + mystring( "'." )).c_str());
    }
}

// -------------------------------------------------------------------------
// ProcessFile: read profile and sample context windows
// -------------------------------------------------------------------------

void ProFreqWriter::ProcessFile( const char* filename, int serial, void* args )
{
    FrequencyMatrix     freq;
    LogOddsMatrix       pssm;
    GapScheme           gaps;
    double              nosperpos;
    int                 nosperfile;
    bool*               mask = NULL;    //mask of used positions
    SimpleVector*       group = NULL;
    CtxtFrequencies*    sample = NULL;
    Serializer          serializer;
    mystring            errstr;
    int                 nopf = GetNoProcessedFiles();
    int                 ctxtsz = GetContextSize();
    int                 step = GetContextStep();
    int                 cextent = ctxtsz + (ctxtsz-1)*(step-1);//context extent
    double      prob = 0.0;
//     const int   centre = GetContextSize() / 2 - ( GetContextSize() % 2 == 0 );
    int         position, n, c, noiters;
    bool        mark;

    if( ctxtsz < 1 )
        throw myruntime_error("ProFreqWriter: ProcessFile: Invalid context length.");

    if( step < 1 || 10 < step || ctxtsz < step )
        throw myruntime_error("ProFreqWriter: ProcessFile: Invalid step.");

    if( !args )
        throw myruntime_error("ProFreqWriter: ProcessFile: Argument expected.");

    group = new SimpleVector( GetDefGroupSize());
    if( !group )
        throw( myruntime_error("ProFreqWriter: ProcessFile: Not enough memory."));

    nosperpos = *( double* )args;
    if( nosperpos <= 0.0 )
        throw myruntime_error("ProFreqWriter: ProcessFile: Wrong argument.");

    if( 1.0 < nosperpos )
        nosperpos = 1.0;

    try {
//         serializer.DeserializeProfile( freq, pssm, gaps, filename );//obsolete
        serializer.ReadProfile( filename, freq, pssm, gaps );

//         freq.CheckForAllZeros();    //necessary to ensure valid scores

        if( pssm.GetEffNoSequences() < GetMinEffThickness())
            return;

        if( pssm.GetColumns() < cextent )
            return;

        nosperfile = ( int )rint( nosperpos * pssm.GetColumns());

        if( nosperfile == 0 ) {
            nosperfile = ( int )rint( 1.0 /( nosperpos * pssm.GetColumns()));
            if( serial % nosperfile )
                return;

            prob = GetNextRndNumber();
            if( prob <= 0.0 || 1.0 <= prob )
                return;

            position = ( int )rint( prob * ( pssm.GetColumns() - 1 ));
//             position -= centre;

            if( pssm.GetColumns() < position + cextent )
                return;

            Sample( filename, freq, pssm, position, &sample );
            if( sample == NULL )
                return;

            group->Push( sample );
            IncActNoSamples();
            GetSamples()->Push( group );
            IncNoProcessedFiles();
        }
        else {
            //initialize mask
            mask = ( bool* )malloc( sizeof( bool ) * pssm.GetColumns());
            if( !mask )
                throw myruntime_error( "ProFreqWriter: Not enough memory.", CRITICAL );
            mark = 0;
            noiters = nosperfile;
            if( 0.5 < nosperpos ) {
                mark = 1;
                noiters = pssm.GetColumns() - nosperfile;
            }
            for( n = 0; n < pssm.GetColumns(); n++ )
                mask[n] = mark;

            // select randomly positions
            for( n = 0; n < noiters; n++ ) {
                prob = GetNextRndNumber();
                if( prob <= 0.0 || 1.0 <= prob )
                    continue;

                position = ( int )rint( prob * ( pssm.GetColumns() - 1 ));
//                 position -= centre;
                mask[position] = !mark;
            }

            // mask the end of sequence to avoid duplicates of samples
//             for( n = pssm.GetColumns() - GetContextSize() + 1; n < pssm.GetColumns(); n+=step )
//                 if( 0 <= n )
//                     mask[n] = 0;

            // make samples for the selected positions
            for( n = 0; n < pssm.GetColumns(); n++ ) {
                if( pssm.GetColumns() < n + cextent )
                    break;

                if( !mask[n] )
                    continue;

                Sample( filename, freq, pssm, n, &sample );
                if( sample == NULL )
                    continue;
                //mark over all context length
                for( c = 0; c < 1/*cextent*/ && n+c < pssm.GetColumns(); c+=step )
                    mask[n+c] = 0;

                group->Push( sample );
                IncActNoSamples();
            }
            if( group->GetSize()) {
                GetSamples()->Push( group );
                IncNoProcessedFiles();
            }
        }

    } catch( myexception const& ex )
    {
        if( ex.eclass() == CRITICAL )
            errstr = ex.what();
        else
            error( mystring( ex.what() + 
                mystring( " Skipping '" ) + filename + mystring( "'." )).c_str());
    }

    if( group && GetNoProcessedFiles() <= nopf ) {
        delete group;
        group = NULL;
    }
    if( mask ) {
        free( mask );
        mask = NULL;
    }
    if( !errstr.empty())
        throw myruntime_error( errstr, CRITICAL );
}

// -------------------------------------------------------------------------
// MixSample: mix sample 
//
void ProFreqWriter::MixSample( CtxtFrequencies** p_sample, bool tonormal )
{
    const char* serr;
    mystring errstr;

    const CtxtCoefficients* cfs = GetMixCoeffs();
    int ctxlen = cfs->GetLength();
    int n, r;

    if( p_sample == NULL || *p_sample == NULL )
        return;

    if( cfs == NULL || cfs->GetCoefficients() == NULL ||
      (*p_sample)->GetLength() != cfs->GetLength())
        throw myruntime_error("ProFreqWriter: Sample: Inconsistent coefficients.");

    double val;
    const double* coeffs = cfs->GetCoefficients();
    CtxtFrequencies* sample = new CtxtFrequencies(1);

    if( !sample )
        throw myruntime_error( "ProFreqWriter: Not enough memory." );

    try {
        for( n = 0; n < ctxlen; n++ ) {
            for( r = 0; r < (*p_sample)->GetEffectiveNoResids(); r++ ) {
                val = (*p_sample)->GetObsFrequency( n, r ) * coeffs[n];
                sample->AddObsFrequency( 0, r, val );
            }
        }
        delete *p_sample;
        *p_sample = NULL;

        if( tonormal )
            if(( serr = SampleToNormal( sample, 0 )) != NULL ) {
                delete sample;
                sample = NULL;
                throw myruntime_error( serr );
                return;
            }
        *p_sample = sample;
        sample = NULL;

    } catch( myexception const& ex ) {
        errstr = ex.what();
    }

    if( sample ) { 
        delete sample; 
        sample = NULL;
    }

    if( !errstr.empty())
        throw myruntime_error( errstr, CRITICAL );
}

// -------------------------------------------------------------------------
// MixSample: mix sample by making a matrix normal variate (A x 3) 
//
void ProFreqWriter::MixSample3( CtxtFrequencies** p_sample, bool tonormal )
{
    const char* serr;
    mystring errstr;

    const CtxtCoefficients* cfs = GetMixCoeffs();
    int ctxlen = cfs->GetLength();
    int cpos = ctxlen >> 1;
    int n, r;

    if( p_sample == NULL || *p_sample == NULL || cpos < 1)
        return;

    if( cfs == NULL || cfs->GetCoefficients() == NULL ||
      (*p_sample)->GetLength() != cfs->GetLength())
        throw myruntime_error("ProFreqWriter: Sample: Inconsistent coefficients.");

    double val, w;
    const double* coeffs = cfs->GetCoefficients();
    const double cwg = coeffs[cpos];
    const int nobs = 3;
    CtxtFrequencies* sample = new CtxtFrequencies( nobs );

    if( !sample )
        throw myruntime_error( "ProFreqWriter: Not enough memory." );

    try {
        w = cwg /( double )( cpos );
        //mix left flang
        for( n = 0; n < cpos; n++ ) {
            for( r = 0; r < (*p_sample)->GetEffectiveNoResids(); r++ ) {
                val = (*p_sample)->GetObsFrequency( n, r ) * ( TIMES2(coeffs[n]) + w );
                sample->AddObsFrequency( 0, r, val );
            }
        }
        //rewrite central vector
        for( r = 0; r < (*p_sample)->GetEffectiveNoResids(); r++ ) {
            val = (*p_sample)->GetObsFrequency( n, r );
            sample->AddObsFrequency( 1, r, val );
        }
        //mix right flang
        for( n++; n < ctxlen; n++ ) {
            for( r = 0; r < (*p_sample)->GetEffectiveNoResids(); r++ ) {
                val = (*p_sample)->GetObsFrequency( n, r ) * ( TIMES2(coeffs[n]) + w );
                sample->AddObsFrequency( 2, r, val );
            }
        }
        delete *p_sample;
        *p_sample = NULL;

        if( tonormal )
            for( n = 0; n < nobs; n++ ) {
                if(( serr = SampleToNormal( sample, n )) != NULL ) {
                    delete sample;
                    sample = NULL;
                    throw myruntime_error( serr );
                    return;
                }
            }
        *p_sample = sample;
        sample = NULL;

    } catch( myexception const& ex ) {
        errstr = ex.what();
    }

    if( sample ) { 
        delete sample; 
        sample = NULL;
    }

    if( !errstr.empty())
        throw myruntime_error( errstr, CRITICAL );
}

// -------------------------------------------------------------------------
// Sample: generate a sample of observed frequencies
//
void ProFreqWriter::Sample( 
    const char* filename,
    FrequencyMatrix& freq, LogOddsMatrix& pssm, int position, 
    CtxtFrequencies** p_sample )
{
    int ctxtsz = GetContextSize();
    bool mix = GetMix();
    bool tonormal = true;//ctxtsz < 2 || !mix;
    char bufstr[KBYTE];
    mystring errstr;

    if( p_sample == NULL )
        return;

    SampleHlp( filename, freq, pssm, position, p_sample, tonormal );

    if( *p_sample == NULL || ctxtsz < 2 || !mix )
        return;

    try {
        MixSample( p_sample, !tonormal );
//         MixSample3( p_sample, !tonormal );

    } catch( myexception const& ex ) {
        sprintf( bufstr, "%s; pos. %d, file %s.",
                ex.what(), position, my_basename( filename ));
        errstr = bufstr;
    }

    if( !errstr.empty())
        throw myruntime_error( errstr, CRITICAL );
}

// -------------------------------------------------------------------------
// SampleHlp: helper method to generate a sample of observed frequencies
//
void ProFreqWriter::SampleHlp( 
    const char* filename,
    FrequencyMatrix& freq, LogOddsMatrix& pssm, int position, 
    CtxtFrequencies** p_sample, bool tonormal )
{
    const bool      bCHACC = false;//check accuracy
    int             ctxtsz = GetContextSize();
    int             step = GetContextStep();
    int             cextent = ctxtsz + (ctxtsz-1)*(step-1);//context extent
    const size_t    cmineffth = GetMinEffThickness();
    const size_t    effnoress = NUMAA;
    char            residue;
    double          f, fsum;
    double          effth;
    int             n, r;
    char            bufstr[KBYTE];
    const char*     serr;
    mystring        errstr;
    const double    accuracy = 1.0e-1;
    const double    highacc = 1.0e-6;
    const double    lambda = LOSCORES.StatisParam( Ungapped, Lambda );

    if( !p_sample )
        return;
    *p_sample = NULL;

    if( pssm.GetColumns() < position + cextent || ctxtsz < 1 || step < 1 || ctxtsz < step )
        return;

    if( pssm.GetColumns() < position + cextent )
        position = pssm.GetColumns() - cextent;
    if( position < 0 )
        position = 0;

    CtxtFrequencies*  sample = new CtxtFrequencies( GetContextSize());

    if( !sample )
        throw myruntime_error( "ProFreqWriter: Not enough memory." );

    try {
        for( n = 0; n < ctxtsz; n++ )
        {
            effth = pssm.GetEffNoSequences();
            residue = pssm.GetResidueAt( position + n*step );
            //avoid sampling if any of the positions in a context holds X
            //assuming that the positions are filtered with SEG
            if( residue == X ) {
                delete sample;
                sample = NULL;
                return;
            }
            //skip sampling if any of the positions has small eff. thicknes
            if( effth < cmineffth ) {
                delete sample;
                sample = NULL;
                return;
            }
            sample->SetResidue( n, residue );
            fsum = 0.0;
            for( r = 0; r < effnoress; r++ ) {
                //obsolete
                //sample->SetObsFrequency( n, r, rint( effth * pssm.GetValueAt( position + n, r )));
                f = pssm.GetValueAt( position + n*step, r );
                f = exp( f * lambda ) * LOSCORES.PROBABility( r );
                if( f < 0.0 ) {
                    sprintf( bufstr, "Negative frequency (%g) at pos. %d, file %s; skipped.",
                            f, position + n, my_basename( filename ));
                    warning( bufstr );
                    delete sample;
                    sample = NULL;
                    return;
                }
                fsum += f;
                sample->SetObsFrequency( n, r, f );
            }
            if( fsum < 1.0 - accuracy || 1.0 + accuracy < fsum ) {
                sprintf( bufstr, "Invalid frequencies (%g) at pos. %d, file %s; skipped.",
                        fsum, position + n, my_basename( filename ));
                warning( bufstr );
                delete sample;
                sample = NULL;
                return;
            }
            //ERROR CORRECTION...
            if( bCHACC ) {
                SampleErrCorrection( sample, n );
                //2nd round CHECK...
                fsum = 0.0;
                for( r = 0; r < effnoress; r++ ) {
                    f = sample->GetObsFrequency( n, r );
                    fsum += f;
                }
                if( fsum < 1.0 - highacc || 1.0 + highacc < fsum ) {
                    sprintf( bufstr, "Error correction failed (%g); pos. %d, file %s; skipped.",
                            fsum, position + n, my_basename( filename ));
                    warning( bufstr );
                    delete sample;
                    sample = NULL;
                    return;
                }
            }
            //TRANSFORM to NORMAL...
            if( tonormal )
                if(( serr = SampleToNormal( sample, n )) != NULL ) {
                    sprintf( bufstr, "%s; pos. %d, file %s; skipped.",
                            serr, position + n, my_basename( filename ));
                    warning( bufstr );
                    delete sample;
                    sample = NULL;
                    return;
                }
        }

        *p_sample = sample;
        sample = NULL;

    } catch( myexception const& ex ) {
        errstr = ex.what();
    }

    if( sample ) {
        delete sample;
        sample = NULL;
    }
    if( !errstr.empty())
        throw myruntime_error( errstr, CRITICAL );
}

// -------------------------------------------------------------------------
// SampleErrCorrection: correct sample error
//
bool ProFreqWriter::SampleErrCorrection( CtxtFrequencies* sample, int n )
{
    const size_t    effnoress = NUMAA;
    double          f, fsum, corr;
    int             r, nopf;
    const double    accuracy = 1.0e-1;
    const double    highacc = 1.0e-6;

    if( sample == NULL )
        return false;

    fsum = 0.0;
    nopf = 0;
    for( r = 0; r < effnoress; r++ ) {
        f = sample->GetObsFrequency( n, r );
        if( f != 0.0 )
            nopf++;
        fsum += f;
    }

    corr = fsum - 1.0;
    if( corr == 0.0 )
        return true;

    if( fabs( corr ) <= highacc ) {
        for( r = 0; r < effnoress; r++ ) {
            f = sample->GetObsFrequency( n, r );
            //correct small overflow error with the first frequency met
            if(( 0.0 < corr && corr <= f ) ||
                ( corr < 0.0 && (( nopf < effnoress && f == 0.0 ) || 
                                ( effnoress <= nopf  && f - corr <= 1.0 )))) {
                f -= corr;
                sample->SetObsFrequency( n, r, f );
                break;
            }
        }
    }
    else if( highacc < fabs( corr ) && nopf ) {
        corr /= ( double )nopf;
        for( r = 0; r < effnoress; r++ ) {
            f = sample->GetObsFrequency( n, r );
            if( f == 0.0 )
                continue;
            f -= corr;
            nopf--;
            //correct overflow/undeflow error
            if( f < 0.0 ) {
                if( nopf )
                    corr -= f /( double )nopf;
                f = 0.0;
            }
            else if( 1.0 < f ) {
                if( nopf )
                    corr += ( 1.0 - f )/( double )nopf;
                f = 1.0;
            }
            sample->SetObsFrequency( n, r, f );
        }
    }
    return true;
}

// -------------------------------------------------------------------------
// SampleToNormal: transform logistic normal sample to normal one
//
const char* ProFreqWriter::SampleToNormal( CtxtFrequencies* sample, int n )
{
    const bool      bCH0 = false;//check 0s
    const bool      bCHACC = false;//check accuracy
    const size_t    effnoress = NUMAA;
    double          f, flast, fsum, corr;
    int             r, nozs, nopf;
    double          fakev = 1.0e-4;//fake value
    const double    accuracy = 1.0e-6;

    if( sample == NULL || effnoress < 1 )
        return NULL;

    fsum = 0.0;
    nozs = 0;
    for( r = 0; r < effnoress; r++ ) {
        f = sample->GetObsFrequency( n, r );
        if( f == 0.0 )
            nozs++;
        fsum += f;
    }
    nopf = effnoress - nozs;

    if( bCHACC )
        if( fsum < 1.0 - accuracy || 1.0 + accuracy < fsum )
            return "Frequencies not conserved";

    if(( bCH0? nozs < 0: 0 < nozs )|| nopf <= 0 )
        return "Invalid frequencies";

    if( bCH0 ) {
        //make sure there are no null frequencies
        fsum = 0.0;
        corr = fakev * ( double )nozs / ( double )nopf;
        for( r = 0; r < effnoress; r++ ) {
            f = sample->GetObsFrequency( n, r );
            if( f == 0.0 ) {
                fsum += fakev;
                //frequency is zero, add fake value
                sample->SetObsFrequency( n, r, fakev );
                continue;
            }
            f -= corr;
            nopf--;
            if( f <= 0.0 ) {
                fsum += fakev;
                //frequency is less than correction value; recalc. `corr'
                sample->SetObsFrequency( n, r, fakev );
                corr += ( fakev - f ) / ( double )nopf;
                continue;
            }
            fsum += f;
            sample->SetObsFrequency( n, r, f );
        }
    }
    flast = sample->GetObsFrequency( n, effnoress - 1 );
    if( flast <= 0.0 )
        return "Null last element";

    if( bCHACC )
        if( fsum < 1.0 - accuracy || 1.0 + accuracy < fsum )
            return "Frequencies not conserved after correction 2";

    //now apply transformation to get normally distributed variable
    for( r = 0; r < effnoress - 1; r++ ) {
        f = sample->GetObsFrequency( n, r );
        if( f <= 0.0 )
            return "Frequencies not conserved after correction 3";
        f = log( f / flast );
        sample->SetObsFrequency( n, r, f );
    }
    sample->SetObsFrequency( n, effnoress - 1, 0.0 );

    return NULL;
}
