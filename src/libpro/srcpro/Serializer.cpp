/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <stdio.h>
#include <unistd.h>

#include "mystring.h"
#include "myexcept.h"

#include "version.h"
#include "Serializer.h"
#include "DistributionMatrix.h"
#include "GapScheme.h"
#include "FrequencyStore.h"



//procomprofile1.0| -0x2a
const char*  Serializer::signature[] = {
    "\nPROFILE FORMAT VERSION ",
    dataversion,
    "\x20\x46\x48\x45\x39\x45\x43\x46\x48\x45\x3c\x3f\x42\x3b\x07\x04\x06\x00",
    NULL
};


// -------------------------------------------------------------------------
// Constructor

Serializer::Serializer()
:
#ifdef SER_ANSI
    fp( NULL )
#else
    fd( -1 )
#endif
{
}

// -------------------------------------------------------------------------
// Destructor

Serializer::~Serializer()
{
}

// -------------------------------------------------------------------------
// SerializeProfile: serializes matrices
// -------------------------------------------------------------------------

void Serializer::SerializeProfile(  const FrequencyMatrix& freq, const LogOddsMatrix& pssm,
                                    const GapScheme& gaps )
{
    PutSignature();
    freq.Serialize( *this );
    pssm.Serialize( *this );
    gaps.Serialize( *this );
}

// -------------------------------------------------------------------------
// SerializeProfile: serializes matrices to file descriptor
// -------------------------------------------------------------------------

void Serializer::SerializeProfile(  const FrequencyMatrix& freq, const LogOddsMatrix& pssm,
                                    const GapScheme& gaps, FILE* ffp )
{
    if( ffp == NULL )
        return;

    if( fp != NULL )
        throw myruntime_error( mystring( "Serializer: Unable to serialize profile(s)." ));

    fp = ffp;

    try {
        SerializeProfile( freq, pssm, gaps );

    } catch( myexception const& ex )
    {
        fp = NULL;
        throw myruntime_error( ex.what(), ex.eclass());
    }
    fp = NULL;
}

// -------------------------------------------------------------------------
// SerializeProfile: serializes matrices to file
// -------------------------------------------------------------------------

void Serializer::SerializeProfile(  const FrequencyMatrix& freq, const LogOddsMatrix& pssm,
                                    const GapScheme& gaps, const char* filename )
{
#ifdef SER_ANSI
    if( fp != NULL )
        throw myruntime_error( mystring( "Serializer: Unable to serialize profile(s)." ));

    fp = fopen( filename, "wb" );
    if( fp == NULL )
        throw myruntime_error(
            mystring( "Serializer: Failed to open file for writing." ));
#else
    if( fd != -1 )
        throw myruntime_error( mystring( "Serializer: Unable to serialize profile(s)." ));

    fd = open( filename, O_WRONLY | O_CREAT | O_TRUNC );
    if( fd == -1 )
        throw myruntime_error(
            mystring( "Serializer: Failed to open file for writing." ));
#endif

    try {
        SerializeProfile( freq, pssm, gaps );

    } catch( myexception const& ex )
    {
#ifdef SER_ANSI
        fclose( fp );
        fp = NULL;
#else
        close( fd );
        fd = -1;
#endif        
        throw myruntime_error( ex.what(), ex.eclass());
    }


#ifdef SER_ANSI
    fclose( fp );
    fp = NULL;
#else
    close( fd );
    fd = -1;
#endif
}

// -------------------------------------------------------------------------
// DeserializeProfile: deserializes profile matrices
// -------------------------------------------------------------------------

void Serializer::DeserializeProfile( FrequencyMatrix& freq, LogOddsMatrix& pssm,
                                     GapScheme& gaps )
{
    GetSignature();
    freq.Deserialize( *this );
    pssm.Deserialize( *this );
    gaps.Deserialize( *this );
}

// -------------------------------------------------------------------------
// DeserializeProfile: deserializes matrices from file descriptor
// -------------------------------------------------------------------------

void Serializer::DeserializeProfile( FrequencyMatrix& freq, LogOddsMatrix& pssm,
                                     GapScheme& gaps, FILE* ffp )
{
    if( ffp == NULL )
        return;

    if( fp != NULL )
        throw myruntime_error( mystring( "Serializer: Unable to deserialize profile(s)." ));

    fp = ffp;

    try {
        DeserializeProfile( freq, pssm, gaps );

    } catch( myexception const& ex )
    {
        fp = NULL;
        throw myruntime_error( ex.what(), ex.eclass());
    }
    fp = NULL;
}

// -------------------------------------------------------------------------
// DeserializeProfile: deserializes profile matrices from file
// -------------------------------------------------------------------------

void Serializer::DeserializeProfile( FrequencyMatrix& freq, LogOddsMatrix& pssm,
                                     GapScheme& gaps, const char* filename )
{
#ifdef SER_ANSI
    if( fp != NULL )
        throw myruntime_error( mystring( "Serializer: Unable to deserialize profile(s)." ));

    fp = fopen( filename, "rb" );
    if( fp == NULL )
        throw myruntime_error(
            mystring( "Serializer: Failed to open file for reading." ));
#else
    if( fd != -1 )
        throw myruntime_error( mystring( "Serializer: Unable to deserialize profile(s)." ));

    fd = open( filename, O_RDONLY );
    if( fd == -1 )
        throw myruntime_error(
            mystring( "Serializer: Failed to open file for reading." ));
#endif

    try {
        DeserializeProfile( freq, pssm, gaps );

    } catch( myexception const& ex )
    {
#ifdef SER_ANSI
        fclose( fp );
        fp = NULL;
#else
        close( fd );
        fd = -1;
#endif
        throw myruntime_error( ex.what(), ex.eclass());
    }


#ifdef SER_ANSI
    fclose( fp );
    fp = NULL;
#else
    close( fd );
    fd = -1;
#endif
}





// =========================================================================
// WriteProfile: write profile to file descriptor
//
void Serializer::WriteProfile(
        FILE* ffp,
        const FrequencyMatrix& frqs,
        const LogOddsMatrix& pssm,
        const GapScheme& gaps,
        int scale )
{
    if( ffp == NULL )
        return;

    if( fp != NULL )
        throw myruntime_error( mystring( "Serializer: Unable to write profile(s)." ));

    fp = ffp;

    try {
        TextWriteProfile( fp, frqs, pssm, gaps, scale );

    } catch( myexception const& ex )
    {
        fp = NULL;
        throw myruntime_error( ex.what(), ex.eclass());
    }
    fp = NULL;
}

// WriteProfile: write profile to file
//
void Serializer::WriteProfile(
        const char* fname,
        const FrequencyMatrix& frqs,
        const LogOddsMatrix& pssm,
        const GapScheme& gaps,
        int scale )
{
#ifdef SER_ANSI
    mystring    errstr;
    int         eclass;

    if( fp != NULL )
        throw myruntime_error( "Serializer: Unable to write profile(s)." );

    fp = fopen( fname, "w" );
    if( fp == NULL )
        throw myruntime_error( "Serializer: Failed to open file for writing." );

    try {
        TextWriteProfile( fp, frqs, pssm, gaps, scale );

    } catch( myexception const& ex )
    {
        errstr = ex.what();
        eclass = ex.eclass();
    }
    fclose( fp );
    fp = NULL;
    if( !errstr.empty())
        throw myruntime_error( errstr, eclass );
#else
    throw myruntime_error( mystring( "Serializer: Unable to write profile(s)." ));
#endif
}

// -------------------------------------------------------------------------
// ReadProfile: read profile from file descriptor
//
void Serializer::ReadProfile( FILE* ffp, FrequencyMatrix& frqs, LogOddsMatrix& pssm, GapScheme& gaps )
{
    if( ffp == NULL )
        return;

    if( fp != NULL )
        throw myruntime_error( mystring( "Serializer: Unable to read profile(s)." ));

    fp = ffp;

    try {
        TextReadProfile( fp, frqs, pssm, gaps );

    } catch( myexception const& ex )
    {
        fp = NULL;
        warning( ex.what());
        throw myruntime_error( ex.what(), ex.eclass());
    }
    fp = NULL;
}

// ReadProfile: read profile from file
//
void Serializer::ReadProfile( const char* fname, FrequencyMatrix& frqs, LogOddsMatrix& pssm, GapScheme& gaps )
{
#ifdef SER_ANSI
    mystring    errstr;
    int         eclass;

    if( fp != NULL )
        throw myruntime_error("Serializer: Unable to read profile(s).");

    fp = fopen( fname, "r" );
    if( fp == NULL )
        throw myruntime_error("Serializer: Failed to open file for reading.");

    try {
        TextReadProfile( fp, frqs, pssm, gaps );

    } catch( myexception const& ex )
    {
        errstr = ex.what();
        eclass = ex.eclass();
        warning( ex.what());
    }
    fclose( fp );
    fp = NULL;
    if( !errstr.empty())
        throw myruntime_error( errstr, eclass );
#else
    throw myruntime_error("Serializer: Unable to read profile(s).");
#endif

}





// =========================================================================
// SerializeFrequencies: serializes frequency vector
// -------------------------------------------------------------------------

void Serializer::SerializeFrequencies( const FrequencyVector& vect )
{
    vect.Serialize( *this );
}

// -------------------------------------------------------------------------
// SerializeFrequencies: serializes frequency vector to file descriptor
// -------------------------------------------------------------------------

void Serializer::SerializeFrequencies( const FrequencyVector& vect, FILE* ffp )
{
    if( ffp == NULL )
        return;

    if( fp != NULL )
        throw myruntime_error( mystring( "Serializer: Unable to serialize frequencies." ));

    fp = ffp;

    try {
        SerializeFrequencies( vect );

    } catch( myexception const& ex )
    {
        fp = NULL;
        throw myruntime_error( ex.what(), ex.eclass());
    }
    fp = NULL;
}

// -------------------------------------------------------------------------
// DeserializeFrequencies: deserializes frequency vector
// -------------------------------------------------------------------------

void Serializer::DeserializeFrequencies( FrequencyVector& vect )
{
    vect.Deserialize( *this );
}

// -------------------------------------------------------------------------
// DeserializeFrequencies: deserializes frequency vector from file
//     descriptor
// -------------------------------------------------------------------------

void Serializer::DeserializeFrequencies( FrequencyVector& vect, FILE* ffp )
{
    if( ffp == NULL )
        return;

    if( fp != NULL )
        throw myruntime_error( mystring( "Serializer: Unable to deserialize frequencies." ));

    fp = ffp;

    try {
        DeserializeFrequencies( vect );

    } catch( myexception const& ex )
    {
        fp = NULL;
        throw myruntime_error( ex.what(), ex.eclass());
    }
    fp = NULL;
}





// =========================================================================
// WriteVector: write vector of frequencies and scores to file descriptor
//
void Serializer::WriteVector( FILE* ffp, const FrequencyVector& vect )
{
    if( ffp == NULL )
        return;

    if( fp != NULL )
        throw myruntime_error( mystring( "WriteVector: Unable to write data vector." ));

    fp = ffp;

    try {
        vect.TextWriteVector( fp );

    } catch( myexception const& ex )
    {
        fp = NULL;
        throw myruntime_error( ex.what(), ex.eclass());
    }
    fp = NULL;
}

// -------------------------------------------------------------------------
// ReadVector: read vector of frequencies and scores from file descriptor
//
void Serializer::ReadVector( FILE* ffp, FrequencyVector& vect )
{
    if( ffp == NULL )
        return;

    if( fp != NULL )
        throw myruntime_error( mystring( "ReadVector: Unable to read data vector." ));

    fp = ffp;

    try {
        vect.TextReadVector( fp );

    } catch( myexception const& ex )
    {
        fp = NULL;
        throw myruntime_error( ex.what(), ex.eclass());
    }
    fp = NULL;
}





// =========================================================================
// Write: static write function
// -------------------------------------------------------------------------

void Serializer::Write( FILE* ffp, char* buffer, size_t size, size_t length )
{
    if( ffp == NULL )
        return;

    size_t  ret;
    ret = fwrite( buffer, size, length, ffp );
    if( ret != length )
        throw myruntime_error( mystring( "Write of data to file failed." ));
}

// -------------------------------------------------------------------------
// Read: static read function
// -------------------------------------------------------------------------

void Serializer::Read( FILE* ffp, char* buffer, size_t size, size_t length )
{
    if( ffp == NULL )
        return;

    size_t  ret;
    ret = fread( buffer, size, length, ffp );
    if( ret != length )
        throw myruntime_error( mystring( "Read of data from file failed." ));
}

// -------------------------------------------------------------------------
// Write: writes information from the buffer to file
// -------------------------------------------------------------------------

void Serializer::Write( char* buffer, size_t size, size_t length )
{
#ifdef SER_ANSI
    if( fp == NULL )
        throw myruntime_error(
            mystring( "Serializer: invalid file descriptor." ));

    size_t  ret;
    ret = fwrite( buffer, size, length, fp );
    if( ret != length )
        throw myruntime_error(
            mystring( "Serializer: failed to write data to file." ));

#else
    if( fd == -1 )
        throw myruntime_error(
            mystring( "Serializer: invalid file descriptor." ));

    ssize_t ret;
    ssize_t written = 0;
    ssize_t towrite = 0;

    written = 0;
    towrite = length * size;

    while( written < towrite ) {
        ret = write( fd, buffer + written, towrite - written );
        if( ret <= 0 )
            throw myruntime_error(
                mystring( "Serializer: failed to write data to file." ));
        written += ret;
    }
#endif
}

// -------------------------------------------------------------------------
// Read: reads information from file to the buffer
// -------------------------------------------------------------------------

void Serializer::Read( char* buffer, size_t size, size_t length )
{
#ifdef SER_ANSI
    if( fp == NULL )
        throw myruntime_error(
            mystring( "Serializer: invalid file descriptor." ));

    size_t  ret;
    ret = fread( buffer, size, length, fp );
    if( ret != length )
        throw myruntime_error(
            mystring( "Serializer: failed to read data from file: wrong format." ));

#else
    if( fd == -1 )
        throw myruntime_error(
            mystring( "Serializer: invalid file descriptor." ));

    ssize_t ret;
    ssize_t haveread = 0;
    ssize_t toread = 0;

    haveread = 0;
    toread = length * size;

    while( haveread < toread ) {
        ret = read( fd, buffer + haveread, toread - haveread );
        if( ret <= 0 )
            throw myruntime_error(
                mystring( "Serializer: failed to read data from file: wrong format." ));
        haveread += ret;
    }
#endif
}

// -------------------------------------------------------------------------
// PutSignature: puts signature in the file opened
// -------------------------------------------------------------------------

void Serializer::PutSignature()
{
    for( int n = 0; signature[n]; n++ )
        Write(( char* )signature[n], 1, strlen( signature[n] ));
}

// -------------------------------------------------------------------------
// GetSignature: gets signature from the file and checks it for validness
// -------------------------------------------------------------------------

void Serializer::GetSignature()
{
    char    loc_signature[CHAR_MAX];
    //
    for( int n = 0; signature[n]; n++ ) {
        Read( loc_signature, 1, strlen( signature[n] ));

        if( strncmp( signature[n], loc_signature, strlen( signature[n] ))) {
            if( n == DATAVER )
                throw myruntime_error(
                    mystring( "Data format version number is invalid." ));
                //warning( "Data formats differ in version number." );
            else
                throw myruntime_error(
                    mystring( "Serializer: data format of the file to be read is invalid." ));
        }
    }
}

