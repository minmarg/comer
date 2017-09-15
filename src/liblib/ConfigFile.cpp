/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "rc.h"
#include "ConfigFile.h"




const char*     ConfigFile::s_equal   = "=";
const char*     ConfigFile::s_bracket = "[]";
const char*     ConfigFile::s_newline = "\n";

//------------------------------------------------------------------------------
// Constructor
//

ConfigFile::ConfigFile( const char* section, const char* filename )
:   section_name( section ),
    filename( filename ),
    errbuffer( NULL ),
    errbuflen( 0 )
{
    errbuflen = KBYTE;
    errbuffer = ( char* )malloc( errbuflen );
    if( !errbuffer )
        throw myruntime_error( mystring( "ConfigFile: Not enough memory." ));
}

//------------------------------------------------------------------------------
// Destructor
//

ConfigFile::~ConfigFile()
{
    if( errbuffer )
        free( errbuffer );
}

//------------------------------------------------------------------------------
// GetInt: gets integer value of the key
//     key,  key name
//     value, address to write key value in
// returns true if the key has been found and its value has been extracted,
//     false otherwise
//

bool ConfigFile::GetInt( const char* key, int* value )
{
    return GetPrivateProfileInt( key, value );
}

//------------------------------------------------------------------------------
// GetDouble: gets double value of the key
//     key,  key name
//     value, address to write key value in
// returns true if the key has been found and its value has been extracted,
//     false otherwise
//

bool ConfigFile::GetDouble( const char* key, double* value )
{
    return GetPrivateProfileDouble( key, value );
}

//------------------------------------------------------------------------------
// GetString: gets string value of the key
//     key, key name
//     value, address of key value
//     size, size of value
// returns true if the key has been found and its value has been extracted,
//     false otherwise
//

bool ConfigFile::GetString(  const char* key, char* value, unsigned size )
{
    return GetPrivateProfileString( key, value, size );
}

//------------------------------------------------------------------------------
// WriteInt: writes integer value of the key to file
//     key, key name
//     value, key value
// returns true if the key has been found; if the key or section haven't been
//     found, they are created and inserted to the file
//

bool ConfigFile::WriteInt( const char* key, int value )
{
    char        strvalue[BUF_MAX];

    sprintf( strvalue, "%d", value );
    return WriteString( key, strvalue );
}

//------------------------------------------------------------------------------
// WriteDouble: writes double value of the key to file
//     key, key name
//     value, key value
// returns true if the key has been found; if the key or section haven't been
//     found, they are created and inserted to the file
//

bool ConfigFile::WriteDouble( const char* key, double value )
{
    char        strvalue[BUF_MAX];

    sprintf( strvalue, "%.6f", value );
    return WriteString( key, strvalue );
}

//------------------------------------------------------------------------------
// WriteString: writes string value of the key to file
//     key, key name
//     value, key value
// returns true if the key has been found; otherwise it'll be inserted
//

bool ConfigFile::WriteString( const char* key, const char* value )
{
    return WritePrivateProfileString( key, value );
}
    
//------------------------------------------------------------------------------
// GetPrivateProfileInt: opens the file to search the key for a value
// returns true if the key has been found and its value has been extracted,
//     false otherwise
//

bool ConfigFile::GetPrivateProfileInt( const char* key, int* value )
{
    char    strvalue[BUF_MAX];
    char*   p;

    if( !value )
        throw myruntime_error( mystring( "ConfigFile: Wrong argument." ));

    if( !GetPrivateProfileString( key, strvalue, BUF_MAX ))
        return false;

    *value = strtol( strvalue, &p, 10 );
    if( errno || *p ) {
        Format_Message( "ConfigFile: Invalid integer value of key  ", key );
        throw myruntime_error( mystring( GetErrBuffer()));
    }

    return true;
}

//------------------------------------------------------------------------------
// GetPrivateProfileDouble: opens the file to search the key for a double value
// returns true if the key has been found and its value has been extracted,
//     false otherwise
//

bool ConfigFile::GetPrivateProfileDouble( const char* key, double* value )
{
    char    strvalue[BUF_MAX];
    char*   p;

    if( !value )
        throw myruntime_error( mystring( "ConfigFile: Wrong argument." ));

    if( !GetPrivateProfileString( key, strvalue, BUF_MAX ))
        return false;

    *value = strtod( strvalue, &p );

    if( errno || *p ) {
        Format_Message( "ConfigFile: Invalid double value of key  ", key );
        throw myruntime_error( mystring( GetErrBuffer()));
    }

    return true;
}

//------------------------------------------------------------------------------
// GetPrivateProfileString: opens the file to search the key for a value of
//     type string
//     key, key name
//     values, address where the results to be written to
//     size, how many bytes <values> can contain
// returns true if the key has been found and its value has been extracted,
//     false otherwise
//

bool ConfigFile::GetPrivateProfileString( const char* key, char* value, unsigned size )
{
    if( !GetSectionName() || !GetFileName())
        throw myruntime_error( mystring( "ConfigFile: No filename or section name given." ));

    FILE*       fd = NULL;
    bool        found = false;

    fd = fopen( GetFileName(), "r" );

    if( fd == NULL ) {
        Format_Message( "ConfigFile: Failed to open file, ", GetFileName());
        throw myruntime_error( mystring( GetErrBuffer()));
    }

    try {
        if( FindFileSection( fd ))
            if( GetSectionKey( fd, key, value, size ))
                found =  true;

    } catch( myexception const& ex )
    {
        fclose( fd );
        throw myruntime_error( ex.what(), ex.eclass());
    }

    fclose( fd );
    return found;
}

//------------------------------------------------------------------------------
// WritePrivateProfileString: opens the file to write the key value
//          sekcijoje duota rakto simboline reiksme
//     key, key name
//     value, symbolic key value
// returns true if the key has been found; if the key or section haven't been
//     found, they are created and inserted to the file
//

bool ConfigFile::WritePrivateProfileString( const char* key, const char* value )
{
    if( !GetSectionName() || !GetFileName())
        throw myruntime_error( mystring( "ConfigFile: No filename or section name given." ));

    FILE*       fd = NULL;
    bool        found = false;

    fd = fopen( GetFileName(), "r+w" );

    if( fd == NULL ) {
        Format_Message( "ConfigFile: Failed to open file, ", GetFileName());
        throw myruntime_error( mystring( GetErrBuffer()));
    }

    try {
        found = FindFileSection( fd, true );
        if( WriteSectionKey( fd, key, value, found ))
            found =  true;

    } catch( myexception const& ex )
    {
        fclose( fd );
        throw myruntime_error( ex.what(), ex.eclass());
    }

    fclose( fd );
    return found;
}

//------------------------------------------------------------------------------
// FindFileSection: finds section in the file by name; if not found, inserts one
//     if needed
//     fd, pointer to file
//     bwrite, true if append section to the file
// returns true if section has been found, false otherwise
//
bool ConfigFile::FindFileSection( FILE* fd, bool bwrite )
{
    if( !GetSectionName())
        throw myruntime_error( mystring( "ConfigFile: No section name given." ));

    if( fd == NULL )
        return false;

    const int   locsize = BUF_MAX;
    char        locbuffer[locsize] = {0};
    char*       p = locbuffer;
    const int   seclen = strlen( GetSectionName());
    const int   nllen = strlen( s_newline );
    size_t      len = nllen;
    size_t      nlns = 0;           //number of newlines until the end of file
    bool        found = false;

    while( !feof( fd )) {

        do {
            p = fgets( locbuffer, locsize, fd );

            len = strlen( locbuffer );

            if( p == NULL )
                break;

            if( len && strcmp( locbuffer + len - nllen, s_newline ) == 0 )
                break;

        } while( 1 );

        if( p == NULL )
            break;

        //omit spaces before checking for newline
        for( p = locbuffer; *p == 32 || *p == 9; p++ );

        if( nllen <= len && strcmp( p, s_newline ) == 0 )
            nlns += nllen + ( p - locbuffer ) - ( p != locbuffer );
        else
            nlns = 0;

        if( *( p = locbuffer ) != s_bracket[0] )
            continue;
        p++;

        if( strncmp( p, GetSectionName(), seclen ) == 0 && p[seclen] == s_bracket[1] ) {
            found = true;
            break;
        }
    }

    if( !bwrite || found )
        return found;

    SeekToPosition( fd, -nlns, SEEK_END );

    if( fprintf( fd, "%s%c%s%c%s", s_newline, s_bracket[0], GetSectionName(), s_bracket[1], s_newline ) < 0 ) {
        Format_Message( "ConfigFile: Unable to write to file ", GetFileName());
        throw myruntime_error( mystring( GetErrBuffer()));
    }

    return found;
}

//------------------------------------------------------------------------------
// GetSectionKey: searches for a key in the file by the name
//     fd, pointer to file
//     key, key name
//     value, key value
//     size, size of value
//     pos, address to contain file position of the key found or the end of the
//          section
// returns true if the key has been found, false otherwise
//

bool ConfigFile::GetSectionKey( FILE* fd, const char* key, char* value, int size, long* pos, long* lastpos )
{
    if( fd == NULL || key == NULL )
        return false;

    const int   locsize = BUF_MAX;
    char        locbuffer[locsize] = {0};
    char*       p = locbuffer;
    char*       lastline = locbuffer;
    int         keylen = 0;
    int         equlen = strlen( s_equal );
    int         nllen = strlen( s_newline );
    size_t      len = 1;
    size_t      plen = 0;
    int         cmpeq;
    bool        found = false;
    bool        nl_assigned = false;


    keylen = strlen( key );

    while( !feof( fd )) {

        if( pos && !nl_assigned )
            *pos = TellPosition( fd );

        lastline = NULL;

        do {
            p = fgets( locbuffer, locsize, fd );

            len = strlen( locbuffer );

            if( p == NULL )
                //if at the end of file
                break;

            if( !lastline && len && *locbuffer == s_bracket[0] )
                //if arrived at new section
                return found;

            if(!lastline ) {
                //omit spaces
                for( lastline = locbuffer; *lastline == 32 || *lastline == 9; lastline++ );
            }

        } while( len && strcmp( locbuffer + len - nllen, s_newline ));

        if( lastline && strcmp( lastline, s_newline ) == 0 )
                    nl_assigned = true;
            else    nl_assigned = false;

        if( lastpos && !nl_assigned )
            *lastpos = TellPosition( fd );

        if( p == NULL )
            break;

        //allow spaces before the key
        for( p = locbuffer; *p == 32 || *p == 9; p++ );

        if( strncmp( key, p, keylen ))
            continue;

        //omit spaces before sign of equality
        for( p += keylen; *p == 32 || *p == 9; p++ );
        if( strncmp( p, s_equal, equlen ))
            continue;

        //omit spaces after equality
        for( p += equlen; *p == 32 || *p == 9; p++ );

        plen = strlen( p );
        if( size < plen )
            plen = size;

        if( value && plen ) {
            strncpy( value, p, plen - 1 );
            value[plen-1] = 0;
            ProcessValue( value );

            found = true;
            break;
        }
    }

    return found;
}

//------------------------------------------------------------------------------
// WriteSectionKey: tries to write a key in the file by searching the key itself
//     at first
//     fd, pointer to file
//     key, key name
//     value, key value
// returns true if the key has been found, false otherwise
//

bool ConfigFile::WriteSectionKey( FILE* fd, const char* key, const char* value, bool section_found )
{
    if( fd == NULL || key == NULL || value == NULL )
        return false;

    char    strvalue[BUF_MAX];
    long    file_pos = 0;
    long    last_pos = 0;
    bool    found = false;

    found = GetSectionKey( fd, key, strvalue, BUF_MAX, &file_pos, &last_pos );
    //we have position to insert key after
    InsertSectionKeyAfter( fd, key, value, file_pos, last_pos, section_found );

    return found;
}

//------------------------------------------------------------------------------
// InsertSectionKeyAfter: tries to insert the key with its value in the file
//     after the position specified
//     fd, pointer to file
//     key, key name
//     value, key value
//     file_pos, file position to insert key after
//

void ConfigFile::InsertSectionKeyAfter(
        FILE* fd,
        const char* key,    const char* value,
        long file_pos,      long last_pos,
        bool section_found )
{
    if( fd == NULL || key == NULL || value == NULL )
        return;

    SeekToPosition( fd, last_pos, SEEK_SET );

    const int   size = KBYTE;
    char*       inner = ( char* )malloc( size );
    char*       p = inner;
    int         cursize = size;
    int         prvsize = 0;
    int         amount = 0;
    size_t      lastread = 0;

    if( inner == NULL )
        throw myruntime_error( mystring( "ConfigFile: Not enough memory." ));

    *inner = 0;

    while( !feof( fd ) && section_found ) {
        lastread = fread( p, 1, size, fd );
        amount += lastread;

        if( ferror( fd )) {
            free( inner );
            Format_Message( "ConfigFile: Read error: ", GetFileName());
            throw myruntime_error( mystring( GetErrBuffer()));
        }

        if( feof( fd ))
            break;

        inner = ( char* )realloc( inner, cursize += size );

        if( inner == NULL )
            throw myruntime_error( mystring( "ConfigFile: Not enough memory." ));

        //each time inner changes its address!!
        p = inner + amount;
    }

    SeekToPosition( fd, file_pos, SEEK_SET );

    if( fprintf( fd, "%s %s %s%s", key, s_equal, value, s_newline ) < 0 ) {
        free( inner );
        Format_Message( "ConfigFile: Unable to insert key in file ", GetFileName());
        throw myruntime_error( mystring( GetErrBuffer()));
    }

    for( p = inner, prvsize = size, cursize = 0; cursize < amount; cursize += size ) {
        if( amount < cursize + size )
            prvsize = amount - cursize;

        fwrite( p, 1, prvsize, fd );

        if( ferror( fd )) {
            free( inner );
            Format_Message( "ConfigFile: Write error: ", GetFileName());
            throw myruntime_error( mystring( GetErrBuffer()));
        }

        p += prvsize;
    }

    free( inner );

    if( section_found )
        FlushToEOF( fd );
}

//------------------------------------------------------------------------------
// FlushToEOF: flushes file from the current position to the end of file
//

void ConfigFile::FlushToEOF( FILE* fd )
{
    if( !fd )
        return;

    char*   temp = NULL;
    long    curpos = TellPosition( fd );    SeekToPosition( fd, 0, SEEK_END );
    long    endpos = TellPosition( fd );
    size_t  nllen = strlen( s_newline );

    long    diff = endpos - curpos;

    if( 0 < diff && diff < nllen )
        diff = nllen;

    if( 0 < diff ) {
        SeekToPosition( fd, curpos, SEEK_SET );
        temp = ( char* )malloc( diff + 1 ); //one for null character
        if( !temp )
            throw myruntime_error( mystring( "ConfigFile: Not enough memory." ));
        memset( temp, 32, diff );
        strcpy( temp + diff - nllen, s_newline );
        fwrite( temp, 1, diff, fd );
        free( temp );
    }
}

//------------------------------------------------------------------------------
// TellPosition: tells current position of the file
//

long ConfigFile::TellPosition( FILE* fd )
{
    long    pos = 0;

    if( !fd )
        return pos;

    pos = ftell( fd );

    if( pos < 0 ) {
        Format_Message( "ConfigFile: Unable to determine position in file ", GetFileName());
        throw myruntime_error( mystring( GetErrBuffer()));
    }

    return pos;
}

//------------------------------------------------------------------------------
// SeekToPosition: seeks to the position in the file
//

void ConfigFile::SeekToPosition( FILE* fd, long pos, int whence )
{
    if( fseek( fd, pos, whence ) < 0 ) {
        Format_Message( "ConfigFile: Unable to seek in file ", GetFileName());
        throw myruntime_error( mystring( GetErrBuffer()));
    }
}

//------------------------------------------------------------------------------
// ProcessValue: processes value of key by removing leading and trailing spaces;
// comments at the end of the value are also flushed;
// NOTE: argument must be null-terminated string
//

void ConfigFile::ProcessValue( char* value )
{
    if( !value )
        return;

    size_t  len = strlen( value );
    bool    text = false;

    for( int n = len - 1; 0 <= n; n-- ) {
        if( value[n] == '#' ) {
            value[n] = 0;
            text = false;
            continue;
        }
        if( !text && ( value[n] == 32 || value[n] == 9 )) {
            value[n] = 0;
            continue;
        }
        text = true;
    }
}

//------------------------------------------------------------------------------
// Format_Message: formats message by appending the filename to the end
//

void ConfigFile::Format_Message( const char* msg, const char* cval )
{
    int             length = GetErrBufLength();
    char*           p = GetErrBuffer();

    if( length <= 0 )
        throw myruntime_error( mystring( "ConfigFile: No allocated space for message." ));

    *p = 0;

    if( !msg || !cval )
        return;

    int             prelen = strlen( msg );
    int             fnmlen = strlen( cval );

    if( length <= prelen )
        prelen  = length - 1;

    memcpy( p, msg, prelen );
    length -= prelen;
    p += prelen;

    if( length < 2 ) {
        *p = 0;
        return;
    }

    if( fnmlen + 1 < length )
        length = fnmlen + 1;

    memcpy( p, cval, length - 1 );
    p[length - 1] = 0;
}

