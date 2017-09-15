/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <stdio.h>
#include <stdarg.h>
#include "myexcept.h"
#include "rc.h"


const char* PROGDIR = NULL;
const char* PROGNAME = NULL;
const char* PROGVERSION = NULL;

int*        __PARGC__ = NULL;
char***     __PARGV__ = NULL;

bool        BEQUIET = false;


// -------------------------------------------------------------------------
// Routines to set the global variables
//
void SetGlobalProgName( const char* name, const char* version )
{
    PROGDIR = my_dirname( name );
    PROGNAME = my_basename( name );
    if( version )
        PROGVERSION = version;
}


void SetArguments( int* pargc, char*** pargv )
{
    __PARGC__ = pargc;
    __PARGV__ = pargv;
}


void SetQuiet( bool quiet )
{
    BEQUIET = quiet;
}

// -------------------------------------------------------------------------
// Printing routines standard to the toolset
//
void error( const char* str, bool putnl )
{
    if( PROGNAME )
        fprintf( stderr, "[%s] ", PROGNAME );
    fprintf( stderr, "Error: %s\n", str );
    if( putnl )
        fprintf( stderr, "\n" );
}


void warning( const char* str, bool putnl )
{
    if( BEQUIET )
        return;

    if( PROGNAME )
        fprintf( stderr, "[%s] ", PROGNAME );
    if( str )
        fprintf( stderr, "Warning: %s\n", str );
    if( putnl )
        fprintf( stderr, "\n" );
}

void message( const char* str, bool putnl )
{
    if( BEQUIET )
        return;

    if( str ) {
        if( PROGNAME )
            fprintf( stderr, "[%s] ", PROGNAME );

        fprintf( stderr, "%s\n", str );
    }

    if( putnl )
        fprintf( stderr, "\n" );
}

void messagec( char sym, bool indent )
{
    if( BEQUIET )
        return;

    if( indent )
        fprintf( stderr, "%20c", sym );
    else
        fprintf( stderr, "%c", sym );
}

void progname_and_version( FILE* fp )
{
    if( PROGNAME )
        fprintf( fp, "%s", PROGNAME );
    if( PROGVERSION )
        fprintf( fp, " %s", PROGVERSION );
    if( PROGNAME || PROGVERSION )
        fprintf( fp, "\n\n" );
}

// -------------------------------------------------------------------------
// print_cmdline: print command line
//
void print_cmdline( TPrintFunction pfunc, void* vpn )
{
    int n;

    if( !vpn )
        return;

    if( __PARGV__ && __PARGC__  && *__PARGV__ ) {
        for( n = 0; n < *__PARGC__; n++ )
            pfunc( vpn, "%s ", (*__PARGV__)[n] );
    }
    pfunc( vpn, "\n" );
}

// -------------------------------------------------------------------------
// Path-processing functions
//
const char* my_basename( const char* name )
{
    const char* bp = NULL;

    if( name )
        for(bp  = name + strlen( name );
            bp != name && bp[-1] != DIRSEP;
            bp-- );
    return  bp;
}

// my_dirname: returns directory name given pathname; if directory is found
//     in the path given, current directory (i.e. '.') is returned
//
const char* my_dirname( const char* name )
{
    static const size_t _size = KBYTE;
    static char dir[_size];
    const char* bp = name;
    const char* rp = name;

    dir[0] = '.'; dir[1] = 0;

    if( !bp )
        return dir;

    for( bp += strlen( name ); bp != name && bp[-1] != DIRSEP; bp-- );

    if( bp != name && *bp == 0 ) {
        //the last character is directory separator
        for( ; bp != name && bp[-1] == DIRSEP; bp-- );

        if( bp != name )
            //find separator before the last one
            for( ; bp != name && bp[-1] != DIRSEP; bp-- );
    }

    if( bp == name )
        //no directory separator has been found
        return dir;

    for( ; bp != name && bp[-1] == DIRSEP; bp-- );

    if( bp == name ) {
        bp++;

        if( *bp == DIRSEP && bp[1] != DIRSEP )
            //there are exactly two separators at the beginning; save them
            bp++;
    }

    if( bp - name < _size ) {
        //if having enough memory allocated
        memcpy( dir, name, bp - name );
        dir[ bp - name ] = 0;
    }

    return dir;
}

// -------------------------------------------------------------------------
// usage: returns <instructions> translated by appropriately inserting
//     program name, version, and date information
//
mystring usage( const char* path, const char* instructions, const char* version, const char* date )
{
    size_t      pos;
    mystring    instr = instructions;
    mystring    prog = path;
    mystring    full;
    bool        first = true;

    if( !path || !instructions )
        return instr;

    if(( pos = prog.rfind( DIRSEP )) != mystring::npos ){
        prog = prog.substr( pos+1 );
    }

    full = prog + mystring( " " );

    if( version && strlen( version ))
        full += version;

    if( date && strlen( date ))
        full += mystring( " (" ) + date + mystring( ")" );

    while(( pos = instr.find( "<>" )) != mystring::npos ){
        instr.erase( pos, 2 );
        instr.insert( pos, first? full: prog );
        if( first ) first = false;
    }
    while(( pos = instr.find( "-*" )) != mystring::npos ){
        instr.erase( pos, 2 );
        for( size_t len = full.length(); len; len-- )
            instr.insert( pos++, "-"  );
    }
    return instr;
}

// -------------------------------------------------------------------------
// file_print: function redirects formatted output string to file; file
//     pointer (vpn) must be intialized earlier and must be valid
//
int file_print( void* vpn, const char* format, ... )
{
    if( !vpn )
        return -1;

    FILE*   fp = ( FILE* )vpn;
    va_list ap;
    int     ret;

    va_start( ap, format );

    ret = vfprintf( fp, format, ap );

    va_end( ap );

    if( ret < 0 )
        throw myruntime_error( mystring( "Formatted print to file failed." ));

    return ret;
}

// -------------------------------------------------------------------------
// string_print: same as file_print except that redirection is to string
//     which is to have enough space to contain information; string
//     pointer (vpn) must be allocated earlier and must be valid
//
int string_print( void* vpn, const char* format, ... )
{
    if( !vpn )
        return -1;

    char*   sp = ( char* )vpn;
    va_list ap;
    int     ret;

    va_start( ap, format );

    ret = vsprintf( sp + strlen( sp ), format, ap );

    va_end( ap );

    if( ret < 0 )
        throw myruntime_error( mystring( "Formatted print to string failed." ));

    return ret;
}

