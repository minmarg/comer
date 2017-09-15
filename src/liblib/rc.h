/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __rc__
#define __rc__

#include <stdio.h>
#include "mystring.h"
#include "myexcept.h"
#include "platform.h"
#include "defines.h"
#include "faccess.h"
#include "pcmath.h"

// typedefs
typedef int ( *TPrintFunction )( void*, const char* format, ... );


// global variables used in the toolset
extern const char*  PROGDIR;
extern const char*  PROGNAME;
extern const char*  PROGVERSION;

extern int*         __PARGC__;
extern char***      __PARGV__;

extern bool         BEQUIET;

// assignment of global variables
void SetQuiet( bool );
void SetGlobalProgName( const char* name, const char* version = NULL );
void SetArguments( int* pargc, char*** pargv );

// printing
void error( const char*, bool = true );
void warning( const char*, bool = true );
void message( const char*, bool = true );
void messagec( char sym, bool indent = false );
void progname_and_version( FILE* );

void print_cmdline( TPrintFunction pfunc, void* vpn );

// path processing
const char* my_basename( const char* );
const char* my_dirname( const char* );

// string translation
mystring usage( const char* progname, const char* instructions, const char* version, const char* date );

// printing-to-stream interface
int file_print( void*, const char* format, ... );
int string_print( void*, const char* format, ... );



#endif//__rc__

