/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __faccess_h__
#define __faccess_h__

#include <stdio.h>

class mystring;

// file access routines
bool file_exists( const char* );
int skip_comments( FILE* fp, char* buffer, size_t bufsize, size_t* readlen, char cc = '#' );
int skip_comments( FILE* fp, mystring& buffer, char cc = '#' );
int read_double( const char* readfrom, size_t readlen, double* membuf, size_t* rbytes );
int read_integer( const char* readfrom, size_t readlen, int* membuf, size_t* rbytes );
int read_llinteger( const char* readfrom, size_t readlen, long long int* membuf, size_t* rbytes );
int read_symbol( const char* readfrom, size_t readlen, char* membuf, size_t* rbytes );

#define ERR_SC_MACC ( 111 )
#define ERR_SC_READ ( 113 )
#define ERR_RD_MACC ( 115 )
#define ERR_RD_NOVL ( 117 )
#define ERR_RD_INVL ( 119 )
#define ERR_RI_MACC ( 131 )
#define ERR_RI_NOVL ( 133 )
#define ERR_RI_INVL ( 135 )
#define ERR_RL_MACC ( 137 )
#define ERR_RL_NOVL ( 139 )
#define ERR_RL_INVL ( 151 )
#define ERR_RS_MACC ( 153 )
#define ERR_RS_INVL ( 155 )

//translation of error messages
inline const char* TranslateReadError( int code )
{
    switch( code ) {
        case 0: return "OK";
        case ERR_SC_MACC: return "skip_comments: Memory access error";
        case ERR_SC_READ: return "skip_comments: Reading error";
        case ERR_RD_MACC: return "read_double: Memory access error";
        case ERR_RD_NOVL: return "read_double: No read double value";
        case ERR_RD_INVL: return "read_double: Invalid double value";
        case ERR_RI_MACC: return "read_integer: Memory access error";
        case ERR_RI_NOVL: return "read_integer: No integer read";
        case ERR_RI_INVL: return "read_integer: Invalid integer value";
        case ERR_RL_MACC: return "read_llinteger: Memory access error";
        case ERR_RL_NOVL: return "read_llinteger: No integer read";
        case ERR_RL_INVL: return "read_llinteger: Invalid integer value";
        case ERR_RS_MACC: return "read_symbol: Memory access error";
        case ERR_RS_INVL: return "read_symbol: Not a single symbol";
    }
    return "Unknown";
}

#endif//__faccess_h__
