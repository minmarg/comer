/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __data_h__
#define __data_h__

#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "defines.h"
#include "pcmath.h"

// effective number of residues (different amino acid types)
#define NUMAA   20

#define NUMALPH ( gNumAlphabet )
#define MINFREQ 0
#define MAXFREQ 100
//count of discrete frequency values (boundary values included)
#define NUMFREQ 101


// Global functions interface
//
char DehashCode( unsigned char );
int HashAlphSymbol( char );
int ResToInt( char );
const int NumAlphabet();
bool IsValidResSym( unsigned char );
//
char DehashSSCode( unsigned char );
int HashSSState( char );
//

//                                012345678901234567890123
static const char*     gAAcids = "ARNDCQEGHILKMFPSTWYVBZX*";

//                                  0123456789012345678901234
static const char*     gAlphabet = "ARNDCQEGHILKMFPSTWYVBZX*-";
static const int       gNumAlphabet = 25;

#define X           22
#define ASTERISK    23
#define GAP         24

//{{SS states
enum SSSTATES {
    SS_C,
    SS_E,
    SS_H,
    SS_NSTATES
};
static const char*  gSSAlphabet = "CEH";
//}}

typedef TScore ( *SCORE_ARRAY )[NUMALPH];
typedef const TScore ( *CONST_SCORE_ARRAY )[NUMALPH];

typedef double ( *DOUBLE_ARRAY )[NUMALPH];
typedef const double ( *CONST_DOUBLE_ARRAY )[NUMALPH];


#endif//__data_h__
