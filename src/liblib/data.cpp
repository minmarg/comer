/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include <ctype.h>
#include "rc.h"
#include "mystring.h"
#include "myexcept.h"
#include "needconfig.h"
#include "data.h"

// -------------------------------------------------------------------------
// size of alphabet of symbols a sequence can consist of
//
const int NumAlphabet()
{
    return strlen( gAlphabet );
}

// -------------------------------------------------------------------------
// DehashCode: returns ascii residue symbol given its code
//
char DehashCode( unsigned char in )
{
    if( in < NUMALPH )
        return gAlphabet[in];

    throw myruntime_error( mystring( "Unrecognized residue hash." ));
}

// -------------------------------------------------------------------------
// HashResidue: get code for alphabet symbol
//
int HashAlphSymbol( char in )
{
    char    upper = toupper( in );

    for( int n = 0; gAlphabet[n]; n++ )
        if( upper == gAlphabet[n] )
            return n;

    char    instr[] = { in, 0 };
    if( upper == 'J' || upper == 'O' || upper == 'U' ) {
        warning(( mystring( "Residue '" ) + instr + mystring( "' replaced by X" )).c_str(), false );
        return X;
    }
    throw myruntime_error( mystring( "Unrecognized residue: '" ) + instr + mystring( "'." ));
}

// -------------------------------------------------------------------------
// ResToInt: converts amino acid letter to the corresponding numeric value
//
int ResToInt( char in )
{
    for( int n = 0; gAAcids[n]; n++ )
        if( in == gAAcids[n])
            return n;

    throw myruntime_error( mystring( "Unrecognized residue." ));
}

// -------------------------------------------------------------------------
// IsValidRes: check if `code' codes for a valid residue
//
bool IsValidResSym( unsigned char code )
{
    if( code < ASTERISK )
        return true;
    return false;
}


// =========================================================================
// DehashSSCode: returns ascii SS state symbol given its code
//
char DehashSSCode( unsigned char in )
{
    if( in < SS_NSTATES )
        return gSSAlphabet[in];
    throw myruntime_error("Unrecognized SS state.");
}

// -------------------------------------------------------------------------
// HashSSState: get code for SS state symbol
//
int HashSSState( char in )
{
    char    upper = toupper( in );

    for( int n = 0; gSSAlphabet[n]; n++ )
        if( upper == gSSAlphabet[n])
            return n;

    char    instr[] = { in, 0 };
    throw myruntime_error( mystring("Unrecognized SS state: '") + instr + mystring("'."));
}

