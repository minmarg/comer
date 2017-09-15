/***************************************************************************
 *   Copyright (C) 2013 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "data.h"
#include "ext/psl.h"
#include "VirtScores.h"


// -------------------------------------------------------------------------
// constructor: initialization
//
VirtScores::VirtScores()
:   scores_( NULL ),
    noplvs_( 0 ),
    notbls_( 0 ),
    prblvs_(),
    levels_(),
    naval_( -9999.0 ),
    card_( 0 )
{
}

// -------------------------------------------------------------------------
// destructor:
//
VirtScores::~VirtScores()
{
    DestroyScores();
}

// -------------------------------------------------------------------------
// ReadScores: read scores from file descriptor
//
void VirtScores::ReadScoresHelper( FILE* fp, Pslvector* scos )
{
    if( fp == NULL )
        return;

    size_t          length, rbts, read;
    const size_t    locsize = TIMES4( KBYTE );
    char            locbuffer[locsize+1] = {0};
    mystring        buffer;
    const mystring  preambl = "VirtScores: ReadScoresHelper: ";
    mystring        errmsg;
    const char*     p, *pp;
    int             emsg;
    const char*     patstrcard = "Cardinality =";
    const char*     patstrNA = "NA";
    const int       lenpatstrNA = strlen( patstrNA );
    int     intval, noelems, card, n, c;
    double  dblval;

    try {
        if( scos == NULL )
            throw myruntime_error("Memory access error.");

        //read cardinality
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        if(( p = strstr( locbuffer, patstrcard )) == NULL )
            throw myruntime_error("Wrong file format.");

        p += strlen( patstrcard );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error("Wrong file format.");

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( intval <= 0.0 || 4000 < intval )
            throw myruntime_error("Invalid value of cardinality.");

        card = intval;
        if(!GetCardinality())
            SetCardinality( card );
        else if( card != GetCardinality())
            throw myruntime_error("Inconsistent score tables.");

        //allocate space for vector
        noelems = card + (( card * ( card-1 )) >> 1 );
        scos->Clear();
        scos->Allocate( noelems );

        //read scores
        for( n = 0; n < card; n++ )
        {
            if(( emsg = skip_comments( fp, buffer )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( !buffer.length())
                throw myruntime_error("Short of scores.");

            p = buffer.c_str();
            pp = p + buffer.length();
            for( c = 0; c <= n; c++, p += rbts )
            {
                if( pp <= p )
                    throw myruntime_error("Short of scores.");

                //check for NA
                for( ; *p == ' ' || *p == '\t'; p++ );
                if( strncmp( p, patstrNA, lenpatstrNA ) == 0 ) {
                    //NA value
                    rbts = lenpatstrNA;
                    scos->Push( GetNAvalue());
                    continue;
                }

                if(( emsg = read_double( p, buffer.length() - size_t( p - buffer.c_str()), &dblval, &rbts )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));

                scos->Push( dblval );
            }
        }
    } catch( myexception const& ex ) {
        errmsg = preambl + ex.what();
    }

    if( !errmsg.empty()) {
        scos->Clear();
        throw myruntime_error( errmsg );
    }
}

