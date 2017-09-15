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
#include "SSSScores.h"


// =========================================================================
//global object definition
SSSScores SSSSCORES;
//

// -------------------------------------------------------------------------
// constructor: initialization
//
SSSScores::SSSScores()
:   VirtScores(),
    sssweight_(0.0)
{
}

// -------------------------------------------------------------------------
// destructor:
//
SSSScores::~SSSScores()
{
}

// -------------------------------------------------------------------------
// ReadScores: read scores from file
//
void SSSScores::ReadScores( const char* filename )
{
    FILE*       fp = NULL;
    mystring    errstr;
    mystring msg = "Reading ";
    msg += filename;

    if( !filename )
        return;

    if(( fp = fopen( filename, "r" )) == NULL )
        throw myruntime_error( mystring("Failed to open file ")+ filename );

    try {
        message( msg.c_str());
        ReadLevels( fp );
        ReadScoreTables( fp );
    } catch( myexception const& ex ) {
        errstr = ex.what();
    }

    fclose( fp );
    if( !errstr.empty())
        throw myruntime_error( errstr );
}

// -------------------------------------------------------------------------
// ReadLevels: read probability levels and those of eff. no. sequences
//
void SSSScores::ReadLevels( FILE* fp )
{
    if( fp == NULL )
        return;

    size_t          length, rbts, read;
    const size_t    locsize = TIMES4( KBYTE );
    char            locbuffer[locsize+1] = {0};
    const mystring  preambl = "SSSScores: ReadLevels: ";
    mystring        errmsg;
    const char*     p, *pp, *ppp;
    int             emsg;
    const char*     patstrssss = "ssss=";
    const char*     patstrpl = "+";
    double  dblval;

    try {
        prblvs_.Clear();
        levels_.Clear();
        //try read levels of eff. no. sequences
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        if(( p = strstr( locbuffer, patstrssss )) == NULL )
            throw myruntime_error("Wrong file format.");

        p += strlen( patstrssss );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error("Wrong file format.");

        if(( pp = strstr( locbuffer, patstrpl )) != NULL ) {
            //read probability levels
            for( ; p < pp; )
            {
                if(( emsg = read_double( p, length - size_t( p - locbuffer ), &dblval, &rbts )) != 0 ) {
                    if( emsg == ERR_RD_NOVL )
                        break;
                    throw myruntime_error( TranslateReadError( emsg ));
                }
                if( dblval < 0.0 || 100.0 <= dblval )
                    throw myruntime_error("SSSScores: ReadLevels: Invalid probability threshold.");
                prblvs_.Push( dblval/100.0 );
                for( p += rbts; p < pp && (*p == ' ' || *p == '\t'); p++);
            }
            p = pp + strlen( patstrpl );
        }

        ppp = locbuffer + length;
        for( ; p < ppp; p += rbts )
        {
            //read levels of eff. no. sequences
            if(( emsg = read_double( p, length - size_t( p - locbuffer ), &dblval, &rbts )) != 0 ) {
                if( emsg == ERR_RD_NOVL )
                    break;
                throw myruntime_error( TranslateReadError( emsg ));
            }
            levels_.Push( dblval );
        }

    } catch( myexception const& ex ) {
        errmsg = preambl + ex.what();
    }

    if( !errmsg.empty()) {
        throw myruntime_error( errmsg );
    }
}

// -------------------------------------------------------------------------
// ReadScoreTables: read score tables from files
//
void SSSScores::ReadScoreTables( FILE* fp )
{
    size_t          length;
    const size_t    locsize = TIMES4( KBYTE );
    char            locbuffer[locsize+1] = {0};
    const mystring  preambl = "SSSScores: ReadScoreTables: ";
    mystring        errmsg;
    const char*     p;
    int             emsg;
    char            levcmb[KBYTE] = {0};
    int szprblvs = SLC_MAX( 1, prblvs_.GetSize());
    int n, m, k1, k2, plvl, lvl, noplvs, notbls;

    if( fp == NULL )
        return;

    if( szprblvs < 1 || 10 < szprblvs )
        throw myruntime_error("SSSScores: ReadScoreTables: "
                              "Invalid number of probability level values.");
    if( levels_.GetSize() < 1 || 10 < levels_.GetSize())
        throw myruntime_error("SSSScores: ReadScoreTables: "
                              "Either too much (>10) or no eff. level values.");

    noplvs = ( szprblvs *( szprblvs + 1 )) >> 1;
    notbls = ( levels_.GetSize() *( levels_.GetSize() + 1 )) >> 1;
    NewScores( noplvs, notbls );

    plvl = 0;
    for( k1 = 0; k1 < szprblvs; k1++ ) {
        for( k2 = k1; k2 < szprblvs; k2++, plvl++ ) {
            lvl = 0;
            for( n = 0; n < levels_.GetSize(); n++ ) {
                for( m = n; m < levels_.GetSize(); m++ ) 
                {
                    if( prblvs_.GetSize())
                        sprintf( levcmb, "%02d%02d %d%d",
                                ( int )rint( prblvs_.GetValueAt(k1) * 100.0),
                                ( int )rint( prblvs_.GetValueAt(k2) * 100.0),
                                ( int )rint( levels_.GetValueAt(n)),
                                ( int )rint( levels_.GetValueAt(m)));
                    else
                        sprintf( levcmb, "0000 %d%d",
                                ( int )rint( levels_.GetValueAt(n)),
                                ( int )rint( levels_.GetValueAt(m)));

                    try {
                        //read next level values waiting scores
                        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                            throw myruntime_error( TranslateReadError( emsg ));

                        if( feof( fp ) || !length )
                            throw myruntime_error("Wrong file format.");

                        if(( strncmp( levcmb, p = locbuffer, strlen(levcmb))) != 0 )
                            throw myruntime_error("Wrong file format.");

                        p += strlen( levcmb );
                        for( ; p && (*p == ' '||*p == '\t'||*p == '\r'||*p == '\n'); p++);

                        if( size_t( p - locbuffer ) < length )
                            throw myruntime_error("Wrong file format.");
                    } catch( myexception const& ex ) {
                        errmsg = preambl + ex.what();
                    }

                    if( !errmsg.empty()) {
                        throw myruntime_error( errmsg );
                    }

                    ReadScoresHelper( fp, scores_[plvl][lvl++]);
                }
            }
        }
    }
}
