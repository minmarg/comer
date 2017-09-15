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
#include "HDPscores.h"


// =========================================================================
//global object for application
HDPscores HDPSCORES;
HDPscores HDPctSCORES;
//

// -------------------------------------------------------------------------
// constructor: initialization
//
HDPscores::HDPscores()
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
HDPscores::~HDPscores()
{
    DestroyScores();
}

// -------------------------------------------------------------------------
// ReadScores: read scores from file
//
void HDPscores::ReadScores( const char* filename )
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
    } catch( myexception const& ex ) {
        errstr = ex.what();
    }

    fclose( fp );
    if( !errstr.empty())
        throw myruntime_error( errstr );

    ReadScoreTables( filename );
}

// -------------------------------------------------------------------------
// ReadScores: read probability levels and those of eff. no. sequences
//
void HDPscores::ReadLevels( FILE* fp )
{
    if( fp == NULL )
        return;

    size_t          length, rbts, read;
    const size_t    locsize = TIMES4( KBYTE );
    char            locbuffer[locsize+1] = {0};
    const mystring  preambl = "HDPscores: ReadLevels: ";
    mystring        errmsg;
    const char*     p, *pp, *ppp;
    int             emsg;
    const char*     patstrhdps = "hdps=";
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

        if(( p = strstr( locbuffer, patstrhdps )) == NULL )
            throw myruntime_error("Wrong file format.");

        p += strlen( patstrhdps );

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
                    throw myruntime_error("HDPscores: ReadLevels: Invalid probability threshold.");
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
void HDPscores::ReadScoreTables( const char* filename )
{
    FILE*       fp = NULL;
    mystring    fname;
    mystring    errstr;
    char        locbuf[KBYTE] = {0};
    int szprblvs = SLC_MAX( 1, prblvs_.GetSize());
    int n, m, k1, k2, plvl, lvl, noplvs, notbls;

    if( !filename )
        return;

    if( szprblvs < 1 || 10 < szprblvs )
        throw myruntime_error("HDPscores: ReadScoreTables: "
                              "Invalid number of hdps probability level values.");
    if( levels_.GetSize() < 1 || 10 < levels_.GetSize())
        throw myruntime_error("HDPscores: ReadScoreTables: "
                              "Either too much (>10) or no hdps level values.");

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
                        sprintf( locbuf, "%d%d%d%d",
                                ( int )rint( prblvs_.GetValueAt(k1) * 100.0),
                                ( int )rint( prblvs_.GetValueAt(k2) * 100.0),
                                ( int )rint( levels_.GetValueAt(n)),
                                ( int )rint( levels_.GetValueAt(m)));
                    else
                        sprintf( locbuf, "%d%d",( int )rint( levels_.GetValueAt(n)),
                                ( int )rint( levels_.GetValueAt(m)));
                    fname = filename + mystring( locbuf );
                    if(( fp = fopen( fname.c_str(), "r" )) == NULL )
                        throw myruntime_error( mystring("Failed to open file ")+ fname );

                    try {
                        ReadScoresHelper( fp, scores_[plvl][lvl++]);
                    } catch( myexception const& ex ) {
                        errstr = ex.what();
                    }

                    fclose( fp );
                    if( !errstr.empty())
                        throw myruntime_error( errstr );
                }
            }
        }
    }
}

// -------------------------------------------------------------------------
// ReadScores: read scores from file descriptor
//
void HDPscores::ReadScoresHelper( FILE* fp, Pslvector* scos )
{
    if( fp == NULL )
        return;

    size_t          length, rbts, read;
    const size_t    locsize = TIMES4( KBYTE );
    char            locbuffer[locsize+1] = {0};
    mystring        buffer;
    const mystring  preambl = "HDPscores: ReadScoresHelper: ";
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

