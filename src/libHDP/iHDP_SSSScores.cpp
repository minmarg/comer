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
#include "iHDP_SSSScores.h"


// =========================================================================
//global object for application
iHDP_SSSScores iHDPSSSSCORES;
//

// -------------------------------------------------------------------------
// constructor: initialization
//
iHDP_SSSScores::iHDP_SSSScores()
:   scores_( NULL ),
    weight_( 0.0 ),
    noeffs_( 0 ),
    noplvs_( 0 ),
    nosssplvs_( 0 ),
    levels_(),
    prblvs_(),
    sssprblvs_(),
    naval_( -9999.0 ),
    card_( 0 ),
    ssscard_( 0 )
{
}

// -------------------------------------------------------------------------
// destructor:
//
iHDP_SSSScores::~iHDP_SSSScores()
{
    DestroyScores();
}

// -------------------------------------------------------------------------
// ReadScores: read scores from file
//
void iHDP_SSSScores::ReadScores( const char* filename )
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
void iHDP_SSSScores::ReadLevels( FILE* fp )
{
    if( fp == NULL )
        return;

    size_t          length, rbts;
    const size_t    locsize = TIMES4( KBYTE );
    char            locbuffer[locsize+1] = {0};
    const mystring  preambl = "iHDP_SSSScores: ReadLevels: ";
    mystring        errmsg;
    const char*     p, *ppp;
    int             emsg;
    const char*     patstrihdpssss = "ihdpssss=";
    double  dblval;

    try {
        levels_.Clear();
        //try read levels of eff. no. sequences
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        if(( p = strstr( locbuffer, patstrihdpssss )) == NULL )
            throw myruntime_error("Wrong file format.");

        p += strlen( patstrihdpssss );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error("Wrong file format.");

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
void iHDP_SSSScores::ReadScoreTables( const char* filename )
{
    FILE*       fp = NULL;
    mystring    fname;
    mystring    errstr;
    char        locbuf[KBYTE] = {0};
    int n;

    if( !filename )
        return;

    if( levels_.GetSize() < 1 || 10 < levels_.GetSize())
        throw myruntime_error("iHDP_SSSScores: ReadScoreTables: "
                              "Either too much (>10) or no eff. level values.");

    noeffs_ = levels_.GetSize();

    for( n = 0; n < levels_.GetSize(); n++ ) {
        sprintf( locbuf, "%d",( int )rint( levels_.GetValueAt(n)));
        fname = filename + mystring( locbuf );
        if(( fp = fopen( fname.c_str(), "r" )) == NULL )
            throw myruntime_error( mystring("Failed to open file ")+ fname );

        try {
            ReadScoresHelper( fp, n );
        } catch( myexception const& ex ) {
            errstr = ex.what();
        }

        fclose( fp );
        if( !errstr.empty())
            throw myruntime_error( errstr );
    }
}

// -------------------------------------------------------------------------
// ReadScores: read scores from file descriptor
//
void iHDP_SSSScores::ReadScoresHelper( FILE* fp, int eff )
{
    if( fp == NULL )
        return;

    size_t          length, rbts;
    const size_t    locsize = TIMES4( KBYTE );
    char            locbuffer[locsize+1] = {0};
    const mystring  preambl = "iHDP_SSSScores: ReadScoresHelper: ";
    mystring        errmsg;
    const char*     p, *pp, *ppp;
    int             emsg, i, pi, spi;
    int             plvl = 0, sssplvl = 0, sssi = 0;
    const char*     patstrihdpssss = "indhdpssss";
    const char*     patstrpl = "+";
    const char*     patstrcard = "Cardinality =";
    const char*     patstrx = "x";
    char            patstrbuf[KBYTE];
    const bool      init = eff < 1;
    double  dblval;
    int     intval;

    if( eff < 0 || levels_.GetSize() <= eff )
        throw myruntime_error("iHDP_SSSScores: ReadScoresHelper: "
                              "Invalid indices of eff. levels.");
    sprintf( patstrbuf, "%s%d=", patstrihdpssss,(int)rint(levels_.GetValueAt(eff)));

    try {
        if( init ) {
            prblvs_.Clear();
            sssprblvs_.Clear();
        }
        //try read probability levels
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        if(( p = strstr( locbuffer, patstrbuf )) == NULL )
            throw myruntime_error("Wrong file format.");

        p += strlen( patstrbuf );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error("Wrong file format.");

        if(( pp = strstr( locbuffer, patstrpl )) == NULL )
            throw myruntime_error("Wrong file format.");

        //read probability levels
        for( i = 0; p < pp; i++ )
        {
            if(( emsg = read_double( p, length - size_t( p - locbuffer ), &dblval, &rbts )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( dblval < 0.0 || 100.0 <= dblval )
                throw myruntime_error("Invalid probability threshold.");
            if( init )
                prblvs_.Push( dblval/100.0 );
            else {
                if( prblvs_.GetValueAt(i) != dblval/100.0 )
                    throw myruntime_error("Inconsistent probability threshold.");
            }
            for( p += rbts; p < pp && (*p == ' ' || *p == '\t'); p++);
        }
        p = pp + strlen( patstrpl );

        //read SS state probability levels
        ppp = locbuffer + length;
        for( i = 0; p < ppp; p += rbts, i++ )
        {
            if(( emsg = read_double( p, length - size_t( p - locbuffer ), &dblval, &rbts )) != 0 ) {
                if( emsg == ERR_RD_NOVL )
                    break;
                throw myruntime_error( TranslateReadError( emsg ));
            }
            if( dblval < 0.0 || 100.0 <= dblval )
                throw myruntime_error("Invalid SS state probability threshold.");
            if( init )
                sssprblvs_.Push( dblval/100.0 );
            else {
                if( sssprblvs_.GetValueAt(i) != dblval/100.0 )
                    throw myruntime_error("Inconsistent SS state probability threshold.");
            }
        }

        if( init ) {
            noplvs_ = prblvs_.GetSize();
            nosssplvs_ = sssprblvs_.GetSize();
        }

        //iterate over probability levels
        for( pi = 0; pi < prblvs_.GetSize(); pi++ ) {
            for( spi = 0; spi < sssprblvs_.GetSize(); spi++ )
            {
                sprintf( patstrbuf, "%02d %02d",
                  ( int )rint( prblvs_.GetValueAt(pi)*100.), 
                  ( int )rint( sssprblvs_.GetValueAt(spi)*100.));
                ;
                if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));

                if( feof( fp ) || !length )
                    throw myruntime_error("Wrong file format.");

                if(( p = strstr( locbuffer, patstrbuf )) == NULL )
                    throw myruntime_error("Wrong file format.");

                p += strlen( patstrbuf );

                if( length <= size_t( p - locbuffer ))
                    throw myruntime_error("Wrong file format.");

                //Cardinality
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

                if( intval <= 0 || 30 < intval )
                    throw myruntime_error("Invalid value of SS state cardinality.");

                if( init && !GetSSSCardinality())
                    SetSSSCardinality( intval );
                else if( intval != GetSSSCardinality())
                    throw myruntime_error("Inconsistent SS state cardinality.");

                if(( p = strstr( p + rbts, patstrx )) == NULL )
                    throw myruntime_error("Wrong file format: Cardinality.");

                p += strlen( patstrx );

                if( length <= size_t( p - locbuffer ))
                    throw myruntime_error("Wrong file format: Cardinality.");

                if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));

                if( intval <= 0 || 4000 < intval )
                    throw myruntime_error("Invalid value of cardinality.");

                if( init && !GetCardinality())
                    SetCardinality( intval );
                else if( intval != GetCardinality())
                    throw myruntime_error("Inconsistent cardinality.");

                if( init && scores_ == NULL )
                    NewScores( noeffs_, noplvs_, nosssplvs_ );

                //read matrix
                ReadMatrix( fp, GetScores( eff, pi, spi ));
            }
        }

    } catch( myexception const& ex ) {
        errmsg = preambl + ex.what();
    }

    if( !errmsg.empty()) {
        throw myruntime_error( errmsg );
    }
}

// -------------------------------------------------------------------------
// ReadMatrix: read a score matrix from file descriptor
//
void iHDP_SSSScores::ReadMatrix( FILE* fp, Pslvector* scos )
{
    if( fp == NULL )
        return;

    size_t          rbts;
    mystring        buffer;
    const mystring  preambl = "iHDP_SSSScores: ReadMatrix: ";
    mystring        errmsg;
    const char*     p, *pp;
    int             emsg;
    const char*     patstrNA = "NA";
    const int       lenpatstrNA = strlen( patstrNA );
    int     noelems, n, c;
    double  dblval;

    try {
        if( scos == NULL )
            throw myruntime_error("Memory access error.");

        if( ssscard_ < 1 )
            throw myruntime_error("Invalid SS state cardinality.");

        if( card_ < 1 )
            throw myruntime_error("Invalid cardinality.");

        //allocate space for vector
        noelems = ssscard_ * card_;
        scos->Clear();
        scos->Allocate( noelems );

        //read scores
        for( n = 0; n < ssscard_; n++ )
        {
            if(( emsg = skip_comments( fp, buffer )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( !buffer.length())
                throw myruntime_error("Short of scores.");

            p = buffer.c_str();
            pp = p + buffer.length();

            for( c = 0; c < card_; c++, p += rbts )
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

