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
#include "HDP_SSSScores.h"


// =========================================================================
//global object for application
HDP_SSSScores HDPSSSSCORES;
//

// -------------------------------------------------------------------------
// constructor: initialization
//
HDP_SSSScores::HDP_SSSScores()
:   scores_( NULL ),
    weight_( 0.0 ),
    noeffs_( 0 ),
    noplvs_( 0 ),
    nosssplvs_( 0 ),
    nosss_( 0 ),
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
HDP_SSSScores::~HDP_SSSScores()
{
    DestroyScores();
}

// -------------------------------------------------------------------------
// ReadScores: read scores from file
//
void HDP_SSSScores::ReadScores( const char* filename )
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
void HDP_SSSScores::ReadLevels( FILE* fp )
{
    if( fp == NULL )
        return;

    size_t          length, rbts;
    const size_t    locsize = TIMES4( KBYTE );
    char            locbuffer[locsize+1] = {0};
    const mystring  preambl = "HDP_SSSScores: ReadLevels: ";
    mystring        errmsg;
    const char*     p, *ppp;
    int             emsg;
    const char*     patstrhdpssss = "hdpssss=";
    double  dblval;

    try {
        levels_.Clear();
        //try read levels of eff. no. sequences
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        if(( p = strstr( locbuffer, patstrhdpssss )) == NULL )
            throw myruntime_error("Wrong file format.");

        p += strlen( patstrhdpssss );

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
void HDP_SSSScores::ReadScoreTables( const char* filename )
{
    FILE*       fp = NULL;
    mystring    fname;
    mystring    errstr;
    char        locbuf[KBYTE] = {0};
    int n, m, lvl, noeffs;

    if( !filename )
        return;

    if( levels_.GetSize() < 1 || 10 < levels_.GetSize())
        throw myruntime_error("HDP_SSSScores: ReadScoreTables: "
                              "Either too much (>10) or no eff. level values.");

    noeffs = ( levels_.GetSize() *( levels_.GetSize() + 1 )) >> 1;
    noeffs_ = noeffs;

    lvl = 0;
    for( n = 0; n < levels_.GetSize(); n++ ) {
        for( m = n; m < levels_.GetSize(); m++, lvl++ ) 
        {
            sprintf( locbuf, "%d%d",( int )rint( levels_.GetValueAt(n)),
                    ( int )rint( levels_.GetValueAt(m)));
            fname = filename + mystring( locbuf );
            if(( fp = fopen( fname.c_str(), "r" )) == NULL )
                throw myruntime_error( mystring("Failed to open file ")+ fname );

            try {
                ReadScoresHelper( fp, n, m, lvl );
            } catch( myexception const& ex ) {
                errstr = ex.what();
            }

            fclose( fp );
            if( !errstr.empty())
                throw myruntime_error( errstr );
        }
    }
}

// -------------------------------------------------------------------------
// ReadScoresHelper: read scores from file descriptor
//
void HDP_SSSScores::ReadScoresHelper( FILE* fp, int eff1, int eff2, int lvl )
{
    if( fp == NULL )
        return;

    size_t          length, rbts;
    const size_t    locsize = TIMES4( KBYTE );
    char            locbuffer[locsize+1] = {0};
    const mystring  preambl = "HDP_SSSScores: ReadScoresHelper: ";
    mystring        errmsg;
    const char*     p, *pp, *ppp;
    int             emsg, i, pi, ppi, spi, sppi, si, ssi;
    int             plvl = 0, sssplvl = 0, sssi = 0;
    const char*     patstrhdpssss = "hdpssss";
    const char*     patstrpl = "+";
    const char*     patstrcard = "Cardinality =";
    const char*     patstrx = "x";
    char            patstrbuf[KBYTE];
    const bool      init = eff1 < 1 && eff2 < 1 ;
    double  dblval;
    int     intval;

    if( eff1 < 0 || levels_.GetSize() <= eff1 || eff2 < 0 || levels_.GetSize() <= eff2 )
        throw myruntime_error("HDP_SSSScores: ReadScoresHelper: "
                              "Invalid indices of eff. levels.");
    sprintf( patstrbuf, "%s%d%d=", patstrhdpssss, 
             (int)rint(levels_.GetValueAt(eff1)), (int)rint(levels_.GetValueAt(eff2)));

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
            noplvs_ = ( prblvs_.GetSize() *( prblvs_.GetSize() + 1 )) >> 1;
            nosssplvs_ = ( sssprblvs_.GetSize() *( sssprblvs_.GetSize() + 1 )) >> 1;
        }

        //iterate over probability levels
        for( pi = 0; pi < prblvs_.GetSize(); pi++ ) {
            for( ppi = pi; ppi < prblvs_.GetSize(); ppi++, plvl++ ) 
            {
                for( spi = 0; spi < sssprblvs_.GetSize(); spi++ ) {
                    for( sppi = spi; sppi < sssprblvs_.GetSize(); sppi++, sssplvl++ ) 
                    {
                        sprintf( patstrbuf, "%02d%02d %02d%02d",
                          ( int )rint( prblvs_.GetValueAt(pi)*100.), 
                          ( int )rint( prblvs_.GetValueAt(ppi)*100.),
                          ( int )rint( sssprblvs_.GetValueAt(spi)*100.), 
                          ( int )rint( sssprblvs_.GetValueAt(sppi)*100.));
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

                        if( init && !GetSSSCardinality()) {
                            SetSSSCardinality( intval );
                            nosss_ = ( intval *( intval + 1 )) >> 1;
                        }
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
                            NewScores( noeffs_, noplvs_, nosssplvs_, nosss_ );

                        //iterate over SS states
                        for( si = 0; si < GetSSSCardinality(); si++ ) {
                            for( ssi = 0; ssi <= si; ssi++ ) {
                                ReadSubmatrix( fp, si, ssi, GetScores( lvl, plvl, sssplvl, sssi++ ));
                            }
                        }
                    }
                }
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
// ReadSubmatrix: read a score submatrix from file descriptor
//
void HDP_SSSScores::ReadSubmatrix( FILE* fp, int sss1, int sss2, Pslvector* scos )
{
    if( fp == NULL )
        return;

    size_t          length, rbts, read;
    const size_t    locsize = TIMES4( KBYTE );
    char            locbuffer[locsize+1] = {0};
    mystring        buffer;
    const mystring  preambl = "HDP_SSSScores: ReadSubmatrix: ";
    mystring        errmsg;
    const char*     p, *pp;
    int             emsg;
    const char*     patstrsss1 = "sss1=";
    const char*     patstrsss2 = "sss2=";
    const char*     patstrNA = "NA";
    const int       lenpatstrNA = strlen( patstrNA );
    int     intval, noelems, n, c, cend;
    double  dblval;

    try {
        if( scos == NULL )
            throw myruntime_error("Memory access error.");

        if( card_ < 1 )
            throw myruntime_error("Invalid cardinality.");

        //read SS state indices
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        if(( p = strstr( locbuffer, patstrsss1 )) == NULL )
            throw myruntime_error("Wrong file format.");

        p += strlen( patstrsss1 );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error("Wrong file format.");

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( intval < 0 || ssscard_ <= intval )
            throw myruntime_error("Invalid SS state index.");

        if( intval != sss1 )
            throw myruntime_error("Inconsistent SS state index.");

        if(( p = strstr( locbuffer, patstrsss2 )) == NULL )
            throw myruntime_error("Wrong file format.");

        p += strlen( patstrsss2 );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error("Wrong file format.");

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( intval < 0 || ssscard_ <= intval )
            throw myruntime_error("Invalid SS state index.");

        if( intval != sss2 )
            throw myruntime_error("Inconsistent SS state index.");

        //allocate space for vector
        if( sss1 == sss2 )
            noelems = ( card_ * ( card_+1 )) >> 1;
        else
            noelems = card_ * card_;
        scos->Clear();
        scos->Allocate( noelems );

        //read scores
        for( n = 0; n < card_; n++ )
        {
            if(( emsg = skip_comments( fp, buffer )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( !buffer.length())
                throw myruntime_error("Short of scores.");

            p = buffer.c_str();
            pp = p + buffer.length();
            if( sss1 == sss2 )
                cend = n + 1;
            else
                cend = card_;

            for( c = 0; c < cend; c++, p += rbts )
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

