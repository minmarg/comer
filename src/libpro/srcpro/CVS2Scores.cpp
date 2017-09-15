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
#include "CVS2Scores.h"


// =========================================================================
//global object definition
CVS2Scores CVS2SCORES;
//

// -------------------------------------------------------------------------
// constructor: initialization
//
CVS2Scores::CVS2Scores()
:   cvsweight_(0.0),
    keys_( NULL ),
    scores_( NULL ),
    shifts_( NULL ),
    levels_(),
    notbls_( 0 ),
    cards_( NULL ),
    scale_( NULL ),
    step_( 2 ),
    naval_( -9999.0 ),
    plusval_( 9999.0 )
{
}

// -------------------------------------------------------------------------
// destructor:
//
CVS2Scores::~CVS2Scores()
{
    DestroyScores();
}

// -------------------------------------------------------------------------
// ReadScores: read scores from file
//
void CVS2Scores::ReadScores( const char* filename )
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
// ReadLevels: read levels of eff. no. sequences
//
void CVS2Scores::ReadLevels( FILE* fp )
{
    if( fp == NULL )
        return;

    size_t          length, rbts;
    const size_t    locsize = TIMES4( KBYTE );
    char            locbuffer[locsize+1] = {0};
    const mystring  preambl = "CVS2Scores: ReadLevels: ";
    mystring        errmsg;
    const char*     p, *ppp;
    int             emsg;
    const char*     patstrcvs2s = "cvs2s=";
    double  dblval;

    try {
        levels_.Clear();
        //try read levels of eff. no. sequences
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        if(( p = strstr( locbuffer, patstrcvs2s )) == NULL )
            throw myruntime_error("Wrong file format.");

        p += strlen( patstrcvs2s );

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
void CVS2Scores::ReadScoreTables( FILE* fp )
{
    size_t          length;
    const size_t    locsize = TIMES4( KBYTE );
    char            locbuffer[locsize+1] = {0};
    const mystring  preambl = "CVS2Scores: ReadScoreTables: ";
    mystring        errmsg;
    const char*     p;
    int             emsg;
    char            levcmb[KBYTE] = {0};
    int n, m, k1, k2, lvl, notbls;

    if( fp == NULL )
        return;

    if( levels_.GetSize() < 1 || 10 < levels_.GetSize())
        throw myruntime_error("CVS2Scores: ReadScoreTables: "
                              "Either too much (>10) or no eff. level values.");

    notbls = ( levels_.GetSize() *( levels_.GetSize() + 1 )) >> 1;
    NewScores( notbls );

    lvl = 0;
    for( n = 0; n < levels_.GetSize(); n++ ) {
        for( m = n; m < levels_.GetSize(); m++ ) 
        {
            sprintf( levcmb, "%d%d",
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

            ReadScoresHelper( fp, keys_[lvl], scores_[lvl], shifts_+lvl, cards_+lvl, scale_+lvl );
            lvl++;
        }
    }
}

// -------------------------------------------------------------------------
// ReadScoresHelper: read scores from file descriptor
//
void CVS2Scores::ReadScoresHelper( FILE* fp, Pslvector* keys, Pslvector* scos, 
                                   int* shift, int* card, double* scale )
{
    if( fp == NULL )
        return;

    size_t          length, rbts1, rbts2;
    const size_t    locsize = KBYTE;
    char            locbuffer[locsize+1] = {0};
    mystring        buffer1, buffer2;
    const mystring  preambl = "CVS2Scores: ReadScoresHelper: ";
    mystring        errmsg;
    const char*     p;
    const char*     p1, *pp1;
    const char*     p2, *pp2;
    int             emsg;
    const char*     patstrscal = "Scale =";
    const char*     patstrdimn = "Dimensions =";
    const char*     patstrcard = "Cardinality =";
    const char*     patstrplus = "+";
    const char*     patstrNA = "NA";
    const int       lenpatstrNA = strlen( patstrNA );
    const int       lenpatstrplus = strlen( patstrplus );
    int     intval, noelems, dimn, rcard, c;
    double  dbl1, dbl2;

    try {
        if( keys == NULL || scos == NULL || shift == NULL || card == NULL || scale == NULL )
            throw myruntime_error("Memory access error: Null parameters.");

        //read scale
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        if(( p = strstr( locbuffer, patstrscal )) != NULL ) {
            p += strlen( patstrscal );

            if( length <= size_t( p - locbuffer ))
                throw myruntime_error("Wrong file format.");

            if(( emsg = read_double( p, length - size_t( p - locbuffer ), &dbl1, &rbts1 )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( dbl1 <= 0.0 || 1000.0 < dbl1 )
                throw myruntime_error("Invalid scale in file.");
            
            *scale = dbl1;

        //read dimensions
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw myruntime_error("Wrong file format.");
        }

        if(( p = strstr( locbuffer, patstrdimn )) == NULL )
            throw myruntime_error("Wrong file format.");

        p += strlen( patstrdimn );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error("Wrong file format.");

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts1 )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( intval < NUMAA-1 || 100*(NUMAA-1) < intval )
            throw myruntime_error("Invalid dimensions in file.");

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

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts1 )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( intval <= 0.0 || 4000 < intval )
            throw myruntime_error("Invalid value of cardinality.");

        rcard = intval;

        //allocate space for vector
        noelems = rcard + 10;//add extra space for missing values
        keys->Clear(); keys->Allocate( noelems );
        scos->Clear(); scos->Allocate( noelems );

        //read keys and scores
        if(( emsg = skip_comments( fp, buffer1 )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( !buffer1.length())
            throw myruntime_error("No keys.");

        if(( emsg = skip_comments( fp, buffer2 )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( !buffer2.length())
            throw myruntime_error("No scores.");

        p1 = buffer1.c_str(); pp1 = p1 + buffer1.length();
        p2 = buffer2.c_str(); pp2 = p2 + buffer2.length();

        for( c = 0; c < rcard; c++, p1 += rbts1, p2 += rbts2 )
        {
            if( pp1 <= p1 )
                throw myruntime_error("Short of keys.");
            if( pp2 <= p2 )
                throw myruntime_error("Short of scores.");

            if(( emsg = read_double( p1, buffer1.length() - size_t(p1-buffer1.c_str()), &dbl1, &rbts1 )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            //check for '+'
            for( ; *p2 == ' ' || *p2 == '\t'; p2++ );
            if( strncmp( p2, patstrplus, lenpatstrplus ) == 0 ) {
                //value '+'
                rbts2 = lenpatstrplus;
                dbl2 = GetValuePlus();
            }
            else
                if(( emsg = read_double( p2, buffer2.length() - size_t(p2-buffer2.c_str()), &dbl2, &rbts2 )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));

            keys->Push( dbl1 );
            scos->Push( dbl2 );
        }
        SmoothScores( keys, scos, shift, card );

    } catch( myexception const& ex ) {
        errmsg = preambl + ex.what();
    }

    if( !errmsg.empty()) {
        keys->Clear();
        scos->Clear();
        throw myruntime_error( errmsg );
    }
}

// -------------------------------------------------------------------------
// SmoothScores: smooth scores and add missing keys
//
void CVS2Scores::SmoothScores( Pslvector* keys, Pslvector* scos, int* shift, int* card )
{
    const mystring  preambl = "CVS2Scores: SmoothScores: ";
    mystring        errmsg;

    if( GetStep() < 1 )
        throw myruntime_error("CVS2Scores: SmoothScores: Invalid step.");

    double  incr = 1./(double)GetStep();//increment
    double  pkk, kk, sc, sm, rsum, avg;
    const int   nSWIN = 5;//window length for calculating average
    const int   nSLAST = 9999;//17;//70;//7;//number of last elements taking smoothed values
    const int   nSWINh = nSWIN >> 1;
    double  rwin[nSWIN];
    int n, nv, rn, p, nn;

    try {
        if( keys == NULL || scos == NULL || shift == NULL || card == NULL )
            throw myruntime_error("Memory access error: Null parameters.");
        if( keys->GetSize() != scos->GetSize())
            throw myruntime_error("Inconsistent keys and scores.");
        if( keys->GetSize() < 1 )
            return;

        Pslvector   smth( scos->GetSize());

        //smooth scores
        memset( rwin, 0, nSWIN * sizeof( double ));
        rsum = 0.0;
        for( n = nv = rn = nn = 0; n < smth.GetSize(); n++ ) {
            sc = scos->GetValueAt(n);
            if( ScoreNA(sc))
                throw myruntime_error("Score NA.");
            if( nv < nSWINh )
                nn++;
            if( !ScorePlus(sc)) {
                if( nSWIN <= rn )
                    rn = 0;
                rsum -= rwin[rn];
                rwin[rn++] = sc;
                rsum += sc;
                avg = rsum/(double)nSWIN;
                nv++;
            }
            if( nSWIN <= nv ) {
                if( !ScorePlus(sc)) {
                    if( nSWIN == nv )
                        for( p = 0; p < nn; p++ )
                            smth.SetValueAt( p, avg );
                    smth.SetValueAt( nn, avg );
                    for( ++nn; ScorePlus(scos->GetValueAt(nn)) && nn < n; nn++ )
                        smth.SetValueAt( nn, avg );
                }
                if( smth.GetSize() == n+1 )
                    for( p = nn; p < smth.GetSize(); p++ )
                        smth.SetValueAt( p, avg );
            }
        }//for(n<smth.GetSize())

        //smooth end scores
        for( n = smth.GetSize()-1; 0 <= n && smth.GetSize()-nSLAST <= n; n-- ) {
            sm = smth.GetValueAt(n);
            scos->SetValueAt( n, sm );
        }

        //insert missing keys along with smoothed scores
        for( n = 0; n < keys->GetSize(); n++ ) {
            kk = keys->GetValueAt(n);
            sc = scos->GetValueAt(n);
            if( n && kk <= pkk )
                throw myruntime_error("Keys not sorted.");
            if( n && incr < kk-pkk ) {
                pkk += incr;
                keys->InsertAt( n, pkk );
                scos->InsertAt( n, sm );
                smth.InsertAt( n, sm );
                continue;
            }
            else
                pkk = kk;
            sm = smth.GetValueAt(n);
            if( ScorePlus(sc))
                scos->SetValueAt( n, sm );
        }

        *card = scos->GetSize();
        *shift = 0;
        for( n = 0; n < keys->GetSize(); n++ ) {
            if( 0.0 <= keys->GetValueAt(n))
                break;
            (*shift)++;
        }

    } catch( myexception const& ex ) {
        errmsg = preambl + ex.what();
    }

    if( !errmsg.empty())
        throw myruntime_error( errmsg );
}

