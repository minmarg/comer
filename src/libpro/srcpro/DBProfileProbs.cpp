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
#include "DBProfileProbs.h"


// -------------------------------------------------------------------------
// constructor: initialization
//
DBProfileProbs::DBProfileProbs()
:   probs_( NULL ),
    card_( 0 )
{
    double  dv;
    int  n, nv;
    int  card;

    efflvs_.Allocate(10);
    lenlvs_.Allocate(10);
    midefflvs_.Allocate(10);
    midlenlvs_.Allocate(10);

    for( dv = 2.; dv <= 14.; dv += 2. ) efflvs_.Push(dv);
    for( nv = 50; nv <= 400; nv *= 2 ) lenlvs_.Push(nv);
    lenlvs_.Push(600);
    lenlvs_.Push(800);

    for( n = 0; n < efflvs_.GetSize()-1; n++ ) {
        dv = (efflvs_.GetValueAt(n) + efflvs_.GetValueAt(n+1)) * .5;
        midefflvs_.Push(dv);
    }
    for( n = 0; n < lenlvs_.GetSize()-1; n++ ) {
        nv = (lenlvs_.GetValueAt(n) + lenlvs_.GetValueAt(n+1)) / 2;
        midlenlvs_.Push(nv);
    }

    card = efflvs_.GetSize() * lenlvs_.GetSize();

    margs_.Reserve( card );

    SetCardinality( card );
}

// -------------------------------------------------------------------------
// destructor:
//
DBProfileProbs::~DBProfileProbs()
{
    DestroyProbs();
}



// -------------------------------------------------------------------------
// ReadProbs: read probabilities from file
//
void DBProfileProbs::ReadProbs( const char* filename )
{
    FILE*       fp = NULL;
    mystring    errstr;
    mystring msg = "Reading probabilities...";
//     msg += filename;

    if( !filename )
        return;

    if(( fp = fopen( filename, "r" )) == NULL ) {
//         throw myruntime_error( mystring("Failed to open file ")+ filename );
        warning("Update of profile pair probabilities ignored.");
        return;
    }

    try {
        message( msg.c_str());
        NewProbs();
        ReadProbsHelper( fp, GetProbs());
    } catch( myexception const& ex ) {
        errstr = ex.what();
    }

    fclose( fp );
    if( !errstr.empty())
        throw myruntime_error( errstr );
}

// -------------------------------------------------------------------------
// ReadProbsHelper: helper method for probabilities read from file 
//  descriptor
//
void DBProfileProbs::ReadProbsHelper( FILE* fp, Pslvector* probs )
{
    if( fp == NULL )
        return;

    size_t          length, rbts;
    const size_t    locsize = TIMES4( KBYTE );
    char            locbuffer[locsize+1] = {0};
    mystring        buffer;
    const mystring  preambl = "DBProfileProbs: ReadProbsHelper: ";
    mystring        errmsg;
    const char*     p, *pp;
    int             emsg;
    const char*     patstrcard = "Cardinality =";
    const char*     patstrNA = "NA";
    const int       lenpatstrNA = strlen( patstrNA );
    int     intval, noelems, card, n, c;
    double  dblval;

    try {
        if( probs == NULL )
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

        if( intval <= 0.0 || 100 < intval )
            throw myruntime_error("Invalid value of cardinality.");

        if( intval != efflvs_.GetSize() * lenlvs_.GetSize())
            throw myruntime_error("Inconsistent cardinality.");

        card = intval;
        if(!GetCardinality())
            SetCardinality( card );
        else if( card != GetCardinality())
            throw myruntime_error("Inconsistent probabilities table.");

        //allocate space for vector
        noelems = card + (( card * ( card-1 )) >> 1 );
        probs->Clear();
        probs->Allocate( noelems );

        //read probabilities
        for( n = 0; n < card; n++ )
        {
            if(( emsg = skip_comments( fp, buffer )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( !buffer.length())
                throw myruntime_error("Short of probabilities.");

            p = buffer.c_str();
            pp = p + buffer.length();
            for( c = 0; c <= n; c++, p += rbts )
            {
                if( pp <= p )
                    throw myruntime_error("Short of probabilities.");

                //check for NA
                for( ; *p == ' ' || *p == '\t'; p++ );
                if( strncmp( p, patstrNA, lenpatstrNA ) == 0 ) {
                    //NA value
                    rbts = lenpatstrNA;
                    probs->Push(0.);
                    continue;
                }

                if(( emsg = read_double( p, buffer.length() - size_t( p - buffer.c_str()), &dblval, &rbts )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));

                probs->Push( dblval );
            }
        }

        VerifyProbs( probs, 1.e-3 );

    } catch( myexception const& ex ) {
        errmsg = preambl + ex.what();
    }

    if( !errmsg.empty()) {
        probs->Clear();
        throw myruntime_error( errmsg );
    }
}



// -------------------------------------------------------------------------
// UpdateMargs: update marginal counts of profile pairs
//
void DBProfileProbs::UpdateMargs( double E1, int L1 )
{
    const mystring  preambl = "DBProfileProbs: UpdateMargs: ";
    int     e1, l1;//level indices
    int     row;

    if( efflvs_.GetSize() < 1 || lenlvs_.GetSize() < 1 )
        throw myruntime_error( preambl + "No levels set.");

    for( e1 = 0; e1 < midefflvs_.GetSize(); e1++ )
        if( E1 < midefflvs_.GetValueAt(e1))
            break;
    for( l1 = 0; l1 < midlenlvs_.GetSize(); l1++ )
        if( L1 < midlenlvs_.GetValueAt(l1))
            break;

    row = e1 * lenlvs_.GetSize() + l1;
    if( margs_.GetSize() <= row )
        throw myruntime_error( preambl + "Memory access error.");

    margs_.AddValueAt( row, 1 );
}

// -------------------------------------------------------------------------
// CalcProbs: calculate probabilities given marginal pair counts
//
void DBProfileProbs::CalcProbs()
{
    const mystring  preambl = "DBProfileProbs: CalcProbs: ";
    const int   card = GetCardinality();
    int         noelems, rowcnt, colcnt;
    int         count = 0;
    int         v, n, c;
    double      count2 = 0.;
    double      val;

    if( card < 1 )
        return;

    if( card != margs_.GetSize())
        throw myruntime_error( preambl + "Inconsistent data sizes.");

    noelems = card + (( card * ( card-1 )) >> 1 );
    NewProbs();
    probs_->Allocate( noelems );

    for( n = 0; n < margs_.GetSize(); n++ )
        count += margs_.GetValueAt(n);

    if( count < 1 )
        throw myruntime_error( preambl + "No marginal count data.");

    count2 = (double)count * (double)count;

    for( v = n = 0; n < card; n++ ) {
        rowcnt = margs_.GetValueAt(n);
        for( c = 0; c <= n; c++, v++ )
        {
            colcnt = margs_.GetValueAt(c);
            if( rowcnt < 1 || colcnt < 1 )
                val = 0.;
            else
                val = (double)rowcnt * (double)colcnt / count2;
            if( c < n )
                val *= 2.;
            probs_->Push(val);
        }
    }

    VerifyProbs( probs_ );
}



// -------------------------------------------------------------------------
// VerifyProbs: verify probabilities
//
void DBProfileProbs::VerifyProbs( Pslvector* probs, double acc )
{
    const mystring  preambl = "DBProfileProbs: VerifyProbs: ";
    static char     locbuf[KBYTE];
    double          sum;

    if( probs == NULL )
        throw myruntime_error( preambl + "Null probabilities.");

    sum = probs->Sum() - 1.;
    if( acc < fabs(sum)) {
        sprintf( locbuf, "%s Probabilities not conserved: %g", preambl.c_str(), sum+1.);
        warning( locbuf );
    }
}



// -------------------------------------------------------------------------
// WriteProbs: write probabilities to file
//
void DBProfileProbs::WriteProbs( const char* filename )
{
    FILE*       fp = NULL;
    mystring    errstr;

    if( !filename )
        return;

    if( GetProbs() == NULL || GetProbs()->GetSize() < 1 )
        return;

    if(( fp = fopen( filename, "w" )) == NULL )
        throw myruntime_error( mystring("Failed to open file for writing: ")+ filename );

    try {
        WriteProbsHelper( fp, GetProbs());
    } catch( myexception const& ex ) {
        errstr = ex.what();
    }

    fclose( fp );
    if( !errstr.empty())
        throw myruntime_error( errstr );
}

// -------------------------------------------------------------------------
// WriteProbsHelper: helper method for probabilities write to file 
//  descriptor
//
void DBProfileProbs::WriteProbsHelper( FILE* fp, Pslvector* probs )
{
    if( fp == NULL )
        return;

    size_t          length, rbts;
    const size_t    locsize = TIMES4( KBYTE );
    char            locbuffer[locsize+1] = {0};
    mystring        buffer;
    const mystring  preambl = "DBProfileProbs: WriteProbsHelper: ";
    mystring        errmsg;
    const char*     p, *pp;
    int             emsg;
    const char*     patstrcard = "Cardinality =";
    const char*     patstrNA = "NA";
    const int       lenpatstrNA = strlen( patstrNA );
    const int       card = GetCardinality();
    int     noelems, v, n, c;
    double  val;

    if( probs == NULL )
        throw myruntime_error( preambl + "Memory access error.");

    noelems = card + (( card * ( card-1 )) >> 1 );

    if( noelems != probs->GetSize())
        throw myruntime_error( preambl + "Inconsistent vector size.");

    fprintf( fp, "%s %d\n", patstrcard, card );

    //write probabilities
    for( v = n = 0; n < card; n++ )
    {
        for( c = 0; c <= n; c++, v++ )
        {
            val = probs->GetValueAt(v);
            if( val == 0.0 )
                fprintf( fp, " %8s", patstrNA );
            else {
                if( val < 5.e-7 )
                    fprintf( fp, " %8d", 0 );
                else
                    fprintf( fp, " %8.6f", val );
            }
        }
        fprintf( fp, "\n");
    }
}
