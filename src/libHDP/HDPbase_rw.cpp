/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "data.h"
#include "rc.h"
#include "HDPbase.h"
#include <libpro/srcpro/datapro.h>

//default menu size
const int gcsz_defmenusz = 100;

// =========================================================================
// Read: read set of vectors from file; vectors are supposed to be divided 
//  into groups
//  dids -- dish indices
//
void HDPbase::ReadGroups( const char* filename, Ivector* dids )
{
    FILE*           fp = NULL;
    size_t          length, rbts, t, n, c;
    const size_t    locsize = 10 * KBYTE;
    char            locbuffer[locsize+1] = {0};
    char*           p;
    int             emsg;

    const char*     patstrlpd = "Log Prob of Data =";
    const char*     patstrgss = "Iteration =";
    const char*     patstrnog = "Number of groups =";
    const char*     patstrnds = "Number of dishes =";
    const char*     patstrtns = "Total number of samples =";
    const char*     patstrssz = "Context length =";
    const char*     patstrgrp = "Group";
    const char*     patstrnot = "Number of tables =";
    const char*     patstrtbl = "Table";
    const char*     patstrtid = "Id =";
    const char*     patstrdid = "Dish id =";
    const char*     patstrnos = "Number of samples =";
    const char*     patstrendofsmp = "*";
    const char*     patstrendoftbl = ";;";
    const char*     patstrendofgrp = "////";

    const int       lenstrendofsmp = strlen( patstrendofsmp );
    const int       lenstrendoftbl = strlen( patstrendoftbl );
    const int       lenstrendofgrp = strlen( patstrendofgrp );
    double          lprob;
    int             iter;
    int             ctxtsize = 0, dim = 0;
    int             totsamples;
    int             notables;
    int             nosamples;
    int             nogroups;
    int             nodishes, nocdshs = 0;
    int             dishid;
    int             npos, kpos, dpos;
    mystring        errmsg;
    Pslvector*      frq = NULL;
    Dish*           dish = NULL;
    Table*          tble = NULL;
    Restaurant*     rest = NULL;


    if( dids == NULL )
        throw myruntime_error("Memory access error (dish indices).");


    if( !filename ||
      ( fp = fopen( filename, "r" )) == NULL )
        throw myruntime_error("Failed to open input file.");


    try{
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error( "Wrong file format." );

        //read iteration number
        if(( p = strstr( locbuffer, patstrgss )) != NULL ) {
            p += strlen( patstrgss );

            if( length <= size_t( p - locbuffer ))
                throw myruntime_error( "Wrong file format." );

            if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &iter, &rbts )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( iter < 0 )
                throw myruntime_error( "Iteration number is invalid." );

            SetIterationRead( iter );

            //read next line
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw myruntime_error( "Wrong file format." );
        }

        //read log probability
        if(( p = strstr( locbuffer, patstrlpd )) != NULL ) {
            p += strlen( patstrlpd );

            if( length <= size_t( p - locbuffer ))
                throw myruntime_error( "Wrong file format." );

            if(( emsg = read_double( p, length - size_t( p - locbuffer ), &lprob, &rbts )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            SetMaxLProbData( lprob );

            //read next line
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw myruntime_error( "Wrong file format." );
        }

        //read number of groups
        if(( p = strstr( locbuffer, patstrnog )) == NULL )
            throw myruntime_error( "Wrong file format." );

        p += strlen( patstrnog );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong file format." );

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &nogroups, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( nogroups < 1 )
            throw myruntime_error( "Number of groups is invalid." );

        //read number of dishes
        if(( p = strstr( locbuffer, patstrnds )) == NULL )
            throw myruntime_error( "Wrong file format." );

        p += strlen( patstrnds );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong file format." );

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &nodishes, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( nodishes < 1 )
            throw myruntime_error( "Number of dishes is invalid." );

        //read total number of samples
        if(( p = strstr( locbuffer, patstrtns )) == NULL )
            throw myruntime_error( "Wrong file format." );

        p += strlen( patstrtns );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong file format." );

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &totsamples, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( totsamples < nogroups || totsamples < nodishes )
            throw myruntime_error( "Total number of samples is invalid." );

        //allocate space for dish ids
        dids->Reserve( TIMES2( totsamples ));
        dids->SetAllToValue( -1 );

        //read size of contexts
        if(( p = strstr( locbuffer, patstrssz )) == NULL )
            throw myruntime_error( "Wrong file format." );

        p += strlen( patstrssz );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong file format." );

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &ctxtsize, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( ctxtsize < 1 )
            throw myruntime_error( "Sample size is invalid." );
        SetCtxtSize( ctxtsize );
//
        InitBasin( totsamples );//initialize the main storage of samples
        InitMenu( totsamples );//initially, allocate a separate dish for each of the samples
        InitChain( nogroups );//allocate chain of restaurants
        ResetTotalNoSamples();
//
        while( !feof( fp )) {
            //read next group
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ))
                break;

            if( !length )
                throw myruntime_error( "Wrong file format." );

            if(( p = strstr( locbuffer, patstrgrp )) == NULL )
                throw myruntime_error( "Wrong file format." );

            //read number of tables
            if(( p = strstr( locbuffer, patstrnot )) == NULL )
                throw myruntime_error( "Wrong file format." );

            p += strlen( patstrnot );

            if( length <= size_t( p - locbuffer ))
                throw myruntime_error( "Wrong file format." );

            if(( emsg = 
                read_integer( p, length - size_t( p - locbuffer ), &notables, &rbts )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( notables < 1 )
                throw myruntime_error( "Number of tables is invalid." );
//
            rest = new Restaurant( nosamples );
            if( rest == NULL )
                throw myruntime_error( "Read: Not enough memory." );
            rest->SetDestroy( true );
            GetChain()->NewRestaurant( rest );
//
            for( t = 0; t < notables && !feof( fp ); t++ )
            {
                //read next table
                if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));

                if( !length )
                    throw myruntime_error( "Wrong file format." );

                if(( p = strstr( locbuffer, patstrtbl )) == NULL )
                    throw myruntime_error( "Wrong file format." );

                //read dish id
                if(( p = strstr( locbuffer, patstrdid )) == NULL )
                    throw myruntime_error( "Wrong file format." );

                p += strlen( patstrdid );

                if( length <= size_t( p - locbuffer ))
                    throw myruntime_error( "Wrong file format." );

                if(( emsg = 
                    read_integer( p, length - size_t( p - locbuffer ), &dishid, &rbts )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));

                if( dishid < 0 )
                    throw myruntime_error( "Dish id is invalid." );
                if( TIMES2( totsamples ) <= dishid )
                    throw myruntime_error( "Unexpectedly large value of Dish id." );

                dish = NULL;
                if( !GetCluster4Each()) {
                    //create new dish or use existing
                    if( dids->GetValueAt( dishid ) < 0 ) {
                        dish = new Dish( GetDefDishSize());
                        if( dish == NULL )
                            throw myruntime_error( "Read: Not enough memory." );
                        dish->SetBasin( GetBasin());
                        kpos = GetMenu()->NewDish( dish );
                        dids->SetValueAt( dishid, kpos );
                        nocdshs++;
                    }
                    else {
                        dish = GetMenu()->GetDishAt( kpos = dids->GetValueAt( dishid ));
                    }
                    if( dish == NULL )
                        throw myruntime_error( "Read: Null dish." );
                }
                //read number of samples
                if(( p = strstr( locbuffer, patstrnos )) == NULL )
                    throw myruntime_error( "Wrong file format." );

                p += strlen( patstrnos );

                if( length <= size_t( p - locbuffer ))
                    throw myruntime_error( "Wrong file format." );

                if(( emsg = 
                    read_integer( p, length - size_t( p - locbuffer ), &nosamples, &rbts )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));

                if( nosamples < 1 )
                    throw myruntime_error( "Number of samples is invalid." );

                tble = NULL;
                if( !GetCluster4Each()) {
                    //create new table
                    tble = new Table( GetDefTableSize());
                    if( tble == NULL )
                        throw myruntime_error( "Read: Not enough memory." );
                    tble->SetBasin( GetBasin());
                    tble->SetMenu( GetMenu());
                    tble->SetDishIndex( kpos );//kpos is from dish assignment above
                    rest->NewTable( tble );
                }
//
                for( n = 0; n < nosamples && !feof( fp ); n++ )
                {
                    if( !ReadSample( fp, &frq, dim ))
                        throw myruntime_error( "Failed to read data." );

                    if( frq == NULL )
                        throw myruntime_error( "Failed to read data." );

                    if( !dim ) {
                        dim = frq->GetSize();
                        GetMenu()->SetDim( dim );
                    }
                    else if( dim != frq->GetSize())
                        throw myruntime_error( "Wrong format: Dimensionality of samples." );
//
                    npos = GetBasin()->NewValue( frq );
//
                    if( GetCluster4Each() || dish == NULL ) {
                        //dish for each sample
                        dish = new Dish( GetDefDishSize());
                        if( dish == NULL )
                            throw myruntime_error( "Read: Not enough memory." );
                        dish->SetBasin( GetBasin());
                        kpos = GetMenu()->NewDish( dish );
                        nocdshs++;
                    }
                    dpos = dish->NewVectorNInd( npos );
//
                    if( GetCluster4Each() || tble == NULL ) {
                        //table for each sample
                        tble = new Table( GetDefTableSize());
                        if( tble == NULL )
                            throw myruntime_error( "Read: Not enough memory." );
                        tble->SetBasin( GetBasin());
                        tble->SetMenu( GetMenu());
                        tble->SetDishIndex( kpos );
                        rest->NewTable( tble );
                    }
                    tble->NewVectorNInd( npos, dpos );
//
                    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                        throw myruntime_error( TranslateReadError( emsg ));

                    if( length < lenstrendofsmp )
                        throw myruntime_error( "Wrong format: End of sample." );

                    for( c = 0; c < lenstrendofsmp; c++ )
                        if( locbuffer[c] != patstrendofsmp[c] )
                            throw myruntime_error( "Wrong format: End of sample." );

                    IncTotalNoSamples();
                }//end of one table
//
                if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));

                if( length < lenstrendoftbl )
                    throw myruntime_error( "Wrong format: End of table." );

                for( c = 0; c < lenstrendoftbl; c++ )
                    if( locbuffer[c] != patstrendoftbl[c] )
                        throw myruntime_error( "Wrong format: End of table." );
            }//end of one group

            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( length < lenstrendofgrp )
                throw myruntime_error( "Wrong format: End of group." );

            for( c = 0; c < lenstrendofgrp; c++ )
                if( locbuffer[c] != patstrendofgrp[c] )
                    throw myruntime_error( "Wrong format: End of group." );
        }

    } catch( myexception const& ex ) {
        errmsg = ex.what();
    }

    fclose( fp );

    if( !errmsg.empty())
        throw myruntime_error( errmsg );

    if(( GetCluster4Each() && nocdshs != totsamples )||
      ( !GetCluster4Each() && nocdshs != nodishes ))
        throw myruntime_error( "Inconsistent number of dishes." );
}

// -------------------------------------------------------------------------
// ReadSample: read observed frequencies from file
//
bool HDPbase::ReadSample( FILE* fp, Pslvector** pfv, int dim )
{
    if( fp == NULL || pfv == NULL )
        return false;

    mystring        errstr;
    int             emsg, n;
    size_t          len, rbts, read = 0;
    const size_t    locsize = 10 * KBYTE;
    char            locbuffer[locsize] = {0};
    int             ctxtlen = GetCtxtSize();
    double          val;
    Pslvector*      f;

    *pfv = NULL;

    if( dim <= 0 )
        dim = GetDefSampleDim() * ctxtlen;
    f = new Pslvector();
    if( f == NULL )
        throw myruntime_error( "ReadSample: Not enough memory." );
    f->Allocate( dim );

    try {
        for( n = 0; n < ctxtlen; read = 0, n++ ) {
            if(( emsg = skip_comments( fp, locbuffer, locsize, &len )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ))
                throw myruntime_error("ReadSample: Unexpected end of file.");

            while( emsg == 0 ) {
                //read frequency
                emsg = read_double( locbuffer + read, len - read, &val, &rbts );
                if( emsg == ERR_RD_NOVL )
                    break;
                if( emsg != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));
                read += rbts;

                f->Push( val );
            }
        }
        if( f->GetSize() <= 0 )
            throw myruntime_error( "ReadSample: No sample values read." );
        *pfv = f;

    } catch( myexception const& ex ) {
        errstr = ex.what();
        if( f ) {
            delete f;
            f = NULL;
        }
    }

    if( !errstr.empty())
        throw myruntime_error( errstr );

    return true;
}





// =========================================================================
// ReadParameters: read HDP parameters
//
void HDPbase::ReadParameters( const char* filename, const Ivector* dids )
{
    FILE*       fp = NULL;
    mystring    errstr;

    if( !filename )
        return;

    if(( fp = fopen( filename, "r" )) == NULL )
        throw myruntime_error( mystring("Failed to open file ")+ filename );

    try {
        ReadParameters( fp, dids );
    } catch( myexception const& ex ) {
        errstr = ex.what();
    }

    fclose( fp );
    if( !errstr.empty())
        throw myruntime_error( errstr );
}

// -------------------------------------------------------------------------
// ReadParameters: read HDP parameters given file descriptor
//
void HDPbase::ReadParameters( FILE* fp, const Ivector* dids )
{
    if( fp == NULL )
        return;
    if( GetMenu() == NULL )
        InitMenu( gcsz_defmenusz );

    ReadPriorParams( fp );
    ReadDishParams( fp, dids );
    if( GetCtxtSize() && GetCtxtSize() != GetMenu()->GetCtx())
        throw myruntime_error("ReadParameters: Inconsistent context length.");
    SetCtxtSize( GetMenu()->GetCtx());
    //probability factors
    CalcPriorProbFact();
}





// =========================================================================
// ReadPriorParams: read prior parameters from file
//
void HDPbase::ReadPriorParams( FILE* fp )
{
    if( fp == NULL )
        return;

    size_t          length, rbts, read;
    const size_t    locsize = 10 * KBYTE;
    char            locbuffer[locsize+1] = {0};
    const mystring  preambl = "ReadPriorParams: ";
    mystring        errmsg;
    char*           p;
    int             emsg;

    const char*     patstrtau = "tau =";
    const char*     patstrgam = "gamma =";
    const char*     patstrdim = "dim =";
    const char*     patstrctx = "ctx =";
    const char*     patstrkp0_a = "a_k0 =";
    const char*     patstrkp0_b = "b_k0 =";
    const char*     patstrkp0 = "kappa0 =";
    const char*     patstrnu0 = "nu0 =";
    const char*     patstrmu0 = "mu0 =";
    const char*     patstrs0m = "S0 =";
    const char*     patstrs0i = "S0^{-1} =";
    const char*     patstrds0 = "ln|S0| =";
    const char*     patstrendofdsh = "////";

    Pslvector*      f = NULL;
    SPDmatrix*      s = NULL;
    double  dblval;
    int     intval, dim, ctx;
    int     i, j;

    try {
        if( GetMenu() == NULL )
            throw myruntime_error("ReadPriorParams: Null menu.");

        //read tau
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        if(( p = strstr( locbuffer, patstrtau )) == NULL )
            throw myruntime_error("Wrong file format.");

        p += strlen( patstrtau );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error("Wrong file format.");

        if(( emsg = read_double( p, length - size_t( p - locbuffer ), &dblval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( dblval <= 0.0 )
            throw myruntime_error("Invalid `tau' value.");

        SetDPMTau( dblval );

        //read gamma
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        if(( p = strstr( locbuffer, patstrgam )) == NULL )
            throw myruntime_error("Wrong file format.");

        p += strlen( patstrgam );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error("Wrong file format.");

        if(( emsg = read_double( p, length - size_t( p - locbuffer ), &dblval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( dblval <= 0.0 )
            throw myruntime_error("Invalid `gamma' value.");

        SetDPMGamma( dblval );

        //read dim
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        if(( p = strstr( locbuffer, patstrdim )) == NULL )
            throw myruntime_error("Wrong file format.");

        p += strlen( patstrdim );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error("Wrong file format.");

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &dim, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( dim < 1 || 10000 < dim )
            throw myruntime_error("Value of `dim' is invalid.");

        if( GetMenu()->GetDim() && dim != GetMenu()->GetDim())
            throw myruntime_error("Inconsistent dimensions.");
        GetMenu()->SetDim( dim );

        //read ctx
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        if(( p = strstr( locbuffer, patstrctx )) == NULL )
            throw myruntime_error("Wrong file format.");

        p += strlen( patstrctx );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error("Wrong file format.");

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &ctx, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( ctx < 1 || 100 < ctx )
            throw myruntime_error("Value of `ctx' is invalid.");

        GetMenu()->SetCtx( ctx );

        //read prior parameter a for kappa_0
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        if(( p = strstr( locbuffer, patstrkp0_a )) == NULL )
            throw myruntime_error("Wrong file format.");

        p += strlen( patstrkp0_a );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error("Wrong file format.");

        if(( emsg = read_double( p, length - size_t( p - locbuffer ), &dblval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( dblval <= 0.0 )
            throw myruntime_error("Value of prior parameter `a' for kappa0 is invalid.");

        GetMenu()->SetKappa0_pp_a( dblval );

        //read prior parameter b for kappa_0
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        if(( p = strstr( locbuffer, patstrkp0_b )) == NULL )
            throw myruntime_error("Wrong file format.");

        p += strlen( patstrkp0_b );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error("Wrong file format.");

        if(( emsg = read_double( p, length - size_t( p - locbuffer ), &dblval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( dblval <= 0.0 )
            throw myruntime_error("Value of prior parameter `b' for kappa0 is invalid.");

        GetMenu()->SetKappa0_pp_b( dblval );

        //read kappa_0
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        if(( p = strstr( locbuffer, patstrkp0 )) == NULL )
            throw myruntime_error("Wrong file format.");

        p += strlen( patstrkp0 );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error("Wrong file format.");

        if(( emsg = read_double( p, length - size_t( p - locbuffer ), &dblval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( dblval <= 0.0 )
            throw myruntime_error("Value of `kappa0' is invalid.");

        GetMenu()->SetKappa0( dblval );

        //read nu_0
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        if(( p = strstr( locbuffer, patstrnu0 )) == NULL )
            throw myruntime_error("Wrong file format.");

        p += strlen( patstrnu0 );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error("Wrong file format.");

        if(( emsg = read_double( p, length - size_t( p - locbuffer ), &dblval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( dblval <= 0.0 )
            throw myruntime_error("Value of `nu0' is invalid.");

        GetMenu()->SetNu0( dblval );

        //read mu_0
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        if(( p = strstr( locbuffer, patstrmu0 )) == NULL )
            throw myruntime_error("Wrong file format.");

        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        f = new Pslvector();
        if( f == NULL )
            throw myruntime_error("Not enough memory.");
        f->Allocate( dim );

        for( i = 0, read = 0; i < dim; i++ ) {
            if(( emsg = read_double( locbuffer + read, length - read, &dblval, &rbts )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));
            read += rbts;
            f->Push( dblval );
        }

        GetMenu()->SetMu0( f );
        f = NULL;

        //read S_0
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        if(( p = strstr( locbuffer, patstrs0m )) == NULL )
            throw myruntime_error("Wrong file format.");

        s = new SPDmatrix( dim );
        if( s == NULL )
            throw myruntime_error("Not enough memory.");

        for( i = 0; i < dim; i++ ) {
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw myruntime_error("Wrong file format.");

            for( j = 0, read = 0; j < dim; j++ ) {
                if(( emsg = read_double( locbuffer + read, length - read, &dblval, &rbts )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));
                read += rbts;
                s->SetValueAt( i, j, dblval );
            }
        }

        GetMenu()->SetS0( s );
        s = NULL;

        //read S_0^{-1}
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        if(( p = strstr( locbuffer, patstrs0i )) == NULL )
            throw myruntime_error("Wrong file format.");

        s = new SPDmatrix( dim );
        if( s == NULL )
            throw myruntime_error("Not enough memory.");

        for( i = 0; i < dim; i++ ) {
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw myruntime_error("Wrong file format.");

            for( j = 0, read = 0; j < dim; j++ ) {
                if(( emsg = read_double( locbuffer + read, length - read, &dblval, &rbts )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));
                read += rbts;
                s->SetValueAt( i, j, dblval );
            }
        }

        GetMenu()->SetInvS0( s );
        s = NULL;

        //read ln|S_0|
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        if(( p = strstr( locbuffer, patstrds0 )) == NULL )
            throw myruntime_error("Wrong file format.");

        p += strlen( patstrds0 );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error("Wrong file format.");

        if(( emsg = read_double( p, length - size_t( p - locbuffer ), &dblval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        GetMenu()->SetLDetS0( dblval );

        //read footer
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        if(( p = strstr( locbuffer, patstrendofdsh )) == NULL )
            throw myruntime_error("Wrong file format.");

    } catch( myexception const& ex ) {
        errmsg = preambl + ex.what();
    }

    if( f ) { delete f; f = NULL; }
    if( s ) { delete s; s = NULL; }

    if( !errmsg.empty())
        throw myruntime_error( errmsg );
}

// -------------------------------------------------------------------------
// PrintPriorParams: print prior parameters to file
//
void HDPbase::PrintPriorParams( FILE* fp )
{
    if( fp == NULL )
        return;
    if( GetMenu() == NULL )
        throw myruntime_error("PrintPriorParams: Null Menu.");
        
    const char*     patstrtau = "tau =";
    const char*     patstrgam = "gamma =";
    const char*     patstrdim = "dim =";
    const char*     patstrctx = "ctx =";
    const char*     patstrkp0_a = "a_k0 =";
    const char*     patstrkp0_b = "b_k0 =";
    const char*     patstrkp0 = "kappa0 =";
    const char*     patstrnu0 = "nu0 =";
    const char*     patstrmu0 = "mu0 =";
    const char*     patstrs0m = "S0 =";
    const char*     patstrs0i = "S0^{-1} =";
    const char*     patstrds0 = "ln|S0| =";
    const char*     patstrendofdsh = "////";

    time_t          t_tt;
    struct tm       t_ttm;
    char            t_ttma[KBYTE];

    time( &t_tt );
    localtime_r( &t_tt, &t_ttm );
    asctime_r( &t_ttm, t_ttma );

    fprintf( fp, "## %s", t_ttma );

    fprintf( fp, "## -- Concentration parameters --\n##\n" );
    fprintf( fp, "%s %.21g\n", patstrtau, GetDPMTau());
    fprintf( fp, "%s %.21g\n\n", patstrgam, GetDPMGamma());

    fprintf( fp, "## -- Prior parameters --\n##\n" );
    fprintf( fp, "%s %d\n%s %d\n%s %.5g\n%s %.5g\n%s %.21g\n%s %.21g\n", 
             patstrdim, GetMenu()->GetDim(), patstrctx, GetMenu()->GetCtx(),
             patstrkp0_a, GetMenu()->GetKappa0_pp_a(), patstrkp0_b, GetMenu()->GetKappa0_pp_b(),
             patstrkp0, GetMenu()->GetKappa0(), patstrnu0, GetMenu()->GetNu0());
    fprintf( fp, "mu0 =\n" );
    if( GetMenu()->GetMu0())
        GetMenu()->GetMu0()->Print( fp, " %.21g");
    fprintf( fp, "S0 =\n" );
    if( GetMenu()->GetS0())
        GetMenu()->GetS0()->Print( fp, " %.21g");
    fprintf( fp, "S0^{-1} =\n" );
    if( GetMenu()->GetInvS0())
        GetMenu()->GetInvS0()->Print( fp, " %.21g");
    fprintf( fp, "ln|S0| = %.21g\n", GetMenu()->GetLDetS0());
    fprintf( fp, "%s\n\n", patstrendofdsh );
}





// =========================================================================
// ReadDishParams: read dish parameters from file
//
void HDPbase::ReadDishParams( FILE* fp, const Ivector* dids )
{
    if( fp == NULL )
        return;

    size_t          length, rbts, read;
    const size_t    locsize = 10 * KBYTE;
    char            locbuffer[locsize+1] = {0};
    const mystring  preambl = "ReadDishParams: ";
    mystring        errmsg;
    char*           p;
    int             emsg;

    const char*     patstrgps = "Number of groups =";
    const char*     patstrnod = "Number of dishes =";
    const char*     patstrdsh = "Dish";
    const char*     patstrdid = "Id=";
    const char*     patstrtm = "time=";
    const char*     patstrnk = "nk =";
    const char*     patstrmk = "mk =";
    const char*     patstrlp = "ln p =";
    const char*     patstrmu = "mu =";
    const char*     patstrsm = "S =";
    const char*     patstrsi = "S^{-1} =";
    const char*     patstrds = "ln|S| =";
    const char*     patstrendofdsh = "////";
    Dish*           dish;
    Pslvector*      f = NULL;
    SPDmatrix*      s = NULL;
    long long int   time;
    double  dblval;
    int     intval;
    int     nodshs, nk, mk, dim, ctx;
    int     i, j, d, dishid;

    try {
        if( GetMenu() == NULL )
            throw myruntime_error("Null Menu.");

        if( dids == NULL && ( GetMenu()->GetActualSize() || GetMenu()->GetSize()))
            throw myruntime_error("Menu configuration already exists.");

        dim = GetMenu()->GetDim();
        ctx = GetMenu()->GetCtx();

        if( dim < 1 || ctx < 1 )
            throw myruntime_error("Invalid dimensions and/or context length.");

        //read Number of groups
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        if(( p = strstr( locbuffer, patstrgps )) == NULL )
            throw myruntime_error("Wrong file format.");

        p += strlen( patstrgps );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error("Wrong file format.");

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( intval < 0.0 )
            throw myruntime_error("Invalid number of groups.");

        SetReadNoGroups( intval );

        //read Number of dishes
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( fp ) || !length )
            throw myruntime_error("Wrong file format.");

        if(( p = strstr( locbuffer, patstrnod )) == NULL )
            throw myruntime_error("Wrong file format.");

        p += strlen( patstrnod );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error("Wrong file format.");

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( intval < 0.0 )
            throw myruntime_error("Invalid number of dishes.");

        if( dids && GetMenu()->GetActualSize() != intval )
            throw myruntime_error("Inconsistent number of dishes.");

        nodshs = intval;

        for( d = 0; d < nodshs; d++ )
        {
            //read Dish index
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw myruntime_error("Wrong file format.");

            if(( p = strstr( locbuffer, patstrdsh )) == NULL )
                throw myruntime_error("Wrong file format.");

            p += strlen( patstrdsh );

            if( length <= size_t( p - locbuffer ))
                throw myruntime_error("Wrong file format.");

            if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( intval != d )
                throw myruntime_error("Invalid dish index.");

            //read dish id
            if(( p = strstr( locbuffer, patstrdid )) == NULL )
                throw myruntime_error("Wrong file format.");

            p += strlen( patstrdid );

            if( length <= size_t( p - locbuffer ))
                throw myruntime_error("Wrong file format.");

            if(( emsg = 
                read_integer( p, length - size_t( p - locbuffer ), &dishid, &rbts )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( dishid < 0 )
                throw myruntime_error("Invalid Dish id.");

            if( dids ) {
                if( dids->GetSize() <= dishid )
                    throw myruntime_error("Wrong dish id.");
                dishid = dids->GetValueAt( dishid );
                if(( dish = GetMenu()->GetDishAt( dishid )) == NULL )
                    throw myruntime_error("Wrong dish id.");
            }
            else {
                dish = new Dish( HDPbase::GetDefDishSize());
                if( dish == NULL )
                    throw myruntime_error("Not enough memory.");
                dish->SetBasin( NULL );//no basin
                dishid = GetMenu()->NewDish( dish );
                if( dishid != d )
                    throw myruntime_error("Invalid new-dish index obtained.");
            }

            //read dish time
            if(( p = strstr( locbuffer, patstrtm )) == NULL )
                throw myruntime_error("Wrong file format.");

            p += strlen( patstrtm );

            if( length <= size_t( p - locbuffer ))
                throw myruntime_error("Wrong file format.");

            if(( emsg = 
                read_llinteger( p, length - size_t( p - locbuffer ), &time, &rbts )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( time < 0 )
                throw myruntime_error("Invalid Dish time.");

            dish->SetTime(( time_t )time );

            //read nk
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw myruntime_error("Wrong file format.");

            if(( p = strstr( locbuffer, patstrnk )) == NULL )
                throw myruntime_error("Wrong file format.");

            p += strlen( patstrnk );

            if( length <= size_t( p - locbuffer ))
                throw myruntime_error("Wrong file format.");

            if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &nk, &rbts )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( nk < 1 )
                throw myruntime_error("Invalid value of `nk'.");

            if( dids && nk != dish->GetActualSize())
                throw myruntime_error("Wrong value of `nk'.");

            dish->SetReadSize( nk );

            //read mk
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw myruntime_error("Wrong file format.");

            if(( p = strstr( locbuffer, patstrmk )) == NULL )
                throw myruntime_error("Wrong file format.");

            p += strlen( patstrmk );

            if( length <= size_t( p - locbuffer ))
                throw myruntime_error("Wrong file format.");

            if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &mk, &rbts )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( mk < 1 || nk < mk )
                throw myruntime_error("Invalid value of `mk'.");

            dish->SetReadNoTables( mk );

            //read log probability
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw myruntime_error("Wrong file format.");

            if(( p = strstr( locbuffer, patstrlp )) == NULL )
                throw myruntime_error("Wrong file format.");

            p += strlen( patstrlp );

            if( length <= size_t( p - locbuffer ))
                throw myruntime_error("Wrong file format.");

            if(( emsg = read_double( p, length - size_t( p - locbuffer ), &dblval, &rbts )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            GetMenu()->SetProbAt( dishid, dblval );

            //read mu
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw myruntime_error("Wrong file format.");

            if(( p = strstr( locbuffer, patstrmu )) == NULL )
                throw myruntime_error("Wrong file format.");

            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw myruntime_error("Wrong file format.");

            f = new Pslvector();
            if( f == NULL )
                throw myruntime_error("Not enough memory.");
            f->Allocate( dim );

            for( i = 0, read = 0; i < dim; i++ ) {
                if(( emsg = read_double( locbuffer + read, length - read, &dblval, &rbts )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));
                read += rbts;
                f->Push( dblval );
            }

            GetMenu()->SetMuVectorAt( dishid, f );//set mu
            f = NULL;

            //read S_0
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw myruntime_error("Wrong file format.");

            if(( p = strstr( locbuffer, patstrsm )) == NULL )
                throw myruntime_error("Wrong file format.");

            s = new SPDmatrix( dim );
            if( s == NULL )
                throw myruntime_error("Not enough memory.");

            for( i = 0; i < dim; i++ ) {
                if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));

                if( feof( fp ) || !length )
                    throw myruntime_error("Wrong file format.");

                for( j = 0, read = 0; j < dim; j++ ) {
                    if(( emsg = read_double( locbuffer + read, length - read, &dblval, &rbts )) != 0 )
                        throw myruntime_error( TranslateReadError( emsg ));
                    read += rbts;
                    s->SetValueAt( i, j, dblval );
                }
            }

            GetMenu()->SetSMatrixAt( dishid, s );
            s = NULL;

#ifdef PROBMTXINVSM
            //read S^{-1}
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw myruntime_error("Wrong file format.");

            if(( p = strstr( locbuffer, patstrsi )) == NULL )
                throw myruntime_error("Wrong file format.");

            s = new SPDmatrix( dim );
            if( s == NULL )
                throw myruntime_error("Not enough memory.");

            for( i = 0; i < dim; i++ ) {
                if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                    throw myruntime_error( TranslateReadError( emsg ));

                if( feof( fp ) || !length )
                    throw myruntime_error("Wrong file format.");

                for( j = 0, read = 0; j < dim; j++ ) {
                    if(( emsg = read_double( locbuffer + read, length - read, &dblval, &rbts )) != 0 )
                        throw myruntime_error( TranslateReadError( emsg ));
                    read += rbts;
                    s->SetValueAt( i, j, dblval );
                }
            }

            GetMenu()->SetInvSMatrixAt( dishid, s );
            s = NULL;
#endif

            //read ln|S|
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw myruntime_error("Wrong file format.");

            if(( p = strstr( locbuffer, patstrds )) == NULL )
                throw myruntime_error("Wrong file format.");

            p += strlen( patstrds );

            if( length <= size_t( p - locbuffer ))
                throw myruntime_error("Wrong file format.");

            if(( emsg = read_double( p, length - size_t( p - locbuffer ), &dblval, &rbts )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            GetMenu()->SetLDetSMAt( dishid, dblval );

            //read footer
            if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ) || !length )
                throw myruntime_error("Wrong file format.");

            if(( p = strstr( locbuffer, patstrendofdsh )) == NULL )
                throw myruntime_error("Wrong file format.");
        }

    } catch( myexception const& ex ) {
        errmsg = preambl + ex.what();
    }

    if( f ) { delete f; f = NULL; }
    if( s ) { delete s; s = NULL; }

    if( !errmsg.empty())
        throw myruntime_error( errmsg );
}

// -------------------------------------------------------------------------
// PrintDishParams: print dish parameters to file
//
void HDPbase::PrintDishParams( FILE* fp )
{
    if( fp == NULL )
        return;
    if( GetMenu() == NULL )
        throw myruntime_error("PrintDishParams: Null Menu.");
    if( GetChain() == NULL )
        throw myruntime_error("PrintDishParams: Null Chain.");

    const char*     patstrlpd = "Log Prob of Data =";
    const char*     patstrgps = "Number of groups =";
    const char*     patstrnod = "Number of dishes =";
    const char*     patstrdsh = "Dish";
    const char*     patstrdid = "Id=";
    const char*     patstrtm = "time=";
    const char*     patstrnk = "nk =";
    const char*     patstrmk = "mk =";
    const char*     patstrlp = "ln p =";
    const char*     patstrmu = "mu =";
    const char*     patstrsm = "S =";
    const char*     patstrsi = "S^{-1} =";
    const char*     patstrds = "ln|S| =";
    const char*     patstrendofdsh = "////";
    const Dish*     dish;
    int c, d;

    fprintf( fp, "## -- Dish parameters --\n##\n" );
    fprintf( fp, "## %s %g\n", patstrlpd, GetLastLProbData());//TODO:get last probability
    fprintf( fp, "%s %d\n", patstrgps, GetChain()->GetSize());
    fprintf( fp, "%s %d\n", patstrnod, GetMenu()->GetActualSize());
    for( c = d = 0; d < GetMenu()->GetSize(); d++ ) {
        dish = GetMenu()->GetDishAt( d );
        if( dish == NULL )
            continue;
        fprintf( fp, "%s %d (%s%d %s%lu)\n", patstrdsh, c++, patstrdid, d, 
                 patstrtm, ( unsigned long int )dish->GetTime());
        fprintf( fp, "%s %d\n", patstrnk, dish->GetActualSize());
        fprintf( fp, "%s %d\n", patstrmk, dish->GetNoTables());
        fprintf( fp, "%s %g\n", patstrlp, GetMenu()->GetProbAt(d));
        fprintf( fp, "%s\n", patstrmu );
        if( GetMenu()->GetMuVectorAt( d ))
            GetMenu()->GetMuVectorAt( d )->Print( fp, " %.21g");
        fprintf( fp, "%s\n", patstrsm );
        if( GetMenu()->GetSMatrixAt( d ))
            GetMenu()->GetSMatrixAt( d )->Print( fp, " %.21g");
#ifdef PROBMTXINVSM
        fprintf( fp, "%s\n", patstrsi );
        if( GetMenu()->GetInvSMatrixAt( d ))
            GetMenu()->GetInvSMatrixAt( d )->Print( fp, " %.21g");
#endif
        fprintf( fp, "%s %.21g\n", patstrds, GetMenu()->GetLDetSMAt( d ));
        fprintf( fp, "%s\n\n", patstrendofdsh );
    }
}

// =========================================================================




// =========================================================================
// PrintGroups: print HDP structure
//
void HDPbase::PrintGroups( const char* filename, double* lprob, int* iter )
{
    FILE*       fp = NULL;
    mystring    errstr;

    if( !filename )
        return;

    if(( fp = fopen( filename, "w" )) == NULL )
        throw myruntime_error( "Failed to open output file." );

    try {
        PrintGroups( fp, lprob, iter );
    } catch( myexception const& ex ) {
        errstr = ex.what();
    }

    fclose( fp );
    if( !errstr.empty())
        throw myruntime_error( errstr );
}

// -------------------------------------------------------------------------
// Print: print HDP structure to file
//
void HDPbase::PrintGroups( FILE* fp, double* lprob, int* iter )
{
    if( fp == NULL )
        return;
    if( GetChain() == NULL || GetMenu() == NULL )
        return;

    const char*     patstrlpd = "Log Prob of Data =";
    const char*     patstrgss = "Iteration =";
    const char*     patstrnog = "Number of groups =";
    const char*     patstrnds = "Number of dishes =";
    const char*     patstrtns = "Total number of samples =";
    const char*     patstrssz = "Context length =";
    const char*     patstrgrp = "Group";
    const char*     patstrnot = "Number of tables =";
    const char*     patstrtbl = "Table";
    const char*     patstrtid = "Id =";
    const char*     patstrdid = "Dish id =";
    const char*     patstrnos = "Number of samples =";
    const char*     patstrendofsmp = "*";
    const char*     patstrendoftbl = ";;";
    const char*     patstrendofgrp = "////";
    const int       dim = GetMenu()->GetDim();
    const int       ctx = GetCtxtSize();
    const int       nvals = dim / ctx;
    const Restaurant*   rest;
    const Table*        tble;
    const Pslvector*    frq;
    int r, tt, t, n, v;

    time_t          t_tt;
    struct tm       t_ttm;
    char            t_ttma[KBYTE];

    time( &t_tt );
    localtime_r( &t_tt, &t_ttm );
    asctime_r( &t_ttm, t_ttma );

    fprintf( fp, "## %s", t_ttma );

    fprintf( fp, "## " ); 
    print_cmdline( &file_print, fp );
    fprintf( fp, "##\n" );

    if( GetResType() == TRT_MHUpdate )
        fprintf( fp, "## MH update %d-%d\n", GetGibbsIt(), GetMHUpdateNum());
    else if( GetResType() == TRT_GibbsUpdate )
        fprintf( fp, "## Gibbs sampling update %d\n", GetGibbsIt());
    else
        fprintf( fp, "## N/A update\n" );

    if( iter )
        fprintf( fp, "%s %d\n", patstrgss, *iter );
    if( lprob )
        fprintf( fp, "%s %g\n", patstrlpd, *lprob );
    fprintf( fp, "##\n" ); 
    fprintf( fp, "%s %d; %s %d; %s %d; %s %d\n", 
              patstrnog, GetChain()->GetSize(), 
              patstrnds, GetMenu()->GetActualSize(),
              patstrtns, GetTotalNoSamples(), patstrssz, GetCtxtSize());
    for( r = 0; r < GetChain()->GetSize(); r++ ) {
        rest = GetChain()->GetRestaurantAt( r );
        if( rest == NULL )
            continue;
        fprintf( fp, "%s %d; %s %d\n", patstrgrp, r, patstrnot, rest->GetActualSize());
        for( tt = t = 0; t < rest->GetSize(); t++ )
        {
            tble = rest->GetTableAt( t );
            if( tble == NULL )
                continue;
            fprintf( fp, "%s %d (%s %d; %s %d); %s %d\n", 
                     patstrtbl, tt++, patstrtid, t, patstrdid, tble->GetDishIndex(),
                     patstrnos, tble->GetActualSize());
            for( n = 0; n < tble->GetSize(); n++ )
            {
                if( tble->GetVectorNIndAt( n ) < 0 )
                    continue;
                frq = tble->GetVectorNAt( n );
                if( frq == NULL )
                    throw myruntime_error("PrintGroups: Null vector.");
                for( v = 0; v < frq->GetSize(); v++ ) {
                    fprintf( fp, " %12.6g", frq->GetValueAt( v ));
                    if(( v+1 ) % nvals == 0 && ( v+1 ) < frq->GetSize())
                        fprintf( fp, "\n");
                }
                fprintf( fp, "\n%s\n", patstrendofsmp );
            }
            fprintf( fp, "%s\n", patstrendoftbl );//end of table
        }
        fprintf( fp, "%s\n\n", patstrendofgrp );
    }
}

// =========================================================================
// PrintDishes: print HDP dish structure
//
void HDPbase::PrintDishes( const char* filename, double* lprob, int* iter )
{
    FILE*       fp = NULL;
    mystring    errstr;

    if( !filename )
        return;

    if(( fp = fopen( filename, "w" )) == NULL )
        throw myruntime_error( mystring("Failed to open file ")+ filename );

    try {
        PrintDishes( fp, lprob, iter );
    } catch( myexception const& ex ) {
        errstr = ex.what();
    }

    fclose( fp );
    if( !errstr.empty())
        throw myruntime_error( errstr );
}

// -------------------------------------------------------------------------
// Print: print HDP structure to file
//
void HDPbase::PrintDishes( FILE* fp, double* lprob, int* iter )
{
    if( fp == NULL )
        return;
    if( GetMenu() == NULL )
        return;

    const char*     patstrlpd = "Log Prob of Data =";
    const char*     patstrgss = "Iteration =";
    const char*     patstrnog = "Number of clusters =";
    const char*     patstrtns = "Total number of samples =";
    const char*     patstrssz = "Context length =";
    const char*     patstrgrp = "Cluster";
    const char*     patstrcid = "Id=";
    const char*     patstrtm  = "time=";
    const char*     patstrnot = "Number of samples =";
    const char*     patstrlp  = "Log Prob =";
    const char*     patstrendofsmp = "*";
    const char*     patstrendofgrp = "////";
    const int       dim = GetMenu()->GetDim();
    const int       ctx = GetCtxtSize();
    const int       nvals = dim / ctx;
    const Dish*         dish;
    const Pslvector*    frq;
    int c, d, v, a;

    time_t          t_tt;
    struct tm       t_ttm;
    char            t_ttma[KBYTE];

    time( &t_tt );
    localtime_r( &t_tt, &t_ttm );
    asctime_r( &t_ttm, t_ttma );

    fprintf( fp, "## %s", t_ttma );

    fprintf( fp, "## " ); 
    print_cmdline( &file_print, fp );
    fprintf( fp, "##\n" ); 

    if( GetResType() == TRT_MHUpdate )
        fprintf( fp, "## MH update %d-%d\n", GetGibbsIt(), GetMHUpdateNum());
    else if( GetResType() == TRT_GibbsUpdate )
        fprintf( fp, "## Gibbs sampling update %d\n", GetGibbsIt());
    else
        fprintf( fp, "## N/A update\n" );

    if( iter )
        fprintf( fp, "%s %d\n", patstrgss, *iter );
    if( lprob )
        fprintf( fp, "%s %g\n", patstrlpd, *lprob );
    fprintf( fp, "##\n" ); 
    fprintf( fp, "%s %d; %s %d; %s %d\n", 
              patstrnog, GetMenu()->GetActualSize(), 
              patstrtns, GetTotalNoSamples(), patstrssz, GetCtxtSize());
    for( c = d = 0; d < GetMenu()->GetSize(); d++ ) {
        dish = GetMenu()->GetDishAt( d );
        if( dish == NULL )
            continue;
        fprintf( fp, "%s %d (%s%d %s%lu); %s %d; %s %g\n", 
                 patstrgrp, c++, patstrcid, d, 
                 patstrtm, ( unsigned long int )dish->GetTime(), 
                 patstrnot, dish->GetActualSize(), 
                 patstrlp, GetMenu()->GetProbAt(d));
        for( v = 0; v < dish->GetSize(); v++ )
        {
            if( dish->GetVectorNIndAt( v ) < 0 )
                continue;
            frq = dish->GetVectorNAt( v );
            if( frq == NULL )
                throw myruntime_error("PrintDishes: Null vector.");
            for( a = 0; a < frq->GetSize(); a++ ) {
                fprintf( fp, " %12.6g", frq->GetValueAt( a ));
                if(( a+1 ) % nvals == 0 && ( a+1 ) < frq->GetSize())
                    fprintf( fp, "\n");
            }
            fprintf( fp, "\n%s\n", patstrendofsmp );
        }
        fprintf( fp, "%s\n\n", patstrendofgrp );
    }
}

// -------------------------------------------------------------------------
// PrintParameters: print parameters of the HDP structure
//
void HDPbase::PrintParameters( const char* filename )
{
    FILE*       fp = NULL;
    mystring    errstr;

    if( !filename )
        return;

    if(( fp = fopen( filename, "w" )) == NULL )
        throw myruntime_error( mystring("Failed to open file ")+ filename );

    try {
        PrintParameters( fp );
    } catch( myexception const& ex ) {
        errstr = ex.what();
    }

    fclose( fp );
    if( !errstr.empty())
        throw myruntime_error( errstr );
}

// -------------------------------------------------------------------------
// PrintParameters: print parameters of the HDP structure
//
void HDPbase::PrintParameters( FILE* fp )
{
    if( fp == NULL )
        return;
    PrintPriorParams( fp );
    PrintDishParams( fp );
}





// =========================================================================
// PrintSummary: print summary information at the end of file
//
void HDPbase::PrintSummary( const char* filename, double* lprob, int* iter )
{
    FILE*       fp = NULL;
    mystring    errstr;

    if( !filename )
        return;

    if(( fp = fopen( filename, "a" )) == NULL )
        throw myruntime_error("Failed to open output report file.");

    try {
        PrintSummary( fp, lprob, iter );
    } catch( myexception const& ex ) {
        errstr = ex.what();
    }

    fclose( fp );
    if( !errstr.empty())
        throw myruntime_error( errstr );
}

// -------------------------------------------------------------------------
//
void HDPbase::PrintSummary( FILE* fp, double* lprob, int* iter )
{
    if( fp == NULL )
        return;
    if( GetMenu() == NULL )
        return;

    const char*     patstrlpd = "Log Prob of Data =";
    const char*     patstrgss = "Iteration =";
    const char*     patstrtau = "tau =";
    const char*     patstrgam = "gamma =";
    const char*     patstrkp0_a = "a_k0 =";
    const char*     patstrkp0_b = "b_k0 =";
    const char*     patstrkp0 = "kappa0 =";
    const char*     patstrnu0 = "nu0 =";
    const char*     patstrdim = "dim =";
    const char*     patstrctx = "ctx =";
    const char*     patstrnog = "Number of clusters =";
    const char*     patstrtns = "Total number of samples =";
    const char*     patstrgrp = "Cluster";
    const char*     patstrcid = "Id =";
    const char*     patstrnot = "nk =";
    const char*     patstrmk  = "mk =";
    const char*     patstrlp  = "ln p =";
    const char*     patstrtm  = "time =";
    const char*     patstrendofgrp = "////";
    const Dish*     dish;

    time_t          t_tt;
    struct tm       t_ttm;
    char            t_ttma[KBYTE];
    int c, d;

    time( &t_tt );
    localtime_r( &t_tt, &t_ttm );
    asctime_r( &t_ttm, t_ttma );

    fprintf( fp, "## %s", t_ttma );

    if( iter )
        fprintf( fp, "%s %d %s\n", patstrgss, *iter, GetRestarted()? "restart": "");

    if( GetResType() == TRT_MHUpdate )
        fprintf( fp, "MH update %d-%d\n", GetGibbsIt(), GetMHUpdateNum());
    else if( GetResType() == TRT_GibbsUpdate )
        fprintf( fp, "Gibbs sampling update %d\n", GetGibbsIt());
    else
        fprintf( fp, "N/A update\n" );

    if( lprob )
        fprintf( fp, "%s %g\n", patstrlpd, *lprob );
    fprintf( fp, "%s %.5g\n", patstrtau, GetDPMTau());
    fprintf( fp, "%s %.5g\n", patstrgam, GetDPMGamma());
    fprintf( fp, "%s %.5g\n%s %.5g\n%s %.5g\n%s %.5g\n%s %d\n%s %d\n", 
             patstrkp0_a, GetMenu()->GetKappa0_pp_a(), patstrkp0_b, GetMenu()->GetKappa0_pp_b(), 
             patstrkp0, GetMenu()->GetKappa0(), patstrnu0, GetMenu()->GetNu0(), 
             patstrdim, GetMenu()->GetDim(), patstrctx, GetMenu()->GetCtx());
    fprintf( fp, "%s %d; %s %d\n", 
              patstrnog, GetMenu()->GetActualSize(), 
              patstrtns, GetTotalNoSamples());
    for( c = d = 0; d < GetMenu()->GetSize(); d++ ) {
        dish = GetMenu()->GetDishAt( d );
        if( dish == NULL )
            continue;
        fprintf( fp, "%s %d (%s %d); %s %d; %s %d; %s %g; %s %lu\n", 
                 patstrgrp, c++, patstrcid, d, 
                 patstrnot, dish->GetActualSize(), patstrmk, dish->GetNoTables(), 
                 patstrlp, GetMenu()->GetProbAt(d), patstrtm, ( unsigned long int )dish->GetTime());
    }
    fprintf( fp, "%s\n\n", patstrendofgrp );
}





// =========================================================================
// PrintBasin: print basin of samples
//
void HDPbase::PrintBasin( FILE* fp )
{
    if( fp == NULL )
        return;
    if( GetBasin() == NULL )
        return;

    const char*     patstrtns = "Total number of samples =";
    const char*     patstrssz = "Sample size =";
    const char*     patstrendofsmp = "*";
    const Pslvector*    frq;
    int n, v;

    fprintf( fp, "%s %d; %s %d\n",
              patstrtns, GetTotalNoSamples(), patstrssz, GetCtxtSize());
    for( n = 0; n < GetBasin()->GetSize(); n++ ) {
        frq = GetBasin()->GetValueAt( n );
        if( frq == NULL )
            continue;
        for( v = 0; v < frq->GetSize(); v++ ) {
            fprintf( fp, " %12.6g", frq->GetValueAt( v ));
        }
        fprintf( fp, "\n%s\n", patstrendofsmp );
    }
}
