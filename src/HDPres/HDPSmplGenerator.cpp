/***************************************************************************
 *   Copyright (C) 2012 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "ext/psl.h"
#include "ext/pslcodes.h"
#include "ext/gamma.h"
#include "ext/rv/rvdisc.h"
#include "ext/rv/rvnorm.h"
#include "ext/rv/rmvnorm.h"
#include "liblib/logitnormal.h"
#include "HDPSmplGenerator.h"

// -------------------------------------------------------------------------
// CLASS HDPSmplGenerator
//
// Constructor
//
HDPSmplGenerator::HDPSmplGenerator( const char* patoutput )
:   patoutname_( patoutput ),
    dim_( 0 ),
    novariates_( 0 ),
    ctxtlen_( 0 ),
    noclusters_( 0 ),
    nocsamples_( 0 ),
    lgstnrmsampling_( false ),
    overlap_( true ),
    smplcounter_( 0 ),
    samples_( NULL ),
    olbuff_( NULL )
{
}

// Default constructor
//
HDPSmplGenerator::HDPSmplGenerator()
:   patoutname_( NULL ),
    dim_( 0 ),
    novariates_( 0 ),
    ctxtlen_( 0 ),
    noclusters_( 0 ),
    nocsamples_( 0 ),
    lgstnrmsampling_( false ),
    overlap_( true ),
    smplcounter_( 0 ),
    samples_( NULL ),
    olbuff_( NULL )
{
    throw( myruntime_error("HDPSmplGenerator: Default construction is not allowed."));
}

// Destructor
//
HDPSmplGenerator::~HDPSmplGenerator()
{
    DestroyOLBuffer();
    DestroySamples();
}

// -------------------------------------------------------------------------
// DestroySamples: destroy all samples in vector
//
void HDPSmplGenerator::DestroySamples()
{
    SimpleVector*   sv;
    double*         smp;
    int n, m;
    if( samples_ ) {
        for( n = 0; n < samples_->GetSize(); n++ ) {
            sv = ( SimpleVector* )samples_->GetValueAt( n );
            if( sv == NULL )
                continue;
            for( m = 0; m < sv->GetSize(); m++ ) {
                smp = ( double* )sv->GetValueAt( m );
                if( smp == NULL )
                    continue;
                free( smp );
            }
            delete sv;
            sv = NULL;
        }
        delete samples_;
        samples_ = NULL;
    }
}

// -------------------------------------------------------------------------
// NewSamples: create new structure of samples
//
void HDPSmplGenerator::NewSamples( int size )
{
    DestroySamples();
    samples_ = new SimpleVector( size );
    if( !samples_ )
        throw( myruntime_error("HDPSmplGenerator: Not enough memory."));
}

// -------------------------------------------------------------------------
// Run: start generating samples
//
void HDPSmplGenerator::Run()
{
    char        strbuf[BUF_MAX];
    double      lgdim, lgdmv, lgvar;
    double      err;
    int         code, adim, vars, ncls;
    mystring    errstr;

    try {
        if( GetDimensions() < 1 )
            throw( myruntime_error("HDPSmplGenerator: Invalid dimensions."));

        if( GetNoVars() < 1 || GetDimensions() < GetNoVars())
            throw( myruntime_error("HDPSmplGenerator: Invalid number of variates."));

        if( GetContextLength() < 1 )
            throw( myruntime_error("HDPSmplGenerator: Invalid context length."));

        if( GetNoClusters() < 1 )
            throw( myruntime_error("HDPSmplGenerator: Invalid number of clusters."));

        if( GetNoSamplesC() < 1 )
            throw( myruntime_error("HDPSmplGenerator: Invalid number of samples."));

        adim = GetDimensions()+1;
        vars = GetNoVars();
        ncls = GetNoClusters();

        if(( code = psl_lngamma_e(( double )( adim+1 ), &lgdim, &err )) != PSL_SUCCESS )
            throw myruntime_error( TranslatePSLError( code ));

        if(( code = psl_lngamma_e(( double )( adim-vars+1 ), &lgdmv, &err )) != PSL_SUCCESS )
            throw myruntime_error( TranslatePSLError( code ));

        if(( code = psl_lngamma_e(( double )( vars+1 ), &lgvar, &err )) != PSL_SUCCESS )
            throw myruntime_error( TranslatePSLError( code ));

        lgdim = exp( lgdim - lgdmv - lgvar ) + 1.e-9;
        if( !isfinite( lgdim ) || ( int )lgdim < ncls )
            throw( myruntime_error("HDPSmplGenerator: Number of clusters exceeds "
                                   "maximum number of combinations."));

        NewSamples( GetNoClusters());
        NewOLBuffer( adim );

        message( "Generating samples..." );
        ResetSampleCounter();
        GenerateRec( vars, adim );
        sprintf( strbuf, "Generated %d cluster(s), %d sample(s).",
                GetSamples()->GetSize(), GetSampleCounter());
        message( strbuf );

        message( "Writing samples..." );
        Output();

    } catch( myexception const& ex )
    {
        errstr = ex.what();
    }

    if( !errstr.empty())
        throw myruntime_error( errstr );
}

// -------------------------------------------------------------------------
// GenerateRec: generate combinations of conserved positions
//  nvs -- number of variates
//  adim = dim + 1 -- augmented dimensions; +1 for extra dimension
//  varndxs -- indexes of variates
//
void HDPSmplGenerator::GenerateRec( int nvs, int adim, int* varndxs )
{
    mystring        preamb = "HDPSmplGenerator: Generate: ";
    mystring        errstr;
    int             novars = GetNoVars();
    int             nocls = GetNoClusters();
    int             d, v, v1;

    try {
        if( novars < 1 )
            throw myruntime_error( preamb + "Invalid number of variates.");
        if( GetSamples() == NULL )
            throw myruntime_error( preamb + "Null samples structure.");

        if( varndxs == NULL ) {
            varndxs = ( int* )malloc( novars * sizeof( int ));
            if( varndxs == NULL )
                throw myruntime_error( preamb + "Not enough memory.");
            for( v = 0; v < novars; v++ )
                varndxs[v] = v;
            GenSamplesForCluster( varndxs, novars );
        }

        if( 0 <= ( v = nvs-1 )) {
            for( d = v+1; d < adim; d++ ) {
                varndxs[v] = d;
                if( nocls <= GetSamples()->GetSize())
                    break;
                GenSamplesForCluster( varndxs, novars );
                GenerateRec( nvs - 1, d, varndxs );//NOTE:recursion
                for( v1 = 0; v1 < nvs; v1++ )
                    varndxs[v1] = v1;
            }
        }

    } catch( myexception const& ex ) {
        errstr = ex.what();
    }

    if( nvs == novars && varndxs ) { free( varndxs ); varndxs = NULL; }
    if( !errstr.empty())
        throw myruntime_error( errstr, CRITICAL );
}

// -------------------------------------------------------------------------
// GenSamplesForCluster: generate samples for one cluster determined by 
//  varndxs -- indexes of variates
//  vsize -- size of `varndxs'
//
void HDPSmplGenerator::GenSamplesForCluster( int* varndxs, int vsize )
{
    if( GetLGSTNSampling())
        GenSamplesForClusterLGSTNRM( varndxs, vsize );
    else
        GenSamplesForClusterMeanVar( varndxs, vsize );
}

// -------------------------------------------------------------------------
// GenSamplesForClusterMeanVar: generate samples for one cluster 
//  determined by variates; generate using mean and std.var.;
//  varndxs -- indexes of variates
//  vsize -- size of `varndxs'
//
void HDPSmplGenerator::GenSamplesForClusterMeanVar( int* varndxs, int vsize )
{
    mystring        preamb = "HDPSmplGenerator: GenSamplesForClusterMeanVar: ";
    mystring        errstr;
    int             dim = GetDimensions();
    int             adim = dim + 1;//augmented dimensions
    int             ctxtlen = GetContextLength();
    int             novars = GetNoVars();
    int             nocls = GetNoClusters();
    int             nosmps = GetNoSamplesC();
    int*            olbuf = GetOLBuffer();
    const double    csum = 0.9;
    const double    lovr = 1.0 - csum;
    double          cvalue, uncval, rv;
    double          fsum, diff;
    int             s, n, d, v, a;
    int             eco;
    double*         sample = NULL, *p;
    SimpleVector*   cluster = NULL;
    MTRng           rng;
    RVNorm          rvnorm( rng );
    static unsigned long tc = 0;

    if( dim < 1 || novars < 1 || adim < novars || vsize != novars )
        throw myruntime_error( preamb + "Invalid number of variates.");

    if( varndxs == NULL || olbuf == NULL )
        throw myruntime_error( preamb + "Memory access error.");
    if( GetSamples() == NULL )
        throw myruntime_error( preamb + "Null samples structure.");

    if( !GetDoOverlap()) {
        for( v = 0; v < vsize; v++ ) {
            if( adim <= varndxs[v])
                throw myruntime_error( preamb + "Invalid variate index.");
            if( olbuf[varndxs[v]])
                return;
        }
        for( v = 0; v < vsize; v++ )
            olbuf[varndxs[v]] = 1;
    }

    rng.Set(( unsigned long )time( NULL )+tc++);

    try {
        if( adim == novars ) {
            cvalue = 1.0 /( double )novars;
            uncval = 0.0;
        }
        else {
            //leave fraction lovr to distribute among other variates
            cvalue = csum /( double )novars;
            uncval = lovr /( double )( adim - novars );
            rvnorm.SetMean( uncval );
            rvnorm.SetStd( uncval /( double )adim );
        }

        cluster = new SimpleVector( nosmps );
        if( cluster == NULL )
            throw myruntime_error( preamb + "Not enough memory.");

        for( s = 0; s < nosmps; s++ )
        {
            sample = ( double* )malloc( ctxtlen * adim * sizeof( double ));
            if( !sample )
                throw myruntime_error( preamb + "Not enough memory.");
            for( n = 0, p = sample; n < ctxtlen; n++, p += adim ) {
                fsum = 0.0;
                for( d = 0; d < adim; d++ )
                {
                    if(( eco = rvnorm.Gen( &rv )) !=0 )
                        throw myruntime_error( preamb + TranslatePSLError( eco ));
                    if( rv < 0.0 )
                        throw myruntime_error( preamb + "Negative random number.");
                    p[d] = rv;
                    fsum += rv;
                }
                for( v = 0; v < vsize; v++ )
                {
                    if( adim <= varndxs[v])
                        throw myruntime_error( preamb + "Invalid variate index.");
                    fsum += cvalue - p[varndxs[v]];
                    p[varndxs[v]] = cvalue;
                }
                diff = 1.0 - fsum;
                if( diff ) {
                    a = -1;
                    for( v = 1; v < vsize; v++ ) {
                        if( 1 < varndxs[v] - varndxs[v-1]) {
                            a = varndxs[v-1]+1;
                            break;
                        }
                    }
                    if( a < 0 )
                        for( v = varndxs[vsize-1]+1; v < adim; v++ ) {
                            a = v;
                            break;
                        }
                    if( a < 0 )
                        for( v = 0; v < adim && v < varndxs[0]; v++ ) {
                            a = v;
                            break;
                        }
                    if( a < 0 )
                        throw myruntime_error( preamb + "Failed to find unconserved position.");
                    p[a] += diff;
                    if( p[a] < 0.0 )
                        throw myruntime_error( preamb + "Negative adjusted frequency.");
                }
                //TRANSFORM to multivariate normal
                LogitNormal2Normal( p, ( size_t )adim, 1.e-6, false );
            }
            //add generated sample
            cluster->Push( sample );
            sample = NULL;
            IncSampleCounter();
        }
        GetSamples()->Push( cluster );
        cluster = NULL;

    } catch( myexception const& ex ) {
        errstr = ex.what();
    }

    if( sample ) { free( sample ); sample = NULL; }
    if( cluster ) { delete cluster; cluster = NULL; }
    if( !errstr.empty())
        throw myruntime_error( errstr, CRITICAL );
}

// -------------------------------------------------------------------------
// GenSamplesForClusterLGSTNRM: generate samples for one cluster 
//  determined by variates; generate from logistic-normal distribution;
//  varndxs -- indexes of variates
//  vsize -- size of `varndxs'
//
void HDPSmplGenerator::GenSamplesForClusterLGSTNRM( int* varndxs, int vsize )
{
    mystring        preamb = "HDPSmplGenerator: GenSamplesForClusterLGSTNRM: ";
    mystring        errstr;
    int             dim = GetDimensions();
    int             adim = dim + 1;//augmented dimensions
    int             ctxtlen = GetContextLength();
    int             novars = GetNoVars();
    int             nocls = GetNoClusters();
    int             nosmps = GetNoSamplesC();
    int*            olbuf = GetOLBuffer();
    const double    csum = 0.9;
    const double    lovr = 1.0 - csum;
    double          cvalue, uncval, rv;
    double          fsum, diff;
    int             s, n, d, v, a;
    int             eco;
    double*         sample = NULL, *p, *pp;
    SimpleVector*   cluster = NULL;
    MTRng           rng;
    RMVNorm         rvnorm( rng );
    const double    accr = 1.e-6;
    static unsigned long tc = 0;

    if( ctxtlen < 1 || dim < 1 || novars < 1 || adim < novars || vsize != novars )
        throw myruntime_error( preamb + "Invalid number of variates.");

    if( varndxs == NULL || olbuf == NULL )
        throw myruntime_error( preamb + "Memory access error.");
    if( GetSamples() == NULL )
        throw myruntime_error( preamb + "Null samples structure.");

    if( !GetDoOverlap()) {
        for( v = 0; v < vsize; v++ ) {
            if( adim <= varndxs[v])
                throw myruntime_error( preamb + "Invalid variate index.");
            if( olbuf[varndxs[v]])
                return;
        }
        for( v = 0; v < vsize; v++ )
            olbuf[varndxs[v]] = 1;
    }

    rng.Set(( unsigned long )time( NULL )+tc++);

    Pslvector   amean( ctxtlen * adim );
    Pslvector   mean( ctxtlen * dim );
    Pslvector   rmv( dim );
    SPDmatrix   omega( ctxtlen );//covariance matrix between context positions
    SPDmatrix   sigma( dim );//covariance between dimension variates
    SPDmatrix   covm( ctxtlen * dim );

    try {
        omega.SetIdentity();

        //set covariance matrix:
        //off-diagonal elements
        for( d = 0; d < dim; d++ )
            for( a = 0; a < dim; a++ )
                sigma.SetValueAt( d, a, 0.5 );
        //diagonal elements
        for( d = 0; d < dim; d++ )
            sigma.SetValueAt( d, d, 1.0 );
        //covariated elements
        for( v = 0; v < vsize; v++ )
            for( a = v+1; a < vsize; a++ ) {
                if( dim <= varndxs[v])
                    continue;
                sigma.SetValueAt( varndxs[v], varndxs[a], 0.5 );
                sigma.SetValueAt( varndxs[a], varndxs[v], 0.5 );
            }
        //take Kronecker product
        covm.KronProduct( omega, sigma );

        if( adim == novars ) {
            cvalue = 1.0 /( double )novars;
            uncval = 0.0;
        }
        else {
            //leave fraction lovr to distribute among other variates
            cvalue = csum /( double )novars;
            uncval = lovr /( double )( adim - novars );
        }

        //set mean vector
        p = amean.GetVector();
        pp = mean.GetVector();
        for( n = 0; n < ctxtlen; n++, p += adim, pp += dim ) {
            fsum = 0.0;
            for( d = 0; d < adim; d++ ) {
                p[d] = uncval;
                fsum += uncval;
            }
            for( v = 0; v < vsize; v++ ) {
                if( adim <= varndxs[v])
                    throw myruntime_error( preamb + "Invalid variate index.");
                fsum += cvalue - p[varndxs[v]];
                p[varndxs[v]] = cvalue;
            }
            diff = 1.0 - fsum;
            if( diff < -accr || accr < diff )
                throw myruntime_error( preamb + "Unconserved mean vector.");
            //TRANSFORM to multivariate normal
            LogitNormal2Normal( p, ( size_t )adim, 1.e-6, false );
            //reduce dimenions of mean vector
            for( d = 0; d < dim; d++ )
                pp[d] = p[d];
        }

        //set mean vector and cov. matrix for generator
        rvnorm.SetMu( &mean );
        rvnorm.SetSigma( &covm );

        cluster = new SimpleVector( nosmps );
        if( cluster == NULL )
            throw myruntime_error( preamb + "Not enough memory.");

        for( s = 0; s < nosmps; s++ )
        {
            sample = ( double* )malloc( ctxtlen * adim * sizeof( double ));
            if( !sample )
                throw myruntime_error( preamb + "Not enough memory.");

            if(( eco = rvnorm.Gen( &rmv )) !=0 )
                throw myruntime_error( preamb + TranslatePSLError( eco ));

            //copy random variable to sample
            for( n = 0, p = sample; n < ctxtlen; n++, p += adim ) {
                for( d = 0; d < dim; d++ )
                    p[d] = rmv.GetValueAt(d);
                p[d] = 0.0;
            }

            //add generated sample
            cluster->Push( sample );
            sample = NULL;
            IncSampleCounter();
        }
        GetSamples()->Push( cluster );
        cluster = NULL;

    } catch( myexception const& ex ) {
        errstr = ex.what();
    }

    if( sample ) { free( sample ); sample = NULL; }
    if( cluster ) { delete cluster; cluster = NULL; }
    if( !errstr.empty())
        throw myruntime_error( errstr, CRITICAL );
}





// -------------------------------------------------------------------------
// Output: print samples to file
//
void HDPSmplGenerator::Output()
{
    FILE*   fp = NULL;
    mystring    preamb = "HDPSmplGenerator: Output: ";
    mystring    fname;

    if( !GetOutputPattern() || !strlen( GetOutputPattern()))
        throw myruntime_error( "HDPSmplGenerator: No filename pattern." );

    //-------
    fname = GetOutputPattern();
    fname += _sch_ext_crt_;
    fname += _sch_ext_grp_;

    fp = fopen( fname.c_str(), "w" );
    if( fp == NULL )
        throw myruntime_error( preamb + "Failed to open file '" + fname + "' for writing." );

    PrintCorrect( fp );
    fclose( fp );

    //-------
    fname = GetOutputPattern();
    fname += _sch_ext_mix_;
    fname += _sch_ext_grp_;

    fp = fopen( fname.c_str(), "w" );
    if( fp == NULL )
        throw myruntime_error( preamb + "Failed to open file '" + fname + "' for writing." );

    PrintMixed( fp );
    fclose( fp );

    //-------
    fname = GetOutputPattern();
    fname += _sch_ext_one_;
    fname += _sch_ext_grp_;

    fp = fopen( fname.c_str(), "w" );
    if( fp == NULL )
        throw myruntime_error( preamb + "Failed to open file '" + fname + "' for writing." );

    PrintOne( fp );
    fclose( fp );

    //-------
    fname = GetOutputPattern();
    fname += _sch_ext_two_;
    fname += _sch_ext_grp_;

    fp = fopen( fname.c_str(), "w" );
    if( fp == NULL )
        throw myruntime_error( preamb + "Failed to open file '" + fname + "' for writing." );

    PrintAny( fp, 2 );
    fclose( fp );

    //-------
    fname = GetOutputPattern();
    fname += _sch_ext_th3_;
    fname += _sch_ext_grp_;

    fp = fopen( fname.c_str(), "w" );
    if( fp == NULL )
        throw myruntime_error( preamb + "Failed to open file '" + fname + "' for writing." );

    PrintAny( fp, 3 );
    fclose( fp );
}

// -------------------------------------------------------------------------
// PrintCorrect: print samples clustered correctly to file
//
void HDPSmplGenerator::PrintCorrect( FILE* fp )
{
    if( fp == NULL )
        return;

    if( !GetSamples())
        return;

    int             dim = GetDimensions();
    int             adim = dim + 1;
    int             ctxtlen = GetContextLength();
    int             novars = GetNoVars();
    int             nocls = GetNoClusters();
    int             nosmps = GetNoSamplesC();
    SimpleVector*   cluster;
    double*         smp, *p;

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
    int tt, dd, nogrps, notss;
    int n, m, c, v;

    fprintf( fp, "%s %d; %s %d; %s %d; %s %d\n", 
              patstrnog, GetSamples()->GetSize(), 
              patstrnds, GetSamples()->GetSize(),
              patstrtns, GetSampleCounter(), patstrssz, ctxtlen );

    nogrps = 0;
    for( dd = n = 0; n < ( int )GetSamples()->GetSize(); n++ ) {
        notss = SLC_MIN( nosmps, ( int )GetSamples()->GetSize());
        fprintf( fp, "%s %d; %s %d\n", patstrgrp, nogrps++, patstrnot, notss );
        for( tt = m = 0; m < notss; m++ ) {
            //distinct cluster each iteration
            cluster = ( SimpleVector* )GetSamples()->GetValueAt( m );
            if( cluster == NULL )
                throw myruntime_error("HDPSmplGenerator: PrintCorrect: Null cluster.");
            if( nosmps != ( int )cluster->GetSize())
                throw myruntime_error("HDPSmplGenerator: PrintCorrect: Invalid number of samples in cluster.");
            smp = ( double* )cluster->GetValueAt( n );
            if( smp == NULL )
                throw myruntime_error("HDPSmplGenerator: PrintCorrect: Null sample.");
            fprintf( fp, "%s %d (%s %d; %s %d); %s %d\n", 
                     patstrtbl, tt++, patstrtid, m, patstrdid, m,
                     patstrnos, 1 );
            for( c = 0, p = smp; c < ctxtlen; c++, p += adim ) {
                for( v = 0; v < dim; v++ )
                    fprintf( fp, " %12.6g", p[v]);
                fprintf( fp, "\n" );
            }
            fprintf( fp, "%s\n", patstrendofsmp );
            fprintf( fp, "%s\n", patstrendoftbl );//end of table
        }
        fprintf( fp, "%s\n\n", patstrendofgrp );
    }
}

// -------------------------------------------------------------------------
// PrintMixed: print samples mixed over clusters to file
//
void HDPSmplGenerator::PrintMixed( FILE* fp )
{
    if( fp == NULL )
        return;

    if( !GetSamples())
        return;

    int             dim = GetDimensions();
    int             adim = dim + 1;
    int             ctxtlen = GetContextLength();
    int             novars = GetNoVars();
    int             nocls = GetNoClusters();
    int             nosmps = GetNoSamplesC();
    SimpleVector*   cluster;
    double*         smp, *p;

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
    int tt, dd, nogrps, notss;
    int n, m, c, v;

    fprintf( fp, "%s %d; %s %d; %s %d; %s %d\n", 
              patstrnog, GetSamples()->GetSize(), 
              patstrnds, GetSamples()->GetSize() * nosmps,
              patstrtns, GetSampleCounter(), patstrssz, ctxtlen );

    nogrps = 0;
    for( dd = n = 0; n < ( int )GetSamples()->GetSize(); n++ ) {
        notss = SLC_MIN( nosmps, ( int )GetSamples()->GetSize());
        fprintf( fp, "%s %d; %s %d\n", patstrgrp, nogrps++, patstrnot, notss );
        for( tt = m = 0; m < notss; m++ ) {
            //distinct cluster each iteration
            cluster = ( SimpleVector* )GetSamples()->GetValueAt( m );
            if( cluster == NULL )
                throw myruntime_error("HDPSmplGenerator: PrintMixed: Null cluster.");
            if( nosmps != ( int )cluster->GetSize())
                throw myruntime_error("HDPSmplGenerator: PrintMixed: Invalid number of samples in cluster.");
            smp = ( double* )cluster->GetValueAt( n );
            if( smp == NULL )
                throw myruntime_error("HDPSmplGenerator: PrintMixed: Null sample.");
            fprintf( fp, "%s %d (%s %d; %s %d); %s %d\n", 
                     patstrtbl, tt++, patstrtid, m, patstrdid, dd++/*grow no.clusters*/,
                     patstrnos, 1 );
            for( c = 0, p = smp; c < ctxtlen; c++, p += adim ) {
                for( v = 0; v < dim; v++ )
                    fprintf( fp, " %12.6g", p[v]);
                fprintf( fp, "\n" );
            }
            fprintf( fp, "%s\n", patstrendofsmp );
            fprintf( fp, "%s\n", patstrendoftbl );//end of table
        }
        fprintf( fp, "%s\n\n", patstrendofgrp );
    }
}

// -------------------------------------------------------------------------
// PrintOne: print samples gathered in one cluster to file
//
void HDPSmplGenerator::PrintOne( FILE* fp )
{
    if( fp == NULL )
        return;

    if( !GetSamples())
        return;

    int             dim = GetDimensions();
    int             adim = dim + 1;
    int             ctxtlen = GetContextLength();
    int             novars = GetNoVars();
    int             nocls = GetNoClusters();
    int             nosmps = GetNoSamplesC();
    SimpleVector*   cluster;
    double*         smp, *p;

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
    int nogrps, notss;
    int n, m, c, v;

    fprintf( fp, "%s %d; %s %d; %s %d; %s %d\n", 
              patstrnog, GetSamples()->GetSize(), 
              patstrnds, 1/*single cluster*/,
              patstrtns, GetSampleCounter(), patstrssz, ctxtlen );

    nogrps = 0;
    for( n = 0; n < ( int )GetSamples()->GetSize(); n++ ) {
        fprintf( fp, "%s %d; %s %d\n", patstrgrp, nogrps++, patstrnot, 1/*one table*/ );
        notss = SLC_MIN( nosmps, ( int )GetSamples()->GetSize());
        fprintf( fp, "%s %d (%s %d; %s %d); %s %d\n", 
                  patstrtbl, 0/*single table*/, patstrtid, 0/*tbl.id*/, patstrdid, 0/*single cluster*/,
                  patstrnos, notss );
        for( m = 0; m < notss; m++ ) {
            //distinct cluster each iteration
            cluster = ( SimpleVector* )GetSamples()->GetValueAt( m );
            if( cluster == NULL )
                throw myruntime_error("HDPSmplGenerator: PrintOne: Null cluster.");
            if( nosmps != ( int )cluster->GetSize())
                throw myruntime_error("HDPSmplGenerator: PrintOne: Invalid number of samples in cluster.");
            smp = ( double* )cluster->GetValueAt( n );
            if( smp == NULL )
                throw myruntime_error("HDPSmplGenerator: PrintOne: Null sample.");
            for( c = 0, p = smp; c < ctxtlen; c++, p += adim ) {
                for( v = 0; v < dim; v++ )
                    fprintf( fp, " %12.6g", p[v]);
                fprintf( fp, "\n" );
            }
            fprintf( fp, "%s\n", patstrendofsmp );
        }
        fprintf( fp, "%s\n", patstrendoftbl );//end of table
        fprintf( fp, "%s\n\n", patstrendofgrp );
    }
}

// -------------------------------------------------------------------------
// PrintOne: print samples randomly distributed among the number of 
//  clusters `noclsts' to file
//
void HDPSmplGenerator::PrintAny( FILE* fp, int noclsts )
{
    if( fp == NULL )
        return;

    if( !GetSamples())
        return;

    int             dim = GetDimensions();
    int             adim = dim + 1;
    int             ctxtlen = GetContextLength();
    int             novars = GetNoVars();
    int             nocls = GetNoClusters();//correct no.clusters
    int             nosmps = GetNoSamplesC();
    SimpleVector*   cluster;
    double*         smp, *p;

    mystring preamb = "HDPSmplGenerator: PrintAny: ";
    mystring merror;
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
    time_t      tm;
    MTRng       lrng_di;//random number generator
    RVDisc      rvdis( lrng_di, RVDisc::TRVDsc_Inverse );
    double*     dprobs = NULL;
    int tt, dd, nogrps, notss;
    int n, m, c, v, cerr;

    if( noclsts < 1 || GetSampleCounter() < noclsts )
        throw myruntime_error( preamb + "Invalid no. clusters to produce.");

    dprobs = ( double* )malloc( noclsts * sizeof( double ));
    if( dprobs == NULL )
        throw myruntime_error( preamb + "Not enough memory.");

    time( &tm );
    lrng_di.Set((unsigned long)(size_t)this + (unsigned long)(size_t)(&lrng_di) + (unsigned long)tm );

    rvdis.SetProbs( dprobs, noclsts );

    dprobs[0] = 1.0 / (double)noclsts;
    for( dd = 1; dd < noclsts; dd++ )
        dprobs[dd] = dprobs[0];

    fprintf( fp, "%s %d; %s %d; %s %d; %s %d\n", 
              patstrnog, GetSamples()->GetSize(), 
              patstrnds, noclsts/*clusters*/,
              patstrtns, GetSampleCounter(), patstrssz, ctxtlen );

    try {
        nogrps = 0;
        for( n = 0; n < ( int )GetSamples()->GetSize(); n++ ) {
            notss = SLC_MIN( nosmps, ( int )GetSamples()->GetSize());
            fprintf( fp, "%s %d; %s %d\n", patstrgrp, nogrps++, patstrnot, notss );
            for( tt = m = 0; m < notss; m++ ) {
                //distinct cluster each iteration
                cluster = ( SimpleVector* )GetSamples()->GetValueAt( m );
                if( cluster == NULL )
                    throw myruntime_error( preamb + "Null cluster.");
                if( nosmps != ( int )cluster->GetSize())
                    throw myruntime_error( preamb + "Invalid number of samples in cluster.");
                smp = ( double* )cluster->GetValueAt( n );
                if( smp == NULL )
                    throw myruntime_error( preamb + "Null sample.");

                //sample dish index
                if(( cerr = rvdis.GenI( &dd )) != PSL_OK )
                    throw myruntime_error( preamb + TranslatePSLError( cerr ));
                //print table
                fprintf( fp, "%s %d (%s %d; %s %d); %s %d\n", 
                        patstrtbl, tt++, patstrtid, m, patstrdid, dd, patstrnos, 1 );

                for( c = 0, p = smp; c < ctxtlen; c++, p += adim ) {
                    for( v = 0; v < dim; v++ )
                        fprintf( fp, " %12.6g", p[v]);
                    fprintf( fp, "\n" );
                }
                fprintf( fp, "%s\n", patstrendofsmp );
                fprintf( fp, "%s\n", patstrendoftbl );//end of table
            }
            fprintf( fp, "%s\n\n", patstrendofgrp );
        }
    } catch( myexception const& ex ) {
        merror = ex.what();
    }

    free( dprobs );
    dprobs = NULL;

    if( !merror.empty())
        throw myruntime_error( merror );
}

