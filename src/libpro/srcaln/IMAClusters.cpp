/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "IMAClusters.h"


////////////////////////////////////////////////////////////////////////////
// CLASS IMAClusters
//
// Constructor
//

IMAClusters::IMAClusters( unsigned reservation )
:   InputMultipleAlignment( reservation ),
    imaobscounts_(),
    imaclusterat_( CLUSTERING_THRESHOLD ),
    imagapignore_( 1.0 ),
    imausedposits_( NULL ),
    imaclustersz_( NULL ),
    imaclusterszalloc_( 0 ),
    imamincluster_( 0 ),
    imamaxcluster_( 0 ),
    imanoclusters_( 0 ),
    imactualnoclusters_( 0 ),
    imavgclustersc_( 0.0 ),
    processed_( false )
{
}

// -------------------------------------------------------------------------
// Destructor
// -------------------------------------------------------------------------

IMAClusters::~IMAClusters()
{
    DestroyUsedPosits();
    DestroyClusterSz();
}

// -------------------------------------------------------------------------
// Make: Make clusters and compute observed frequencies
// -------------------------------------------------------------------------

void IMAClusters::Make()
{
    SetProcessed( false );
//     PreprocessAlignment();
    if( !size() )
        throw myruntime_error( mystring( "IMAClusters: No sequences in multiple sequence alignment." ));

    PurgeAtSequenceIdentity();

    InitQueryDescription();

    if( GetUsingSEGFilter())
        RefineWithSEG();

    ResetMinMaxClusterSz();
    ResetCurrentNoClusters();
    ResetActualNoClusters();

    AllocClusterSz( size());

    PurgePositions();
    MakeClusters();
    ComputeFrequencies();
}

// -------------------------------------------------------------------------
// AllocClusterSz: allocates memory for vector of cluster sizes

void IMAClusters::AllocClusterSz( size_t nocells )
{
    DestroyClusterSz();

    imaclustersz_ = ( size_t* )malloc( sizeof( size_t ) * nocells );

    if( imaclustersz_ == NULL )
        throw myruntime_error( mystring( "IMAClusters: Not enough memory." ));

    imaclusterszalloc_ = nocells;

    memset( imaclustersz_, 0, sizeof( size_t ) * nocells );
}

// -------------------------------------------------------------------------
// DestroyClusterSz: deallocates memory of vector of cluster sizes

void IMAClusters::DestroyClusterSz()
{
    if( imaclustersz_ ) {
        free( imaclustersz_ );
        imaclustersz_ = NULL;
        imaclusterszalloc_ = 0;
    }
}

// -------------------------------------------------------------------------
// DestroyUsedPosits: deallocates memory of vector of used positions
// -------------------------------------------------------------------------

void IMAClusters::DestroyUsedPosits()
{
    if( imausedposits_ ) {
        free( imausedposits_ );
        imausedposits_ = NULL;
    }
}

// -------------------------------------------------------------------------
// PurgePositions: Purges positions that contain gap with percentage
//     greater than class's object allows
// -------------------------------------------------------------------------

void IMAClusters::PurgePositions()
{
    PosDescriptionVector*   query = SequenceAt( 0 );

    if( query == NULL )
        throw myruntime_error( mystring( "IMAClusters: No sequences." ));

    size_t                  no_ress;
    size_t                  no_gaps;
    size_t                  usize = query->size();
    PosDescriptionVector*   seq = NULL;

    if( !usize )
        throw myruntime_error( mystring( "IMAClusters: Sequence of length 0." ));

    DestroyUsedPosits();


    imausedposits_ = ( bool* )malloc( sizeof( bool ) * usize );

    if( imausedposits_ == NULL )
        throw myruntime_error( mystring( "IMAClusters: Not enough memory." ));

    memset( imausedposits_, 0, sizeof( bool ) * usize );

    for( size_t p = 0; p < usize; p++ )
    {
        no_ress = no_gaps = 0;

        for( size_t i = 0; i < size(); i++ ) {
            seq = SequenceAt( i );
            if( seq == NULL || !seq->GetUsed())
                continue;

            no_ress++;

            if( seq->ResidueAt( p ) == GAP )
                no_gaps++;
        }
        if( no_ress && GetIgnoreGapPercentage() <= ( double ) no_gaps / no_ress )
            imausedposits_[p] = false;
        else
            imausedposits_[p] = true;
    }
}

// -------------------------------------------------------------------------
// MakeClusters: Make clusters given sequence identity percentage at which
//     clustering to be performed
// -------------------------------------------------------------------------

void IMAClusters::MakeClusters()
{
    size_t  segment;    //interval of comparable pair segments
    size_t  matches;    //number of residue matches in the sequences
    bool            u1, u2;
    unsigned char   r1, r2;
    size_t  i, j;

    PosDescriptionVector*   query = SequenceAt( 0 );

    if( query == NULL )
        throw myruntime_error( mystring( "IMAClusters: No sequences." ));

    if( size() == 0 )
        throw myruntime_error( mystring( "IMAClusters: No sequences." ));


    size_t  usize = query->size();


    for( i = 0; i + 1 < size(); i++ )
    {
        PosDescriptionVector*   first = SequenceAt( i );
        if( first == NULL || !first->GetUsed())
            continue;

        for( j = i + 1; j < size(); j++ )
        {
            PosDescriptionVector*   secnd = SequenceAt( j );
            if( secnd == NULL || !secnd->GetUsed())
                continue;

            segment = 0;
            matches = 0;

            for( size_t p = 0; p < usize; p++ )
            {
                if( ! GetUsedPositionAt( p ))
                    continue;

                u1 = first->IsUsedAt( p );
                u2 = secnd->IsUsedAt( p );
                r1 = first->ResidueAt( p );
                r2 = secnd->ResidueAt( p );

                if( !u1 || !u2 ) {
                    continue;
                }
                if( r1 == X || r2 == X ) {
                    continue;
                }

                if( r1 == GAP || r2 == GAP ) {
                    continue;
                }

                segment++;
                if( r1 == r2 )
                    matches++;
            }
            if( segment && GetClusteringThreshold() <= ( double( matches ) / segment ))
                ClusterSame( first, secnd );
        }
    }

    DetermineClusterSizes();
}

// -------------------------------------------------------------------------
// ClusterSame: Unite pair of sequences in the same cluster
// -------------------------------------------------------------------------

void IMAClusters::ClusterSame( PosDescriptionVector* first, PosDescriptionVector* secnd )
{
#ifdef __DEBUG__
    if( !first || !secnd )
        throw myruntime_error( mystring( "IMAClusters: Wrong pair of sequences to cluster." ));
#endif

    size_t  firstclus;
    size_t  secndclus;

    if( GetCurrentNoClusters() <= first->GetCluster()) {
        if( GetCurrentNoClusters() <= secnd->GetCluster())
        {
            first->SetCluster( GetCurrentNoClusters());
            secnd->SetCluster( GetCurrentNoClusters());
            IncCurrentNoClusters();
        }
        else
            first->SetCluster( secnd->GetCluster());
    }
    else
        if( GetCurrentNoClusters() <= secnd->GetCluster())
            secnd->SetCluster( first->GetCluster());
        else {
            //both sequences have already been assigned to clusters
            firstclus = first->GetCluster();
            secndclus = secnd->GetCluster();
            if( secndclus < firstclus ) {
                secndclus = first->GetCluster();
                firstclus = secnd->GetCluster();
            }
            //merge two clusters
            for( size_t n = 0; n < size(); n++ )
                if( SequenceAt( n ) &&
                    SequenceAt( n )->GetCluster() == secndclus )
                    SequenceAt( n )->SetCluster( firstclus );
        }
}

// -------------------------------------------------------------------------
// DetermineClusterSizes: determine sizes of clusters
// -------------------------------------------------------------------------

void IMAClusters::DetermineClusterSizes()
{
    size_t  no_clusseqs = 0;

    for( size_t n = 0; n < size(); n++ )
    {
        PosDescriptionVector*   seqn = SequenceAt( n );
        if( !seqn || !seqn->GetUsed())
            continue;

        if( GetCurrentNoClusters() <= seqn->GetCluster()) {
            seqn->SetCluster( GetCurrentNoClusters());
            IncCurrentNoClusters();
        }
        no_clusseqs++;
        IncClusterSzAt( seqn->GetCluster());
        if( GetClusterSzAt( seqn->GetCluster()) == 1 )
            IncActualNoClusters();
    }

    if( GetActualNoClusters())
        SetAvgClusterSize(( double )no_clusseqs / GetActualNoClusters());
}

// -------------------------------------------------------------------------
// ComputeFrequencies: Computes frequencies after clusers have been made
// -------------------------------------------------------------------------

void IMAClusters::ComputeFrequencies()
{
    bool            u1, u2;
    unsigned char   r1, r2;
    size_t          firstclusz;
    size_t          secndclusz;
    double          weight1, weight2;
    double          weight;

    PosDescriptionVector*   query = SequenceAt( 0 );

    if( query == NULL )
        throw myruntime_error( mystring( "IMAClusters: No sequences." ));

    if( size() == 0 )
        throw myruntime_error( mystring( "IMAClusters: No sequences." ));

    const size_t    nosqs = size();
    const size_t    usize = query->size();
    bool            smask[nosqs][usize];
    size_t  i, j, p;

    for( i = 0; i < nosqs; i++ )
        for( p = 0; p < usize; p++ )
            smask[i][p] = false;

    for( i = 0; i + 1 < nosqs; i++ )
    {
        PosDescriptionVector*   first = SequenceAt( i );
        if( first == NULL || !first->GetUsed())
            continue;

//         if( GetCurrentNoClusters() <= first->GetCluster())
//             continue;

        firstclusz = GetClusterSzAt( first->GetCluster());

        if( firstclusz == 0 )
            continue;

        if( GetMaxClusterSz() < firstclusz )    SetMaxClusterSz( firstclusz );
        if( firstclusz < GetMinClusterSz())     SetMinClusterSz( firstclusz );

        weight1 = 1.0 / ( double )firstclusz;

        for( j = i + 1; j < nosqs; j++ )
        {
            PosDescriptionVector*   secnd = SequenceAt( j );
            if( secnd == NULL || !secnd->GetUsed())
                continue;

//             if( GetCurrentNoClusters() <= secnd->GetCluster())
//                 continue;

            if( first->GetCluster() == secnd->GetCluster())
                continue;

            secndclusz = GetClusterSzAt( secnd->GetCluster());

            if( secndclusz == 0 )
                continue;

            if( i == 0 && j + 1 == nosqs ) {
                if( GetMaxClusterSz() < secndclusz )    SetMaxClusterSz( secndclusz );
                if( secndclusz < GetMinClusterSz())     SetMinClusterSz( secndclusz );
            }

            weight2 = 1.0 / ( double )secndclusz;
            weight = weight1 * weight2;

            for( p = 0; p < usize; p++ )
            {
                if( ! GetUsedPositionAt( p ))
                    continue;

                u1 = first->IsUsedAt( p );
                u2 = secnd->IsUsedAt( p );
                r1 = first->ResidueAt( p );
                r2 = secnd->ResidueAt( p );

                if( !u1 || !u2 )
                    continue;
//                 if( r1 == X || r1 == GAP || r2 == X || r2 == GAP )
//                     continue;
                //take into account explicit known residues only 
                if( NUMAA <= r1 || NUMAA <= r2 )
                    continue;

                if( !GetProcessed())
                    SetProcessed( true );

                if( !smask[i][p]) {
                    GetObsCounts().IncResidCountBy( r1, weight1 );
                    smask[i][p] = true;
                }
                if( !smask[j][p]) {
                    GetObsCounts().IncResidCountBy( r2, weight2 );
                    smask[j][p] = true;
                }
                GetObsCounts().IncPairCountBy( r1, r2, weight );
            }
        }
    }
}

// -------------------------------------------------------------------------
// OutputCounts: Outputs all information of observed amino acid counts
// -------------------------------------------------------------------------

void IMAClusters::OutputCounts( const char* filename )
{
    FILE*   fp = stdout;

    if( !GetProcessed()) {
        warning("No counts obtained (too few clusters?).");
        return;
    }

    if( filename && strlen( filename ))
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw myruntime_error( mystring( "IMAClusters: Failed to open file for writing." ));

    fprintf( fp, "## Alignment processed: %s\n", GetName()? GetName(): "" );
    if( GetIgnoreGapsInQuery())
        fprintf( fp, "##  (Positions with gaps in query ignored)\n");
    fprintf( fp, "##  (Positions with %d%% gap content ignored)\n", ( int )rint( GetIgnoreGapPercentage() * 100.0 ));
    fprintf( fp, "##  - Total sequences,      %10d\n", size());
    fprintf( fp, "##  - Clusters constructed, %10d (Clustering threshold %3d%%)\n",
                        GetActualNoClusters(), ( int )rint( GetClusteringThreshold() * 100.0 ));
    fprintf( fp, "##  - Minimum cluster size, %10d\n", GetMinClusterSz());
    fprintf( fp, "##  - Maximum cluster size, %10d\n", GetMaxClusterSz());
    fprintf( fp, "##  - Average cluster size, %10.2lf\n##\n", GetAvgClusterSize());

    GetObsCounts().Write( fp );

    if( fp != stdout )
        fclose( fp );
}

