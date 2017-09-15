/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __IMAClusters__
#define __IMAClusters__

#include "compdef.h"
#include "debug.h"
#include "rc.h"
#include "data.h"
#include "mystring.h"
#include "myexcept.h"
#include "InputMultipleAlignment.h"
#include "IMACounts.h"


// _________________________________________________________________________
// CLASS IMAClusters
//
class IMAClusters: public InputMultipleAlignment
{
public:
    IMAClusters( unsigned = ALLOCSEQ );
    virtual ~IMAClusters();

    void    Make();

    void    OutputCounts( const char* filename );

    double  GetClusteringThreshold() const          { return imaclusterat_; }
    void    SetClusteringThreshold( double value )  { imaclusterat_ = value; }

    double  GetIgnoreGapPercentage() const          { return imagapignore_; }
    void    SetIgnoreGapPercentage( double value )  { imagapignore_ = value; }

protected:
    void    MakeClusters();
    void    ClusterSame( PosDescriptionVector*, PosDescriptionVector* );
    void    DetermineClusterSizes();
    void    ComputeFrequencies();

    IMACounts&          GetObsCounts()          { return imaobscounts_; }
    const IMACounts&    GetObsCounts() const    { return imaobscounts_; }

    size_t  GetCurrentNoClusters() const    { return imanoclusters_; }
    void    IncCurrentNoClusters()          { imanoclusters_++; }
    void    ResetCurrentNoClusters()        { imanoclusters_ = 0; }

    size_t  GetActualNoClusters() const     { return imactualnoclusters_; }
    void    IncActualNoClusters()           { imactualnoclusters_++; }
    void    ResetActualNoClusters()         { imactualnoclusters_ = 0; }

    double  GetAvgClusterSize() const           { return imavgclustersc_; }
    void    SetAvgClusterSize( double value )   { imavgclustersc_ = value; }

    void    AllocClusterSz( size_t );
    void    DestroyClusterSz();

    void    PurgePositions();               //...allocation of used positions as well
    void    DestroyUsedPosits();

    bool    GetUsedPositionAt( size_t pos ) const;

    size_t  GetClusterSzAt( size_t pos ) const;
    void    IncClusterSzAt( size_t pos );
    void    SetClusterSzAt( size_t pos, size_t value );

    size_t  GetMinClusterSz() const         { return imamincluster_; }
    void    SetMinClusterSz( size_t value ) { imamincluster_ = value; }

    size_t  GetMaxClusterSz() const         { return imamaxcluster_; }
    void    SetMaxClusterSz( size_t value ) { imamaxcluster_ = value; }

    void    ResetMinMaxClusterSz()          { imamincluster_ = ( size_t )-1; imamaxcluster_ = 0; }

    bool    GetProcessed() const { return processed_; }
    void    SetProcessed( bool value ) { processed_ = value; }

private:
    IMACounts   imaobscounts_;      //observed amino acid and its pair counts
    double      imaclusterat_;      //clustering threshold
    double      imagapignore_;      //threshold to ignore positions with gaps exceeding that percentage
    bool*       imausedposits_;     //vector of used positions
    size_t*     imaclustersz_;      //vector of cluster sizes
    size_t      imaclusterszalloc_; //allocated number of vector cells for imaclustersz_
    size_t      imamincluster_;     //minimum cluster size
    size_t      imamaxcluster_;     //maximum cluster size
    size_t      imanoclusters_;     //current number of clusters
    size_t      imactualnoclusters_;//actual number of clusters
    double      imavgclustersc_;    //average cluster size
    bool        processed_;         //flag of precssing
};

////////////////////////////////////////////////////////////////////////////
// Class IMAClusters inlines
//
// -------------------------------------------------------------------------
// GetIgnorePositionAt: returns flag of whether a position is to be used

inline
bool IMAClusters::GetUsedPositionAt( size_t pos ) const
{
#ifdef __DEBUG__
    if( imausedposits_ == NULL ||
        SequenceAt( 0 ) == NULL ||
        SequenceAt( 0 )->size() <= pos )
        throw myruntime_error( mystring( "IMAClusters: Memory access error." ));
#endif
    return imausedposits_[pos];
}

// -------------------------------------------------------------------------
// GetClusterSzAt: returns size of cluster specified by a position

inline
size_t IMAClusters::GetClusterSzAt( size_t pos ) const
{
#ifdef __DEBUG__
    if( imaclustersz_ == NULL || imaclusterszalloc_ <= pos )
        throw myruntime_error( mystring( "IMAClusters: Memory access error." ));
#endif
    return imaclustersz_[pos];
}

// IncClusterSzAt: increases by one size of cluster at a position
//
inline
void IMAClusters::IncClusterSzAt( size_t pos )
{
#ifdef __DEBUG__
    if( imaclustersz_ == NULL || imaclusterszalloc_ <= pos )
        throw myruntime_error( mystring( "IMAClusters: Memory access error." ));
#endif
    imaclustersz_[pos]++;
}

// SetClusterSzAt: sets size of cluster at a position
//
inline
void IMAClusters::SetClusterSzAt( size_t pos, size_t value )
{
#ifdef __DEBUG__
    if( imaclustersz_ == NULL || imaclusterszalloc_ <= pos )
        throw myruntime_error( mystring( "IMAClusters: Memory access error." ));
#endif
    imaclustersz_[pos] = value;
}



#endif//__IMAClusters__

