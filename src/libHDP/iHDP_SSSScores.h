/***************************************************************************
 *   Copyright (C) 2013 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __iHDP_SSSScores__
#define __iHDP_SSSScores__

#include <stdio.h>
#include <stdlib.h>

#include "ext/psl.h"
#include "ext/pslvector.h"
#include "rc.h"
#include "myexcept.h"

#if !defined( USEiHDPSSSSCORES )
//#define USEiHDPSSSSCORES
#endif

//global object for application
class iHDP_SSSScores;
extern iHDP_SSSScores iHDPSSSSCORES;

// -------------------------------------------------------------------------
// class iHDP_SSSScores: alignment-independent SS state and Hierarchical 
//  Dirichlet Process model-derived scores
//
class iHDP_SSSScores
{
public:
    iHDP_SSSScores();
    ~iHDP_SSSScores();

    double          GetiHDPSSSWeight() const { return weight_; }
    void            SetiHDPSSSWeight( double value ) { weight_ = value; }

    int             GetCardinality() const { return card_; }
    int             GetSSSCardinality() const { return ssscard_; }
    int             GetNoEffLvs() const { return noeffs_; }
    int             GetNoPrbLvs() const { return noplvs_; }
    int             GetNoSSSPrbLvs() const { return nosssplvs_; }

    double          GetScore( int clid, int sssid, double prb, double sssprb, double ens ) const;
    void            ReadScores( const char* filename );

    bool            ScoreNA( double value ) const { return value <= naval_; }

protected:
    void            ReadLevels( FILE* fp );
    void            ReadScoreTables( const char* filename );
    void            ReadScoresHelper( FILE* fp, int eff );
    void            ReadMatrix( FILE* fp, Pslvector* scos );

    void            SetCardinality( int value ) { card_ = value; }
    void            SetSSSCardinality( int value ) { ssscard_ = value; }
    double          GetNAvalue() const { return naval_; }

    Pslvector*      GetScores( int lvl, int plvl, int sssplvl );
    const Pslvector*GetScores( int lvl, int plvl, int sssplvl ) const;
    void            NewScores( int nolevels, int noplvs, int nosssplvs );
    void            DestroyScores();

private:
    Pslvector**** scores_;//scores
    double      weight_;//weight of all score tables
    int         noeffs_;//number of levels for effective no. seqs 
    int         noplvs_;//number of probability levels
    int         nosssplvs_;//number of SS state probability levels
    Pslvector   levels_;//levels of eff. no. sequences
    Pslvector   prblvs_;//probability (of clusters) levels 
    Pslvector   sssprblvs_;//SS state probability levels 
    const double naval_;//value for NA field
    int         card_;//cardinality
    int         ssscard_;//SS state cardinality
};

// -------------------------------------------------------------------------
// INLINES
//
// DestroyScores: destroy all score tables
inline
void iHDP_SSSScores::DestroyScores()
{
    int n, p, ssp;
    if( scores_ ) {
        for( n = 0; n < noeffs_; n++ ) {
            if( scores_[n] == NULL )
                continue;
            for( p = 0; p < noplvs_; p++ ) {
                if( scores_[n][p] == NULL )
                    continue;
                for( ssp = 0; ssp < nosssplvs_; ssp++ ) {
                    if( scores_[n][p][ssp])
                        delete scores_[n][p][ssp];
                    scores_[n][p][ssp] = NULL;
                }
                free( scores_[n][p]);
                scores_[n][p] = NULL;
            }
            free( scores_[n]);
            scores_[n] = NULL;
        }
        free( scores_ );
        scores_ = NULL;
    }
    nosssplvs_ = 0;
    noplvs_ = 0;
    noeffs_ = 0;
}
// NewScores: allocate new score table for each level of eff. no sequences
//  nolevels: no. levels for eff. no. sequences
//  noplvs: no. probability levels
//  nosssplvs: no. SS state probability levels
inline
void iHDP_SSSScores::NewScores( int nolevels, int noplvs, int nosssplvs )
{
    int n, p, ssp;
    DestroyScores();
    if( nolevels < 1 || noplvs < 1 || nosssplvs < 1 )
        return;
    scores_ = ( Pslvector**** )malloc( nolevels * sizeof( Pslvector* ));
    if( scores_ == NULL )
        throw myruntime_error("iHDP_SSSScores: NewScores: Not enough memory.");
    noeffs_ = nolevels;
    noplvs_ = noplvs;
    nosssplvs_ = nosssplvs;

    for( n = 0; n < noeffs_; n++ ) {
        scores_[n] = ( Pslvector*** )malloc( noplvs * sizeof( Pslvector* ));
        if( scores_[n] == NULL )
            throw myruntime_error("iHDP_SSSScores: NewScores: Not enough memory.");
        for( p = 0; p < noplvs_; p++ ) {
            scores_[n][p] = ( Pslvector** )malloc( nosssplvs * sizeof( Pslvector* ));
            if( scores_[n][p] == NULL )
                throw myruntime_error("iHDP_SSSScores: NewScores: Not enough memory.");
            for( ssp = 0; ssp < nosssplvs_; ssp++ ) {
                scores_[n][p][ssp] = new Pslvector();
                if( scores_[n][p][ssp] == NULL )
                    throw myruntime_error("iHDP_SSSScores: NewScores: Not enough memory.");
            }
        }
    }
}

// -------------------------------------------------------------------------
// GetScores: get score table corresponding to the levels of eff. no 
//  seqns, and probabilities
inline
Pslvector* iHDP_SSSScores::GetScores( int lvl, int plvl, int sssplvl )
{
    if( lvl < 0 || noeffs_ <= lvl )
        throw myruntime_error("iHDP_SSSScores: GetScores: Memory access error.");
    if( plvl < 0 || noplvs_ <= plvl )
        throw myruntime_error("iHDP_SSSScores: GetScores: Memory access error.");
    if( sssplvl < 0 || nosssplvs_ <= sssplvl )
        throw myruntime_error("iHDP_SSSScores: GetScores: Memory access error.");
    return scores_[lvl][plvl][sssplvl];
}

// -------------------------------------------------------------------------
// GetScores: get score table corresponding to the levels of eff. no 
//  seqns, and probabilities
inline
const Pslvector* iHDP_SSSScores::GetScores( int lvl, int plvl, int sssplvl ) const
{
    if( lvl < 0 || noeffs_ <= lvl )
        throw myruntime_error("iHDP_SSSScores: GetScores: Memory access error.");
    if( plvl < 0 || noplvs_ <= plvl )
        throw myruntime_error("iHDP_SSSScores: GetScores: Memory access error.");
    if( sssplvl < 0 || nosssplvs_ <= sssplvl )
        throw myruntime_error("iHDP_SSSScores: GetScores: Memory access error.");
    return scores_[lvl][plvl][sssplvl];
}

// -------------------------------------------------------------------------
// GetScore: get score at table position identified by the HDP cluster and 
//  SS state; table is selected according to probability levels `prb' and 
//  `sssprb', and eff. no. sequences given by `ens'
inline
double iHDP_SSSScores::GetScore( int clid, int sssid, double prb, double sssprb, double ens ) const
{
    int     e, pi, spi;//level indices
    int     n, tmp;
    int     nolvs = levels_.GetSize();
    int     noplvs = SLC_MAX( 1, prblvs_.GetSize());
    int     nosssplvs = SLC_MAX( 1, sssprblvs_.GetSize());
    double  dtmp;
    double  pi1, pi2, ppi1, ppi2;
    double  e1, e2, ee1, ee2;
    const Pslvector* scos = NULL;
    static char locbuf[KBYTE];

    //translate given ss state and its probability to joint ss state id
    sssid = ( int )sssid * 10 + int( sssprb * 10.0 );

    if( clid < 0 || card_ <= clid )
        throw myruntime_error("iHDP_SSSScores: GetScore: Memory access error.");

    if( sssid < 0 || ssscard_ <= sssid )
        throw myruntime_error("iHDP_SSSScores: GetScore: Memory access error.");

    if( nolvs < 1 )
        throw myruntime_error("iHDP_SSSScores: GetScore: No levels for eff. no. sqns.");
    if( ens < levels_.GetValueAt(0))
        return 0.0;

    if( noplvs < 1 )
        throw myruntime_error("iHDP_SSSScores: GetScore: No probability levels.");
    if( prblvs_.GetSize())
        if( prb < prblvs_.GetValueAt(0))
            return 0.0;

    if( nosssplvs < 1 )
        throw myruntime_error("iHDP_SSSScores: GetScore: No SS state probability levels.");
    if( sssprblvs_.GetSize())
        if( sssprb < sssprblvs_.GetValueAt(0))
            return 0.0;

    //get corresponding score table
    //
    //index corresponding to level of eff. no. sequences
    for( e = 0; e < nolvs; e++ ) {
        e1 = levels_.GetValueAt(e);
        e2 = 20.01;
        if( e + 1 < nolvs )
            e2 = levels_.GetValueAt(e+1);
        if( e1 <= ens && ens < e2 )
            break;
    }
    if( nolvs <= e ) {
        sprintf( locbuf, "Segment covering eff. no. sequences not found: %g.", ens );
        throw myruntime_error("iHDP_SSSScores: GetScore: " + mystring( locbuf ));
    }

    //index corresponding to probability level
    pi = 0;
    if( prblvs_.GetSize()) {
        for( pi = 0; pi < noplvs; pi++ ) {
            pi1 = prblvs_.GetValueAt(pi);
            pi2 = 1.01;
            if( pi + 1 < noplvs )
                pi2 = prblvs_.GetValueAt(pi+1);
            if( pi1 <= prb && prb < pi2 )
                break;
        }
        if( noplvs <= pi ) {
            sprintf( locbuf, "Segment covering probability level value not found: %g.", prb );
            throw myruntime_error("iHDP_SSSScores: GetScore: " + mystring( locbuf ));
        }
    }

    //index corresponding to SS state probability level
    spi = 0;
    if( sssprblvs_.GetSize()) {
        for( spi = 0; spi < nosssplvs; spi++ ) {
            pi1 = sssprblvs_.GetValueAt(spi);
            pi2 = 1.01;
            if( pi + 1 < nosssplvs )
                pi2 = sssprblvs_.GetValueAt(spi+1);
            if( pi1 <= sssprb && sssprb < pi2 )
                break;
        }
        if( nosssplvs <= spi ) {
            sprintf( locbuf, "Segment covering SS state probability level value not found: %g.", 
                      sssprb );
            throw myruntime_error("iHDP_SSSScores: GetScore: " + mystring( locbuf ));
        }
    }


    //score table
    scos = GetScores( e, pi, spi );
    if( scos == NULL )
        throw myruntime_error("iHDP_SSSScores: GetScore: Null score table.");

    //calculate index within score table
    n = clid;
    if( sssid )
        n += card_ * sssid;

    if( scos->GetSize() <= n )
        throw myruntime_error("iHDP_SSSScores: GetScore: Memory access error.");
    return scos->GetValueAt(n);
}


#endif//__iHDP_SSSScores__
