/***************************************************************************
 *   Copyright (C) 2013 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __HDPscores__
#define __HDPscores__

#include <stdio.h>
#include <stdlib.h>

#include "ext/psl.h"
#include "ext/pslvector.h"
#include "rc.h"
#include "myexcept.h"

//global object for application
class HDPscores;
extern HDPscores HDPSCORES;
extern HDPscores HDPctSCORES;

// -------------------------------------------------------------------------
// class HDPscores: Hierarchical Dirichlet Process model-derived scores
//
class HDPscores
{
public:
    HDPscores();
    ~HDPscores();

    int             GetCardinality() const { return card_; }
    int             GetNoPrbLvs() const { return noplvs_; }
    int             GetNoTables() const { return notbls_; }

    double          GetScore( int row, int col, 
                              double fstprb, double secprb, double fstens, double secens ) const;
    void            ReadScores( const char* filename );

    bool            ScoreNA( double value ) const { return value <= naval_; }

protected:
    void            ReadLevels( FILE* fp );
    void            ReadScoreTables( const char* filename );
    void            ReadScoresHelper( FILE*, Pslvector* scos );

    void            SetCardinality( int value ) { card_ = value; }
    double          GetNAvalue() const { return naval_; }

    Pslvector*      GetScores( int plvl, int lvl );
    const Pslvector*GetScores( int plvl, int lvl ) const;
    void            NewScores( int noplvs, int nolevels );
    void            DestroyScores();

private:
    Pslvector*** scores_;//scores
    int         noplvs_;//number of probability levels
    int         notbls_;//number of score tables for each prob. level
    Pslvector   prblvs_;//probability (of clusters) levels 
    Pslvector   levels_;//levels of eff. no. sequences
    const double naval_;//value for NA field
    int         card_;//cardinality
};

// -------------------------------------------------------------------------
// INLINES
//
// DestroyScores: destroy all score tables
inline
void HDPscores::DestroyScores()
{
    int n, m;
    if( scores_ ) {
        for( m = 0; m < noplvs_; m++ ) {
            if( scores_[m] == NULL )
                continue;
            for( n = 0; n < notbls_; n++ ) {
                if( scores_[m][n])
                    delete scores_[m][n];
                scores_[m][n] = NULL;
            }
            free( scores_[m] );
            scores_[m] = NULL;
        }
        free( scores_ );
        scores_ = NULL;
        notbls_ = 0;
    }
    noplvs_ = 0;
}
// NewScores: allocate new score table for each level of eff. no sequences
//  noplvs: no. probability level values
//  nolevels: no. level values of eff. no. sequences
inline
void HDPscores::NewScores( int noplvs, int nolevels )
{
    int n, m;
    DestroyScores();
    if( noplvs < 1 || nolevels < 1 )
        return;
    scores_ = ( Pslvector*** )malloc( noplvs * sizeof( Pslvector** ));
    if( scores_ == NULL )
        throw myruntime_error("HDPscores: NewScores: Not enough memory.");
    noplvs_ = noplvs;
    notbls_ = nolevels;

    for( m = 0; m < noplvs_; m++ ) {
        scores_[m] = ( Pslvector** )malloc( nolevels * sizeof( Pslvector* ));
        if( scores_[m] == NULL )
            throw myruntime_error("HDPscores: NewScores: Not enough memory.");
        for( n = 0; n < notbls_; n++ ) {
            scores_[m][n] = new Pslvector();
            if( scores_[m][n] == NULL )
                throw myruntime_error("HDPscores: NewScores: Not enough memory.");
        }
    }
}

// -------------------------------------------------------------------------
// GetScores: get score table corresponding to the levels of probability and 
//  eff. no seqns
inline
Pslvector* HDPscores::GetScores( int plvl, int lvl )
{
    if( plvl < 0 || noplvs_ <= plvl )
        throw myruntime_error("HDPscores: GetScores: Memory access error.");
    if( lvl < 0 || notbls_ <= lvl )
        throw myruntime_error("HDPscores: GetScores: Memory access error.");
    return scores_[plvl][lvl];
}

// -------------------------------------------------------------------------
// GetScores: get score table corresponding to the levels of probability and 
//  eff. no seqns
inline
const Pslvector* HDPscores::GetScores( int plvl, int lvl ) const
{
    if( plvl < 0 || noplvs_ <= plvl )
        throw myruntime_error("HDPscores: GetScores: Memory access error.");
    if( lvl < 0 || notbls_ <= lvl )
        throw myruntime_error("HDPscores: GetScores: Memory access error.");
    return scores_[plvl][lvl];
}

// -------------------------------------------------------------------------
// GetScore: get score at table position identified by row, col; table is 
//  selected according to probability levels `fstprb' and `secprb' and eff. 
//  no. sequences given by `fstens' and `secens' 
inline
double HDPscores::GetScore( int row, int col, 
    double fstprb, double secprb, double fstens, double secens ) const
{
    int     e, ee, ln, pi, ppi, plvl;//level indices
    int     n, tmp;
    int     noplvs = SLC_MAX( 1, prblvs_.GetSize());
    int     nolvs = levels_.GetSize();
    double  dtmp;
    double  pi1, pi2, ppi1, ppi2;
    double  e1, e2, ee1, ee2;
    const Pslvector* scos = NULL;
    static char locbuf[KBYTE];

    if( row < 0 || col < 0 || card_ <= row || card_ <= col )
        throw myruntime_error("HDPscores: GetScore: Memory access error.");

    if( noplvs < 1 )
        throw myruntime_error("HDPscores: GetScore: No probability levels.");
    if( prblvs_.GetSize())
        if( fstprb < prblvs_.GetValueAt(0) || secprb < prblvs_.GetValueAt(0))
            return 0.0;
    if( secprb < fstprb ) {
        dtmp = secprb; secprb = fstprb; fstprb = dtmp; 
    }

    if( nolvs < 1 )
        throw myruntime_error("HDPscores: GetScore: No levels of eff. no. sqns.");
    if( fstens < levels_.GetValueAt(0) || secens < levels_.GetValueAt(0))
        return 0.0;
    if( secens < fstens ) {
        dtmp = secens; secens = fstens; fstens = dtmp; 
    }

    //get corresponding score table
    //index corresponding to probability level
    plvl = 0;
    if( prblvs_.GetSize()) {
        for( pi = 0; pi < noplvs; pi++ ) {
            pi1 = prblvs_.GetValueAt(pi);
            pi2 = 1.01;
            if( pi + 1 < noplvs )
                pi2 = prblvs_.GetValueAt(pi+1);
            if( pi1 <= fstprb && fstprb < pi2 )
                break;
        }
        for( ppi = pi; ppi < noplvs; ppi++ ) {
            ppi1 = prblvs_.GetValueAt(ppi);
            ppi2 = 1.01;
            if( ppi + 1 < noplvs )
                ppi2 = prblvs_.GetValueAt(ppi+1);
            if( ppi1 <= secprb && secprb < ppi2 )
                break;
        }
        if( noplvs <= pi || noplvs <= ppi || ppi < pi ) {
            sprintf( locbuf, "Segment covering probability level value not found: %g vs %g.", 
                     fstprb, secprb );
            throw myruntime_error("HDPscores: GetScore: " + mystring( locbuf ));
        }
        //index calculated as sum of arithm. series
        plvl = (( pi*( TIMES2(noplvs)-pi+1 ))>>1 ) + ppi-pi;
    }
    //index corresponding to level of eff. no. sequences
    for( e = 0; e < nolvs; e++ ) {
        e1 = levels_.GetValueAt(e);
        e2 = 20.01;
        if( e + 1 < nolvs )
            e2 = levels_.GetValueAt(e+1);
        if( e1 <= fstens && fstens < e2 )
            break;
    }
    for( ee = e; ee < nolvs; ee++ ) {
        ee1 = levels_.GetValueAt(ee);
        ee2 = 20.01;
        if( ee + 1 < nolvs )
            ee2 = levels_.GetValueAt(ee+1);
        if( ee1 <= secens && secens < ee2 )
            break;
    }
    if( nolvs <= e || nolvs <= ee || ee < e ) {
        sprintf( locbuf, "Segment covering eff. no. sequences not found: %g vs %g.", fstens, secens );
        throw myruntime_error("HDPscores: GetScore: " + mystring( locbuf ));
    }
    //index calculated as sum of arithm. series
    ln = (( e*( TIMES2(nolvs)-e+1 ))>>1 ) + ee-e;

    scos = GetScores( plvl, ln );
    if( scos == NULL )
        throw myruntime_error("HDPscores: GetScore: Null score table.");

    //calculate index within score table
    if( row < col ) {
        tmp = row; row = col; col = tmp; 
    }
    n = col;
    if( row )
        n += row + (( row * ( row-1 )) >> 1 );
    if( scos->GetSize() <= n )
        throw myruntime_error("HDPscores: GetScore: Memory access error.");
    return scos->GetValueAt(n);
}


#endif//__HDPscores__
