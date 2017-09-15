/***************************************************************************
 *   Copyright (C) 2013 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __HDP_SSSScores__
#define __HDP_SSSScores__

#include <stdio.h>
#include <stdlib.h>

#include "ext/psl.h"
#include "ext/pslvector.h"
#include "rc.h"
#include "myexcept.h"


//global object declaration
class HDP_SSSScores;
extern HDP_SSSScores HDPSSSSCORES;

// -------------------------------------------------------------------------
// class HDP_SSSScores: Hierarchical Dirichlet Process model-derived scores
//
class HDP_SSSScores
{
public:
    HDP_SSSScores();
    ~HDP_SSSScores();

    double          GetHDPSSSWeight() const { return weight_; }
    void            SetHDPSSSWeight( double value ) { weight_ = value; }

    int             GetCardinality() const { return card_; }
    int             GetSSSCardinality() const { return ssscard_; }
    int             GetNoEffLvs() const { return noeffs_; }
    int             GetNoPrbLvs() const { return noplvs_; }
    int             GetNoSSSPrbLvs() const { return nosssplvs_; }

    double          GetScore( 
                        int row, int col, 
                        int sssrow, int ssscol, 
                        double fstprb, double secprb, double fstsssprb, double secsssprb, 
                        double fstens, double secens ) const;
    void            ReadScores( const char* filename );

    bool            ScoreNA( double value ) const { return value <= naval_; }

protected:
    void            ReadLevels( FILE* fp );
    void            ReadScoreTables( const char* filename );
    void            ReadScoresHelper( FILE* fp, int eff1, int eff2, int lvl );
    void            ReadSubmatrix( FILE* fp, int sss1, int sss2, Pslvector* scos );

    void            SetCardinality( int value ) { card_ = value; }
    void            SetSSSCardinality( int value ) { ssscard_ = value; }
    double          GetNAvalue() const { return naval_; }

    Pslvector*      GetScores( int lvl, int plvl, int sssplvl, int sssi );
    const Pslvector*GetScores( int lvl, int plvl, int sssplvl, int sssi ) const;
    void            NewScores( int nolevels, int noplvs, int nosssplvs, int nosss );
    void            DestroyScores();

private:
    Pslvector***** scores_;//scores
    double      weight_;//weight of all score tables
    int         noeffs_;//number of levels for effective no. seqs 
    int         noplvs_;//number of probability levels
    int         nosssplvs_;//number of SS state probability levels
    int         nosss_;//number of pair SS states
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
void HDP_SSSScores::DestroyScores()
{
    int n, p, ssp, sss;
    if( scores_ ) {
        for( n = 0; n < noeffs_; n++ ) {
            if( scores_[n] == NULL )
                continue;
            for( p = 0; p < noplvs_; p++ ) {
                if( scores_[n][p] == NULL )
                    continue;
                for( ssp = 0; ssp < nosssplvs_; ssp++ ) {
                    if( scores_[n][p][ssp] == NULL )
                        continue;
                    for( sss = 0; sss < nosss_; sss++ ) {
                        if( scores_[n][p][ssp][sss])
                            delete scores_[n][p][ssp][sss];
                        scores_[n][p][ssp][sss] = NULL;
                    }
                    free( scores_[n][p][ssp]);
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
    nosss_ = 0;
    nosssplvs_ = 0;
    noplvs_ = 0;
    noeffs_ = 0;
}
// NewScores: allocate new score table for each level of eff. no sequences
//  nolevels: no. levels for eff. no. sequences
//  noplvs: no. probability levels
//  nosssplvs: no. SS state probability levels
//  nosss: no. SS states
inline
void HDP_SSSScores::NewScores( int nolevels, int noplvs, int nosssplvs, int nosss )
{
    int n, p, ssp, sss;
    DestroyScores();
    if( nolevels < 1 || noplvs < 1 || nosssplvs < 1 || nosss < 1 )
        return;
    scores_ = ( Pslvector***** )malloc( nolevels * sizeof( Pslvector* ));
    if( scores_ == NULL )
        throw myruntime_error("HDP_SSSScores: NewScores: Not enough memory.");
    noeffs_ = nolevels;
    noplvs_ = noplvs;
    nosssplvs_ = nosssplvs;
    nosss_ = nosss;

    for( n = 0; n < noeffs_; n++ ) {
        scores_[n] = ( Pslvector**** )malloc( noplvs * sizeof( Pslvector* ));
        if( scores_[n] == NULL )
            throw myruntime_error("HDP_SSSScores: NewScores: Not enough memory.");
        for( p = 0; p < noplvs_; p++ ) {
            scores_[n][p] = ( Pslvector*** )malloc( nosssplvs * sizeof( Pslvector* ));
            if( scores_[n][p] == NULL )
                throw myruntime_error("HDP_SSSScores: NewScores: Not enough memory.");
            for( ssp = 0; ssp < nosssplvs_; ssp++ ) {
                scores_[n][p][ssp] = ( Pslvector** )malloc( nosss * sizeof( Pslvector* ));
                if( scores_[n][p][ssp] == NULL )
                    throw myruntime_error("HDP_SSSScores: NewScores: Not enough memory.");
                for( sss = 0; sss < nosss_; sss++ ) {
                    scores_[n][p][ssp][sss] = new Pslvector();
                    if( scores_[n][p][ssp][sss] == NULL )
                        throw myruntime_error("HDP_SSSScores: NewScores: Not enough memory.");
                }
            }
        }
    }
}

// -------------------------------------------------------------------------
// GetScores: get score table corresponding to the levels of eff. no seqns, 
//  probabilities and pair SS state index
inline
Pslvector* HDP_SSSScores::GetScores( int lvl, int plvl, int sssplvl, int sssi )
{
    if( lvl < 0 || noeffs_ <= lvl )
        throw myruntime_error("HDP_SSSScores: GetScores: Memory access error.");
    if( plvl < 0 || noplvs_ <= plvl )
        throw myruntime_error("HDP_SSSScores: GetScores: Memory access error.");
    if( sssplvl < 0 || nosssplvs_ <= sssplvl )
        throw myruntime_error("HDP_SSSScores: GetScores: Memory access error.");
    if( sssi < 0 || nosss_ <= sssi )
        throw myruntime_error("HDP_SSSScores: GetScores: Memory access error.");
    return scores_[lvl][plvl][sssplvl][sssi];
}

// -------------------------------------------------------------------------
// GetScores: get score table corresponding to the levels of eff. no seqns, 
//  probabilities and pair SS state index
inline
const Pslvector* HDP_SSSScores::GetScores( int lvl, int plvl, int sssplvl, int sssi ) const
{
    if( lvl < 0 || noeffs_ <= lvl )
        throw myruntime_error("HDP_SSSScores: GetScores: Memory access error.");
    if( plvl < 0 || noplvs_ <= plvl )
        throw myruntime_error("HDP_SSSScores: GetScores: Memory access error.");
    if( sssplvl < 0 || nosssplvs_ <= sssplvl )
        throw myruntime_error("HDP_SSSScores: GetScores: Memory access error.");
    if( sssi < 0 || nosss_ <= sssi )
        throw myruntime_error("HDP_SSSScores: GetScores: Memory access error.");
    return scores_[lvl][plvl][sssplvl][sssi];
}

// -------------------------------------------------------------------------
// GetScore: get score at table position identified by HDP clusters and SS 
//  states; table is selected according to probability levels `fstprb' and 
//  `secprb', `fstsssprb' and `secsssprb', and eff. no. sequences given by 
//  `fstens' and `secens' 
inline
double HDP_SSSScores::GetScore( 
    int row, int col, 
    int sssrow, int ssscol, 
    double fstprb, double secprb, double fstsssprb, double secsssprb, 
    double fstens, double secens ) const
{
    int     e, ee, ln, pi, ppi, plvl, sssplvl, sssi;//level indices
    int     n, tmp;
    int     nolvs = levels_.GetSize();
    int     noplvs = SLC_MAX( 1, prblvs_.GetSize());
    int     nosssplvs = SLC_MAX( 1, sssprblvs_.GetSize());
    int     nosss = nosss_;
    double  dtmp;
    double  pi1, pi2, ppi1, ppi2;
    double  e1, e2, ee1, ee2;
    const Pslvector* scos = NULL;
    static char locbuf[KBYTE];

    if( row < 0 || col < 0 || card_ <= row || card_ <= col )
        throw myruntime_error("HDP_SSSScores: GetScore: Memory access error.");

    if( sssrow < 0 || ssscol < 0 || ssscard_ <= sssrow || ssscard_ <= ssscol )
        throw myruntime_error("HDP_SSSScores: GetScore: Memory access error.");

    if( nolvs < 1 )
        throw myruntime_error("HDP_SSSScores: GetScore: No levels for eff. no. sqns.");
    if( fstens < levels_.GetValueAt(0) || secens < levels_.GetValueAt(0))
        return 0.0;
    if( secens < fstens ) {
        dtmp = secens; secens = fstens; fstens = dtmp; 
    }

    if( noplvs < 1 )
        throw myruntime_error("HDP_SSSScores: GetScore: No probability levels.");
    if( prblvs_.GetSize())
        if( fstprb < prblvs_.GetValueAt(0) || secprb < prblvs_.GetValueAt(0))
            return 0.0;
    if( secprb < fstprb ) {
        dtmp = secprb; secprb = fstprb; fstprb = dtmp; 
    }

    if( nosssplvs < 1 )
        throw myruntime_error("HDP_SSSScores: GetScore: No SS state probability levels.");
    if( sssprblvs_.GetSize())
        if( fstsssprb < sssprblvs_.GetValueAt(0) || secsssprb < sssprblvs_.GetValueAt(0))
            return 0.0;
    if( secsssprb < fstsssprb ) {
        dtmp = secsssprb; secsssprb = fstsssprb; fstsssprb = dtmp; 
    }

    //get corresponding score table
    //
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
        throw myruntime_error("HDP_SSSScores: GetScore: " + mystring( locbuf ));
    }
    //index calculated as sum of arithm. series
    ln = (( e*( TIMES2(nolvs)-e+1 ))>>1 ) + ee-e;

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
            throw myruntime_error("HDP_SSSScores: GetScore: " + mystring( locbuf ));
        }
        //index calculated as sum of arithm. series
        plvl = (( pi*( TIMES2(noplvs)-pi+1 ))>>1 ) + ppi-pi;
    }

    //index corresponding to SS state probability level
    sssplvl = 0;
    if( sssprblvs_.GetSize()) {
        for( pi = 0; pi < nosssplvs; pi++ ) {
            pi1 = sssprblvs_.GetValueAt(pi);
            pi2 = 1.01;
            if( pi + 1 < nosssplvs )
                pi2 = sssprblvs_.GetValueAt(pi+1);
            if( pi1 <= fstsssprb && fstsssprb < pi2 )
                break;
        }
        for( ppi = pi; ppi < nosssplvs; ppi++ ) {
            ppi1 = sssprblvs_.GetValueAt(ppi);
            ppi2 = 1.01;
            if( ppi + 1 < nosssplvs )
                ppi2 = sssprblvs_.GetValueAt(ppi+1);
            if( ppi1 <= secsssprb && secsssprb < ppi2 )
                break;
        }
        if( nosssplvs <= pi || nosssplvs <= ppi || ppi < pi ) {
            sprintf( locbuf, "Segment covering SS state probability level value not found: %g vs %g.", 
                     fstsssprb, secsssprb );
            throw myruntime_error("HDP_SSSScores: GetScore: " + mystring( locbuf ));
        }
        //index calculated as sum of arithm. series
        sssplvl = (( pi*( TIMES2(nosssplvs)-pi+1 ))>>1 ) + ppi-pi;
    }

    //index of submatrix determined by SS state values
    if( sssrow < ssscol ) {
        tmp = sssrow; sssrow = ssscol; ssscol = tmp; 
    }
    sssi = ssscol;
    if( sssrow )
        sssi += sssrow + (( sssrow * ( sssrow-1 )) >> 1 );
    if( nosss <= sssi ) {
        sprintf( locbuf, "Invalid pair SS state index: %g vs %g.", fstsssprb, secsssprb );
        throw myruntime_error("HDP_SSSScores: GetScore: " + mystring( locbuf ));
    }


    //score table
    scos = GetScores( ln, plvl, sssplvl, sssi );
    if( scos == NULL )
        throw myruntime_error("HDP_SSSScores: GetScore: Null score table.");

    //calculate index within score table
    if( sssrow == ssscol ) {
        if( row < col ) {
            tmp = row; row = col; col = tmp; 
        }
        n = col;
        if( row )
            n += row + (( row * ( row-1 )) >> 1 );
    }
    else {
        n = col;
        if( row )
            n += card_ * row;
    }
    if( scos->GetSize() <= n )
        throw myruntime_error("HDP_SSSScores: GetScore: Memory access error.");
    return scos->GetValueAt(n);
}


#endif//__HDP_SSSScores__
