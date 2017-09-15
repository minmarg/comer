/***************************************************************************
 *   Copyright (C) 2013 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __CVS2Scores__
#define __CVS2Scores__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ext/psl.h"
#include "ext/pslvector.h"
#include "rc.h"
#include "data.h"
#include "myexcept.h"

//global object declaration
class CVS2Scores;
extern CVS2Scores CVS2SCORES;

// -------------------------------------------------------------------------
// class CVS2Scores: implements the translation of context normal vector 
//  scores to alignment scores
//
class CVS2Scores
{
public:
    CVS2Scores();
    virtual ~CVS2Scores();

    double          GetScore( double cvscore, double fstens, double secens ) const;
    void            ReadScores( const char* filename );

    double          GetCVSWeight() const { return cvsweight_; }
    void            SetCVSWeight( double value ) { cvsweight_ = value; }

    int             GetNoTables() const { return notbls_; }

    bool            ScoreNA( double value ) const { return value <= naval_; }
    bool            ScorePlus( double value ) const { return plusval_ <= value; }

protected:
    void            ReadLevels( FILE* fp );
    void            ReadScoreTables( FILE* fp );
    void            ReadScoresHelper( FILE* fp, Pslvector* keys, Pslvector* scos, 
                                      int* shift, int* card, double* scale );
    void            SmoothScores( Pslvector* keys, Pslvector* scos, int* shift, int* card );

    int             GetStep() const { return step_; }
    double          GetNAvalue() const { return naval_; }
    double          GetValuePlus() const { return plusval_; }

    Pslvector*      GetKeys( int lvl );
    const Pslvector*GetKeys( int lvl ) const;
    Pslvector*      GetScores( int lvl );
    const Pslvector*GetScores( int lvl ) const;
    int             GetShift( int lvl ) const;
    double          GetScale( int lvl ) const;
    void            NewScores( int nolevels );
    void            DestroyScores();

private:
    double  cvsweight_;//weight of translated vector scores
    Pslvector** keys_;//scores of normal vectors
    Pslvector** scores_;//translated scores
    int*        shifts_;//number of negative keys for each table
    Pslvector   levels_;//levels of eff. no. sequences
    int         notbls_;//number of score tables 
    int*        cards_;//cardinalities for each table
    double*     scale_;//scale factors for each table
    const int   step_;//number of intermediate values +1 between two adjacent integers
    const double naval_;//NA value
    const double plusval_;//value for '+'
};

// -------------------------------------------------------------------------
// INLINES
//
// DestroyScores: destroy all score tables
inline
void CVS2Scores::DestroyScores()
{
    int m;
    if( shifts_ ) {
        free( shifts_ );
        shifts_ = NULL;
    }
    if( cards_ ) {
        free( cards_ );
        cards_ = NULL;
    }
    if( scale_ ) {
        free( scale_ );
        scale_ = NULL;
    }
    if( keys_ ) {
        for( m = 0; m < notbls_; m++ ) {
            if( keys_[m] == NULL )
                continue;
            delete keys_[m];
            keys_[m] = NULL;
        }
        free( keys_ );
        keys_ = NULL;
    }
    if( scores_ ) {
        for( m = 0; m < notbls_; m++ ) {
            if( scores_[m] == NULL )
                continue;
            delete scores_[m];
            scores_[m] = NULL;
        }
        free( scores_ );
        scores_ = NULL;
    }
    notbls_ = 0;
}

// NewScores: allocate new score table for each level of eff. no sequences
//  nolevels: no. level values of eff. no. sequences
inline
void CVS2Scores::NewScores( int nolevels )
{
    int m;
    DestroyScores();
    if( nolevels < 1 )
        return;
    keys_ = ( Pslvector** )malloc( nolevels * sizeof( Pslvector* ));
    scores_ = ( Pslvector** )malloc( nolevels * sizeof( Pslvector* ));
    shifts_ = ( int* )malloc( nolevels * sizeof( int ));
    cards_ = ( int* )malloc( nolevels * sizeof( int ));
    scale_ = ( double* )malloc( nolevels * sizeof( double ));
    if( keys_ == NULL || scores_ == NULL || shifts_ == NULL || cards_ == NULL || scale_ == NULL )
        throw myruntime_error("CVS2Scores: NewScores: Not enough memory.");
    memset( shifts_, 0, nolevels * sizeof( int ));
    memset( cards_, 0, nolevels * sizeof( int ));
    notbls_ = nolevels;

    for( m = 0; m < notbls_; m++ )
        scale_[m] = 1.0;

    for( m = 0; m < notbls_; m++ ) {
        keys_[m] = new Pslvector();
        scores_[m] = new Pslvector();
        if( keys_[m] == NULL || scores_[m] == NULL )
            throw myruntime_error("CVS2Scores: NewScores: Not enough memory.");
    }
}

// -------------------------------------------------------------------------
// GetKeys: get keys corresponding to the level of eff. no seqns
inline
Pslvector* CVS2Scores::GetKeys( int lvl )
{
    if( lvl < 0 || notbls_ <= lvl )
        throw myruntime_error("CVS2Scores: GetKeys: Memory access error.");
    return keys_[lvl];
}

// GetKeys: get keys corresponding to the level of eff. no seqns
inline
const Pslvector* CVS2Scores::GetKeys( int lvl ) const
{
    if( lvl < 0 || notbls_ <= lvl )
        throw myruntime_error("CVS2Scores: GetKeys: Memory access error.");
    return keys_[lvl];
}

// -------------------------------------------------------------------------
// GetScores: get score table corresponding to the level of eff. no seqns
inline
Pslvector* CVS2Scores::GetScores( int lvl )
{
    if( lvl < 0 || notbls_ <= lvl )
        throw myruntime_error("CVS2Scores: GetScores: Memory access error.");
    return scores_[lvl];
}

// GetScores: get score table corresponding to the level of eff. no seqns
inline
const Pslvector* CVS2Scores::GetScores( int lvl ) const
{
    if( lvl < 0 || notbls_ <= lvl )
        throw myruntime_error("CVS2Scores: GetScores: Memory access error.");
    return scores_[lvl];
}

// -------------------------------------------------------------------------
// GetShift: get shift number for the table corresponding to the level of 
//  eff. no seqns
inline
int CVS2Scores::GetShift( int lvl ) const
{
    if( lvl < 0 || notbls_ <= lvl )
        throw myruntime_error("CVS2Scores: GetShift: Memory access error.");
    return shifts_[lvl];
}

// -------------------------------------------------------------------------
// GetScale: get scale factor for keys corresponding to the level of 
//  eff. no seqns
inline
double CVS2Scores::GetScale( int lvl ) const
{
    if( lvl < 0 || notbls_ <= lvl )
        throw myruntime_error("CVS2Scores: GetScale: Memory access error.");
    return scale_[lvl];
}

// -------------------------------------------------------------------------
// GetScore: get score at table position identified by log-odds score 
//  `cvscore' of normal vectors
inline
double CVS2Scores::GetScore( double cvscore, double fstens, double secens ) const
{
    int     e, ee, ln;//level indices
    int     n, tmp, shft, size;
    int     step = GetStep();
    int     nolvs = levels_.GetSize();
    double  dtmp, scl;
    double  e1, e2, ee1, ee2;
    const Pslvector* keys = NULL;
    const Pslvector* scos = NULL;
    static char locbuf[KBYTE];

    if( step < 1 )
        throw myruntime_error("CVS2Scores: GetScore: Invalid step.");
    if( nolvs < 1 )
        throw myruntime_error("CVS2Scores: GetScore: No levels of eff. no. sqns.");
    if( fstens < levels_.GetValueAt(0) || secens < levels_.GetValueAt(0))
        return 0.0;
    if( secens < fstens ) {
        dtmp = secens; secens = fstens; fstens = dtmp; 
    }

    //get corresponding score table index
    // corresponding to level of eff. no. sequences
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
        throw myruntime_error("CVS2Scores: GetScore: " + mystring( locbuf ));
    }
    //index calculated as sum of arithm. series
    ln = (( e*( TIMES2(nolvs)-e+1 ))>>1 ) + ee-e;

    keys = GetKeys(ln);
    scos = GetScores(ln);
    if( scos == NULL || keys == NULL )
        throw myruntime_error("CVS2Scores: GetScore: Null score table.");
    shft = GetShift(ln);
    if( shft < 0 )
        throw myruntime_error("CVS2Scores: GetScore: Negative shift.");
    scl = GetScale(ln);
    if( scl <= 0.0 )
        throw myruntime_error("CVS2Scores: GetScore: Invalid scale.");

    size = scos->GetSize();

    if( size < 1 )
        return 0.0;

    if( scl != 1.0 )
        cvscore *= scl;

    if( cvscore <= keys->GetValueAt(0))
        return scos->GetValueAt(0);
    if( keys->GetValueAt(size-1) <= cvscore )
        return scos->GetValueAt(size-1);

    //calculate index within score table
    if( step == 2 )
        n = shft + (int)rint(TIMES2(cvscore));
    else if( step == 1 )
        n = shft + (int)rint(cvscore);
    else
        n = shft + (int)rint((double)step*cvscore);

    if( size <= n || n < 0 )
        throw myruntime_error("CVS2Scores: GetScore: Memory access error.");
    return scos->GetValueAt(n);
}


#endif//__CVS2Scores__
