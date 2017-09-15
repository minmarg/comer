/***************************************************************************
 *   Copyright (C) 2013 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __DBProfileProbs__
#define __DBProfileProbs__

#include <stdio.h>
#include <stdlib.h>

#include "ext/psl.h"
#include "ext/pslvector.h"
#include "ext/ivector.h"
#include "rc.h"
#include "myexcept.h"


// -------------------------------------------------------------------------
// class DBProfileProbs: Probabilities for database profiles 
//
class DBProfileProbs
{
public:
    DBProfileProbs();
    ~DBProfileProbs();

    double          GetProbAt( int n ) const;
    double          GetProbAt( int row, int col ) const;
    double          GetProbAt( int e1, int l1, int e2, int l2 ) const;
    double          GetProb( double E1, int L1, double E2, int L2 ) const;
    bool            ProbsReady() const { return probs_ && probs_->GetSize(); }

    const Pslvector&GetEffLvs() const { return efflvs_; }
    const Ivector&  GetLenLvs() const { return lenlvs_; }

    int             GetIndex( double E1, int L1 ) const;
    void            GetLevIndices( double E1, int L1, int* e1, int* l1 ) const;

    int             GetCardinality() const { return card_; }
    bool            Prob0( double value ) const { return value <= 0.; }

    void            ReadProbs( const char* filename );
    void            WriteProbs( const char* filename );

    void            ClearMargs();
    void            UpdateMargs( double E1, int L1 );
    void            CalcProbs();

    const Pslvector*GetProbs() const;

protected:
    void            ReadProbsHelper( FILE*, Pslvector* probs );
    void            WriteProbsHelper( FILE*, Pslvector* probs );
    void            SetCardinality( int value ) { card_ = value; }

    Pslvector*      GetProbs();
    void            NewProbs();
    void            DestroyProbs();
    void            VerifyProbs( Pslvector* probs, double acc = 1.e-6 );

protected:
    Pslvector*  probs_;//probabilities
    Ivector     margs_;//marginal pair counts
    Pslvector   efflvs_;//levels of eff. no. sequences
    Ivector     lenlvs_;//distinct values of profile lengths
    Pslvector   midefflvs_;//intermediate levels of eff. no. sequences
    Ivector     midlenlvs_;//intermediate values of profile lengths
    int         card_;//cardinality
};

// -------------------------------------------------------------------------
// INLINES
//
// ClearMargs: clear marginal counts
inline
void DBProfileProbs::ClearMargs()
{
    margs_.Clear();
    margs_.Reserve( GetCardinality());
}
// DestroyProbs: destroy probabilities
inline
void DBProfileProbs::DestroyProbs()
{
    if( probs_ ) {
        delete probs_;
        probs_ = NULL;
    }
}
// NewProbs: allocate new probabilities table
inline
void DBProfileProbs::NewProbs()
{
    DestroyProbs();
    if( efflvs_.GetSize() < 1 || lenlvs_.GetSize() < 1 )
        return;
    probs_ = new Pslvector();
    if( probs_ == NULL )
        throw myruntime_error("DBProfileProbs: NewProbs: Not enough memory.");
}

// -------------------------------------------------------------------------
// GetProbs: get probabilities table
inline
Pslvector* DBProfileProbs::GetProbs()
{
    return probs_;
}

// -------------------------------------------------------------------------
// GetProbs: get probabilities table
inline
const Pslvector* DBProfileProbs::GetProbs() const
{
    return probs_;
}

// -------------------------------------------------------------------------
// GetIndex: get nearest level indices for the given E1 and L1
inline
void DBProfileProbs::GetLevIndices( double E1, int L1, int* e1, int* l1 ) const
{
    if( e1 == NULL || l1 == NULL )
        throw myruntime_error("DBProfileProbs: GetLevIndices: Null parameters.");

    for( *e1 = 0; *e1 < midefflvs_.GetSize(); *e1++ )
        if( E1 < midefflvs_.GetValueAt(*e1))
            break;
    for( *l1 = 0; *l1 < midlenlvs_.GetSize(); *l1++ )
        if( L1 < midlenlvs_.GetValueAt(*l1))
            break;
}

// -------------------------------------------------------------------------
// GetIndex: get index (matrix row/col) given E1 and L1
inline
int DBProfileProbs::GetIndex( double E1, int L1 ) const
{
    int     e1, l1;//level indices
    int     row;

    for( e1 = 0; e1 < midefflvs_.GetSize(); e1++ )
        if( E1 < midefflvs_.GetValueAt(e1))
            break;
    for( l1 = 0; l1 < midlenlvs_.GetSize(); l1++ )
        if( L1 < midlenlvs_.GetValueAt(l1))
            break;

    row = e1 * lenlvs_.GetSize() + l1;
    return row;
}

// -------------------------------------------------------------------------
// GetProbAt: get probability at the specified location 
inline
double DBProfileProbs::GetProbAt( int n ) const
{
    if( probs_ == NULL || n < 0 || probs_->GetSize() <= n )
        throw myruntime_error("DBProfileProbs: GetProbAt: Memory access error.");
    return probs_->GetValueAt(n);
}

// -------------------------------------------------------------------------
// GetProb: get probability at table position given by row and col
inline
double DBProfileProbs::GetProbAt( int row, int col ) const
{
    int nn;

    if( probs_ == NULL || probs_->GetSize() < 1 )
        return 0.;

    if( row < col )
        throw myruntime_error("DBProfileProbs: GetProbAt: Invalid indices.");
    if( row < 0 )
        throw myruntime_error("DBProfileProbs: GetProbAt: Negative indices.");

    nn = col;
    if( row )
        nn += row + (( row * ( row-1 )) >> 1 );
    if( probs_->GetSize() <= nn )
        throw myruntime_error("DBProfileProbs: GetProbAt: Memory access error.");
    return probs_->GetValueAt(nn);
}

// -------------------------------------------------------------------------
// GetProbAt: get probability at table position identified by level 
//  indices e1,l1,e2,l2; 
inline
double DBProfileProbs::GetProbAt( int e1, int l1, int e2, int l2 ) const
{
    int row, col;
    row = e1 * lenlvs_.GetSize() + l1;
    col = e2 * lenlvs_.GetSize() + l2;
    return GetProbAt( row, col );
}

// -------------------------------------------------------------------------
// GetProb: get probability at table position identified by E1,L1,E2,L2; 
inline
double DBProfileProbs::GetProb( double E1, int L1, double E2, int L2 ) const
{
    int     row, col;
    int     tmp;
//     static char locbuf[KBYTE];

    if( probs_ == NULL || probs_->GetSize() < 1 )
        return 0.;

    if( efflvs_.GetSize() < 1 || lenlvs_.GetSize() < 1 )
        throw myruntime_error("DBProfileProbs: GetProb: No probabilities data.");

    row = GetIndex( E1, L1 );
    col = GetIndex( E2, L2 );
    if( row < col ) {
        tmp = col; col = row; row = tmp;
    }

    return GetProbAt( row, col );
}


#endif//__DBProfileProbs__
