/***************************************************************************
 *   Copyright (C) 2013 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __SSSScores__
#define __SSSScores__

#include <stdio.h>
#include <stdlib.h>

#include "ext/psl.h"
#include "ext/pslvector.h"
#include "rc.h"
#include "data.h"
#include "myexcept.h"
#include "VirtScores.h"

//global object declaration
class SSSScores;
extern SSSScores SSSSCORES;

// -------------------------------------------------------------------------
// class SSSScores: Secondary structure state scores
//
class SSSScores: public VirtScores
{
public:
    SSSScores();
    virtual ~SSSScores();

    double          GetScore( char fstsss, char secsss, 
                              double fstprb, double secprb, double fstens, double secens ) const;
    void            ReadScores( const char* filename );

    double          GetSSSWeight() const { return sssweight_; }
    void            SetSSSWeight( double value ) { sssweight_ = value; }

protected:
    void            ReadLevels( FILE* fp );
    void            ReadScoreTables( FILE* fp );

private:
    double  sssweight_;//weight of SS state scores
};

// -------------------------------------------------------------------------
// INLINES
//
// -------------------------------------------------------------------------
// GetScore: get score at table position identified by SS states 
//  `fstsss' and `secsss' and probabilities `fstprb' and `secprb'
inline
double SSSScores::GetScore( char fstsss, char secsss, 
    double fstprb, double secprb, double fstens, double secens ) const
{
    int row = ( int )fstsss * 10 + int( fstprb * 10.0 );
    int col = ( int )secsss * 10 + int( secprb * 10.0 );
    return VirtScores::GetScore( row, col, fstprb, secprb, fstens, secens );
}


#endif//__SSSScores__
