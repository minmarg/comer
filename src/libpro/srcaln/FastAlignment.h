/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __FastAlignment__
#define __FastAlignment__

#include "rc.h"
#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "mystring.h"
#include "myexcept.h"
#include "libpro/srcsco/AbstractScoreMatrix.h"


// _________________________________________________________________________
// Class FastAlignment
//
// Fast ungapped alignments and final pseudo score production
//

class FastAlignment
{
    enum ACluster {     //clusters in the dynamic programming matrix
        aDValues,       //accumulated scores
        aMaxvals,       //maximum values in the upper-left area from the current position
        szACluster
    };

    enum AIndices {     //Constants for indices
        aRow,           //row index
        aColumn,        //column index
        aMarker,        //alignment break marker
        aLength,        //alignment length
        szIndices,
        aFirst = aRow,
        aSecnd = aColumn,
        szDims = 2
    };

public:
	FastAlignment( const AbstractScoreMatrix*  scsystem, const char* queryres = NULL, const char* sbjctres = NULL );
	~FastAlignment();

    void                        Run();

    int                         GetQuerySize() const            { return querySize; }
    int                         GetSubjectSize() const          { return subjectSize; }

    int                         GetAlnSteps() const             { return alnsteps; }
    double                      GetScore() const                { return alnscore; }

    void                        Print( char* sp );                  //print alignment with statistical significance
    void                        Print( FILE* fp );                  //print alignment with statistical significance
                                                                    //universal print method
    void                        Print( TPrintFunction, void* vpn );
    void                        PrintDPMatrix( FILE* fp ) const;

protected:
    void                        Align();
    void                        MakeAlignmentPath();

    void                        SetAlnSteps( int steps )        { alnsteps = steps; }
    void                        SetScore( double value )        { alnscore = value; }

    const AbstractScoreMatrix*  GetScoreMatrix() const;

    const char*                 GetQueryResidues() const        { return queryresidues; }
    const char*                 GetSubjectResidues() const      { return sbjctresidues; }

    void                        PrintIsolatedSegment(
                                    TPrintFunction  print_func,
                                    void*           vpn,
                                    mystring&       query,
                                    mystring&       sbjct,
                                    mystring&       match,
                                    int             path_step_from,
                                    int             path_step_to
                                );

private:
    double      ( **F )[szACluster];        //augmented dynamic programming matrix
    int         ( **IND )[szIndices];       //matrix of indices of the maximum values

    int         querySize;                  //length of query
    int         subjectSize;                //length of subject

    const AbstractScoreMatrix*  scoreSystem;//score system

    const char* queryresidues;              //query residues
    const char* sbjctresidues;              //subject residues

    int         ( *PATH )[ szDims ];        //alignment path containing indices of two profiles in each position

    double      alnscore;                   //alignment score
    int         alnsteps;                   //length of slignment
};


////////////////////////////////////////////////////////////////////////////
// INLINES
//
// -------------------------------------------------------------------------
// GetScoreMatrix: obtains one of the two possible scoring matrices
// -------------------------------------------------------------------------

inline
const AbstractScoreMatrix* FastAlignment::GetScoreMatrix() const
{
    return scoreSystem;
}


#endif//__FastAlignment__
