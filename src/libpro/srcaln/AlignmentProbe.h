/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __AlignmentProbe__
#define __AlignmentProbe__

#include "rc.h"
#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "mystring.h"
#include "myexcept.h"


// _________________________________________________________________________
// Class AlignmentProbe
//
// Simple gapped alignment probe for statistical significance
//

class AlignmentProbe
{
    enum ACluster {
        aDValues,//accumulated scores
        szACluster
    };
    enum AIndices {
        aDir,
        szIndices
    };
    enum ADims {
        aFirst,
        aSecnd,
        szDims
    };

public:
    AlignmentProbe( const int querylen = 0, const char* queryres = NULL, 
                    const int sbjctlen = 0, const char* sbjctres = NULL
    );
    ~AlignmentProbe();

    void                        Run();

    int         GetQuerySize() const { return querySize_; }
    int         GetSubjectSize() const { return subjectSize_; }

    int         GetAlnSteps() const { return alnsteps_; }
    double      GetScore() const { return alnscore_; }
    int         GetIdentities() const { return identities_; }

    void        Print( char* sp );                  //print alignment
    void        Print( FILE* fp );                  //print alignment
    void        Print( TPrintFunction, void* vpn );//print method
    void        PrintDPMatrix( FILE* fp ) const;

protected:
    void        Align();
    void        MakeAlignmentPath();

    void        SetAlnSteps( int steps ) { alnsteps_ = steps; }
    void        SetScore( double value ) { alnscore_ = value; }

    const char* GetQueryResidues() const { return queryresidues_; }
    const char* GetSubjectResidues() const { return sbjctresidues_; }

private:
    double      (**F_)[szACluster];//augmented dynamic programming matrix
    int         (**IND_)[szIndices];//matrix of indices of maximum values

    int         querySize_;//length of query
    int         subjectSize_;//length of subject

    const char* queryresidues_;//query residues
    const char* sbjctresidues_;//subject residues

    int         (*PATH_)[szDims];//alignment path

    int         maxsrowind_;//row index of the maximum score
    int         maxscolind_;//row index of the maximum score
    double      alnscore_;//alignment score
    int         alnsteps_;//length of alignment
    int         identities_;//#identities found in alignment
};


////////////////////////////////////////////////////////////////////////////
// INLINES
//


#endif//__AlignmentProbe__
