/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __AlignmentSimulation__
#define __AlignmentSimulation__

#include "debug.h"
#include "rc.h"
#include "data.h"
#include "InputMultipleAlignment.h"

// _________________________________________________________________________
// CLASS AlignmentSimulation
//
class AlignmentSimulation
{
public:
    // typedefs ...
    //
    AlignmentSimulation(
                const char* output, int noalns, int length, int thickness,
                double gap_freq, bool allowgaps, long seed = -1 );
    ~AlignmentSimulation();

    long        Run();

    const char* GetPattern() const      { return pattern; }
    int         GetNoAlns() const       { return no_alns; }
    int         GetAlnLength() const    { return aln_length; }
    int         GetAlnThickness() const { return aln_thickness; }
    bool        GetAllowGaps() const    { return allow_gaps; }

protected:
                                            //initialize interval of probabilities
    void            InitProbs( double* probs, double gap_freq = 0.0 );
    void            RegenerateAlignment();  //generation of alignments
    void            GenerateVector( PosDescriptionVector*, bool first );
                                            //generation of a residue code
    unsigned char   GenerateResidue( bool gaps_allowed = true );

private:
    const char*     pattern;            //name of file to write results to
    int             no_alns;            //number of different alignments to generate
    int             aln_length;         //length of ailgnment to be produced
    int             aln_thickness;      //thickness of alignments
    bool            allow_gaps;         //whether to allow gaps in the first sequence of the alignment

    double          wogprobs[NUMALPH+1];//interval of probabilities when no gap frequency is given
    double          resprobs[NUMALPH+1];//interval of probabilities given gap frequency

    const double    gap_frequency;      //gap frequency: how often (in percents) gap occurs among other residues
    long            rng_seed;           //seed for random number generator

    InputMultipleAlignment* alignment;  //input alignmemt information
};



////////////////////////////////////////////////////////////////////////////
// Class AlignmentSimulation inlines
//
// -------------------------------------------------------------------------


#ifdef __DEBUG__
#endif


#endif//__AlignmentSimulation__
