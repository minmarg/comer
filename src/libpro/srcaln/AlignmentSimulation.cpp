/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <stdlib.h>
#include <string.h>

#include "mystring.h"
#include "myexcept.h"

#include "AlignmentSimulation.h"



////////////////////////////////////////////////////////////////////////////
// CLASS AlignmentSimulation
//
// Constructor
//
AlignmentSimulation::AlignmentSimulation(
        const char* output,
        int noalns, int length, int thickness,
        double gap_freq, bool allowgaps,
        long seed )
:
    pattern( output ),
    no_alns( noalns ),
    aln_length( length ),
    aln_thickness( thickness ),
    allow_gaps( allowgaps ),
    gap_frequency( gap_freq ),
    rng_seed( seed )
{
    InitProbs( wogprobs );
    InitProbs( resprobs, gap_frequency );

    if( rng_seed == -1 )
        rng_seed = -rand();
    if( 0 < rng_seed )
        rng_seed = -rng_seed;

    alignment = new InputMultipleAlignment();
}

// Destructor
//
AlignmentSimulation::~AlignmentSimulation()
{
    if( alignment )
        delete alignment;
}

// -------------------------------------------------------------------------
// InitProbs: initialize interval of probabilities used to generate residues
// -------------------------------------------------------------------------

void AlignmentSimulation::InitProbs( double* probs, double gap_freq )
{
    if( !probs )
        return;

    double  res_freq = 1.0 - gap_freq;  //frequency for all other residues
    double  recip_fr = 1.0 / res_freq;  //inverse of res_freq;
    int     n;

    probs[0] = 0.0;

    for( n = 1; n < NUMALPH + 1 && probs[n-1] < ( res_freq - 0.0001 ); n++ )
        probs[n] = ( probs[n-1] * recip_fr + LOSCORES.PROBABility( n-1 )) * res_freq;

    for( ; n < NUMALPH + 1; n++ )
        probs[n] = 1.0;
}

// -------------------------------------------------------------------------
// Run: start generating multiple sequence alignments
// -------------------------------------------------------------------------

long AlignmentSimulation::Run()
{
    if( !alignment || !GetPattern())
        throw myruntime_error( mystring( "Unable to start generating alignments." ));

    const int   size = KBYTE;
    char        name[size];
    int         num = 0;
    int         len = strlen( GetPattern());
    const char* ext = GetPattern() + len - 1;

    if( size < len + 12 )
        throw myruntime_error( mystring( "Filename pattern given is too long." ));

    for( ; ext != GetPattern() && *ext != DIRSEP && *ext != '.'; ext-- );
    if( ext == GetPattern() || *ext == DIRSEP )
        ext =  GetPattern() + len;

    for( int n = 0; n < GetNoAlns(); n++ ) {
        RegenerateAlignment();

        size_t  base = ext - GetPattern();
        strncpy( name, GetPattern(), base );
        sprintf( name + base, "_%d%s", num++, ext );

        alignment->PutAlignment( name );
    }

    return rng_seed;
}

// -------------------------------------------------------------------------
// RegenerateAlignment: regenerates multiple sequence alignment
// -------------------------------------------------------------------------

void AlignmentSimulation::RegenerateAlignment()
{
    alignment->clear();
    alignment->SetTitle( "Query" );

    for( int n = 0; n < GetAlnThickness(); n++ ) {
        PosDescriptionVector*   one = alignment->NewPositionVector();
        GenerateVector( one, n == 0 );
        alignment->push( one );
    }
}

// -------------------------------------------------------------------------
// GenerateVector: generate residue sequence
// -------------------------------------------------------------------------

void AlignmentSimulation::GenerateVector( PosDescriptionVector* one, bool first )
{
    if( !one )
        throw myruntime_error( mystring( "Unable to generate residue sequence." ));

    bool    gaps_allowed = GetAllowGaps()? true: !first;

    for( int n = 0; n < GetAlnLength(); n++ ) {
        one->push( GenerateResidue( gaps_allowed ));
    }
}

// -------------------------------------------------------------------------
// GenerateResidue: generate single residue
// -------------------------------------------------------------------------

unsigned char AlignmentSimulation::GenerateResidue( bool gaps_allowed )
{
    int             n;
    double          prob = LecuyerRand( &rng_seed );
    double*         prob_intervals = gaps_allowed? resprobs: wogprobs;

    for( n = 1; n < NUMALPH + 1; n++ )
        if( prob_intervals[n-1] <= prob && prob < prob_intervals[n] )
            break;

    if( GAP - 5 < n - 1 )
        //if a number generated does not correspond to any non-abstract residue assign it to gap
        return GAP;

    return n - 1;
}
