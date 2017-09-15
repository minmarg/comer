/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __ProfileShuffler__
#define __ProfileShuffler__

#include "debug.h"
#include "rc.h"
#include "data.h"
#include "mystring.h"
#include "myexcept.h"
#include "libpro/srcpro/Database.h"
#include "libpro/srcpro/GapScheme.h"
#include "libpro/srcpro/DistributionMatrix.h"
#include "InputMultipleAlignment.h"


// _________________________________________________________________________
// CLASS ProfileShuffler
//
class PositionProbabilities
{
public:
    PositionProbabilities();
    ~PositionProbabilities();

    void    Allocate( size_t size );        //allocate space for the probability and index entries
    size_t  GetIndex( double probability ); //get an index of position vector given value of probability distribution function
    void    Push( double prob );            //push values given the probability of the positional vector and its index
    bool    Check();                        //check if probability distribution function is valid

protected:
    size_t      GetNoEntries() const    { return no_entries; }

    double      GetDistFunctionAt( size_t pos ) const;
    void        SetDistFunctionAt( size_t pos, double value );

    size_t      GetCounter() const      { return counter; }
    void        IncCounter()            { counter++; }

    void        Destroy();                  //destroy entries

private:
    double*     distFunction;   //probability distribution function
    size_t      no_entries;     //number of entries
    size_t      counter;        //counter for entries
};


// _________________________________________________________________________
// CLASS ProfileShuffler
//
class ProfileShuffler
{
public:
    // typedefs ...
    //
    ProfileShuffler(
                const char* database,
                const char* output,
                int noprofs,
                int length,
                int thickness,
                bool textformat,
                long seed = -1
    );
    ~ProfileShuffler();

    long        Run();

    const char* GetDatabase() const         { return database_name; }
    const char* GetPattern() const          { return pattern; }
    int         GetNoProfiles() const       { return no_profiles; }
    int         GetProfileLength() const    { return profile_length; }
    int         GetMinThickness() const     { return min_thickness; }
    bool        GetPrintText() const        { return print_text; }

protected:
    explicit ProfileShuffler();
                                                                            //simulate profiles
    void            SimulateProfile( FrequencyMatrix&, LogOddsMatrix&, GapScheme& );

    void            InitProbs( double* probs, double gap_freq = 0.0 );      //initialize interval of probabilities
    unsigned char   GenerateResidue();                                      //generate single residue

    void            FillPositionProbabilities();                            //fill probability distribution function

    size_t          GetNoSequences() const      { return no_db_sequences; } //number of sequences in database
    Uint64          GetDbSize() const           { return db_size; }         //size of database

    void            SetNoSequences( size_t value )  { no_db_sequences = value; }
    void            SetDbSize( Uint64 value )       { db_size = value; }

    const double* const GetResidueProbs() const     { return residueprobs; }
    const char*     GetResidues() const         { return residues; }
    char*           GetResidues()               { return residues; }

    long            GetCurrentSeedValue() const { return rng_seed; }

private:
    const char*     database_name;      //profile database name
    const char*     pattern;            //name of file to write results to
    int             no_profiles;        //number of different alignments to generate
    int             profile_length;     //length of ailgnment to be produced
    int             min_thickness;      //thickness of alignments
    bool            print_text;         //print profiles in text in addition
                                        //interval of probabilities when no gap frequency is given
    double          residueprobs[NUMALPH+1];
    char*           residues;           //amino acid sequence of length profile_length

    long            rng_seed;           //seed for random number generator

    Database        profile_db;         //profile database
    size_t          no_db_sequences;    //number of sequences within the database
    Uint64          db_size;            //size of the database expressed in number of positions

    static const char*      s_text_ext; //file extension for profiles in the text format
    PositionProbabilities   positionProbs;  //probabilities of position vectors
};



////////////////////////////////////////////////////////////////////////////
// Class PositionProbabilities inlines
//
inline
double PositionProbabilities::GetDistFunctionAt( size_t pos ) const
{
#ifdef __DEBUG__
    if( GetNoEntries() <= pos )
        throw myruntime_error( mystring( "PositionProbabilities: Memory access error." ));
#endif
    return distFunction[pos];
}

// -------------------------------------------------------------------------
//
inline
void PositionProbabilities::SetDistFunctionAt( size_t pos, double value )
{
#ifdef __DEBUG__
    if( GetNoEntries() <= pos )
        throw myruntime_error( mystring( "PositionProbabilities: Memory access error." ));
#endif
    distFunction[pos] = value;
}

////////////////////////////////////////////////////////////////////////////
// Class ProfileShuffler inlines
//
// -------------------------------------------------------------------------


#ifdef __DEBUG__
#endif


#endif//__ProfileShuffler__
