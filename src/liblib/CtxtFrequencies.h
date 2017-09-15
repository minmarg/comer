/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __CtxtFrequencies__
#define __CtxtFrequencies__

#include "compdef.h"
#include "debug.h"
#include "rc.h"
#include "data.h"

#include "mystring.h"
#include "myexcept.h"


// _________________________________________________________________________
// CLASS CtxtFrequencies
//
class CtxtFrequencies
{
public:
    CtxtFrequencies();
    CtxtFrequencies( size_t size );

    ~CtxtFrequencies();

    size_t          GetLength() const               { return length_; }
    size_t          GetEffectiveNoResids() const    { return NUMAA; }

    char            GetResidue( size_t n ) const;
    void            SetResidue( size_t n, char value );
    char*           GetResidueAddrAt( size_t n );

    double          GetObsFrequency( size_t n, size_t res ) const;
    void            SetObsFrequency( size_t n, size_t res, double value );
    void            AddObsFrequency( size_t n, size_t res, double value );    //add value
    const double ( *GetObsFreqsAddrAt( size_t n ) const )[NUMALPH];
    double*         GetObsFreqsAddrAt( size_t n );

    double          GetMinValue() const             { return minval_; }
    double          GetMaxValue() const             { return maxval_; }

    void            FindMinMaxValues( double* min = NULL, double* max = NULL );

    virtual bool    Read( FILE* );
    virtual void    Write( FILE* ) const;
    void            WriteDoubles( FILE* ) const;
    void            WriteDoubles1( FILE*, bool title ) const;

protected:
    virtual void    clear();                                //clear all information
    virtual void    destroy();                              //deallocate memory
    virtual void    reallocate( size_t howmuch );           //memory allocation

    void            SetMinValue( double value )     { minval_ = value; }
    void            SetMaxValue( double value )     { maxval_ = value; }

    size_t          GetAllocated() const            { return allocated_; }
    void            IncLength()                     { length_++; }

private:
    void    SetAllocated( size_t poss ) { allocated_ = poss; }
    void    SetLength( size_t len )     { length_ = len; }

private:
    size_t  allocated_;                 //how many positions is allocated
    size_t  length_;                    //length of profile
    char*   sequence_;                  //amino acid sequence
    double  ( *obsfreqs_ )[NUMALPH];    //observed frequencies embedded in profile
    double  minval_;                    //min value
    double  maxval_;                    //max value
};

// _________________________________________________________________________
// CLASS CtxtProScores
//
class CtxtProScores: public CtxtFrequencies
{
public:
    CtxtProScores();
    ~CtxtProScores();

    virtual void    Read( const char* filename );
    virtual void    Write( FILE* ) const;

protected:
    void            ReadThickness( const char* readfrom, size_t readlen, int* membuf );
    size_t          ReadResidue( const char* readfrom, size_t readlen, char* membuf );
    size_t          ReadValues( const char* readfrom, size_t readlen, double* membuf, bool allownegs = true );

    virtual void    clear();                                //clear all information
    virtual void    destroy();                              //deallocate memory
    virtual void    reallocate( size_t howmuch );           //memory allocation

    int*    GetEffThicknessAddr()           { return &effthickness_; }
    int     GetEffThickness() const         { return effthickness_; }
    void    SetEffThickness( int value )    { effthickness_ = value; }

    double          GetScore( size_t n, size_t res ) const;
    void            SetScore( size_t n, size_t res, double value );
    const double ( *GetScoresAddrAt( size_t n ) const )[NUMALPH];
    double*         GetScoresAddrAt( size_t n );

private:
    int     effthickness_;              //effectivenes thickness of profile
    char*   sequence_;                  //amino acid sequence
    double  ( *scores_ )[NUMALPH];      //profile scores
};

////////////////////////////////////////////////////////////////////////////
// Class CtxtFrequencies inlines
//
// -------------------------------------------------------------------------
// GetResidue: get residue at the position
//
inline
char CtxtFrequencies::GetResidue( size_t n ) const
{
#ifdef __DEBUG__
    if( GetLength() <= n )
        throw myruntime_error( mystring( "CtxtFrequencies: Memory access error." ));
#endif
    return sequence_[n];
}

// SetResidue: set residue
//
inline
void CtxtFrequencies::SetResidue( size_t n, char value )
{
#ifdef __DEBUG__
    if( GetAllocated() <= n )
        throw myruntime_error( mystring( "CtxtFrequencies: Memory access error." ));
#endif
    if( GetLength() <= n )
        SetLength( n + 1 );
    sequence_[n] = value;
}

// GetResidueAddrAt: get address of residue at the position
//
inline
char* CtxtFrequencies::GetResidueAddrAt( size_t n )
{
#ifdef __DEBUG__
    if( GetAllocated() <= n )
        throw myruntime_error( mystring( "CtxtFrequencies: Memory access error." ));
#endif
    return sequence_ + n;
}

// -------------------------------------------------------------------------
// GetObsFrequency: get observed frequency for residue at the position
//
inline
double CtxtFrequencies::GetObsFrequency( size_t n, size_t res ) const
{
#ifdef __DEBUG__
    if( GetLength() <= n )
        throw myruntime_error( mystring( "CtxtFrequencies: Memory access error." ));

    if( NUMALPH <= res )
        throw myruntime_error( mystring( "CtxtFrequencies: Memory access error." ));
#endif
    return obsfreqs_[n][res];
}

// SetObsFrequency: sets observed frequency
//
inline
void CtxtFrequencies::SetObsFrequency( size_t n, size_t res, double value )
{
#ifdef __DEBUG__
    if( GetAllocated() <= n || NUMALPH <= res )
        throw myruntime_error( mystring( "CtxtFrequencies: Memory access error." ));
#endif
    if( GetLength() <= n )
        SetLength( n + 1 );
    obsfreqs_[n][res] = value;
}

// AddObsFrequency: add observed frequency to existing value
//
inline
void CtxtFrequencies::AddObsFrequency( size_t n, size_t res, double value )
{
#ifdef __DEBUG__
    if( GetAllocated() <= n || NUMALPH <= res )
        throw myruntime_error( mystring( "CtxtFrequencies: Memory access error." ));
#endif
    if( GetLength() <= n )
        SetLength( n + 1 );
    obsfreqs_[n][res] += value;
}

// GetObsFreqsAddrAt: get address of observed frequencies at the position
//
inline
const double ( *CtxtFrequencies::GetObsFreqsAddrAt( size_t n ) const )[NUMALPH]
{
#ifdef __DEBUG__
    if( GetLength() <= n )
        throw myruntime_error( mystring( "CtxtFrequencies: Memory access error." ));
#endif
    return obsfreqs_ + n;
}

// GetObsFreqsAddrAt: get address of observed frequencies at the position
//
inline
double* CtxtFrequencies::GetObsFreqsAddrAt( size_t n )
{
#ifdef __DEBUG__
    if( GetAllocated() <= n )
        throw myruntime_error( mystring( "CtxtFrequencies: Memory access error." ));
#endif
    return obsfreqs_[n];
}

////////////////////////////////////////////////////////////////////////////
// Class CtxtProScores inlines
//
// -------------------------------------------------------------------------
// GetScore: get score for residue at the position
//
inline
double CtxtProScores::GetScore( size_t n, size_t res ) const
{
#ifdef __DEBUG__
    if( GetLength() <= n )
        throw myruntime_error( mystring( "CtxtFrequencies: Memory access error." ));

    if( NUMALPH <= res )
        throw myruntime_error( mystring( "CtxtFrequencies: Memory access error." ));
#endif
    return scores_[n][res];
}

// SetScore: set score
//
inline
void CtxtProScores::SetScore( size_t n, size_t res, double value )
{
#ifdef __DEBUG__
    if( GetLength() <= n || NUMALPH <= res )
        throw myruntime_error( mystring( "CtxtFrequencies: Memory access error." ));
#endif
    scores_[n][res] = value;
}

// GetScoresAddrAt: get address of scores at the position
//
inline
const double ( *CtxtProScores::GetScoresAddrAt( size_t n ) const )[NUMALPH]
{
#ifdef __DEBUG__
    if( GetLength() <= n )
        throw myruntime_error( mystring( "CtxtFrequencies: Memory access error." ));
#endif
    return scores_ + n;
}

// GetScoresAddrAt: get address of scores at the position
//
inline
double* CtxtProScores::GetScoresAddrAt( size_t n )
{
#ifdef __DEBUG__
    if( GetAllocated() <= n )
        throw myruntime_error( mystring( "CtxtFrequencies: Memory access error." ));
#endif
    return scores_[n];
}


#endif//__CtxtFrequencies__
