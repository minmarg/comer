/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __segdata__
#define __segdata__

#include <stdio.h>

#define MAX_LOGGAMMA    10000

//STATIC DATA

static const size_t s_SZPRSQR =  101;
static const size_t s_SZPRLOG =  200;
static const size_t s_SZPRENT =  200;
static const size_t s_SZPRLGM = 1000;

// _________________________________________________________________________
// Abstract structure of precomputed values
//
struct _TPRECOMPUTED
{
    _TPRECOMPUTED( size_t );
    virtual ~_TPRECOMPUTED();

    void    Init();
    double  GetValueOf( size_t v ) const
    {
        if( v < GetSize())
            return data[v];
        return Compute( v );
    }

    size_t  GetSize() const  { return no_vals; }

    virtual double  Compute( size_t v ) const = 0;

private:
    const size_t    no_vals;
    double*         data;
};

// _________________________________________________________________________
// Precomputed squared values
//
struct _TSQUARES: public _TPRECOMPUTED
{
    _TSQUARES( size_t = s_SZPRSQR );
    virtual ~_TSQUARES();

    virtual double  Compute( size_t v ) const;
};


// _________________________________________________________________________
// Precomputed logarithm values
//
struct _TLOGARITHMS: public _TPRECOMPUTED
{
    _TLOGARITHMS( size_t = s_SZPRLOG );
    virtual ~_TLOGARITHMS();

    virtual double  Compute( size_t v ) const;
};


// _________________________________________________________________________
// Precomputed entropy values
//
struct _TPARTIAL_ENTROPIES: public _TPRECOMPUTED
{
    _TPARTIAL_ENTROPIES( size_t = s_SZPRENT );
    virtual ~_TPARTIAL_ENTROPIES();

    virtual double  Compute( size_t v ) const;
};


// _________________________________________________________________________
// Precomputed log-gamma function values
//
struct _TLOG_GAMMA
{
    _TLOG_GAMMA();
    ~_TLOG_GAMMA();

    double  GetValueOf( size_t v );
    size_t  GetSize() const  { return no_vals; }

private:
    void    Precompute( size_t newsize );
    void    SetSize( size_t newsize )   { no_vals = newsize; }

    size_t          no_vals;
    double*         data;
};

// =========================================================================
// Declarations of extern constant variables
//

extern const struct _TSQUARES               PRESQUARES;
extern const struct _TLOGARITHMS            LOGARITHMS;
extern const struct _TPARTIAL_ENTROPIES     PRT_ENTROPIES;
extern struct _TLOG_GAMMA                   LOG_GAMMA;


#endif//__segdata__
