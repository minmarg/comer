/***************************************************************************
 *   Copyright (C) 2009 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __CtxtCoefficients__
#define __CtxtCoefficients__

#include <math.h>
#include "compdef.h"
#include "debug.h"
#include "rc.h"
#include "data.h"

#include "mystring.h"
#include "myexcept.h"


// _________________________________________________________________________
// CLASS CtxtCoefficients
//
class CtxtCoefficients
{
public:
    CtxtCoefficients( size_t size );
    CtxtCoefficients( size_t size, double weight );

    ~CtxtCoefficients();

    size_t          GetLength() const { return length_; }

    double          GetCentralWeight() const { return cweight_; }
    void            SetCentralWeight( double value ) { cweight_ = value; }
    void            PutCentralWeight();

    double          GetMultiplier() const { return multip_; }

    const double*   GetCoefficients() const { return coeffs_; }
    double          GetCoefficientAt( size_t n ) const;
    void            SetCoefficientAt( size_t n, double value );

    const double*   GetLogCoefficients() const { return logcoeffs_; }
    double          GetLogCoefficientAt( size_t n ) const;
    void            SetLogCoefficientAt( size_t n );

    void            FindCoefficients();

//     bool            Read( FILE* );
    void            Write( FILE* ) const;
    void            WriteLogs( FILE* ) const;

protected:
    explicit CtxtCoefficients();

    static void     fdfunction( double x, double* f, double* df, void* );
    void            SetMultiplier( double value ) { multip_ = value; }
    double          ExpandCoefficients();

private:
    void            Init( size_t size );
    void            Destroy();
    void            SetLength( size_t value ) { length_ = value; }

private:
    size_t  length_;        //length of coefficient vector
    double  cweight_;       //central weight
    double  multip_;        //multiplier to be found
    double* coeffs_;        //coefficients
    double* logcoeffs_;     //precomputed logs of coefficients
};


////////////////////////////////////////////////////////////////////////////
// Class CtxtCoefficients inlines
//
// -------------------------------------------------------------------------
// GetCoefficientAt: get coefficient at the position
//
inline
double CtxtCoefficients::GetCoefficientAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !coeffs_ || GetLength() <= n )
        throw myruntime_error( mystring( "CtxtCoefficients: Memory access error." ));
#endif
    return coeffs_[n];
}

// SetCoefficientAt: set coefficient value at the position
//
inline
void CtxtCoefficients::SetCoefficientAt( size_t n, double value )
{
#ifdef __DEBUG__
    if( !coeffs_ || GetLength() <= n )
        throw myruntime_error( mystring( "CtxtCoefficients: Memory access error." ));
#endif
    coeffs_[n] = value;
}

// PutCoefficientAt: put coefficient value at the central position
//
inline
void CtxtCoefficients::PutCentralWeight()
{
    const double    value = GetCentralWeight();
    const double    precd = 0.0001;
    const size_t    cpost = GetLength() / 2;

    if( value <= 0.0 + precd || 1.0 - precd <= value )
        throw myruntime_error( mystring( "CtxtCoefficients: Invalid central coefficient." ));
    if( GetLength() < 1 || !( GetLength() % 2 ))
        throw myruntime_error( mystring( "CtxtCoefficients: Number of coefficients should be odd." ));

    coeffs_[cpost] = value;
}

// -------------------------------------------------------------------------
// GetLogCoefficientAt: get log of coefficient at the position
//
inline
double CtxtCoefficients::GetLogCoefficientAt( size_t n ) const
{
#ifdef __DEBUG__
    if( !logcoeffs_ || GetLength() <= n )
        throw myruntime_error( mystring( "CtxtCoefficients: Memory access error." ));
#endif
    return logcoeffs_[n];
}

// SetLogCoefficientAt: set log of coefficient at the position
//
inline
void CtxtCoefficients::SetLogCoefficientAt( size_t n )
{
#ifdef __DEBUG__
    if( !logcoeffs_ || GetLength() <= n )
        throw myruntime_error( mystring( "CtxtCoefficients: Memory access error." ));
#endif
    double  value = GetCoefficientAt( n );
    if( 0.0 < value )
        logcoeffs_[n] = log( value );
    else if( 0.0 == value )
        logcoeffs_[n] = -9999.;
    else
        throw myruntime_error( mystring( "CtxtCoefficients: Log of non-positive number." ));
}

#endif//__CtxtCoefficients__
