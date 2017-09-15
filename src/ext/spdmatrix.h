/***************************************************************************
 *   Copyright (C) 2010 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __spdmatrix__
#define __spdmatrix__

#include "nsmatrix.h"

#ifndef SPDmatrixTESTPRINT
// #define SPDmatrixTESTPRINT
#endif

// -------------------------------------------------------------------------
// Symmetric, Positive Definite square matrix
//
class SPDmatrix: public NSmatrix
{
public:
    SPDmatrix( int nrc );
    SPDmatrix( const SPDmatrix& );
    virtual ~SPDmatrix();

    virtual SPDmatrix&  operator=( const SPDmatrix& );
    virtual SPDmatrix&  operator=( const Pslmatrix& );
    virtual void        priverror( const char* ) const;

    int     CholeskyDecompose();
    int     CDedSolve( Pslvector& xb ) const;
    int     CDedInvert();
    int     CDedDet( double* ) const;
    int     CDedLogDet( double* ) const;

protected:
    explicit    SPDmatrix();
};


// -------------------------------------------------------------------------
//

#endif//__spdmatrix__
