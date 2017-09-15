/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __root_h__
#define __root_h__

// definition of a function type
typedef int ( *TRFunction )( double, double*, double*, void* );

// root finding by the Newton-Raphson method
int rootByNRandBisection(
        TRFunction fdfunction,
        double x1,
        double x2,
        double* root,
        double xacc,
        int maxit,
        void* params );


#endif//__root_h__

