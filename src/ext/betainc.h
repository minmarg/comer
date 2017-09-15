/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __psl_betainc__
#define __psl_betainc__


// Regularized incomplete Beta function 
int psl_betainc_e( const double a, const double b, const double x, 
                   double* result, double* err );


#endif//__psl_betainc__
