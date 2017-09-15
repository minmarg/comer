/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __psl_gammainc__
#define __psl_gammainc__


// Normalized incomplete Gamma Function Q (normalized upper inc. Gamma function)
int psl_gammainc_Q_e( const double a, const double x, double* result, double* err );

// Normalized incomplete Gamma Function P (normalized lower inc. Gamma function, complementary to Q: 1-Q)
int psl_gammainc_P_e( const double a, const double x, double* result, double* err );

// Unnormalized incomplete Gamma Function (unnormalized upper inc. Gamma function)
int psl_gammainc_e( const double a, const double x, double* result, double* err );


#endif//__psl_gammainc__
