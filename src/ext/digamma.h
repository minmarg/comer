/***************************************************************************
 *   Copyright (C) 2009 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __psl_digamma__
#define __psl_digamma__


int psl_psi_int_e( const int n, double* result, double* err );
int psl_psi_e( const double x, double* result, double* err );

// Trigamma function, psi'(n), n>0
int psl_psi_1_int_e( const int n, double* result, double* err );

// Trigamma function, psi'(n) for general x
int psl_psi_1_e( const double x, double* result, double* err );

// Polygamma function, psi^(n)(x) for n>=0, x>0
int psl_psi_n_e( const int n, const double x, double* result, double* err );


#endif//__psl_digamma__
