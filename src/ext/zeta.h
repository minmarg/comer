/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __psl_zeta__
#define __psl_zeta__


// The Hurwitz zeta function, SUM (k+q)^(-s)
int psl_hzeta_e( const double s, const double q, double* result, double* err );


// The Riemann zeta function, SUM_k (k)^(-s) for s<>1
int psl_zeta_e( const double s, double* result, double* err );

// The Riemann zeta function, SUM_k (k)^(-n) for integer n<>1
int psl_zeta_int_e( const int n, double* result, double* err );


// The Riemann zeta function minus 1 for s<>1;
int psl_zetam1_e( const double s, double* result, double* err );

// The Riemann zeta function minus 1 for integer n<>1;
int psl_zetam1_int_e( const int n, double* result, double* err );


// The eta function, (1-2^(1-n))zeta(n) for integer n 
int psl_eta_int_e( int n, double* result, double* err );

// The eta function, (1-2^(1-s))zeta(s)
int psl_eta_e( const double s, double* result, double* err );




#endif//__psl_zeta__
