/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __psl_exp__
#define __psl_exp__


// exp(x)
int psl_exp_e( const double x, double* result, double* err );

// exp(x) with add. error estimate
int psl_exp_err_e( const double x, const double dx, double* result, double* err );


// y exp(x)
int psl_exp_mult_e( const double x, const double y, double* result, double* err );

// y exp(x) with add. error estimate
int psl_exp_mult_err_e( const double x, const double dx, const double y, const double dy, double* result, double* err );


// exp(x)-1; accurate for small x
int psl_expm1_e( const double x, double* result, double* err );


// (exp(x)-1)/x; accurate for small x
int psl_exprel_e( const double x, double* result, double* err );

// (exp(x)-1)/x; using continued fraction representation
int psl_exprel_n_CF_e( const double N, const double x, double* result, double* err );

// 2(exp(x)-1-x)/x^2; accurate for small x
int psl_exprel_2_e( double x, double* result, double* err );

// N-relative exponential, n-th generalization of  exprel and exprel_2
int psl_exprel_n_e( const int N, const double x, double* result, double* err );



#endif//__psl_exp__
