/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __psl_expint__
#define __psl_expint__


// The exponential integral E1, Re[ INT (exp(-x*t)/t) dt ]
int psl_expint_E1_e( const double x, double* result, double* err );

// Implementation of the exponential integral E1, allowing for scaling by exp(x)
int psl_expint_E1_scaled_e( const double x, double* result, double* err );


// The second-order exponential integral E2, Re[ INT (exp(-x*t)/t^2) dt ]
int psl_expint_E2_e( const double x, double* result, double* err );

// Implementation of the second-order exponential integral E2
int psl_expint_E2_scaled_e( const double x, double* result, double* err );


// The exponential integral of order n, En, Re[ INT (exp(-x*t)/t^n) dt ]
int psl_expint_En_e( const int n, const double x, double* result, double* err );

// Implementation of the exponential integral of order n, En
int psl_expint_En_scaled_e( const int n, const double x, double* result, double* err );


// The exponential integral Ei, -PV[ INT_-x (exp(-t)/t) dt ]
int psl_expint_Ei_e( const double x, double* result, double* err );

// Implementation of the exponential integral Ei, allowing for scaling by exp(x)
int psl_expint_Ei_scaled_e( const double x, double* result, double* err );



#endif//__psl_expint__
