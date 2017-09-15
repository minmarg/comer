/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __psl_erf__
#define __psl_erf__


// error function
int psl_erf_e( double x, double* result, double* err );


// Complementary error function, 1-erf(x)
int psl_erfc_e( double x, double* result, double* err  );

// log of the complementary error function, log(erfc(x))
int psl_log_erfc_e( double x, double* result, double* err );


// Standard Gaussian probability density function
int psl_erf_Z_e( double x, double* result, double* err );

// The complement of standard Gaussian probability function
int psl_erf_Q_e( double x, double* result, double* err );

// The hazard function of the standard normal distribution (inverse Mill's ratio)
int psl_hazard_e( double x, double* result, double* err );




#endif//__psl_erf__
