/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __psl_gamma__
#define __psl_gamma__

// * The maximum x such that gamma(x) is not considered an overflow
#define SLC_G_GAMMA_XMAX        ( 171.0 )

//* The maximum n such that psl_fact(n) does not give an overflow
#define SLC_G_FACT_NMAX         ( 170 )

//* The maximum n such that psl_doublefact(n) does not give an overflow
#define SLC_G_DOUBLEFACT_NMAX   ( 297 )


// Log(Gamma(x)); for x<0 the real part of the result is returned
int psl_lngamma_e( double x, double* result, double* err );

// Log(Gamma(x)); value of Gamma can be reconstructed by sgn*exp(result)
int psl_lngamma_sgn_e( double x, double* result, double* err, double* sgn );

// Gamma(x); the maximum value of x is SLC_G_GAMMA_XMAX
int psl_gamma_e( const double x, double* result, double* err );

// Regulated Gamma function, Gamma*(x) for x>0
int psl_gammastar_e( const double x, double* result, double* err );

// Reciprocal of the Gamma function, 1/Gamma(x)
int psl_gammainv_e( const double x, double* result, double* err );


// Factorial, n! = Gamma(n+1); the maximum value of n is SLC_G_FACT_NMAX
int psl_fact_e( const unsigned int n, double* result, double* err );

// Double factorial, n!!; the maximum value of n is SLC_G_DOUBLEFACT_NMAX
int psl_doublefact_e( const unsigned int n, double* result, double* err );

// Log of factorial, log(n!) = log(Gamma(n+1))
int psl_lnfact_e( const unsigned int n, double* result, double* err );

// Log of double factorial, log(n!!)
int psl_lndoublefact_e( const unsigned int n, double* result, double* err );


// Combinatorial factor n choose m, n!/(m!-(n-m)!)
int psl_choose_e( unsigned int n, unsigned int m, double* result, double* err );

// Log of n choose m, log(n!) - log(m!) - log((n-m)!)
int psl_lnchoose_e( unsigned int n, unsigned int m, double* result, double* err );

// Taylor coefficient, x^n/n! for x>=0, n>=0
int psl_taylorcoeff_e( const int n, const double x, double* result, double* err );


#endif//__psl_gamma__
