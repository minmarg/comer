/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __psl_beta__
#define __psl_beta__


// Log(Beta(x,y));
int psl_lnbeta_e( const double x, const double y, double* result, double* err );

// Log(Beta(x,y)); value of Beta can be reconstructed by sgn*exp(result)
int psl_lnbeta_sgn_e( const double x, const double y, 
                      double* result, double* err, double* sgn );

// Beta(x,y);
int psl_beta_e( const double x, const double y, double* result, double* err );


#endif//__psl_beta__
