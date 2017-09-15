/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __psl_log__
#define __psl_log__


// log(x) for x>0
int psl_log_e( const double x, double* result, double* err );

// log(|x|) for x!=0
int psl_log_abs_e( const double x, double* result, double* err );

// log(1+x) for x>-1; accurate for small x
int psl_log_1plusx_e( const double x, double* result, double* err );

// log(1+x)-x for x>-1; accurate for small x
int psl_log_1plusx_mx_e( const double x, double* result, double* err );



#endif//__psl_log__
