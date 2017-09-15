/***************************************************************************
 *   Copyright (C) 2009 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __psl_cheb__
#define __psl_cheb__


// data type for a Chebyshev series over a given interval
//
typedef struct Tcheb_series_struct {
    double*     c;          // coefficients
    int         order;      // order of expansion
    double      a;          // lower interval point
    double      b;          // upper interval point
    int         order_sp;   // effective single precision order
} Tcheb_series;

int cheb_eval_e( const Tcheb_series* cs, const double x, double* result, double* err );

int psl_multiply_e( const double x, const double y, double* result, double* err );
int psl_multiply_err_e( const double x, const double dx, const double y, const double dy, double* result, double* err );


#endif//__psl_cheb__
