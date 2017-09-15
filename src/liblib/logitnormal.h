/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __logitnormal_h__
#define __logitnormal_h__

// interface
void LogitNormalErrCorr( double* f, size_t fsz, double tol );
void LogitNormal2Normal( double* f, size_t fsz, double tol, bool correct = true );
void Normal2LogitNormal( const double* nv, size_t nvsz, double* lnv, size_t lnvsz );

#endif//__logitnormal_h__
