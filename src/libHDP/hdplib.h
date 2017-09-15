/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __hdplib__
#define __hdplib__

#include <float.h>
#include "ext/psl.h"

#ifdef DBL_MIN
#   define LOC_DBL_MIN DBL_MIN
#else
#   define LOC_DBL_MIN SLC_DBL_MIN
#endif

#ifdef DBL_MAX
#   define LOC_DBL_MAX DBL_MAX
#else
#   define LOC_DBL_MAX SLC_DBL_MAX
#endif

//test HDP structure inline
#define HDPTESTINL

//include calculations of/with inverse scale matrix
// #define PROBMTXINVSM

#endif//__hdplib__
