/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __types_h__
#define __types_h__



#include <sys/types.h>


#ifdef __GNUC__

#   include <inttypes.h>
#   include <sys/param.h>

#   define Int32   int32_t
#   define Uint32  uint32_t

#   define Int64   int64_t
#   define Uint64  uint64_t


#elif defined( OS_NT )

#   include <stddef.h>

#   define Int32   __int32
#   define Uint32  unsigned __int32

#   define Int64   __int64
#   define Uint64  unsigned __int64


#else

#   error "OTHER SYSTEMS THAN   GNUC (LINUX/UNIX) and WINDOWS   ARE NOT SUPPORTED!"

#endif



#endif//__types_h__
