/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __hsort_h__
#define __hsort_h__

#include <stdlib.h>

//definition of function types
typedef int ( *TSCompFunction )( const void*, size_t, size_t );
typedef int ( *TSSwapFunction )( void*, size_t, size_t );

//sorting of data with heap sort
void HeapSort( void* vector, size_t size, TSCompFunction cmfunc, TSSwapFunction swfunc );

//sorting of data with heap sort placing indices in array `inds' of size `size'
void HeapSortInd( size_t* inds, const void* vector, size_t size, TSCompFunction cmfunc );


#endif//__hsort_h__

