/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __RJHashing__
#define __RJHashing__

#include "types.h"

#define hashsize(n) (( Uint32 )1 << ( n ))
#define hashmask(n) ( hashsize( n ) - 1 )


Uint32 hashlittle( const void *key, size_t length, Uint32 initval );


#endif//__RJHashing__

