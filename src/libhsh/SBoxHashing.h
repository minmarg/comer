/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __SBoxHashing__
#define __SBoxHashing__

#include "types.h"

Uint32 sboxhash( const void *key, size_t length );

extern const Uint32     sbox[256];


#endif//__SBoxHashing__

