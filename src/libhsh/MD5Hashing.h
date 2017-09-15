/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __MD5Hashing__
#define __MD5Hashing__

#include "types.h"

#define SZCIPHER 4
#define SZINPUT  16

void md5hashing( Uint32 cipher[SZCIPHER], const void *key, size_t length );



#endif//__MD5Hashing__

