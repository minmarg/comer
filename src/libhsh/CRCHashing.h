/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __CRCHashing__
#define __CRCHashing__

#include "types.h"


void build_table();
Uint32 crchash( const void *key, size_t length );

extern const Uint32     gencrctab[256];



#endif//__CRCHashing__

