/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __platform_h__
#define __platform_h__




#ifdef __GNUC__

#   define DIRSEP       '/'
#   define DIRSEPSTR    "/"
#   define NEWLINE      "\n"

#elif defined( OS_NT )

#   define DIRSEP       '\\'
#   define DIRSEPSTR    "\\"
#   define NEWLINE      "\n\r"

#else

#   error "OTHER SYSTEMS THAN   GNUC (LINUX/UNIX) and WINDOWS   ARE NOT SUPPORTED!"

#endif


#define UPDIR ".."



#endif//__platform_h__
