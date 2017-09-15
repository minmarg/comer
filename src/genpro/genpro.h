/***************************************************************************
 *   Copyright (C) 2014 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __genpro_h__
#define __genpro_h__

#define DEFAULT_NO_PROFILES ( 1 )
#define DEFAULT_FRAGLEN ( 11 )

// Version history:
//

static const char*  version = "1.01";
static const char*  verdate = "";

static const char*  makeinst = "\n\
<>\n\
\n\
Generate random profiles using a library of profiles with SS predictions.\n\
(C)2014 Mindaugas Margelevicius,Institute of Biotechnology,Vilnius\n\
\n\
\n\
Usage:\n\
<> -o <dirname> -n <pattern> ( -d <dirname> | <profile1> <profile2> ... ) [options]\n\
\n\
Parameters (with default values in parenthesis):\n\
\n\
-o <dirname>    [Dirname]   Name of directory to write generated profiles in.\n\
-n <name_pat>   [Filename]  Output filename pattern of generated profiles.\n\
-d <dirname>    [Dirname]   Name of directory to read profiles from.\n\
\n\
Profile construction options:\n\
-l <length>     [Integer]   Generate profiles of about this length.\n\
-N <count>      [Integer]   Number of profiles to generate.                   (1)\n\
-f <length>     [Integer]   Length of indivisible profile fragment.          (11)\n\
-s                          Operate with secondary structure elements.\n\
\n\
Startup:\n\
-v                          Enable warnings.\n\
-h                          This text.\n\
\n\
Examples:\n\
<> -o mydir -n prolen200 -l 200 -N 1000 -d ./my_profiles\n\
<> -o mydir -n prolen200 -l 200 -N 10 d2zfga_.pro d1xioa_.pro\n\
<> -o mydir -n prolen200 -l 200 -N 1000 *.pro\n\
\n\
";

#endif//__genpro_h__
