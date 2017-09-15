/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __makedb_h__
#define __makedb_h__

#define DEFAULT_DISTRIBUTION_TYPE   ( DISCRETE )

// Version history:
//
// 1.03 . inclusion of profile posterior probabilities


static const char*  version = "1.04";
static const char*  verdate = "";

static const char*  makeinst = "\n\
<>\n\
\n\
Make profile database.\n\
(C)2015 Mindaugas Margelevicius,VU Institute of Biotechnology,Vilnius\n\
\n\
\n\
Usage:\n\
<> -o <output> ( -d <directory> | <profile1> <profile2> ... ) [options]\n\
\n\
Parameters (with default values in parenthesis):\n\
\n\
-o <output>     [Filename]  Name of database to be created.\n\
-d <directory>  [Dirname]   Name of directory to read profiles from.\n\
\n\
Database construction options:\n\
-t <distrib>    (simple |   Type of distribution according to which      (simple)\n\
                 multin |   vector probabilities are computed:\n\
                 profile)    simple,  simple discrete distribution,\n\
                             multin,  multinomial distribution,\n\
                             profile, distribution of profile vectors.\n\
\n\
SEG options:\n\
-U                          Invoke low-complexity filtering of profiles.\n\
-w <window>     [Integer]   Window length.                                  ( 12)\n\
-f <low>        [Real]      Low entropy threshold.                          (2.2)\n\
-F <high>       [Real]      High entropy threshold.                         (2.5)\n\
-D <distance>   [Real]      Distance of equivalence between profile       (12.96)\n\
                            vectors.\n\
\n\
Startup:\n\
-v                          Enable warnings.\n\
-h                          This text.\n\
\n\
Examples:\n\
<> -o my_db -t uniform -d ./my_profiles\n\
<> -o my_db -t profile d_70_1_2.pro b_119_1_1.pro c_17_1_1.pro c_69_1_5.pro\n\
<> -o my_db *.pro\n\
\n\
";

#endif//__makedb_h__
