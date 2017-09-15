/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __pscores_h__
#define __pscores_h__

// Version history:
//
// ----
//
// 1.01 . .


static const char*  version = "1.01";
static const char*  verdate = "";

static const char*  instructions = "\n\
<>\n\
\n\
Derive substitution scores from a set of observed frequency files.\n\
(C)2008 Mindaugas Margelevicius,VU Institute of Biotechnology,Vilnius\n\
\n\
\n\
Usage:\n\
<> -o <output> ( -d <directory> | <freqfile1> <freqfile2>... ) [options]\n\
\n\
Parameters (with default values in parenthesis):\n\
\n\
-d <directory>  [Dirname]   Directory to read frequencies from.\n\
-o <output>     [Filename]  Output file of scores.                       (stdout)\n\
\n\
Scale options:\n\
-n <scale>      [1-100]     Scale for scores (1/n bits).\n\
                            By default, entropy determines it.\n\
\n\
Startup:\n\
-v                          Enable warnings.\n\
-h                          This text.\n\
\n\
Examples:\n\
<> -o my_scores.txt -d ./my_frequencies\n\
<> -o my_scores.txt d_70_1_2.cnt b_119_1_1.cnt c_17_1_1.cnt c_69_1_5.cnt\n\
<> -o my_scores.txt *.cnt\n\
\n\
";

#endif//__pscores_h__
