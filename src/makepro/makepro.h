/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __makepro_h__
#define __makepro_h__

// Version history:
// 1.03, . profile posterior probabilities included;
//         enriched SS format


static const char*  version = "1.04";
static const char*  verdate = "";

static const char*  makeinst = "\n\
<>\n\
\n\
Make profile from multiple sequence alignment in FASTA or STOCKHOLM.\n\
(C)2015 Mindaugas Margelevicius,VU Institute of Biotechnology,Vilnius\n\
\n\
\n\
Usage:\n\
<> -i <input> -o <output> [-f <outfile>] [-p <options>]\n\
\n\
Parameters:\n\
\n\
-i <input>      [Filename]  Input multiple alignment file in fasta.\n\
-o <output>     [Filename]  Output profile.\n\
-f <file>       [Filename]  Output profile in human-readable format.\n\
-p <options>    [Filename]  Input file of options;\n\
                            By default, the file in installation\n\
                            directory of this package is searched.\n\
\n\
-v                          Enable warnings.\n\
-h                          This text.\n\
\n\
";

#endif//__makepro_h__
