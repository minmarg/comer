/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __pcluster_h__
#define __pcluster_h__

#define DEFAULT_SEQIDENTITY_PERCENTAGE ( 0.94 )
#define DEFAULT_POSGAPIGNORE_PERCENTAGE ( 1.00 )


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
Clustering of sequences within multiple sequence alignment in FASTA.\n\
(C)2008 Mindaugas Margelevicius,VU Institute of Biotechnology,Vilnius\n\
\n\
\n\
Usage:\n\
<> -i <input> -o <output> [options]\n\
\n\
Parameters (with default values in parenthesis):\n\
\n\
-i <input>      [Filename]  Input multiple alignment (MA) file in fasta.\n\
-o <output>     [Filename]  Output file of observed frequencies.         (stdout)\n\
\n\
Clustering options:\n\
-c <percentage> [1-100]     Cluster sequences at this level of               (62)\n\
                            sequence identity.\n\
\n\
Processing options:\n\
-t <percentage> [0-100]     Filter MA at this level of sequence identity.    (94)\n\
-g                          Ignore MA positions with gaps in query (1st).\n\
-p <percentage> [0-100]     Ignore MA positions in which percentage of      (100)\n\
                            gaps is greater than specified value.\n\
\n\
Filtering options:\n\
-U                          Invoke high-complexity (column) filtering.\n\
-w <window>     [Integer]   Window length.                                   (30)\n\
-f <low>        [Real]      Low entropy threshold.                          (3.3)\n\
-F <high>       [Real]      High entropy threshold.                         (3.5)\n\
\n\
Startup:\n\
-v                          Enable warnings.\n\
-h                          This text.\n\
\n\
";

#endif//__pcluster_h__
