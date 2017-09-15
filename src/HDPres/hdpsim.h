/***************************************************************************
 *   Copyright (C) 2012 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __hdpsim_h__
#define __hdpsim_h__

#define DEFAULT_DIM ( 19 )
#define DEFAULT_NUMBER_OF_VARS ( 3 )
#define DEFAULT_CTXT_LENGTH ( 1 )
#define DEFAULT_NUMBER_OF_CLUSTERS ( 200 )
#define DEFAULT_NUMBER_OF_SAMPLES ( 200 )


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
Simulate data (multivariate-normals) for HDP-clustering.\n\
(C)2012 Mindaugas Margelevicius,VU Institute of Biotechnology,Vilnius\n\
\n\
\n\
Usage:\n\
<> -o <output> [options]\n\
\n\
Parameters (with default values in parenthesis):\n\
\n\
-o <output>     [Filename]  Pattern of output files of grouped vectors;\n\
                            <output>.crt.grp -- correctly clustered vectors,\n\
                            <output>.mix.grp -- vectors mixed across groups.\n\
                            <output>.one.grp -- vectors in one cluster.\n\
\n\
Sampling options:\n\
-d <number>     [Integer]   Dimensions of vectors.                           (19)\n\
-a <number>     [Integer]   Number of variates determining one cluster.       (3)\n\
-c <length>     [1-100]     Context length.                                   (1)\n\
-k <number>     [Integer]   Number of clusters.                             (200)\n\
-n <number>     [Integer]   Number of vectors to sample in each cluster;    (200)\n\
                            (must be <= k).\n\
-p                          Sample from logistic-normal distribution.\n\
-l                          Do not overlap positions of support variates.\n\
\n\
Startup:\n\
-v                          Enable warnings.\n\
-h                          This text.\n\
\n\
Examples:\n\
<> -o clusters.out -o my_sim\n\
\n\
";

#endif//__hdpsim_h__
