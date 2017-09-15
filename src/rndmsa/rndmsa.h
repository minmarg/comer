/***************************************************************************
 *   Copyright (C) 2014 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __rndmsa_h__
#define __rndmsa_h__

#define DEFAULT_MSA_LEN ( 200 )
#define DEFAULT_NO_MSAs ( 1 )
#define DEFAULT_MSA_WDT ( 1 )
#define DEFAULT_MSA_NOISE ( 0.1 )

// Version history:
//

static const char*  version = "1.01";
static const char*  verdate = "";

static const char*  makeinst = "\n\
<>\n\
\n\
Generate random multiple sequence alignments.\n\
(C)2014 Mindaugas Margelevicius,Institute of Biotechnology,Vilnius\n\
\n\
\n\
Usage:\n\
<> -o <dirname> -n <pattern> [options]\n\
\n\
Parameters (with default values in parenthesis):\n\
\n\
-o <dirname>    [Dirname]   Name of directory to write generated alignments in.\n\
-n <name_pat>   [Filename]  Output filename pattern of generated alignments.\n\
\n\
General options:\n\
-l <length>     [Integer]   Length of the first sequence in alignment.      (200)\n\
-N <count>      [Integer]   Number of alignments to generate.                 (1)\n\
\n\
Construction options:\n\
-t <probabs>    [List]      Comma-sep. list of transition probabilities    (Opt.)\n\
                            (MM,MI,MD,IM,II,ID,DM,DI,DD).\n\
-p <probabs>    [List]      Comma-sep. list of residue probabilities       (Opt.)\n\
                            (A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V).\n\
-m <profile>    [Filename]  Sample sequences according to the model of     (Opt.)\n\
                            profile (mutex to -t and -p).\n\
-r <noise>      [0.0-1.0]   Distort probabilities at this level of noise.   (0.1)\n\
\n\
Number of sequences:\n\
-e <number>     [Integer]   Desired eff. number of sequences [1-20] of alignment.\n\
-b <number>     [Integer]   Num. of sequences within alignment (mutex to -e). (1)\n\
\n\
Startup:\n\
-v                          Enable warnings.\n\
-h                          This text.\n\
\n\
Examples:\n\
<> -o mydir -n msalen200 -l 200 -N 1000 -b 100\n\
<> -o mydir -n msalen200 -l 200 -N 10   -e 12\n\
<> -o mydir -n msalen200 -l 200 -N 1000 -e 12 -t 0.97,0.01,0.02,0.31,0.69,0,0.28,0,0.72\n\
\n\
";

#endif//__rndmsa_h__
