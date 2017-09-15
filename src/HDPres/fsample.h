/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __fsample_h__
#define __fsample_h__

#define DEFAULT_MIN_EFFECTIVE_THICKNESS ( 3 )
#define DEFAULT_NUMBER_OF_SAMPLES ( 1000000 )
#define DEFAULT_CTXT_LENGTH ( 1 )
#define DEFAULT_CTXT_STEP ( 1 )
#define DEFAULT_CTXT_CWGHT ( 1.0 )


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
Sample profile column frequencies.\n\
(C)2012 Mindaugas Margelevicius,VU Institute of Biotechnology,Vilnius\n\
\n\
\n\
Usage:\n\
<> -o <output> ( -d <directory> | <profile1> <profile2>... ) [options]\n\
\n\
Parameters (with default values in parenthesis):\n\
\n\
-d <directory>  [Dirname]   Directory to read profiles from.\n\
-o <output>     [Filename]  Output file of grouped frequencies.\n\
\n\
Profile options:\n\
-t <number>     [Integer]   Ignore profiles with effective number of          (3)\n\
                            sequences less than this number.\n\
\n\
Sampling options:\n\
-n <number>     [Integer]   Number of frequency vectors to sample.          (1e6)\n\
-c <length>     [1-100]     Context length.                                   (1)\n\
-p <step>       [1-10]      Step between adjacent context positions.          (1)\n\
-m                          Mix all context positions.\n\
-w <weight>     (0.0-1.0]   Weight of central context position.               (1)\n\
-a                          Arrange samples in clusters one per\n\
                            group (profile)\n\
\n\
Startup:\n\
-s <number>     [Natural]   Seed for random number generator               (Opt.)\n\
-v                          Enable warnings.\n\
-h                          This text.\n\
\n\
Examples:\n\
<> -o clusters.out -d ./my_profiles\n\
<> -o clusters.out 15925557.pro 15597784.pro 22538090.pro\n\
<> -o clusters.out *.pro -s 988102063\n\
\n\
";

#endif//__fsample_h__
