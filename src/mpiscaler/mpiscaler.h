/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __mpiscaler_h__
#define __mpiscaler_h__

#define DEFAULT_INFORMATION_CONTENT     (   0.0 )
#define DEFAULT_PRECISION               ( AbstractScoreMatrix::NoScaling )
#define DEFAULT_MASKING                 ( Unmasked )
#define DEFAULT_SCORING_METHOD          ( AbstractScoreMatrix::ParallelUniversal )

// Version history:
//


static const char*  version = "1.01";
static const char*  verdate = "";



static const char*  instructions = "\n\
<>\n\
\n\
Parallel scaler of profile score matrix. The program computes statistical\n\
parameter lambda and optionally scales matrix to attain specific value of\n\
lambda.\n\
IMPORTANT: The program changes configuration (var/comer.par).\n\
NOTE: The program runs via MPI; be aware of having MPICH2 installed.\n\
(C)2008 Mindaugas Margelevicius,VU Institute of Biotechnology,Vilnius\n\
\n\
Usage:\n\
<> -d <database> [options]\n\
\n\
Parameters (with default values in parenthesis):\n\
\n\
-d <database>   [Filename]  Name of database made by makedb.\n\
-o <output>     [Filename]  File to write computed values to.            (stdout)\n\
\n\
Masking options:\n\
-I <info>       [Real]      Ignore positions with information content       (0.0)\n\
                            less than specified value.\n\
\n\
Scaling options:\n\
-l <lambda>     [(0.0-1.0)] Scale matrix to have specified value of        (Opt.)\n\
                            statistical parameter lambda.\n\
                            If not used, the program solves for statistical\n\
                            parameters and exits.\n\
-s (pro|LSO|HDP)            Scoring method to use.                          (pro)\n\
\n\
Startup:\n\
-h                          This text.\n\
\n\
";

#endif//__mpiscaler_h__
