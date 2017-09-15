/***************************************************************************
 *   Copyright (C) 2014 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __genmsa_h__
#define __genmsa_h__

#define DEFAULT_NO_MSAS ( 1 )
#define DEFAULT_FRAGLEN ( 11 )

// Version history:
//

static const char*  version = "1.01";
static const char*  verdate = "";

static const char*  makeinst = "\n\
<>\n\
\n\
Generate random multiple alignments using a library of source MSAs.\n\
(C)2014 Mindaugas Margelevicius,Institute of Biotechnology,Vilnius\n\
\n\
\n\
Usage:\n\
<> -o <dirname> -n <pattern> ( -d <dirname> | <msa1> <msa2> ... ) [options]\n\
\n\
Parameters (with default values in parenthesis):\n\
\n\
-o <dirname>    [Dirname]   Name of directory to write generated alignments in.\n\
-n <name_pat>   [Filename]  Output filename pattern of generated alignments.\n\
-d <dirname>    [Dirname]   Name of directory to read alignments in FASTA from.\n\
\n\
Profile construction options:\n\
-l <length>     [Integer]   Generate MSAs of about this length.\n\
-N <count>      [Integer]   Number of MSAs to generate.                       (1)\n\
-f <length>     [Integer]   Length of indivisible alignment fragment.        (11)\n\
\n\
Startup:\n\
-v                          Enable warnings.\n\
-h                          This text.\n\
\n\
Examples:\n\
<> -o mydir -n msalen200 -l 200 -N 1000 -d ./my_msas\n\
<> -o mydir -n msalen200 -l 200 -N 10 d2zfga_.fa d1xioa_.fa\n\
<> -o mydir -n msalen200 -l 200 -N 1000 *.fa\n\
\n\
";

#endif//__genmsa_h__
