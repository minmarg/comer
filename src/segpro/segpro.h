/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __segpro_h__
#define __segpro_h__

// Version history:
//


static const char*  version = "1.01";
static const char*  verdate = "";

static const char*  instructions = "\n\
<>\n\
\n\
SEG filter for profiles or multiple alignments given in FASTA or STOCKHOLM.\n\
(C)2008 Mindaugas Margelevicius,VU Institute of Biotechnology,Vilnius\n\
\n\
\n\
Usage:\n\
<> -i <input> -o <output> [options]\n\
\n\
Parameters (with default values in parenthesis):\n\
\n\
-i <input>      [Filename]  Either multiple alignment file in fasta or\n\
                            profile made by makepro.\n\
-o <output>     [Filename]  Output file.\n\
-p <file>       [Filename]  Profile in human-readable format.              (Opt.)\n\
-s <file>       [Filename]  Filtered profile sequence in fasta.            (Opt.)\n\
\n\
Output options:\n\
-L <width>      [Integer]   Line width for fasta sequences in output.        (60)\n\
                            -1 -- Do not wrap sequences.\n\
\n\
Profile construction options:\n\
-t <identity>   [1-100]     Ignore sequences in alignment file with          (94)\n\
                            this or higher level of sequence identity.\n\
\n\
SEG options:\n\
-w <window>     [Integer]   Window length.                                   (12)\n\
-f <low>        [Real]      Low entropy threshold.                          (2.2)\n\
-F <high>       [Real]      High entropy threshold.                         (2.5)\n\
-D <distance>   [Real]      Distance of equivalence between profile       (12.96)\n\
                            vectors.\n\
\n\
Startup:\n\
-v                          Enable warnings.\n\
-h                          This text.\n\
\n\
";

#endif//__segpro_h__
