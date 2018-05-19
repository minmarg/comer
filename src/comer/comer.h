/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __comer_h__
#define __comer_h__


#define DEFAULT_SCORING_SCHEME          ( AbstractScoreMatrix::ProfileSpecific )
#define DEFAULT_STATISTICAL_BEHAVIOUR   ( AbstractScoreMatrix::ComputeStatistics )
#define DEFAULT_PRECISION               ( AbstractScoreMatrix::AutoScalling )
#define DEFAULT_MASKING                 ( MaskToIgnore )
#define DEFAULT_ALNALGORITHM            ( ProfileAlignment::Optimizational )

// Version history:
//
// 1.03 . Improved scoring; improved SS processing


static const char*  version = "1.05";
static const char*  verdate = "";



static const char*  instructions = "\n\
<>\n\
\n\
Protein remote homology search tool.\n\
(C)2015 Mindaugas Margelevicius,VU Institute of Biotechnology,Vilnius\n\
\n\
Usage:\n\
<> -i <query> -d <database> [-o <output>] [-p <options>]\n\
\n\
Parameters:\n\
\n\
-i <query>      [Filename]  Either multiple alignment file in fasta or\n\
                            profile made by makepro.\n\
-d <database>   [Filename]  Either name of database made by makedb or\n\
                            another multiple alignment file.\n\
-o <output>     [Filename]  Output file of alignments.\n\
-p <options>    [Filename]  Input file of options;\n\
                            By default, the file in installation\n\
                            directory of this package is searched.\n\
\n\
-v                          Enable warnings.\n\
-h                          This text.\n\
\n\
";

#endif//__comer_h__
