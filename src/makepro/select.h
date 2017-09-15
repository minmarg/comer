/***************************************************************************
 *   Copyright (C) 2013 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __select_h__
#define __select_h__

// Version history:
//

static const char*  version = "1.10";
static const char*  verdate = "";

static const char*  selinst = "\n\
<>\n\
\n\
Select non-redundant set of sequences from multiple alignment.\n\
(C)2013 Mindaugas Margelevicius,Institute of Biotechnology,Vilnius\n\
\n\
\n\
Usage:\n\
<> -i <input> -o <output> [-p <options>] [-t]\n\
\n\
Parameters:\n\
\n\
-i <input>      [Filename]  Input multiple alignment file either in\n\
                            FASTA or in STOCKHOLM.\n\
-o <output>     [Filename]  Output file of multiple alignment.\n\
-p <options>    [Filename]  Input file of options;\n\
                            By default, the file in installation\n\
                            directory of this package is searched.\n\
-t                          Do not apply the selection of sequences;\n\
                            just convert the input.\n\
\n\
-v                          Enable warnings.\n\
-h                          This text.\n\
\n\
";

#endif//__select_h__
