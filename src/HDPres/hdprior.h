/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __hdprior_h__
#define __hdprior_h__

#define DEFAULT_UNINF_SCALEFAC  ( 1.0 )
#define DEFAULT_DEGF_ADJUSTMENT ( -1.0 )

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
Calculate prior parameters under Hierarchical Dirichlet Process.\n\
(C)2012 Mindaugas Margelevicius,VU Institute of Biotechnology,Vilnius\n\
\n\
\n\
Usage:\n\
<> -i <input> -o <output> [options]\n\
\n\
Parameters (with default values in parenthesis):\n\
\n\
-i <input>      [Filename]  Input file of grouped frequencies.\n\
-o <output>     [Filename]  Name of Output file.\n\
\n\
Clustering options:\n\
-n <number>     [Real]      Adjustment to the degrees of freedom in terms of (-1)\n\
                            dimensions over 2. (Valid values from -1).\n\
-s <scale>      [Real]      Uninformative prior parameters with the           (1)\n\
                            given scale factor.\n\
-p                          Calculate prior parameters from file.\n\
\n\
Startup:\n\
-v                          Enable warnings.\n\
-h                          This text.\n\
\n\
";

#endif//__hdprior_h__
