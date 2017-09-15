/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __hdpclust_h__
#define __hdpclust_h__

#define DEFAULT_PARALLEL_CLUST ( 3 )
#define DEFAULT_NUMBER_OF_SCANS ( 1000 )
#define DEFAULT_NUMBER_OF_IRGS_SCANS ( -1 )
#define DEFAULT_NUMBER_OF_MH_UPDATES ( 1 )
#define DEFAULT_TAU ( -1.0 )
#define DEFAULT_GAMMA ( -1.0 )
#define DEFAULT_KAPPA0 ( -1.0 )
#define DEFAULT_NU0 ( -1.0 )
#define DEFAULT_UNINF_SCALEFAC ( 1.0 )
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
Cluster profile column frequencies under Hierarchical Dirichlet\n\
Process.\n\
NOTE: The program runs via MPI; requires MPICH2 installed.\n\
(C)2012 Mindaugas Margelevicius,VU Institute of Biotechnology,Vilnius\n\
\n\
\n\
Usage:\n\
<> -r <input> -o <output> [options]\n\
\n\
Parameters (with default values in parenthesis):\n\
\n\
-r <input>      [Filename]  Name pattern of restart-with files;\n\
                            <input>.grp -- grouped frequencies,\n\
                            <input>.par -- file of parameters.\n\
-o <output>     [Filename]  Name pattern of output files.\n\
\n\
Clustering options:\n\
-e                          A cluster for each initial sample vector.\n\
-a [1|2|3]                  Clustering of parallel processing in turn of      (3)\n\
                            1 -- single vector/table,\n\
                            2 -- single group,\n\
                            3 -- all groups.\n\
-t                          Do not apply whole table sampling in a scan.\n\
-k <number>     [Integer]   Number of sampling scans (iterations).         (1000)\n\
-u <number>     [Integer]   Number of M-H updates in a single scan.           (1)\n\
-m <number>     [Integer]   Number of M-H interm. restricted sampling scans; (-1)\n\
                            -1 means likelihood sampling instead.\n\
-l [1|2]                    1/2 -- Do split/merge proposals only in M-H updates.\n\
-j                          Use modified JN MCMC algorithm in M-H updates.\n\
-d                          Select clusters uniformly in M-H updates.\n\
-f                          Select cluster vectors uniformly in M-H updates.\n\
-y <number>     [Real]      Prior value of concentration parameter tau.      (-1)\n\
                            -1 means to sample parameter tau.\n\
-g <number>     [Real]      Prior value of concentration parameter gamma.    (-1)\n\
                            -1 means to sample parameter gamma.\n\
-b <number>     [Real]      Prior value of scale factor kappa_0.             (-1)\n\
                            -1 means to sample scale factor.\n\
-c <number>     [Real]      Prior value of degrees of freedom.               (-1)\n\
                            -1 means to sample degrees of freedom.\n\
-n <number>     [Real]      Adjustment to the degrees of freedom in terms of (-1)\n\
                            dimensions over 2. (Valid values from -1).\n\
-s <scale>      [Real]      Scale factor for uninformative prior parameters.  (1)\n\
-p                          Calculate prior parameters from file of grouped\n\
                            frequencies.\n\
\n\
Startup:\n\
-v                          Enable warnings.\n\
-h                          This text.\n\
\n\
";

#endif//__hdpclust_h__
