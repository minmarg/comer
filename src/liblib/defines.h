/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __defines_h__
#define __defines_h__


// default filename of scores to be read 
#define DEFAULT_PSCORES_FILENAME    "PSCORES"

// default filename of mixture parameters to be read
#define DEFAULT_MIXTURES_FILENAME   "MIXTURES"

// progressive number of memory units corresponding to the profile
//      positions to allocate at once
#define ALLOCPOS    1000
// more than average number of sequences supposed to be in multiple alignment;
//      reserve memory with blocks of this size
#define ALLOCSEQ    100
// initial hit reserve
#define ALLOCHITS   1024

// constant definitions for free memory reservation
#define BUF_MAX     256

// constant definitions for free memory reservation
#define KBYTE       1024

// Sequence identity level at which sequences to be purged
#define IDENTITY_LEVEL  0.94

// Sequence identity percentage at which sequences to be clustered
#define CLUSTERING_THRESHOLD    0.62


// default cost for opening gap
// in case of position-specific costs are not in use, otherwise
// gap opening costs are position-specific and they are computed on the fly
#define DEFAULTGAPOPENCOST          ( -7 )

// default cost for extending a gap
#define DEFAULTEXTENDCOST           ( -1 )

// deletion probability weight
#define ALN_DELETION_COEFFICIENT    ( 0.6 )


//Evalue threshold for a pair of profiles above which gap prob. factor takes effect
#define ALN_GAPPROBABFACTOR_EVALUE          ( 1.0e-5 )

//Argument weight in expression of computing gap prob. factor
#define ALN_GAPPROBABFACTOR_ARGWEIGHT       ( 0.4 )

//Argument shift in expression of computing gap prob. factor
#define ALN_GAPPROBABFACTOR_ARGSHIFT        ( 0.0 )


// Numerator of expression to compute 1st-pass autocorrection
#define ALN_AUTOCORRECTION_NUMERATOR_1ST        ( 5.0 )

// Numerator of expression to compute 1st-pass positional autocorrections
#define ALN_POSIT_AUTOCORRECTION_NUMERATOR_1ST  ( 3.2 )

// Numerator of expression to compute 2nd-pass upper bound for autocorrection
#define ALN_AUTOCORRECTION_NUMERATOR_2ND    ( 4.0 )

// Logarithmic scale to compute 2nd-pass autocorrection
// #define ALN_AUTOCORRECTION_LOG_SCALE        ( 6.0 )
#define ALN_AUTOCORRECTION_LOG_SCALE        ( 14.0 )

// Denominator scale to compute 2nd-pass autocorrection
#define ALN_AUTOCORRECTION_DENOM_SCALE      ( 0.12 )


// Upper bound of information content threshold used in 2nd-pass computations
#define ALN_AUTOINFORMATION_UPPERBOUND      ( 0.30 )

// Numerator of expression to compute 2nd-pass information content threshold ( 1.602059991327962285 )
// #define ALN_AUTOINFORMATION_NUMERATOR2      ( 1.6 )
#define ALN_AUTOINFORMATION_NUMERATOR2      ( 4.4 )

// Logarithmic scale to compute 2nd-pass information content threshold
#define ALN_AUTOINFORMATION_SCALE2          ( 4.0 )

// Numerator of alternative expression to compute 2nd-pass information content threshold
#define ALN_AUTOINFORMATION_NUMERATOR2ALT   ( 1.0 )

// Logarithmic scale to alternatively compute 2nd-pass information content threshold
#define ALN_AUTOINFORMATION_SCALE2ALT       ( 3.0 )


// maximum number of columns a sequence or profile can have;
// if exceeded the program terminates
#define MAXCOLUMNS  50000

// alignment output width per line
#define OUTPUTWIDTH 60

// output indent
#define OUTPUTINDENT 13

// maximum description length in the program's output
#define MAX_DESCRIPTION_LENGTH 4096

// maximum of iterations before terminating scaling of scores
// to obtain specific statistical parameter lambda
#define MAX_SCALE_ITERATIONS 10

// whether to use binary search in scaling of scores to
// obtain lambda
#define SCALE_BINARY_SEARCH

// maximum number of iterations allowed to perform before
// the root of conservation function is found
#define MAXIT 100

// accuracy to which lambda parameter should be calculated
#define LAMBDA_ACCURACY ( 1.e-5 )

// accuracy of constant in relative entropy expression which is
// used to adjust scores in order to obtain desired relative entropy
#define ACCURACY_OF_CONSTANT_FOR_H ( 1.e-5 )

// small number that is appended to the lower bound of lambda
// domain (positive)
#define EPSILON ( 1.e-6 )

// small number that is appended to the lower bound of constant c
#define EPSILON_C ( 1.e-6 )


// tolerance for optimization of target frequencies by the Newton's method 
#define TARGET_FREQ_OPT_TOLERANCE   ( 1.e-6 )

// maximum number of Newton's method iterations in optimizing the target frequencies
#define TARGET_FREQ_OPT_MAXITERATIONS   ( 1000 )


// number of COLUMN residues to take into window 
#define MAC_SEGSEQ_WIN_LENGTH       ( 30 )
// log2( 30/2 ) *.85, *.90
#define MAC_SEGSEQ_LOW_ENTROPY      ( 3.315 )
#define MAC_SEGSEQ_HIGH_ENTROPY     ( 3.510 )

// minimum required window size in positions of extent
#define MAC_EXTENT_MINWIN           ( 7 )
// minimum required sequence percentage an extent must cover
#define MAC_EXTENT_SEQ_COVER        ( 0.05 )

// relative weight for residue pseudocount frequencies, beta: 7, 9, 13
#define MAC_WEIGHT_PSEUDO_COUNTS    ( 7 )


// value to which frequencies always sum
#define FREQUENCY_SUM 100

// scaling constant used to scale profile scores before
// they are processed for statistical significance parameters
#define SCALE_CONSTANT 200

// scaling constant used to scale information content value
// before it is translated into integer
#define INFO_SCALE_CONSTANT 100

// scaling constant used to scale profile-profile scores before
// using them to align two profiles
#define PP_SCALE_CONSTANT 32

// maximum range of scores allowed
#define MAX_RANGE_OF_SCORES 2000

// maximum number of iterations allowed to perform while
// computing infinite sum of the parameter K
#define PARAMETER_K_MAXIT 100

// accuracy to compute alignment probability expressed in
// infinite sum to; the probability is used to compute K
#define PARAMETER_K_ACCURACY ( 1.e-4 )

// maximum number of times to iterate searching for fixed
// value of length adjustment expression
#define LENGTH_ADJUSTMENT_MAXIT 20



#endif//__defines_h__

