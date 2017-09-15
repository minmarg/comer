/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __compdef_h__
#define __compdef_h__


//Apply multiple alignment thickness constraints to profiles while deriving alignment between them
#define APPLY_THICKNESS_CONSTRAINTS

//Computation of score probabilities using multinomial distribution
#define USE_MULTINOMIAL_PROBS

//Universal scoring scheme based on distribution of frequency vectors
#define UNIVERSALSCORES


#if 0   //#ifdef UNIVERSALSCORES
    //scale profile scores before making final profiles
#   undef   SCALE_PROFILES
#else
#   define  SCALE_PROFILES
#endif

//Use profile background probabilities
// #define USEPROFBACKPROBS


//flag deciding the type of scores to be used;
// if not double, then type int will be used
//#define USE_DOUBLE_AS_SCORE_TYPE


typedef int     TScore;




#endif//__compdef_h__
