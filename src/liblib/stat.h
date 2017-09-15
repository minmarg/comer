/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __stat_h__
#define __stat_h__

#include "debug.h"

#define LARGE_EVALUE  1.0e+7

class StatModel {
public:
    enum Type {
        whole_DS,           //the whole dataset
        align7DS,           //the dataset composed to encompass the alignments of length >7
        whole_DistHom,      //the whole dataset of the distant homologous instances
        PC_DIST_HOM_0_2,    //estimates were derived using version 0.2
        no_types            //number of different types available
    };
    StatModel( Type, int = 0, int = 0 );

    double ComputeLambda( int len_first, int len_secnd );   //compute lambda statistical parameter
    double ComputeMu( int len_first, int len_secnd );       //compute mu statistical parameter

    double ComputeAdjustedLambda( int, int, double scale_lambda, double entropy );  //compute entropy-adjusted parameter lambda
    double ComputeAdjustedMu( int, int, double scale_lambda, double entropy );      //compute entropy-adjusted parameter mu

protected:
    StatModel(){};

    double  GetScaleLambda() const  { return scale_lambda_; }
    double  GetEntropy() const      { return entropy_; }

    void    SetScaleLambda( double value )  { scale_lambda_ = value; }
    void    SetEntropy( double value )      { entropy_ = value; }

    void    SaveValues( int len1, int len2 );                                       //saves values
    void    SaveValues( int len1, int len2, double scale_lambda, double entropy );  //saves values overloaded

private:
    // estimates of the statistical parameters
    // parameters for lambda (decay) estimate
    static const double   lambda_a;
    static const double   lambda_b;

    // parameters for mu (location) estimate
    static const double   mu_a;
    static const double   mu_b;

    // ----------
    // parameters for lambda (decay) estimate; model align7
    static const double   lambda7_a;
    static const double   lambda7_b;

    // parameters for mu (location) estimate; model align7
    static const double   mu7_a;
    static const double   mu7_b;

    // ----------
    // parameters for lambda (decay) estimate; distant homologies (20% identity)
    static const double   lambdaDH_a;
    static const double   lambdaDH_b;

    // parameters for mu (location) estimate; distant homologies (20% identity)
    static const double   muDH_a;
    static const double   muDH_b;

    // ----------
    // parameters for lambda (decay) estimate; distant homologies (20% identity) PC_DIST_HOM_0_2
    static const double   lambdaPCDH_02_a[];
    static const double   lambdaPCDH_02_b[];

    // parameters for mu (location) estimate; distant homologies (20% identity) PC_DIST_HOM_0_2
    static const double   muPCDH_02_a[];
    static const double   muPCDH_02_b[];

    //type describing a model
    Type    type;

    double  lambda;     //precomputed value for lambda
    double  mu;         //computed value for mu

    double  adjusted_lambda;    //precomputed value for adjusted lambda
    double  adjusted_mu;        //computed value for adjusted mu

    int     firstSize;  //length of the first profile
    int     secndSize;  //length of the second profile
    double  product;    //product of the lengths

    double  effective_first;    //effective length of the first sequence profile
    double  effective_secnd;    //effective length of the second sequence profile
    double  effective_product;  //product of the effective lengths
    double  scale_lambda_;      //scale parameter lambda which was used to compute effective lengths
    double  entropy_;           //entropy value which was used to compute effective lengths
};


#endif//__stat_h__
