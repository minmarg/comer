/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mystring.h"
#include "myexcept.h"

#include "stat.h"

// estimates of the statistical parameters
// parameters for the lambda (decay) estimate
const double StatModel::lambda_a =  2.2235;
const double StatModel::lambda_b = -0.2399;

// parameters for the mu (location) estimate
const double StatModel::mu_a =  15.3629;
const double StatModel::mu_b = -31.1988;

// ----------
// parameters for the lambda (decay) estimate; model align7
const double StatModel::lambda7_a =  1.4476;
const double StatModel::lambda7_b = -0.2157;

// parameters for the mu (location) estimate; model align7
const double StatModel::mu7_a =   5.5684;
const double StatModel::mu7_b = -25.0678;
//
// ----------
// estimates of the statistical parameters for distant homologies (20% identity)
// parameters for the lambda (decay) estimate
const double StatModel::lambdaDH_a =  1.5584;
const double StatModel::lambdaDH_b = -0.2010;

// parameters for the mu (location) estimate; distant homologies (20% identity)
const double StatModel::muDH_a =   5.2981;
const double StatModel::muDH_b = -24.6283;
//
// ----------
// parameters for the lambda (decay) estimate; distant homologies (20% identity) PC_DIST_HOM_0_2
const double StatModel::lambdaPCDH_02_a[] = {  1.2090,  1.1382 };
const double StatModel::lambdaPCDH_02_b[] = { -0.1810, -0.1767 };

// parameters for the mu (location) estimate; distant homologies (20% identity) PC_DIST_HOM_0_2
const double StatModel::muPCDH_02_a[] = {   5.7805,   5.5371 };
const double StatModel::muPCDH_02_b[] = { -28.4889, -25.1996 };
// ----------


// -------------------------------------------------------------------------
// constructor: 
// -------------------------------------------------------------------------

StatModel::StatModel( Type t, int sz1, int sz2 )
:   type( t ),
    lambda( 0.0 ),
    mu( 0.0 ),

    adjusted_lambda( 0.0 ),
    adjusted_mu( 0.0 ),

    firstSize( sz1 ),
    secndSize( sz2 ),

    effective_first( 0.0 ),
    effective_secnd( 0.0 ),

    effective_product( 0.0 ),
    scale_lambda_( 0.0 ),
    entropy_( 0.0 )
{
    product = firstSize * secndSize;
}

// -------------------------------------------------------------------------
// ComputeLambda: computes statistical significance parameter lambda
//
//                b ln(m n )
//   lambda = a e                  , where m and n are profile lengths
//
// -------------------------------------------------------------------------

double StatModel::ComputeLambda( int len_first, int len_secnd )
{
    SaveValues( len_first, len_secnd );

    switch( type ){
        case whole_DS:
                lambda = lambda_a * ::pow( product, lambda_b );
            break;

        case align7DS:
                lambda = lambda7_a * ::pow( product, lambda7_b );
            break;

        case whole_DistHom:
                lambda = lambdaDH_a * ::pow( product, lambdaDH_b );
            break;

        case PC_DIST_HOM_0_2:
                lambda = lambdaPCDH_02_a[0] * ::pow( product, lambdaPCDH_02_b[0] );
            break;

        default:
            throw myruntime_error( mystring( "Request for unknown statistical model." ));
    };
    return lambda;
}

// -------------------------------------------------------------------------
// ComputeMu: computes statistical significance parameter mu
//
//   mu = a_ ln( ln( m n ) / lambda ) + b_   , model == whole
//   mu = a_ ln( m n ) + b_                  , model == align7
//   mu = a_ ln( m n ) + b_                  , model == distant homologies
//
// -------------------------------------------------------------------------

double StatModel::ComputeMu( int len_first, int len_secnd )
{
    SaveValues( len_first, len_secnd );

    switch( type ){
        case whole_DS:
                // assumed lambda was computed before
                if( lambda == 0.0 )
                    ComputeLambda( len_first, len_secnd );
                if( lambda == 0.0 )
                    throw myruntime_error( mystring( "Wrong statistical parameter lambda returned." ));
                mu = mu_a * log( log( product ) / lambda ) + mu_b;
            break;

        case align7DS:
                mu = mu7_a * log( product ) + mu7_b;
            break;

        case whole_DistHom:
                mu = muDH_a * log( product ) + muDH_b;
            break;

        case PC_DIST_HOM_0_2:
                mu = muPCDH_02_a[0] * log( product ) + muPCDH_02_b[0];
            break;

        default:
            throw myruntime_error( mystring( "Request for unknown statistical model." ));
    };
    return mu;
}

// -------------------------------------------------------------------------
// ComputeAdjustedLambda: computes entropy-adjusted parameter lambda;
//     lambda is computed in dependece on the product of effective lengths
// -------------------------------------------------------------------------

double StatModel::ComputeAdjustedLambda( int len_first, int len_secnd, double scale_lambda, double entropy )
{
    SaveValues( len_first, len_secnd, scale_lambda, entropy );

    if( effective_first <= 0 || effective_secnd <= 0 )
        return LARGE_EVALUE;

    if( GetEntropy() == 0.0 )
        adjusted_lambda = ComputeLambda( len_first, len_secnd );
    else
        adjusted_lambda = lambdaPCDH_02_a[1] * ::pow( effective_product, lambdaPCDH_02_b[1] );

    return adjusted_lambda;
}

// -------------------------------------------------------------------------
// ComputeAdjustedMu: computes entropy-adjusted parameter mu;
//     mu is computed in dependece on the product of effective lengths;
//     effective length is expressed as
//          length - scale_lambda * avg_mu / entropy,
//     where avg_mu is average score for unrelated sequence pairs given
//     lengths
// -------------------------------------------------------------------------

double StatModel::ComputeAdjustedMu( int len_first, int len_secnd, double scale_lambda, double entropy )
{
    SaveValues( len_first, len_secnd, scale_lambda, entropy );

    if( effective_first <= 0 || effective_secnd <= 0 )
        return LARGE_EVALUE;

    if( GetEntropy() == 0.0 )
        adjusted_mu = ComputeMu( len_first, len_secnd );
    else
        adjusted_mu = muPCDH_02_a[1] * log( effective_product ) + muPCDH_02_b[1];

    return adjusted_mu;
}

// -------------------------------------------------------------------------
// SaveValues: saves values
// -------------------------------------------------------------------------

void StatModel::SaveValues( int len_first, int len_secnd )
{
#ifdef __DEBUG__
    if( len_first <= 0 || len_secnd <= 0 )
        throw myruntime_error( mystring( "Unable to compute statistical parameters." ));
#endif

    if( len_first != firstSize || len_secnd != secndSize )
    {
        firstSize = len_first;
        secndSize = len_secnd;
        product = firstSize * secndSize;
    }
}

// -------------------------------------------------------------------------
// SaveValues: saves values
// -------------------------------------------------------------------------

void StatModel::SaveValues( int len_first, int len_secnd, double scale_lambda, double entropy )
{
    if( len_first == firstSize && len_secnd == secndSize &&
        GetScaleLambda() == scale_lambda && GetEntropy() == entropy )
        return;

    if( scale_lambda < 0.0 || entropy < 0.0 )
        throw myruntime_error( mystring( "Unable to compute statistical parameters." ));


    SaveValues( len_first, len_secnd );

    SetScaleLambda( scale_lambda );
    SetEntropy( entropy );

    if( GetEntropy() == 0.0 ) {
        effective_first = len_first;
        effective_secnd = len_secnd;
        effective_product = product;
        return;
    }

    //average expected alignment length
    double  avg_align_length = ComputeMu( len_first, len_secnd ) * GetScaleLambda() / GetEntropy();

    effective_first = len_first - avg_align_length;
    effective_secnd = len_secnd - avg_align_length;

    effective_product = effective_first * effective_secnd;
}
