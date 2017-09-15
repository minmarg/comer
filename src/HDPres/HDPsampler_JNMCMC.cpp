/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "data.h"
#include "ext/psl.h"
#include "ext/gamma.h"
#include "HDPsampler.h"


//random number generators
MTRng   RNGjl;
MTRng   RNGja;

// Mettropolis-Hastings update calculations
// Implementation of a split-merge MCMC procedure proposed by Jain & Neal
//
// =========================================================================
// MHupdate: Mettropolis-Hastings update
//  d1, d2 -- dish indices of vectors v1 and v2
//  v1, v2 -- vector indices in dishes d1 and d2, respectively
//
bool HDPsampler::MHupdate( int d1, int v1, int d2, int v2 )
{
    mystring errstr;
    Dish*   dsh = NULL;
    double  ltprat;//log ratio of state transition probabilities
    double  lpprat;//log of prior probability ratio
    double  llhrat;//log ratio of likelihoods
    double  accprob;//M-H acceptance probability
    double  urn = 0.0;
    int dn1 = -1, dn2 = -1;//indices of new auxiliary dishes
    int vn1 = -1, vn2 = -1;//new indices of vectors in auxiliary dishes
    bool    accept = false;
    bool    rok = true;

    try {
        if( GetNoRstrGSScans() < 0 ) {
            if( rok )
                rok = MHInitLaunchStateThroughSampling( d1, v1, d2, v2, &dn1, &dn2, &vn1, &vn2 );
        } else {
            if( rok )
                rok = MHInitLaunchState( d1, v1, d2, v2, &dn1, &dn2, &vn1, &vn2 );
            if( rok )
                rok = MHIntmRstrGSScans( dn1, dn2, vn1, vn2 );
        }
        if( rok )
            rok = MHTransProbRatio( d1, v1, d2, v2, dn1, dn2, vn1, vn2, &ltprat );
        if( rok )
            rok = MHPriorProbRatio( d1, v1, d2, v2, dn1, dn2, vn1, vn2, &lpprat );
        if( rok )
            rok = MHLikelihoodRatio( d1, v1, d2, v2, dn1, dn2, vn1, vn2, &llhrat );

        if( rok ) {
            accprob = ltprat + lpprat + llhrat;
            if( accprob < cdp_LOG_DBL_MIN )
                accprob = 0.0;
            else if( accprob < cdp_LOG_DBL_MAX )
                accprob = exp( accprob );
            urn = RNGja.GetDouble();

            if( urn <= accprob ) {
                accept = true;
                if( !MHSaveSplitMergeInfo( d1, v1, d2, v2, dn1, dn2, vn1, vn2 ))
                    accept = false;
            }
        }
    } catch( myexception const& ex ) {
        errstr = ex.what();
    }

    //remove auxiliary dishes: first, the latest one
    if( 0 <= dn2 && dn2 < GetMenu()->GetSize()) {
        dsh = GetMenu()->GetDishAt( dn2 );
        GetMenu()->RemTmpDishAt( dn2, dsh );
    }
    if( 0 <= dn1 && dn1 < GetMenu()->GetSize()) {
        dsh = GetMenu()->GetDishAt( dn1 );
        GetMenu()->RemTmpDishAt( dn1, dsh );
    }

    if( !errstr.empty())
        throw myruntime_error( errstr );

    return accept;
}





// -------------------------------------------------------------------------
// MHInitLaunchState: Initialize launch state
//  d1, d2 -- indices of dishes having vectors v1 and v2
//  v1, v2 -- vector indices in dishes d1 and d2, respectively
//  pd1, pd2 -- indices of new auxiliary dishes
//  pv1, pv2 -- new vector indices in dishes pd1 and pd2, respectively
//
bool HDPsampler::MHInitLaunchState( int d1, int v1, int d2, int v2,
                int* pd1, int* pd2, int* pv1, int* pv2 )
{
    Dish*   dish1 = NULL;
    Dish*   dish2 = NULL;
    const Pslvector* vec1;
    const Pslvector* vec2;
    Dish*   mhdish1 = NULL; //temporary dish 1 for M-H update
    Dish*   mhdish2 = NULL; //temporary dish 2 for M-H update
    int     mhdd1, mhvv1;
    int     mhdd2, mhvv2;
    Dish*   dsh = NULL;
    double  urn = 0.0;
    int     novs = 0;
    int n, b;

    if( pd1 == NULL || pd2 == NULL ||
        pv1 == NULL || pv2 == NULL )
        throw myruntime_error( "MHInitLaunchState: Memory access error." );

    if( GetMenu() == NULL )
        throw myruntime_error( "MHInitLaunchState: Null menu." );

    if( d1 < 0 || GetMenu()->GetSize() <= d1 ||
        d2 < 0 || GetMenu()->GetSize() <= d2 )
        throw myruntime_error( "MHInitLaunchState: Invalid dish indices." );

    if(( dish1 = GetMenu()->GetDishAt( d1 )) == NULL ||
       ( dish2 = GetMenu()->GetDishAt( d2 )) == NULL )
        throw myruntime_error( "MHInitLaunchState: Null dishes at indices." );

    if( v1 < 0 || dish1->GetSize() <= v1 || 
        v2 < 0 || dish2->GetSize() <= v2 )
        throw myruntime_error( "MHInitLaunchState: Invalid vector indices." );

    if(( vec1 = dish1->GetVectorNAt( v1 )) == NULL ||
       ( vec2 = dish2->GetVectorNAt( v2 )) == NULL )
        throw myruntime_error( "MHInitLaunchState: Null vectors at dish indices." );

    novs = dish1->GetActualSize();
    if( dish1 != dish2 )
        novs += dish2->GetActualSize();

    if( novs < 3 )
        //M-H proposal will not be accepted if no.vectors<3
        return false;

    //create two auxiliary dishes
    mhdish1 = new Dish( GetDefDishSize());
    mhdish2 = new Dish( GetDefDishSize());
    if( mhdish1 == NULL || mhdish2 == NULL )
        throw myruntime_error( "MHInitLaunchState: Not enough memory." );
    mhdish1->SetBasin( GetBasin());
    mhdish2->SetBasin( GetBasin());
    mhdd1 = GetMenu()->NewTmpDish( mhdish1 );
    mhdd2 = GetMenu()->NewTmpDish( mhdish2 );

    //copy v1 and v2 to distinct dishes
    mhvv1 = mhdish1->NewVectorNInd( dish1->GetVectorNIndAt( v1 ));
    mhvv2 = mhdish2->NewVectorNInd( dish2->GetVectorNIndAt( v2 ));

    //distribute other vectors from d1 and d2 randomly among two new dishes
    for( dsh = dish1;; dsh = dish2 ) {
        for( n = 0; n < dsh->GetSize(); n++ ) {
            if(( dsh == dish1 && n == v1 )||( dsh == dish2 && n == v2 ))
                //omit those already assigned
                continue;
            b = dsh->GetVectorNIndAt( n );
            if( b < 0 )
                continue;
            urn = RNGjl.GetDouble();
            if( urn <= 0.5 )
                mhdish1->NewVectorNInd( b );
            else
                mhdish2->NewVectorNInd( b );
        }
        if( dish1 == dish2 || dsh == dish2 )
            //if two dishes are the same or both have been processed
            break;
    }

    if( novs != mhdish1->GetActualSize() + mhdish2->GetActualSize())
        throw myruntime_error( "MHInitLaunchState: Inconsistent no. vectors in new dishes." );

    //calculate parameters of new auxiliary dishes
    RecalcDishParams( mhdd1 );
    RecalcDishParams( mhdd2 );

    *pd1 = mhdd1; *pv1 = mhvv1;
    *pd2 = mhdd2; *pv2 = mhvv2;
    return true;
}





// -------------------------------------------------------------------------
// MHInitLaunchStateThroughSampling: Initialize launch state through 
//  sampling according to likelihoods of dishes
//  d1, d2 -- indices of dishes having vectors v1 and v2
//  v1, v2 -- vector indices in dishes d1 and d2, respectively
//  pd1, pd2 -- indices of new auxiliary dishes
//  pv1, pv2 -- new vector indices in dishes pd1 and pd2, respectively
//
bool HDPsampler::MHInitLaunchStateThroughSampling( int d1, int v1, int d2, int v2,
                int* pd1, int* pd2, int* pv1, int* pv2 )
{
    const int   prd = 5;
    mystring preamb = "HDPsampler: MHInitLaunchStateThroughSampling: ";
    mystring errstr;

    if( GetMenu() == NULL )
        throw myruntime_error( preamb + "Null menu." );

    Dish*   dish1 = NULL;
    Dish*   dish2 = NULL;
    const Pslvector* vec1;
    const Pslvector* vec2;
    const Pslvector* vec;
    const double nu0 = GetMenu()->GetNu0();
    double* values = GetResMBufferValues();
    double  lprob1, lprob2;
    Dish*   mhdish1 = NULL; //temporary dish 1 for M-H update
    Dish*   mhdish2 = NULL; //temporary dish 2 for M-H update
    int     mhdd1, mhvv1;
    int     mhdd2, mhvv2;
    int     fact1, fact2;
    int     fnewt, newk;//fake index of table and index of new dish
    Dish*   dsh = NULL;
    int     novs = 0;
    int n, vloc, b, pp;

    if( pd1 == NULL || pd2 == NULL ||
        pv1 == NULL || pv2 == NULL )
        throw myruntime_error( preamb + "Memory access error." );

    if( values == NULL )
        throw myruntime_error( preamb + "Null buffer values.");

    if( d1 < 0 || GetMenu()->GetSize() <= d1 ||
        d2 < 0 || GetMenu()->GetSize() <= d2 )
        throw myruntime_error( preamb + "Invalid dish indices." );

    if(( dish1 = GetMenu()->GetDishAt( d1 )) == NULL ||
       ( dish2 = GetMenu()->GetDishAt( d2 )) == NULL )
        throw myruntime_error( preamb + "Null dishes at indices." );

    if( v1 < 0 || dish1->GetSize() <= v1 || 
        v2 < 0 || dish2->GetSize() <= v2 )
        throw myruntime_error( preamb + "Invalid vector indices." );

    if(( vec1 = dish1->GetVectorNAt( v1 )) == NULL ||
       ( vec2 = dish2->GetVectorNAt( v2 )) == NULL )
        throw myruntime_error( preamb + "Null vectors at dish indices." );

    novs = dish1->GetActualSize();
    if( dish1 != dish2 )
        novs += dish2->GetActualSize();

    if( novs < 3 )
        //M-H proposal will not be accepted if no.vectors<3
        return false;

    //create two auxiliary dishes
    mhdish1 = new Dish( GetDefDishSize());
    mhdish2 = new Dish( GetDefDishSize());
    if( mhdish1 == NULL || mhdish2 == NULL )
        throw myruntime_error( preamb + "Not enough memory." );
    mhdish1->SetBasin( GetBasin());
    mhdish2->SetBasin( GetBasin());
    mhdd1 = GetMenu()->NewTmpDish( mhdish1 );
    mhdd2 = GetMenu()->NewTmpDish( mhdish2 );

    //copy v1 and v2 to distinct dishes
    mhvv1 = mhdish1->NewVectorNInd( dish1->GetVectorNIndAt( v1 ));
    mhvv2 = mhdish2->NewVectorNInd( dish2->GetVectorNIndAt( v2 ));

    //calculate parameters of new auxiliary dishes
    RecalcDishParams( mhdd1 );
    RecalcDishParams( mhdd2 );

    //reset nu0; nu0 is used just to calculate probabilities;
    //it is not used in dish parameter calculations
//     if( 1.0 < nu0 )
//         GetMenu()->SetNu0( 1.0 );

    try {
        //distribute other vectors from d1 and d2 randomly among two new dishes
        for( dsh = dish1;; dsh = dish2 ) {
            for( n = 0; n < dsh->GetSize(); n++ ) {
                if(( dsh == dish1 && n == v1 )||( dsh == dish2 && n == v2 ))
                    //omit those already assigned
                    continue;
                b = dsh->GetVectorNIndAt( n );
                if( b < 0 )
                    continue;
                vec = dsh->GetVectorNAt( n );
                if( vec == NULL )
                    throw myruntime_error( "Null vector." );

                //sample between two dishes
                pp = 0;
                fact1 = mhdish1->GetActualSize();
                fact2 = mhdish2->GetActualSize();
                if( fact1 < 1 || fact2 < 1 )
                    throw myruntime_error( "Invalid sizes of aux. dishes." );
                ProbVecOfDish( vec, mhdd1, &lprob1 );
                ProbVecOfDish( vec, mhdd2, &lprob2 );

                values[pp++] = -1.0;//fake value (table index)
                values[pp++] = ( double )mhdd1;//dish index
                values[pp++] = ( double )1;//fact1;//factor of dish 1; //NOTE:1
                values[pp++] = lprob1;
                values[pp++] = lprob1;//for normalized prob. value
                values[pp++] = -1.0;//fake value (table index)
                values[pp++] = ( double )mhdd2;//dish index
                values[pp++] = ( double )1;//fact2;//factor of dish 2; //NOTE:1
                values[pp++] = lprob2;
                values[pp++] = lprob2;//for normalized prob. value

                //sample
                if( !SlaveProcessSampleFromProbs( values, pp, &fnewt, &newk, true ))
                    return false;

                if( newk != mhdd1 && newk != mhdd2 )
                    throw myruntime_error( "Invalid index of aux. dish." );

                if( values[3] < 0.0 || 1.0 < values[3])
                    throw myruntime_error( "Prob. range error." );
                if( values[3+prd] < 0.0 ||( newk == mhdd2 && 1.0 < values[3+prd]))
                    throw myruntime_error( "Prob. range error." );

                //IMPORTANT: add vector to dish AFTER updating dish params
                AdjustDishParams( newk, vec, true/*add*/);
                if( newk == mhdd1 )
                    vloc = mhdish1->NewVectorNInd( b );
                else
                    vloc = mhdish2->NewVectorNInd( b );
            }
            if( dish1 == dish2 || dsh == dish2 )
                //if two dishes are the same or both have been processed
                break;
        }
    } catch( myexception const& ex ) {
        errstr = ex.what();
    }

//     //set back nu0
//     GetMenu()->SetNu0( nu0 );

    if( !errstr.empty())
        throw myruntime_error( preamb + errstr );

    if( novs != mhdish1->GetActualSize() + mhdish2->GetActualSize())
        throw myruntime_error( preamb + "Inconsistent no. vectors in new dishes." );

    *pd1 = mhdd1; *pv1 = mhvv1;
    *pd2 = mhdd2; *pv2 = mhvv2;
    return true;
}





// -------------------------------------------------------------------------
// MHIntmRstrGSScans: perform intermediate restricted Gibbs sampling 
//  scans from the initial launch state
//  dn1, dn2 -- indices of new auxiliary dishes
//  vn1, vn2 -- new vector indices in auxiliary dishes
//
bool HDPsampler::MHIntmRstrGSScans( int dn1, int dn2, int vn1, int vn2 )
{
    int i;
    for( i = 0; i < GetNoRstrGSScans(); i++ ) {
        if( !MHIntmRstrGSSingle( dn1, dn2, vn1, vn2 ))
            return false;
    }
    return true;
}

// -------------------------------------------------------------------------
// MHIntmRstrGSSingle: perform single intermediate restricted Gibbs 
//  sampling scan from the initial launch state
//  dn1, dn2 -- indices of new auxiliary dishes
//  vn1, vn2 -- new vector indices in auxiliary dishes
//  varray -- array of vector indices assigned to one of two original dishes
//  ltprob -- log of transition probability from the current to new state
// NOTE: `varray' is supposed to have allocated memory of size that of basin
//
bool HDPsampler::MHIntmRstrGSSingle( int dn1, int dn2, int vn1, int vn2, 
                                      int* varray, double* ltprob )
{
    const int   prd = 5;
    Dish*   mhdish1 = NULL; //temporary dish 1
    Dish*   mhdish2 = NULL; //temporary dish 2
    Dish*   dsh = NULL;
    const Pslvector* vec1, *vec2;
    const Pslvector* vec;
    double* values = GetResMBufferValues();
    double  lprob1, lprob2;
    int     fact1, fact2;
    int     fnewt, newk;//fake index of table and index of new dish
    int n, vloc, bn, dn, pp;

    if( GetMenu() == NULL )
        throw myruntime_error( "MHInterRstrGSSingle: Null menu." );

    if( values == NULL )
        throw myruntime_error("MHInterRstrGSSingle: Null buffer values.");

    if( dn1 < 0 || GetMenu()->GetSize() <= dn1 ||
        dn2 < 0 || GetMenu()->GetSize() <= dn2 )
        throw myruntime_error( "MHInterRstrGSSingle: Invalid dish indices." );

    if(( mhdish1 = GetMenu()->GetDishAt( dn1 )) == NULL ||
       ( mhdish2 = GetMenu()->GetDishAt( dn2 )) == NULL )
        throw myruntime_error( "MHInterRstrGSSingle: Null dishes at indices." );

    if( vn1 < 0 || mhdish1->GetSize() <= vn1 || 
        vn2 < 0 || mhdish2->GetSize() <= vn2 )
        throw myruntime_error( "MHInterRstrGSSingle: Invalid vector indices." );

    if(( vec1 = mhdish1->GetVectorNAt( vn1 )) == NULL ||
       ( vec2 = mhdish2->GetVectorNAt( vn2 )) == NULL )
        throw myruntime_error( "MHInterRstrGSSingle: Null vectors at dish indices." );

    //set the processed flag of all vectors to false
    for( n = 0; n < mhdish1->GetSize(); n++ ) mhdish1->SetProcessedAt( n, false );
    for( n = 0; n < mhdish2->GetSize(); n++ ) mhdish2->SetProcessedAt( n, false );

    mhdish1->SetProcessedAt( vn1, true );//preset support vectors as processed
    mhdish2->SetProcessedAt( vn2, true );

    if( ltprob )
        *ltprob = 0.0;

    //iterate over all vectors from the two new dishes
    for( dsh = mhdish1, dn = dn1;; dsh = mhdish2, dn = dn2 ) {
        for( n = 0; n < dsh->GetSize(); n++ ) {
            if(( dsh == mhdish1 && n == vn1 )||( dsh == mhdish2 && n == vn2 ))
                //omit support vectors
                continue;
            bn = dsh->GetVectorNIndAt( n );
            if( bn < 0 )
                continue;
            vec = dsh->GetVectorNAt( n );
            if( vec == NULL )
                throw myruntime_error( "MHInterRstrGSSingle: Null vector." );

            if( dsh->GetProcessedAt( n ))
                continue;//already processed

            //sample between two dishes
            pp = 0;
            fact1 = mhdish1->GetActualSize();
            fact2 = mhdish2->GetActualSize();
            if( fact1 < 1 || fact2 < 1 )
                throw myruntime_error( "MHInterRstrGSSingle: Invalid sizes of aux. dishes." );
            if( dsh == mhdish1 ) {
                ProbVecOfDishExc( vec, dn1, &lprob1 );
                ProbVecOfDish( vec, dn2, &lprob2 );
                fact1--;
            }
            else {
                ProbVecOfDish( vec, dn1, &lprob1 );
                ProbVecOfDishExc( vec, dn2, &lprob2 );
                fact2--;
            }
            values[pp++] = -1.0;//fake value (table index)
            values[pp++] = ( double )dn1;//dish index
            values[pp++] = ( double )fact1;//factor of dish 1
            values[pp++] = lprob1;
            values[pp++] = lprob1;//for normalized prob. value
            values[pp++] = -1.0;//fake value (table index)
            values[pp++] = ( double )dn2;//dish index
            values[pp++] = ( double )fact2;//factor of dish 2
            values[pp++] = lprob2;
            values[pp++] = lprob2;//for normalized prob. value
            //if not given varray (vector original positions), sample;
            //otherwise, it will be needed all normalized probability values
            if( !SlaveProcessSampleFromProbs( values, pp, &fnewt, &newk, varray? false: true ))
                return false;

            if( varray ) {
                //if given array of vector original indices,
                //explicitly sample vectors back to their original state
                if( varray[bn] != 1 && varray[bn] != 2 )
                    throw myruntime_error( "MHInterRstrGSSingle: Invalid value from temp. array." );
                if( varray[bn] == 1 )
                      newk = dn1;
                else  newk = dn2;
            }

            if( newk != dn1 && newk != dn2 )
                throw myruntime_error( "MHInterRstrGSSingle: Invalid index of aux. dish." );

            //set processed
            dsh->SetProcessedAt( n, true );

            if( values[3] < 0.0 || 1.0 < values[3])
                throw myruntime_error( "MHInterRstrGSSingle: Prob. range error." );
            if( values[3+prd] < 0.0 ||( newk == dn2 && 1.0 < values[3+prd]))
                throw myruntime_error( "MHInterRstrGSSingle: Prob. range error." );
            //update log product of probs: on output `values' contain normalized probabilities
            if( ltprob ) {
              if( newk == dn1 ) *ltprob += ( values[3] <= 0.0 )? cdp_LOG_DBL_MIN: log( values[3]);
                else    *ltprob += ( values[3+prd] <= 0.0 )? cdp_LOG_DBL_MIN: log( values[3+prd]);
            }

            if( newk == dn )
                //if to be moved to the same dish
                continue;

            if( dsh->GetActualSize() < 2 )
                throw myruntime_error( "MHInterRstrGSSingle: Invalid aux. dish size." );

            //IMPORTANT: remove vector from dish AFTER updating dish params
            AdjustDishParams( dn, vec, false/*subtract*/);
            dsh->RemValueAt( n, vec );
            //IMPORTANT: add vector to dish AFTER updating dish params
            AdjustDishParams( newk, vec, true/*add*/);
            if( dsh == mhdish1 ) {
                vloc = mhdish2->NewVectorNInd( bn );
                mhdish2->SetProcessedAt( vloc, true );
            }
            else {
                vloc = mhdish1->NewVectorNInd( bn );
                mhdish1->SetProcessedAt( vloc, true );
            }
        }
        if( dsh == mhdish2 )
            //if both have been processed
            break;
    }

    //check whether all vectors have been processed
    for( dsh = mhdish1;; dsh = mhdish2 ) {
        for( n = 0; n < dsh->GetSize(); n++ ) {
            if( dsh->GetVectorNIndAt( n ) < 0 )
                continue;
            if( !dsh->GetProcessedAt( n ))
              throw myruntime_error("MHInterRstrGSSingle: Not all vectors processed.");
        }
        if( dsh == mhdish2 )
            break;
    }
    return true;
}





// -------------------------------------------------------------------------
// MHTransProbRatio: calculate log ratio of transition probabilities 
//  d1, d2 -- indices of dishes having vectors v1 and v2
//  v1, v2 -- vector indices in dishes d1 and d2, respectively
//  dn1, dn2 -- indices of new auxiliary dishes
//  vn1, vn2 -- new vector indices in auxiliary dishes
//  ltprob -- log of transition probability from the current to the proposal
//      state
//
bool HDPsampler::MHTransProbRatio( int d1, int v1, int d2, int v2, 
                                   int dn1, int dn2, int vn1, int vn2, double* ltprat )
{
    double  ltprob;

    if( ltprat == NULL )
        throw myruntime_error("MHTransProbRatio: Memory access error.");

    if( !MHTransitionProb( d1, v1, d2, v2, dn1, dn2, vn1, vn2, &ltprob ))
        return false;

    //if proposal state k* is merge state, then q(k*|k)=1;
    //MHTransitionProb calculates q(k|k*) properly (merge-to original state transition);
    //log ratio of q(k|k*)/q(k*|k) becomes log q(k|k*)
    *ltprat = ltprob;
    if( d1 == d2 )
        //if proposal state k* is split state, then q(k|k*)=1;
        //so that log ratio of q(k|k*)/q(k*|k) becomes -log q(k*|k)
        *ltprat = -*ltprat;

    return true;
}

// -------------------------------------------------------------------------
// MHTransitionProb: calculate transition probability from the original 
//  state to the proposal state
//  d1, d2 -- indices of dishes having vectors v1 and v2
//  v1, v2 -- vector indices in dishes d1 and d2, respectively
//  dn1, dn2 -- indices of new auxiliary dishes
//  vn1, vn2 -- new vector indices in auxiliary dishes
//  ltprob -- log of transition probability from the current to the proposal
//      state
//
bool HDPsampler::MHTransitionProb( int d1, int v1, int d2, int v2, 
                                   int dn1, int dn2, int vn1, int vn2, double* ltprob )
{
    if( ltprob == NULL )
        throw myruntime_error("MHTransitionProb: Memory access error.");

    if( d1 == d2 ) {
        //if the original dishes of two vectors are equal,
        //  conduct ONE final restricted Gibbs sampling scan and
        //  calculate transition probability from the 
        //  launch state to the final split state
        if( !MHIntmRstrGSSingle( dn1, dn2, vn1, vn2, NULL, ltprob ))
            return false;
        return true;
    }

    //if the original dishes of two vectors are distinct,
    //  calculate transition probability from the
    //  launch state to the original split state;
    //(in this case, it'll be transition prob. from proposal to the original state)

    mystring    errstr;
    bool        code = false;
    Dish*   dish1 = NULL;   //original dish 1
    Dish*   dish2 = NULL;   //original dish 2
    const Pslvector* vec1, *vec2;
    Dish*   dsh = NULL;
    int*    varray = NULL;
    int n, bn, dn;

    if( GetBasin() == NULL )
        throw myruntime_error( "MHTransitionProb: Null menu." );

    if( GetMenu() == NULL )
        throw myruntime_error( "MHTransitionProb: Null menu." );

    if( d1 < 0 || GetMenu()->GetSize() <= d1 ||
        d2 < 0 || GetMenu()->GetSize() <= d2 )
        throw myruntime_error( "MHTransitionProb: Invalid dish indices." );

    if(( dish1 = GetMenu()->GetDishAt( d1 )) == NULL ||
       ( dish2 = GetMenu()->GetDishAt( d2 )) == NULL )
        throw myruntime_error( "MHTransitionProb: Null dishes at indices." );

    if( v1 < 0 || dish1->GetSize() <= v1 || 
        v2 < 0 || dish2->GetSize() <= v2 )
        throw myruntime_error( "MHInitLaunchState: Invalid vector indices." );

    if(( vec1 = dish1->GetVectorNAt( v1 )) == NULL ||
       ( vec2 = dish2->GetVectorNAt( v2 )) == NULL )
        throw myruntime_error( "MHTransitionProb: Null vectors at dish indices." );

    varray = ( int* )malloc( GetBasin()->GetSize() * sizeof( int ));
    if( varray == NULL )
        throw myruntime_error( "MHTransitionProb: Not enough memory.");
    memset( varray, 0, GetBasin()->GetSize() * sizeof( int ));

    try {
        //assign vectors from the dishes 1 and 2 to original positions
        for( dsh = dish1, dn = 1;; dsh = dish2, dn = 2 ) {
            for( n = 0; n < dsh->GetSize(); n++ ) {
                if(( dsh == dish1 && n == v1 )||( dsh == dish2 && n == v2 ))
                    //omit support vectors
                    continue;
                bn = dsh->GetVectorNIndAt( n );
                if( bn < 0 )
                    continue;
                if( GetBasin()->GetSize() <= bn )
                    throw myruntime_error("MHTransitionProb: Invalid vector basin index.");
                varray[bn] = dn;
            }
            if( dsh == dish2 )
                //if both have been processed
                break;
        }

        //explicitly sample back from the launch state to the original state
        code = MHIntmRstrGSSingle( dn1, dn2, vn1, vn2, varray, ltprob );

    } catch( myexception const& ex ) {
        errstr = ex.what();
    }

    if( varray ) {
        free( varray );
        varray = NULL;
    }
    if( !errstr.empty())
        throw myruntime_error( errstr );

    return code;
}





// -------------------------------------------------------------------------
// MHCalcNoTablesOfDishes: calculate number of tables the given dishes have
//  d1, d2 -- indices of dishes having original couple of vectors (v1 and v2)
//  v1, v2 -- vector indices in dishes d1 and d2, respectively
//  dn1, dn2 -- indices of new auxiliary dishes
//  vn1, vn2 -- new vector indices in auxiliary dishes
//  md1, md2, mm -- number of tables to calculate for either 
//      dishes dn1, dn2, and d1==d2 (split state) or
//      dishes d1, d2, and merged dish d1<-d2 (merge state)
//  nots -- total number of tables
//
bool HDPsampler::MHCalcNoTablesOfDishes( int d1, int v1, int d2, int v2, 
                                    int dn1, int dn2, int vn1, int vn2, 
                                    int* md1, int* md2, int* mm, int* nots )
{
    if( GetMenu() == NULL || GetChain() == NULL )
        throw myruntime_error("MHCalcNoTablesOfDishes: Null HDP structures.");

    if( md1 == NULL || md2 == NULL || mm == NULL || nots == NULL )
        throw myruntime_error("MHCalcNoTablesOfDishes: Memory access error.");

    int*  dishnminvecs = GetLocIntBuf();//indicator of dish vectors
    const Restaurant*   rest = NULL;
    const Table*        tbl = NULL;
    const Dish*         dsh = NULL;
    const Pslvector*    frq = NULL;
    mystring    preamb = "MHCalcNoTablesOfDishes: ";
    const Dish*         dishn1 = NULL;//new dish 1
    const Dish*         dishn2 = NULL;//new dish 2
    const Dish*         dishnmin = NULL;//dish smaller of the two
    bool    split = ( d1 == d2 );//split state
    int r, t, k, v, vn, b, bn;
    int nk1, nk2;
    int ns1, ns2;//numbers of tables of dishes in split state
    int nmvmin;//number of vectors moved to min dish
    int szt;//table size

    if( d1 < 0 || GetMenu()->GetSize() <= d1 ||
        d2 < 0 || GetMenu()->GetSize() <= d2 )
        throw myruntime_error( preamb + "Invalid dish indices.");

    if( dn1 < 0 || GetMenu()->GetSize() <= dn1 ||
        dn2 < 0 || GetMenu()->GetSize() <= dn2 )
        throw myruntime_error( preamb + "Invalid dish indices.");

    if(( dishn1 = GetMenu()->GetDishAt( dn1 )) == NULL ||
      ( dishn2 = GetMenu()->GetDishAt( dn2 )) == NULL )
        throw myruntime_error( preamb + "Null dishes at indices.");

    dishnmin = dishn1;
    if( dishn2->GetActualSize() < dishn1->GetActualSize())
        dishnmin = dishn2;

    if( dishnminvecs == NULL )
        throw myruntime_error( preamb + "Null local buffer of integers.");
    if( GetSizeLocIntBuf() < GetBasin()->GetSize())
        throw myruntime_error( preamb + "Short of size of local buffer of integers.");
    ResetLocIntBuf( 0 );

    for( vn = 0; vn < dishnmin->GetSize(); vn++ ) {
        bn = dishnmin->GetVectorNIndAt( vn );
        if( bn < 0 )
            continue;
        dishnminvecs[bn] = 1;//dish has vector
    }

    *md1 = 0; *md2 = 0; *mm = 0;
    *nots = 0;

    for( r = 0; r < GetChain()->GetSize(); r++ ) {
        rest = GetChain()->GetRestaurantAt( r );
        if( rest == NULL )
            throw myruntime_error( preamb + "Null restaurant.");
        nk1 = nk2 = 0;
        ns1 = ns2 = 0;
        for( t = 0; t < rest->GetSize(); t++ ) {
            tbl = rest->GetTableAt( t );
            if( tbl == NULL )
                continue;
            k = tbl->GetDishIndex();
            if( k < 0 )
                throw myruntime_error( preamb + "Invalid table's dish index.");
            (*nots)++;
            if( k == d1 )
                nk1++;
            else if( k == d2 )
                //for merge state calculate tables for dish 2
                nk2++;
            if( !split || k != d1 )
                //if not split state or table does not belong to dish
                continue;

            //find which vectors from table are to be moved to new dish (table)
            nmvmin = szt = 0;
            for( v = 0; v < tbl->GetSize(); v++ ) {
                b = tbl->GetVectorNIndAt( v );
                if( b < 0 )
                    continue;
                szt++;
                if( dishnminvecs[b])
                    nmvmin++;
                //{{NOTE:to be removed below
//                 for( vn = 0; vn < dishnmin->GetSize(); vn++ ) {
//                     bn = dishnmin->GetVectorNIndAt( vn );
//                     if( bn < 0 )
//                         continue;
//                     if( bn == b ) {
//                         nmvmin++;
//                         break;
//                     }
//                 }
                //}}
            }
            if(( nmvmin <= 0 && dishnmin == dishn1 )||( nmvmin == szt && dishnmin != dishn1 )) {
                //all vectors moved to new dish (table)
//                 //NOTE:obsolete
//                 if( ns2 <= 0 )
//                     ns2++;//increase if table not yet created
                ns2++;//increase no.tables
            }
            else if(( nmvmin <= 0 && dishnmin != dishn1 )||( nmvmin == szt && dishnmin == dishn1 )) {
                //all vectors retained in original dish
                ns1++;//increase no.tables
            }
            else {
                //some vectors are moved to new dish (table)
                ns1++;//increase no.tables
//                 //NOTE:obsolete
//                 if( ns2 <= 0 )
//                     ns2++;//increase if table not yet created
                ns2++;//increase no.tables
            }
        }
        if( split ) {
//             //NOTE:obsolete
//             if( 1 < ns2 )
//                 throw myruntime_error( preamb + "Invalid no. tables for new dish.");
            *mm += nk1;
            *md1 += ns1;
            //when splitting, the old dish will have smaller or equal number of tables;
            //for new dish just one or zero table per restautant will be assigned
            *md2 += ns2;
        }
        else {
            *md1 += nk1;
            *md2 += nk2;
            //no actual tables are removed when dishes are merged;
            //instead, tables from the second dish are relabelled to have dish 1
            *mm += nk1 + nk2;
        }
    }

    return true;
}

// -------------------------------------------------------------------------
// MHPriorProbRatio: calculate log ratio of state prior probabilities
//  d1, d2 -- indices of dishes having vectors v1 and v2
//  v1, v2 -- vector indices in dishes d1 and d2, respectively
//  dn1, dn2 -- indices of new auxiliary dishes
//  vn1, vn2 -- new vector indices in auxiliary dishes
//  lpprat -- log of prior probability ratio (proposal to the original state)
//
bool HDPsampler::MHPriorProbRatio( int d1, int v1, int d2, int v2, 
                                   int dn1, int dn2, int vn1, int vn2, double* lpprat )
{
    Dish*   dish1 = NULL;//dish 1
    Dish*   dish2 = NULL;//dish 2
    double  gamma = GetDPMGamma();
    double  lgmd1, lgmd2, lgmm;//log gamma values
    double  err;
    mystring    preamb = "MHPriorProbRatio: ";
    bool        split = ( d1 == d2 );//split state
//  md1, md2, mm -- number of tables to calculate for either 
//      dishes dn1, dn2, and d1==d2 (split state) or
//      dishes d1, d2, and merged dish d1<-d2 (merge state)
//  nots -- total number of tables
    int     md1, md2, mm, nots;
    int     t, code;

    if( lpprat == NULL )
        throw myruntime_error( preamb + "Memory access error.");

    //dish-related concentration parameter
    if( gamma <= 0.0 )
        throw myruntime_error( preamb + "Invalid conc. parameter gamma.");

    if( GetMenu() == NULL )
        throw myruntime_error( preamb + "Null menu.");

    if( d1 < 0 || GetMenu()->GetSize() <= d1 ||
        d2 < 0 || GetMenu()->GetSize() <= d2 )
        throw myruntime_error( preamb + "Invalid dish indices." );

    if(( dish1 = GetMenu()->GetDishAt( d1 )) == NULL ||
      ( dish2 = GetMenu()->GetDishAt( d2 )) == NULL )
        throw myruntime_error( preamb + "Null dishes at indices.");

    if( !MHCalcNoTablesOfDishes( d1, v1, d2, v2,  dn1, dn2, vn1, vn2, 
                                    &md1, &md2, &mm, &nots ))
        return false;

    if( md1 < 1 || md2 < 1 || 
        mm < md1 || mm < md2 || 
        nots < md1 || nots < md2 || nots < mm )
        throw myruntime_error( preamb + "Invalid table amount values obtained.");

    if( split ) {
        //proposal state is split state;
        //calculate split to original ratio of prior probs
        if( mm != dish1->GetNoTables())
            throw myruntime_error( preamb + "Inconsistent number of tables for dish.");

        if( md1 + md2 < mm )
            throw myruntime_error( preamb + "Invalid numbers of tables in split state.");

        if(( code = psl_lngamma_e( md2, &lgmd2, &err )) != PSL_SUCCESS )
            throw myruntime_error( TranslatePSLError( code ));

        *lpprat = log( gamma ) + lgmd2;

        if( mm != md1 ) {
            if(( code = psl_lngamma_e( md1, &lgmd1, &err )) != PSL_SUCCESS )
                throw myruntime_error( TranslatePSLError( code ));
            if(( code = psl_lngamma_e( mm, &lgmm, &err )) != PSL_SUCCESS )
                throw myruntime_error( TranslatePSLError( code ));
            *lpprat += lgmd1 - lgmm;
        }

        //NOTE: APPROXIMATION: ignore a factor below to ensure valid Markov chain!!
        //in split state, total no. tables will be nots + md1 + md2 
        //for( t = nots + 1; t <= nots + md1 + md2 - mm; t++ )
        //    *lpprat -= log( gamma + t - 1.0 );
    }
    else {
        //proposal state is merge state;
        //calculate merge to original (that is split) ratio of prior probs
        if( md1 != dish1->GetNoTables() || md2 != dish2->GetNoTables())
            throw myruntime_error( preamb + "Inconsistent number of tables for dishes.");

        if( mm != md1 + md2 )
            throw myruntime_error( preamb + "Unmatched numbers of tables in merge state.");

        if(( code = psl_lngamma_e( md1, &lgmd1, &err )) != PSL_SUCCESS )
            throw myruntime_error( TranslatePSLError( code ));
        if(( code = psl_lngamma_e( md2, &lgmd2, &err )) != PSL_SUCCESS )
            throw myruntime_error( TranslatePSLError( code ));
        if(( code = psl_lngamma_e( mm, &lgmm, &err )) != PSL_SUCCESS )
            throw myruntime_error( TranslatePSLError( code ));

        *lpprat = lgmm - lgmd1 - lgmd2 - log( gamma );
    }

    return true;
}





// -------------------------------------------------------------------------
// MHLikelihoodRatio: calculate ratio of log likelihood of proposal-to-
//  original states
//  d1, d2 -- indices of dishes having vectors v1 and v2
//  v1, v2 -- vector indices in dishes d1 and d2, respectively
//  dn1, dn2 -- indices of new auxiliary dishes
//  vn1, vn2 -- new vector indices in auxiliary dishes
//  llhrat -- log ratio of likelihoods (proposal-to-original states)
//
bool HDPsampler::MHLikelihoodRatio( int d1, int v1, int d2, int v2, 
                                    int dn1, int dn2, int vn1, int vn2, double* llhrat )
{
    double  llhd1, llhd2;//likelihood of separate dishes
    double  llhd;//likelihood of merged dish

    if( llhrat == NULL )
        throw myruntime_error("MHLikelihoodRatio: Memory access error.");

    if( d1 == d2 ) {
        //proposal state is split state;
        //calculate split to original ratio of likelihoods
        if( !MHLogLikelihood( dn1, dn1, &llhd1 )) //lhood of one dish
            return false;
        if( !MHLogLikelihood( dn2, dn2, &llhd2 )) //lhood of one dish
            return false;
        if( !MHLogLikelihood( d1, d1, &llhd )) //lhood of one (merged) original dish (d1==d2)
            return false;
    }
    else {
        //proposal state is merge state;
        //calculate merge to original (ie. split) ratio of likelihoods
        if( !MHLogLikelihood( d1, d1, &llhd1 )) //lhood of one dish
            return false;
        if( !MHLogLikelihood( d2, d2, &llhd2 )) //lhood of one dish
            return false;
        if( !MHLogLikelihood( d1, d2, &llhd )) //lhood of two merged dishes
            return false;
    }

    *llhrat = llhd1 + llhd2 - llhd;
    if( d1 != d2 )
        //change sign if proposal state is merge state
        *llhrat = -*llhrat;

    return true;
}

// -------------------------------------------------------------------------
// MHLogLikelihoodObs: calculate log likelihood of observations in given dish 
//  d1, d2 -- indices of dishes to calculate likelihood for
//  lvlhood -- log likelihood of observations in dishes d1, d2
//
bool HDPsampler::MHLogLikelihoodObs( int d1, int d2, double* lvlhood )
{
    mystring errstr;
    Dish*   dish1 = NULL;
    Dish*   dish2 = NULL;
    Dish*   dsh = NULL;
    Dish*   tmpdsh = NULL;
    int     tmpdd;
    const Pslvector* vec;
    double  lprob;
    int     novs;
    int n, vloc, bn;

    if( GetMenu() == NULL )
        throw myruntime_error( "MHLogLikelihood: Null menu." );

    if( d1 < 0 || GetMenu()->GetSize() <= d1 ||
        d2 < 0 || GetMenu()->GetSize() <= d2 )
        throw myruntime_error( "MHLogLikelihood: Invalid dish indices." );

    if(( dish1 = GetMenu()->GetDishAt( d1 )) == NULL ||
       ( dish2 = GetMenu()->GetDishAt( d2 )) == NULL )
        throw myruntime_error( "MHLogLikelihood: Null dishes at index locations." );

    if( lvlhood )
        *lvlhood = 0.0;

    //create temporary dish
    tmpdsh = new Dish( GetDefDishSize());
    if( tmpdsh == NULL )
        throw myruntime_error( "MHLogLikelihood: Not enough memory." );
    tmpdsh->SetBasin( GetBasin());
    tmpdd = GetMenu()->NewTmpDish( tmpdsh );

    try {
        //calculate likelihood by sequentually adding vectors to tmp. dish
        for( novs = 0, dsh = dish1;; dsh = dish2 ) {
            novs += dsh->GetActualSize();
            for( n = 0; n < dsh->GetSize(); n++ ) {
                //include all vectors in processing
                bn = dsh->GetVectorNIndAt( n );
                if( bn < 0 )
                    continue;
                vec = dsh->GetVectorNAt( n );
                if( vec == NULL )
                    throw myruntime_error( "MHLogLikelihood: Null vector." );

                //probability of vector to belong to tmp. dish
                if( tmpdsh->GetSize() <= 0 ) PriorProbVec( vec, &lprob );
                else ProbVecOfDish( vec, tmpdd, &lprob );

                if( lvlhood )
                    *lvlhood += lprob;

                //IMPORTANT: add vector to tmp dish AFTER updating dish params
                AdjustDishParams( tmpdd, vec, true/*add*/);
                vloc = tmpdsh->NewVectorNInd( bn );
            }
            if( dish1 == dish2 || dsh == dish2 )
                //if two dishes are the same or both have been processed
                break;
        }

        if( novs != tmpdsh->GetActualSize())
            throw myruntime_error( "MHLogLikelihood: Inconsistent no. vectors in dishes." );

    } catch( myexception const& ex ) {
        errstr = ex.what();
    }

    GetMenu()->RemTmpDishAt( tmpdd, tmpdsh );

    if( !errstr.empty())
        throw myruntime_error( errstr );

    return true;
}

// -------------------------------------------------------------------------
// MHLogLikelihood: calculate log likelihood of observations in given 
//  dish
//  d1, d2 -- indices of dishes to calculate likelihood for
//  lvlhood -- log likelihood of observations in dishes d1, d2
//
bool HDPsampler::MHLogLikelihood( int d1, int d2, double* lvlhood )
{
    mystring errstr;
    Dish*   dish1 = NULL;
    Dish*   dish2 = NULL;
    Dish*   dsh = NULL;
    const Pslvector* vec;
    double  lprob;
    int     novs;
    int n, bn;

    if( GetMenu() == NULL )
        throw myruntime_error( "MHLogLikelihood: Null menu." );

    if( d1 < 0 || GetMenu()->GetSize() <= d1 ||
        d2 < 0 || GetMenu()->GetSize() <= d2 )
        throw myruntime_error( "MHLogLikelihood: Invalid dish indices." );

    if(( dish1 = GetMenu()->GetDishAt( d1 )) == NULL ||
       ( dish2 = GetMenu()->GetDishAt( d2 )) == NULL )
        throw myruntime_error( "MHLogLikelihood: Null dishes at index locations." );

    if( lvlhood )
        *lvlhood = 0.0;

    novs = dish1->GetActualSize();
    if( d1 != d2 )
        novs += dish2->GetActualSize();

    //create temporary table
    Table   tmptbl( novs );
    tmptbl.SetBasin( GetBasin());
    tmptbl.SetMenu( GetMenu());

    try {
        //soft copy vectors to table
        for( novs = 0, dsh = dish1;; dsh = dish2 ) {
            novs += dsh->GetActualSize();
            for( n = 0; n < dsh->GetSize(); n++ ) {
                //include all vectors
                bn = dsh->GetVectorNIndAt( n );
                if( bn < 0 )
                    continue;
                vec = dsh->GetVectorNAt( n );
                if( vec == NULL )
                    throw myruntime_error( "MHLogLikelihood: Null vector." );

                //dish index of vector 0 won't be in use in tmp. table
                tmptbl.NewVectorNInd( bn, 0 );
            }
            if( dish1 == dish2 || dsh == dish2 )
                //if two dishes are the same or both have been processed
                break;
        }

        if( novs != tmptbl.GetActualSize())
            throw myruntime_error( "MHLogLikelihood: Inconsistent no. vectors in dishes." );

        //prior probability of all vectors
        PriorProbMtx( &tmptbl, &lprob );

        if( lvlhood )
            *lvlhood = lprob;

    } catch( myexception const& ex ) {
        errstr = ex.what();
    }

    if( !errstr.empty())
        throw myruntime_error( errstr );

    return true;
}





// =========================================================================
// MHSaveSplitMergeInfo: buffer listing of vectors obtained from the 
//  split-merge procedure
//  d1, d2 -- indices of dishes having vectors v1 and v2
//  v1, v2 -- vector indices in dishes d1 and d2, respectively
//  dn1, dn2 -- indices of new auxiliary dishes
//  vn1, vn2 -- new vector indices in auxiliary dishes
//
bool HDPsampler::MHSaveSplitMergeInfo( int d1, int v1, int d2, int v2,
                                       int dn1, int dn2, int vn1, int vn2 )
{
    Dish*   dish1 = NULL; //dish 1
    Dish*   dish2 = NULL; //dish 2
    Dish*   mhdish1 = NULL; //auxiliary dish 1
    Dish*   mhdish2 = NULL; //auxiliary dish 2
    Dish*   dsh = NULL;
    double* values = GetResMBufferValues();
    int n, bn;
    int pp = 0;

    if( GetMenu() == NULL )
        throw myruntime_error( "MHSaveSplitMergeInfo: Null menu." );

    if( values == NULL )
        throw myruntime_error("MHSaveSplitMergeInfo: Null buffer values.");

    if( d1 < 0 || GetMenu()->GetSize() <= d1 ||
        d2 < 0 || GetMenu()->GetSize() <= d2 )
        throw myruntime_error( "MHSaveSplitMergeInfo: Invalid dish indices." );

    if(( dish1 = GetMenu()->GetDishAt( d1 )) == NULL ||
       ( dish2 = GetMenu()->GetDishAt( d2 )) == NULL )
        throw myruntime_error( "MHSaveSplitMergeInfo: Null dishes at index locations." );

    //mhdish1 and mhdish2 contain final split state if d1==d2;
    //otherwise, mhdish1 and mhdish2 have been split to the 
    // original split state d1 and d2.
    //Hence, it is required to move vectors from, say, mhdish2
    // either to new dish (case d1==d2), or to dish d1 
    // (since mhdish2 holds the same set of vectors as the original dish2)

    if( d1 == d2 ) {
        //to split state;
        //launch state dishes are required just in split proposal
        if( dn1 < 0 || GetMenu()->GetSize() <= dn1 ||
            dn2 < 0 || GetMenu()->GetSize() <= dn2 )
            throw myruntime_error( "MHSaveSplitMergeInfo: Invalid auxiliary dish indices." );

        if(( mhdish1 = GetMenu()->GetDishAt( dn1 )) == NULL ||
          ( mhdish2 = GetMenu()->GetDishAt( dn2 )) == NULL )
            throw myruntime_error( "MHSaveSplitMergeInfo: Null auxiliary dishes." );

        values[pp++] = ( double )d2;//from common dish (d2==d1)
        values[pp++] = ( double )-1;//to new dish
        dsh = mhdish2;
    }
    else {
        //to merge state
        values[pp++] = ( double )d2;//from dish d2
        values[pp++] = ( double )d1;//to dish d1
        dsh = dish2;
    }

    for( n = 0; n < dsh->GetSize(); n++ ) {
        //include all vectors
        bn = dsh->GetVectorNIndAt( n );
        if( bn < 0 )
            continue;
        values[pp++] = ( double )bn;//dish index
    }

    SetResMBufferNoVals( pp );
    return true;
}
