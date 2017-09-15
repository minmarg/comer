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


// Mettropolis-Hastings update calculations
// Implementation of a simplified split-merge MCMC procedure
//
// =========================================================================
// MHupdate: Mettropolis-Hastings update
//  d1, d2 -- dish indices of vectors v1 and v2
//  v1, v2 -- vector indices in dishes d1 and d2, respectively
//
bool HDPsampler::MHupdateS( int d1, int v1, int d2, int v2 )
{
    mystring errstr;
    Dish*   dsh = NULL;
    double  llhrat;//log ratio of likelihoods
    double  accprob;//M-H acceptance probability
    double  urn = 0.0;
    int dn1 = -1, dn2 = -1;//indices of new auxiliary dishes
    int vn1 = -1, vn2 = -1;//new indices of vectors in auxiliary dishes
    bool    accept = false;
    bool    rok = true;

    try {
        if( d1 == d2 ) {
            //launch state is required just for split proposal
            if( GetNoRstrGSScans() < 0 ) {
                if( rok )
                    rok = MHInitLaunchStateThroughSampling( d1, v1, d2, v2, &dn1, &dn2, &vn1, &vn2 );
            } else {
                if( rok )
                    rok = MHInitLaunchState( d1, v1, d2, v2, &dn1, &dn2, &vn1, &vn2 );
                if( rok )
                    rok = MHIntmRstrGSScans( dn1, dn2, vn1, vn2 );
            }
        }
        if( rok )
            rok = MHLikelihoodRatio( d1, v1, d2, v2, dn1, dn2, vn1, vn2, &llhrat );

        if( rok ) {
            accprob = llhrat;
            if( accprob < cdp_LOG_DBL_MIN )
                accprob = 0.0;
            else if( accprob < cdp_LOG_DBL_MAX )
                accprob = exp( accprob );
            urn = RNGja.GetDouble();

// fprintf(stderr,"\n\n d1=%d(%d) v1=%d d2=%d(%d) v2=%d\n",d1,GetMenu()->GetDishAt(d1)->GetActualSize(),v1,d2,GetMenu()->GetDishAt(d2)->GetActualSize(),v2);
// GetMenu()->GetDishAt(d1)->GetVectorNAt(v1)->Print(stderr," %12.6g");
// GetMenu()->GetDishAt(d2)->GetVectorNAt(v2)->Print(stderr," %12.6g");
// fprintf(stderr," accprob=%g(%g) urn=%g\n",accprob,llhrat,urn);

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
