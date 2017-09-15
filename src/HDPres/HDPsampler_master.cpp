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
#include <time.h>

#include "data.h"
#include "ext/psl.h"
#include "ext/gamma.h"
#include "ext/digamma.h"
#include "ext/rv/rvdisc.h"
#include "smpl/cars.h"
#include "smpl/slicesmpl.h"
#include "root/rcodes.h"
#include "root/root.h"
#include "HDPsampler.h"

//random number generators
MTRng   RNGi;
MTRng   RNGc;

// =========================================================================
// MasterProcess: master process
//
bool HDPsampler::MasterProcess()
{
    if( !mMPIIAmMaster())
        return false;

    if( mMPIMyRank() < 0 || mMPIRingSize() < 1 )
        throw myruntime_error( "MasterProcess: Invalid MPI rank or ring size." );

    if( GetBasin() == NULL || GetMenu() == NULL || GetChain() == NULL )
        return false;

    double      maxlprob = GetMaxLProbData();
    int         itrsread = GetIterationRead();
    int         maxiters = GetNoIterations();
    char        strbuf[BUF_MAX];
    const int   nopus = mMPIRingSize() - 1;//number of processing units
    int         norests = GetChain()->GetActualSize();
    int         notbls;
    int         noups, err;
    int         ms, n, nn, j;
    int         it = 0;

    time_t          t_tt;
    struct tm       t_ttm;
    char            t_ttma[KBYTE];

    if( 0 <= maxlprob )
        maxlprob = -1.e9;
    if( 0 < itrsread ) {
        it = itrsread + 1;
        maxiters += itrsread + 1;
    }

    sprintf( strbuf, "%2cNo.rests.=%d No.dishes=%d No.vectors=%d", 32,
             norests, GetMenu()->GetActualSize(), GetBasin()->GetActualSize());
    MasterMessage( strbuf );
    SetRestarted( true );


// Table tb1( GetDefTableSize()), tb2( GetDefTableSize());
// Dish* dish1 = new Dish( GetDefDishSize());
// Dish* dish2 = new Dish( GetDefDishSize());
// tb1.SetBasin( GetBasin());
// tb2.SetBasin( GetBasin());
// dish1->SetBasin( GetBasin());
// dish2->SetBasin( GetBasin());
// int dd1 = GetMenu()->NewTmpDish( dish1 );
// int dd2 = GetMenu()->NewTmpDish( dish2 );
// double lprob1, lprob2;
// const Pslvector* vec1, *vec2;
// vec1 = GetBasin()->GetValueAt( 0 );
// for( int b = 0; b < GetBasin()->GetSize(); b++ ) {
//     vec2 = GetBasin()->GetValueAt( b );
//     tb1.NewVectorNInd( 0, 0 );
//     tb2.NewVectorNInd( b, 0 );
// 
// //     if( dish1->GetSize() <= 0 ) PriorProbVec( vec1, &lprob1 ); else ProbVecOfDish( vec1, dd1, &lprob1 );
// //     if( dish2->GetSize() <= 0 ) PriorProbVec( vec2, &lprob2 ); else ProbVecOfDish( vec2, dd2, &lprob2 );
// 
// //     if( dish1->GetSize() <= 0 ) PriorProbMtx( &tb1, &lprob1 ); else ProbMtxOfDish( &tb1, dd1, &lprob1 );
// //     if( dish2->GetSize() <= 0 ) PriorProbMtx( &tb2, &lprob2 ); else ProbMtxOfDish( &tb2, dd2, &lprob2 );
// 
//     PriorProbMtx( &tb1, &lprob1 );
//     PriorProbMtx( &tb2, &lprob2 );
// 
//     fprintf(stderr,"dish1: size=%d lprob1=%g; dish2: size=%d lprob2=%g\n",
//         dish1->GetActualSize(),lprob1,dish2->GetActualSize(),lprob2);
//     AdjustDishParams( dd1, vec1, true/*add*/); dish1->NewVectorNInd( 0 );
//     AdjustDishParams( dd2, vec2, true/*add*/); dish2->NewVectorNInd( b );
//     ;
// }
// GetMenu()->RemTmpDishAt( dd2, dish2 );
// GetMenu()->RemTmpDishAt( dd1, dish1 );
// return true;




    for( ; it < maxiters ||( !GetNoIterations() && it == maxiters ); it++ )
    {
        time( &t_tt );
        localtime_r( &t_tt, &t_ttm );
        asctime_r( &t_ttm, t_ttma );
        if( strlen(t_ttma) < KBYTE && 1 < strlen(t_ttma))
            t_ttma[strlen(t_ttma)-1] = 0;

        sprintf( strbuf, "%2cITERATION %d (%s) %s", 32, it, t_ttma, GetRestarted()? "restart": "");
        MasterMessage( strbuf );
        SetGibbsIt( it );

        if( it && ( GetNoMHUpdates() <= 0 || !GetModJNMCMC()))
            //HYPERPARAMETERS
            if( !MasterSampleHyperparameters())
                return false;

        //{{MCMC PROCEDURE
        for( ms = 0; ms < GetNoMHUpdates(); ms++ ) {
            if( it || ms )
                if( GetModJNMCMC())
                    //HYPERPARAMETERS are required only in modified JN MCMC sampling
                    if( !MasterSampleHyperparameters())
                        return false;

            SetResType( TRT_MHUpdate );
            SetMHUpdateNum( ms );
            sprintf( strbuf, "%4cM-H upd. %d (merge rej.`-' acc.`+'; split rej.`=' acc.`*')", 32, ms );
            MasterMessage( strbuf );
            if( !MasterDistMCMCindices( &noups ))//broadcasts indices
                return false;
            sprintf( strbuf, "%8cExpected simultaneous tests for updates = %d (%d dish(es))", 32, 
                    noups, GetMenu()->GetActualSize());
            MasterMessage( strbuf );
            if( !MasterCmdMCMCaccpd())//broadcasts complete listing
                return false;
            if( !ProcessCmdMCMCcomplete( &noups ))
                return false;
            sprintf( strbuf, "%8cM-H updates accomplished = %d (%d dish(es))", 32, 
                    noups, GetMenu()->GetActualSize());
            MasterMessage( strbuf );
//             if(( err = BlockUntilAllReady( false )) != MOK ) {
//                 sprintf( strbuf, "MasterProcess: Barrier failed: Error code %d.", err );
//                 error( strbuf );
//                 return false;
//             }
            SaveBestAndLast( &maxlprob, it, maxiters <= it );
            SetRestarted( false );
        }
        //}}

        //{{GIBBS SCAN
        if( it < maxiters ) {
            SetResType( TRT_GibbsUpdate );
            sprintf( strbuf, "%4cGIBBS SAMPLING SCAN", 32 );
            MasterMessage( strbuf );
            sprintf( strbuf, "%4cProcessing VECTORS", 32 );
            MasterMessage( strbuf );

            if( GetParallelProc() == 3 ) {
                GetMenu()->CalcNoTables();
                n = 0;
                nn = norests - 1;
                sprintf( strbuf, "%6cProcessing restaurants from %d to %d (%d dish(es), %d table(s))",
                          32, n, nn, GetMenu()->GetActualSize(), GetMenu()->GetNoTables());
                MasterMessage( strbuf );
                if( !MasterProcessRestVectors( n, nn ))
                    return false;
            }
            else {
                for( n = 0, nn = nopus - 1; n < norests; n += nopus, nn += nopus ) {
                    if( norests <= nn )
                        nn = norests - 1;

                    for( j = n, notbls = 0; j <= nn; j++ ) 
                        notbls += GetChain()->GetRestaurantAt( j )->GetActualSize();
                    sprintf( strbuf, "%6cProcessing restaurants from %d to %d (%d dish(es), %d table(s))",
                              32, n, nn, GetMenu()->GetActualSize(), notbls );
                    MasterMessage( strbuf );
                    if( !MasterProcessRestVectors( n, nn ))
                        return false;
                }
            }

            if( !GetNoTableSampling()) {
                sprintf( strbuf, "%4cProcessing TABLES", 32 );
                MasterMessage( strbuf );

                if( GetParallelProc() == 3 ) {
                    GetMenu()->CalcNoTables();
                    n = 0;
                    nn = norests - 1;
                    sprintf( strbuf, "%6cProcessing restaurants from %d to %d (%d dish(es), %d table(s))",
                              32, n, nn, GetMenu()->GetActualSize(), GetMenu()->GetNoTables());
                    MasterMessage( strbuf );
                    if( !MasterProcessRestTables( n, nn ))
                        return false;
                }
                else {
                    for( n = 0, nn = nopus - 1; n < norests; n += nopus, nn += nopus ) {
                        if( norests <= nn )
                            nn = norests - 1;

                        for( j = n, notbls = 0; j <= nn; j++ ) 
                            notbls += GetChain()->GetRestaurantAt( j )->GetActualSize();
                        sprintf( strbuf, "%6cProcessing restaurants from %d to %d (%d dish(es), %d table(s))",
                                  32, n, nn, GetMenu()->GetActualSize(), notbls );
                        MasterMessage( strbuf );
                        if( !MasterProcessRestTables( n, nn ))
                            return false;
                    }
                }
            }
            SaveBestAndLast( &maxlprob, it, ((it+1)%10 == 0 )||( maxiters <= it+1 ));
            SetRestarted( false );
        }
        //}}

        if( GetRestarted())
            SaveBestAndLast( &maxlprob, it, true );

        //{{PARALLEL TeST over all computing nodes
        if(( (it+1)%10 == 0 )||( maxiters <= it+1 )) {
            sprintf( strbuf, "%4cPARALLEL TEST...", 32 );
            MasterMessage( strbuf );
            if( !MasterCmdInitTest())//initialize testing
                return false;
            if( !MasterCmdValsTest())//compare received data with master's data
                return false;
            sprintf( strbuf, "%6cdone.", 32 );
            MasterMessage( strbuf );
        }
        //}}
    }

    return true;
}

// =========================================================================
// SaveSaveBestAndLast: save best and last configuration
//
bool HDPsampler::MasterSampleHyperparameters()
{
    char        strbuf[BUF_MAX];
    bool        sampletau = GetTauPreset() <= 0.0;
    bool        samplegamma = GetGammaPreset() <= 0.0;

    if( !mMPIIAmMaster())
        return false;
    if( GetBasin() == NULL || GetMenu() == NULL || GetChain() == NULL )
        return false;

    GetMenu()->CalcNoTables();
    GetChain()->CalcNoVectors();

    if( GetKappa0Preset() <= 0.0 ) {
        //kappa0
        sprintf( strbuf, "%4cSampling kappa0", 32 );
        MasterMessage( strbuf );
        if( !MasterprocessSampleKappa0())
            return false;
        sprintf( strbuf, "%6cKAPPA0=%.6g (no.dishes=%d)" , 32, 
                GetMenu()->GetKappa0(), GetMenu()->GetActualSize());
        MasterMessage( strbuf );
        sprintf( strbuf, "%6cDelivering and adjusting according to new value of kappa0", 32 );
        MasterMessage( strbuf );
        if( !MasterprocessBcastKappa0())
            return false;
    }

    if( GetNu0Preset() <= 0.0 ) {
        //nu0
        sprintf( strbuf, "%4cSampling nu0", 32 );
        MasterMessage( strbuf );
        if( !MasterprocessSampleNu0())
            return false;
        sprintf( strbuf, "%6cNU0=%.6g (no.dishes=%d)" , 32, 
                GetMenu()->GetNu0(), GetMenu()->GetActualSize());
        MasterMessage( strbuf );
        sprintf( strbuf, "%6cDelivering new value of nu0", 32 );
        MasterMessage( strbuf );
        if( !MasterprocessBcastNu0())
            return false;
    }

    //concentration parameters
    if( sampletau || samplegamma ) {
        sprintf( strbuf, "%4cSampling concentration parameters", 32 );
        MasterMessage( strbuf );
    }
    if( sampletau ) {
        sprintf( strbuf, "%6cSampling tau", 32 );
        MasterMessage( strbuf );
        if( !MasterprocessSampleTau())
            return false;
        sprintf( strbuf, "%6cTAU=%.6g (no.tables=%d, no.rests.=%d)" , 32, 
                GetDPMTau(), GetMenu()->GetNoTables(), GetChain()->GetActualSize());
        MasterMessage( strbuf );
    }

    if( samplegamma ) {
        sprintf( strbuf, "%6cSampling gamma", 32 );
        MasterMessage( strbuf );
        if( !MasterprocessSampleGamma())
            return false;
        sprintf( strbuf, "%6cGAMMA=%.6g (no.tables=%d, no.dishes=%d)" , 32, 
                GetDPMGamma(), GetMenu()->GetNoTables(), GetMenu()->GetActualSize());
        MasterMessage( strbuf );
    }

    if( sampletau || samplegamma ) {
        sprintf( strbuf, "%6cDelivering concentration parameters", 32 );
        MasterMessage( strbuf );
        if( !MasterprocessBcastConcParams())
            return false;
    }

#ifdef HDPTESTINL
    if( !TestMenu())
        return false;
#endif
    return true;
}





// =========================================================================
// SaveSaveBestAndLast: save best and last configuration
//
void HDPsampler::SaveBestAndLast( double* maxlprob, int iter, bool last )
{
    if( !GetOutputFile())
        return;

    double      lprob;
    char        strbuf[BUF_MAX];
    mystring    fbest = GetOutputFile();
    mystring    flast = GetOutputFile();

    fbest += _sch_ext_best_;
    flast += _sch_ext_last_;

    sprintf( strbuf, "%4cSaving current configuration...", 32 );
    MasterMessage( strbuf );

    CalcProbOfData( &lprob );
    sprintf( strbuf, "%4cLog Probability of Data = %g", 32, lprob );
    MasterMessage( strbuf );
    PrintSummary(( flast + _sch_ext_sum_ ).c_str(), &lprob, &iter );
    if( last )
        Save( flast, &lprob, &iter );

    if( maxlprob && *maxlprob < lprob ) {
        *maxlprob = lprob;
        SetMaxLProbData( lprob );
        Save( fbest, &lprob, &iter );
    }
}

// -------------------------------------------------------------------------
// Save: save current configuration
//
void HDPsampler::Save( mystring& fname, double* lprob, int* iter )
{
    PrintGroups(( fname + _sch_ext_grp_ ).c_str(), lprob, iter );
    PrintDishes(( fname + _sch_ext_dsh_ ).c_str(), lprob, iter );
    PrintParameters(( fname + _sch_ext_par_ ).c_str());
}





// =========================================================================
// MasterDistMCMCindices: generate and broadcast MCMC indices of dishes
//
bool HDPsampler::MasterDistMCMCindices( int* nupdates )
{
    if( !mMPIIAmMaster())
        return false;
    if( GetBasin() == NULL || GetMenu() == NULL )
        return false;

    const int   nopus = mMPIRingSize() - 1;//number of processing units
    int szbs = GetBasin()->GetSize();//size of basin
    int szmn = GetMenu()->GetSize();//size of menu
    int nods = GetMenu()->GetActualSize();//number of dishes
    int noups = SLC_MIN( nopus, nods );//number of MCMC updates in single broadcast

    const double pacc = 1.e-4;
    double* data = GetMBufferValues();
    double* dprobs = NULL;
    double  prb, consv;
    Dish*   dsh;
    RVDisc  rvdis( RNGi, RVDisc::TRVDsc_Inverse );
    int     szad;//actual dish size
    int v, vv;//vector indices in dishes
    int d, dd;//dish indices
    int pp, u, n;
    int len, code;
    mystring preamb = "MasterDistMCMCindices: ";
    mystring merror;

    if( data == NULL )
        throw myruntime_error( preamb + "Null data buffer.");
    if( nopus < 1 )
        throw myruntime_error( preamb + "Invalid MPI ring size.");

    dprobs = ( double* )malloc( szmn * sizeof( double ));
    if( dprobs == NULL )
        throw myruntime_error( preamb + "Not enough memory.");

    rvdis.SetProbs( dprobs, szmn );

    for( d = 0, consv = 0.0; d < szmn; d++ ) {
        dprobs[d] = 0.0;
        dsh = GetMenu()->GetDishAt( d );
        if( dsh == NULL )
            continue;
        dsh->SetProcessed( false );
        szad = dsh->GetActualSize();
        if( GetMHSampleDishUnf())
            dprobs[d] = 1.0 /( double )nods;
        else
            dprobs[d] = ( double )szad /( double )szbs;
        consv += dprobs[d];
    }

    try {
        if( consv < 1.0 - pacc || 1.0 + pacc < consv )
            throw myruntime_error( preamb + "Probabilities not conserved.");

        for( pp = 0, u = 0; u < noups; ) {
            //generate first dish's index
            if(( code = rvdis.GenI( &d )) != PSL_OK )
                throw myruntime_error( preamb + TranslatePSLError( code ));

            if( szmn <= d )
                break;

            dsh = GetMenu()->GetDishAt( d );
            if( dsh == NULL || dsh->GetProcessed())
                throw myruntime_error( preamb + "Dish 1 already processed.");
            dsh->SetProcessed( true );
            szad = dsh->GetActualSize();

            //generate second dish's index
            if(( code = rvdis.GenI( &dd )) != PSL_OK )
                throw myruntime_error( preamb + TranslatePSLError( code ));

            if( szmn <= dd || GetSplitProposalsOnly())
                dd = d;

            if( dd == d && szad < 2 && !GetSplitProposalsOnly()) {
                if( dprobs[d] < 0.51 /*&& !u */) {
                    dsh->SetProcessed( false );
                    continue;
                }
                break;
            }

            if( dd != d ) {
                dsh = GetMenu()->GetDishAt( dd );
                if( dsh == NULL || dsh->GetProcessed())
                    throw myruntime_error( preamb + "Dish 2 already processed.");
                dsh->SetProcessed( true );
            }

            if( GetMergeProposalsOnly() && dd == d ) {
                dsh->SetProcessed( false );
                continue;
            }

            if( dd != d || 2 <= szad ) {
                //generate indices of vectors from the dishes
                if( !MasterGenVecIndex( d, -1, &v ))
                    throw myruntime_error( preamb + "Failed to generate 1st vector index.");
                if( !MasterGenVecIndex( dd, ( d == dd )? v: -1, &vv ))
                    throw myruntime_error( preamb + "Failed to generate 2nd vector index.");
                ;

                data[pp++] = ( double )d;
                data[pp++] = ( double )v;
                data[pp++] = ( double )dd;
                data[pp++] = ( double )vv;
                u++;
            }
            //adjust properly probabilities
            prb = 1.0 - dprobs[d];
            dprobs[d] = 0.0;
            if( dd != d ) {
                prb -= dprobs[dd];
                dprobs[dd] = 0.0;
            }

            if( prb <= 0.0 )
                break;

            for( n = 0; n < szmn; n++ )
                if( dprobs[n])
                    dprobs[n] /= prb;
        }
    } catch( myexception const& ex ) {
        merror = ex.what();
    }

    free( dprobs );
    dprobs = NULL;

    if( !merror.empty()) {
        throw myruntime_error( merror );
        return false;
    }

    if( nupdates )
        *nupdates = pp >> 2;//divide by 4

//     if( pp <= 0 )
//         return true;

    SetMBufferHead( THeadDat );
    SetMBufferCmd( TComInitMCMC );
    SetMBufferNoVals( pp );
    len = SetMBufferCRC();
    if( GetMaxSizeOfResMBuffer() <= len || GetMaxSizeOfDataMessage() <= len ) {
        error("MasterDistMCMCindices: Short of size in buffers. Increase it and recompile.");
        return false;
    }
//     if(( code = BcastMPIMessage( GetMBuffer(), /*len*/GetMaxSizeOfDataMessage(), 
//                   false/*throw*/)) != MOK )
//         return false;
    if(( code = BcastlenMPIMessage( GetMBuffer(), len, GetMaxSizeOfDataMessage(), 
                  false/*throw*/)) != MOK )
        return false;
 
    return true;
}

// =========================================================================
// MasterDistMCMCindices: generate and broadcast MCMC indices of dishes
//
bool HDPsampler::MasterDistMCMCindicesObs( int* nupdates )
{
    if( !mMPIIAmMaster())
        return false;
    if( GetBasin() == NULL || GetMenu() == NULL )
        return false;

    const int   nopus = mMPIRingSize() - 1;//number of processing units
    int szbs = GetBasin()->GetSize();//size of basin
    int szmn = GetMenu()->GetSize();//size of menu
    int nods = GetMenu()->GetActualSize();//number of dishes
    int noups = SLC_MIN( nopus, nods );//number of MCMC updates in single broadcast

    double* data = GetMBufferValues();
    Dish*   dsh;
    double  urv;
    bool    round = false;
    int     szad, szd;//dish actual size and size
    int v, vv;//vector indices in dishes
    int d, dd;//dish indices
    int pp, len, err, u;

    if( data == NULL )
        throw myruntime_error("MasterDistMCMCindices: Null address of data buffer.");
    if( nopus < 1 )
        throw myruntime_error( "MasterDistMCMCindices: Invalid MPI ring size." );

    for( d = 0; d < szmn; d++ ) {
        dsh = GetMenu()->GetDishAt( d );
        if( dsh == NULL )
            continue;
        dsh->SetProcessed( false );
    }

    //generate index of the first dish
    urv = RNGi.GetDouble();
    dd = d = ( int )rint( urv *( double )( szmn-1 ));
    if( d <= 0 )
        round = true;

    for( pp = 0, u = 0; u < noups && ( !round || d < szmn ); dd = ++d ) {
        if( !round && szmn <= d ) {
            round = true;
            dd = d = 0;
        }
        dsh = GetMenu()->GetDishAt( d );
        if( dsh == NULL || dsh->GetProcessed())
            continue;
        dsh->SetProcessed( true );
        szad = dsh->GetActualSize();
        szd = dsh->GetSize();
        //generate index of the first vector in dish d
        if( !MasterGenVecIndex( d, -1, &v ))
            throw myruntime_error("MasterDistMCMCindices: Failed to generate vector index.");

        //dicide whether the same dish for the second index will be used in MCMC
        urv = RNGi.GetDouble();
        if( szbs < szad )
            throw myruntime_error("MasterDistMCMCindices: Dish size greater than basin's one.");
        if(( double )szad /( double )szbs < urv )
        {
            //generate index of the second distinct dish
            urv = RNGi.GetDouble();
            dd = d + 1 + ( int )rint( urv *( double )( szmn-2-d ));
            for( ; dd < GetMenu()->GetSize() && 
              (!( dsh = GetMenu()->GetDishAt( dd )) || dsh->GetProcessed()); dd++ );
            if( dd < GetMenu()->GetSize())
                dsh->SetProcessed( true );
            else
                dd = d;
        }
        //reduce basin size by the sizes of the processed dish and next processed dish
        szbs -= szad;
        if( d != dd && dd < GetMenu()->GetSize())
            szbs -= dsh->GetActualSize();
        //generate index of the second vector in dish dd
        if( !MasterGenVecIndex( dd, ( d == dd )? v: -1, &vv )) {
            if( d == dd )
                continue;
            else
                throw myruntime_error("MasterDistMCMCindices: Failed to generate 2nd index.");
        }
        data[pp++] = ( double )d;
        data[pp++] = ( double )v;
        data[pp++] = ( double )dd;
        data[pp++] = ( double )vv;
        u++;
    }

    if( nupdates )
        *nupdates = pp >> 2;//divide by 4

//     if( pp <= 0 )
//         return true;

    SetMBufferHead( THeadDat );
    SetMBufferCmd( TComInitMCMC );
    SetMBufferNoVals( pp );
    len = SetMBufferCRC();
    if( GetMaxSizeOfResMBuffer() <= len || GetMaxSizeOfDataMessage() <= len ) {
        error("MasterDistMCMCindices: Short of size in buffers. Increase it and recompile.");
        return false;
    }
//     if(( err = BcastMPIMessage( GetMBuffer(), /*len*/GetMaxSizeOfDataMessage(), 
//                   false/*throw*/)) != MOK )
//         return false;
    if(( err = BcastlenMPIMessage( GetMBuffer(), len, GetMaxSizeOfDataMessage(), 
                  false/*throw*/)) != MOK )
        return false;
 
    return true;
}

// -------------------------------------------------------------------------
// MasterGenDishIndex: generate vector index in dish given by index
//  dd -- dish index
//  vv -- vector index in dish
//  vvv -- vector index to be returned
//
bool HDPsampler::MasterGenVecIndex( int dd, int vv, int* vvv )
{
    if( vv < 0 || GetMHSampleVectUnf())
        return MasterGenVecIndexS( dd, vv, vvv );
    return MasterGenVecIndexD( dd, vv, vvv );
}

// -------------------------------------------------------------------------
// MasterGenDishIndexS: generate vector index in dish given by index 
//  (simple algorithm)
//  dd -- dish index
//  vv -- vector index in dish
//  vvv -- vector index to be returned
//
bool HDPsampler::MasterGenVecIndexS( int dd, int vv, int* vvv )
{
    mystring preamb = "HDPsampler: MasterGenVecIndexS: ";
    Dish*   dsh;
    double  urv = RNGi.GetDouble();
    int     szd;//dish size
    int     pv, v;

    if( vvv == NULL )
        return false;
    if( GetBasin() == NULL || GetMenu() == NULL )
        return false;
    if( dd < 0 || GetMenu()->GetSize() <= dd )
        throw myruntime_error( preamb + "Invalid dish index.");
    if(( dsh = GetMenu()->GetDishAt( dd )) == NULL )
        throw myruntime_error( preamb + "Null dish at given location.");

    szd = dsh->GetSize();
    pv = ( int )rint( urv *( double )( szd - 1 ));
    for( v = pv; v < szd &&( dsh->GetVectorNIndAt( v ) < 0 || v == vv ); v++ );
    if( szd <= v ) {
        for( v = 0; v < szd && v < pv &&( dsh->GetVectorNIndAt( v ) < 0 || v == vv ); v++ );
        if( szd <= v || pv <= v )
            return false;
    }
    *vvv = v;
    return true;
}

// -------------------------------------------------------------------------
// MasterGenDishIndexD: generate vector index in dish given by index 
//  (euclidean-distance algorithm)
//  dd -- dish index
//  vv -- vector index in dish
//  vvv -- vector index to be returned
//
bool HDPsampler::MasterGenVecIndexD( int dd, int vv, int* vvv )
{
    mystring    preamb = "HDPsampler: MasterGenVecIndexD: ";
    //buffer data is used at this stage
    double*     values = GetResMBufferValues();
    if( values == NULL )
        throw myruntime_error( preamb + "Null buffer.");

    if( vv < 0 )
        return false;

    RVDisc      rvdsc( RNGi, RVDisc::TRVDsc_Inverse );
    Dish*       dsh;
    Pslvector*  vec, vvec, vdlt;
    double  sum, val;
    int     szd, aszd;//dish size
    int     v, vndx;
    int     ecd;

    if( vvv == NULL )
        throw myruntime_error( preamb + "Null vector address.");
    if( GetBasin() == NULL || GetMenu() == NULL )
        throw myruntime_error( preamb + "Null menu structures.");
    if( dd < 0 || GetMenu()->GetSize() <= dd )
        throw myruntime_error( preamb + "Invalid dish index.");
    if(( dsh = GetMenu()->GetDishAt( dd )) == NULL )
        throw myruntime_error( preamb + "Null dish at given location.");
    if( dsh->GetVectorNIndAt(vv) < 0 )
        throw myruntime_error( preamb + "Invalid vector index.");
    if( dsh->GetActualSize() < 2 )
        throw myruntime_error( preamb + "Invalid dish size.");
    if(( vec = dsh->GetVectorNAt(vv)) == NULL )
        throw myruntime_error( preamb + "Null support vector.");

    sum = 0.0;
    vvec = *vec;
    szd = dsh->GetSize();
    aszd = dsh->GetActualSize();
    for( v = 0; v < szd; v++ ) {
        vndx = dsh->GetVectorNIndAt(v);
        if( vndx < 0 || v == vv ) {
            values[v] = 0.0;
            continue;
        }
        vdlt = vvec;
        vec = dsh->GetVectorNAt(v);
        if( vec == NULL )
            throw myruntime_error( preamb + "Null vector.");
        if(( ecd = vdlt.Superposition( -1.0, *vec )) != 0 )
            throw myruntime_error( preamb + TranslatePSLError(ecd));
        values[v] = vdlt.Norm2();
        sum += values[v];
    }

    //normalize
    if( sum <= 0.0 ) {
        //all vectors in dish are equal
        sum = 0.0;
        val = 1.0 /(double)(aszd-1);
        for( v = 0; v < szd; v++ ) {
            vndx = dsh->GetVectorNIndAt(v);
            if( vndx < 0 || v == vv )
                continue;
            values[v] = val;
            sum += val;
        }
    }
    else {
        for( v = 0; v < szd; v++ )
            if( values[v])
                values[v] /= sum;
    }
    if( sum <= 0.0 )
        throw myruntime_error( preamb + "Zero normalizing sum.");

    rvdsc.SetProbs( values, szd );

    //generate index of second support vector 
    if(( ecd = rvdsc.GenI( vvv )) != PSL_OK )
        throw myruntime_error( preamb + TranslatePSLError(ecd));

    if( szd <= *vvv )
        return false;

    return true;
}





// =========================================================================
// MasterCmdMCMCaccpd: process command sent by slaves, related to MH 
//  updates
//
bool HDPsampler::MasterCmdMCMCaccpd()
{
    if( !mMPIIAmMaster())
        return false;

    char        strbuf[BUF_MAX];
    const int   nopus = mMPIRingSize() - 1;//number of processing units
    double*     data = GetMBufferValues();
    double*     values = GetResMBufferValues();
    mystring    merror;
    int head, cmd, novals;
    int len, err;
    int p, pp, n;

    if( mMPIMyRank() < 0 )
        throw myruntime_error( "MasterCmdMCMCaccpd: Invalid MPI rank." );
    if( nopus < 1 )
        throw myruntime_error( "MasterCmdMCMCaccpd: Invalid ring size." );
    if( data == NULL || values == NULL )
        throw myruntime_error("MasterCmdMCMCaccpd: Null buffers.");

    //wait for data (dish and vector indices)
    for( p = 0, pp = 0; p < nopus; p++ ) {
        //all processes should reply
        err = ReceiveMPIMessage( GetResMBuffer(), GetMaxSizeOfResMBuffer(), 
                                  false/*throw*/);
        if( err != MOK ) {
            sprintf( strbuf, "MasterCmdMCMCaccpd: Receive failed: Error code %d.", err );
            merror = strbuf;
            continue;
        }
        if( !merror.empty())
            continue;
        try {
            head = GetResMBufferHead();
            cmd = GetResMBufferCmd();
            novals = GetResMBufferNoVals();
            if( !CheckResMBufferCRC())
                throw myruntime_error("MasterCmdMCMCaccpd: Invalid CRC.");

            if( head != THeadDat )
                throw myruntime_error("MasterCmdMCMCaccpd: Invalid header.");
            if( cmd != TComMCMCaccpd )
                throw myruntime_error("MasterCmdMCMCaccpd: Invalid command.");

            //rewrite all data to another data buffer
            data[pp++] = ( double )novals;
            for( n = 0; n < novals; n++ )
                data[pp++] = GetResMBufferValueAt( n );

            if( 1 < novals ) {
                if( 2 < novals ) {//accepted
                    if(( int )data[pp-novals+1] == -1 )
                          hdmessage( 4 );//split
                    else  hdmessage( 2 );//merge
                }
                else if(( int )data[pp-novals] == ( int )data[pp-novals+1])
                      hdmessage( 5 );//split
                else  hdmessage( 3 );//merge
            }
            else
                hdmessage( 1 );

        } catch( myexception const& ex ) {
            merror = ex.what();
            continue;
        }
    }

    hdmessage( -1 );
    message( NULL );

    if( !merror.empty()) {
        error( merror.c_str());
        return false;
    }

    //broadcast complete message of indices of vectors to be moved to other dishes
    SetMBufferHead( THeadDat );
    SetMBufferCmd( TComMCMCcomplete );
    SetMBufferNoVals( pp );
    len = SetMBufferCRC();
    if( GetMaxSizeOfResMBuffer() <= len || GetMaxSizeOfDataMessage() <= len ) {
        error("MasterCmdMCMCaccpd: Short of size in buffers. Increase it and recompile.");
        return false;
    }
//     if(( err = BcastMPIMessage( GetMBuffer(), /*len*/GetMaxSizeOfDataMessage(), 
//                   false/*throw*/)) != MOK )
//         return false;
    if(( err = BcastlenMPIMessage( GetMBuffer(), len, GetMaxSizeOfDataMessage(), 
                  false/*throw*/)) != MOK )
        return false;

    if( !merror.empty()) {
        error( merror.c_str());
        return false;
    }

    return true;
}





// =========================================================================
// MasterProcessRestVectors: process vectors at all tables in given 
//     restaurants
//
bool HDPsampler::MasterProcessRestVectors( int n, int nn )
{
    if( !mMPIIAmMaster())
        return false;

    if( mMPIMyRank() < 0 || mMPIRingSize() < 1 )
        throw myruntime_error( "MasterProcessRestVectors: Invalid MPI rank or ring size." );

    if( GetBasin() == NULL || GetMenu() == NULL || GetChain() == NULL )
        return false;

    mystring    preamb = "MasterProcessRestVectors: ";
    char        strbuf[BUF_MAX];
    const int   nopus = mMPIRingSize() - 1;//number of processing units
    int         norests = GetChain()->GetActualSize();
    double*     values = GetResMBufferValues();
    int*        rests = NULL;
    mystring    merror;
    int     len, err;
    int     p, pp, cc = 0;

    if( n < 0 || norests <= n || nn < n || norests <= nn )
        throw myruntime_error( preamb + "Invalid restaurant indices.");
    if( values == NULL )
        throw myruntime_error( preamb + "Null address of buffer.");

    rests = ( int* )malloc(( nn-n+1 )*sizeof( int ));
    if( rests == NULL )
        throw myruntime_error( preamb + "Not enough memory.");
    memset( rests, 0, ( nn-n+1 )*sizeof( int ));

    SetMBufferHead( THeadMsg );
    SetMBufferCmd( TComInitVec );
    SetMBufferNoVals( 2 );
    SetMBufferValueAt( 0, ( double )n );
    SetMBufferValueAt( 1, ( double )nn );
    len = SetMBufferCRC();
    if( GetMaxSizeOfResMBuffer() <= len || GetMaxSizeOfDataMessage() <= len ) {
        error("MasterProcessRestVectors: Short of size in buffers. Increase it and recompile.");
        return false;
    }
//     if(( err = BcastMPIMessage( GetMBuffer(), /*len*/GetMaxSizeOfDataMessage(), 
//             false/*throw*/)) != MOK )
//         return false;
    if(( err = BcastlenMPIMessage( GetMBuffer(), len, GetMaxSizeOfDataMessage(), 
            false/*throw*/)) != MOK )
        return false;

    try {
        //iterate until slaves send empty messages
        for( ;; ) {
            //wait for data (table and dish indices)
            for( p = 0, pp = 0; p < nopus; p++ ) {
                //all processes should reply
                err = ReceiveMPIMessage( GetMBuffer(), GetMaxSizeOfDataMessage(), 
                                          false/*throw*/);
                if( err != MOK ) {
                    sprintf( strbuf, "MasterProcessRestVectors: Receive failed: Error code %d.", err );
                    merror = strbuf;
                    continue;
                }
                if( !merror.empty())
                    continue;
                try {
                    MasterProcessCmdValsVec( n, nn, values, &pp, rests, ( nn-n+1 ));
                } catch( myexception const& ex ) {
                    merror = ex.what();
                    continue;
                }
            }
            if( !merror.empty()) {
                error( merror.c_str());
                break;
            }
            hdmessage( 1 );
            if( pp <= 0 ) {
                //move on to next restaurants
                hdmessage( -1 );
                message( NULL );
                break;
            }

            //broadcast message of the vectors have migrated
            SetResMBufferHead( THeadDat );
            SetResMBufferCmd( TComMigrVec );
            SetResMBufferNoVals( pp );
            len = SetResMBufferCRC();
            if( GetMaxSizeOfResMBuffer() <= len || GetMaxSizeOfDataMessage() <= len )
                throw myruntime_error( preamb + "Short of size in buffers. "
                                      "Increase it and recompile.");
//             if(( err = BcastMPIMessage( GetResMBuffer(), /*len*/GetMaxSizeOfDataMessage(), 
//                           false/*throw*/)) != MOK )
//                 throw myruntime_error( preamb + "Broadcast failed.");
            if(( err = BcastlenMPIMessage( GetResMBuffer(), len, GetMaxSizeOfDataMessage(), 
                          false/*throw*/)) != MOK )
                throw myruntime_error( preamb + "Broadcast failed.");

            if( !merror.empty()) {
                error( merror.c_str());
                break;
            }
        }
    } catch( myexception const& ex ) {
        merror = ex.what();
    }

    if( merror.empty()) {
        for( p = 0; p < nn - n; p++ )
            if( rests[p] == 0 ) {
                merror = preamb + "Not all restaurants processed.";
                break;
            }
    }
    if( rests ) { free( rests ); rests = NULL; }
    if( !merror.empty()) {
        error( merror.c_str());
        return false;
    }

    return true;
}





// =========================================================================
// MasterProcessCmdValsVec: process command of sampling from multivariate 
//  probabilities
//
void HDPsampler::MasterProcessCmdValsVec( int rb, int re, double* values, int* pp, int* rests, int rsz )
{
    if( !mMPIIAmMaster())
        return;

    if( mMPIMyRank() < 0 || mMPIRingSize() < 1 )
        throw myruntime_error( "MasterProcessCmdValsVec: Invalid MPI rank or ring size." );
    if( values == NULL )
        throw myruntime_error("MasterProcessCmdValsVec: Null address of buffer.");
    if( pp == NULL || *pp < 0  )
        throw myruntime_error( "MasterProcessCmdValsVec: Memory access error." );

    if( GetBasin() == NULL || GetMenu() == NULL || GetChain() == NULL )
        return;

    int*        dnxmap = GetLocIntBuf();//map of dish indices
    Restaurant* rest = NULL;
    Table*      tbl = NULL;
    Table*      newtbl = NULL;
    Dish*       newdsh = NULL;
    Pslvector*  vec = NULL;
    int head, cmd, novals;
    int pr, r, t, v, b, k;
    int newt, newk, svnt, svnk;
    int mm, offset;

    if( dnxmap == NULL )
        throw myruntime_error("MasterProcessCmdValsVec: Null local buffer of integers.");
    if( GetSizeLocIntBuf() < GetBasin()->GetSize())
        throw myruntime_error("MasterProcessCmdValsVec: Short of size of local buffer of integers.");
    ResetLocIntBuf( -1 );

    head = GetMBufferHead();
    cmd = GetMBufferCmd();
    novals = GetMBufferNoVals();
    if( !CheckMBufferCRC())
        throw myruntime_error("MasterProcessCmdValsVec: Invalid CRC.");

    if( head != THeadDat )
        throw myruntime_error("MasterProcessCmdValsVec: Invalid header.");
    if( cmd != TComValsVec )
        throw myruntime_error("MasterProcessCmdValsVec: Invalid command.");

    if( novals <= 0 )
        return;

    for( mm = 0, offset = 0, pr = -1; mm < novals; ) {
        newtbl = NULL;
        newdsh = NULL;
        r = ( int )GetMBufferValueAt( mm++ );//restaurant index
        t = ( int )GetMBufferValueAt( mm++ );//table index
        v = ( int )GetMBufferValueAt( mm++ );//vector index
        b = ( int )GetMBufferValueAt( mm++ );//vector index in basin
        svnt = newt = ( int )GetMBufferValueAt( mm++ );//table for vector to migrate to
        svnk = newk = ( int )GetMBufferValueAt( mm++ );//dish for vector to have at table

        if( novals < mm )
            throw myruntime_error("MasterProcessCmdValsVec: Unexpected number of values.");

        if( rest && pr != r )
            offset += rest->GetSize();
        pr = r;

        //chain stays static in process
        if( r < rb || re < r || GetChain()->GetActualSize() <= r ||
          ( rest = GetChain()->GetRestaurantAt( r )) == NULL )
            throw myruntime_error("MasterProcessCmdValsVec: Invalid restaurant index received.");

        if( GetSizeLocIntBuf() <= offset + ( rest->GetSize()))
            throw myruntime_error("MasterProcessCmdValsVec: Short of size of local buffer of integers.");

        if( rests && r-rb < rsz )
            rests[r-rb] = 1;

        if( t < 0 || rest->GetSize() <= t ||
        ( tbl = rest->GetTableAt( t )) == NULL )
            throw myruntime_error("MasterProcessCmdValsVec: Invalid table index received.");

        if(( k = tbl->GetDishIndex()) < 0 || GetMenu()->GetSize() <= k )
            throw myruntime_error("MasterProcessCmdValsVec: Invalid table's dish index.");

        if( v < 0 || tbl->GetSize() <= v ||
        ( vec = tbl->GetVectorNAt( v )) == NULL )
            throw myruntime_error("MasterProcessCmdValsVec: Invalid vector index received.");

        if( b < 0 || GetBasin()->GetSize() <= b ||
            vec != GetBasin()->GetValueAt( b ))
            throw myruntime_error("MasterProcessCmdValsVec: Inconsistent vector index.");

        if( 0 <= newt && ( rest->GetSize() <= newt ||
          ( newtbl = rest->GetTableAt( newt )) == NULL ))
            throw myruntime_error("MasterProcessCmdValsVec: Invalid table index received.");

        //slave and master dish indices may differ when creating new table with new dish
        if( 0 <= newt && 0 <= dnxmap[newt+offset])
            svnk = newk = dnxmap[newt+offset];

        if( 0 <= newk && ( GetMenu()->GetSize() <= newk ||
          ( newdsh = GetMenu()->GetDishAt( newk )) == NULL )) {
              MasterMessage( "New dish instead of vanished one." );
              newtbl = NULL;
              newdsh = NULL;
              newt = newk = -1;
              svnt = svnk = -1;
        }

        if(!( newt == t && newk == k )) {
            //values of newt and newk will be changed if they're <0
            if( !MoveVector( rest, t, k, tbl, v, vec, newt, newk, newtbl ))
                throw myruntime_error("MasterProcessCmdValsVec: Moving of vector failed.");

            values[(*pp)++] = ( double )rb;//begin index of restaurants
            values[(*pp)++] = ( double )re;//end index of restaurants
            values[(*pp)++] = ( double )r;//restaurant index
            values[(*pp)++] = ( double )t;//table index
            values[(*pp)++] = ( double )v;//vector index
            values[(*pp)++] = ( double )b;//vector index in basin
            values[(*pp)++] = ( double )svnt;//table of migration of vector
            values[(*pp)++] = ( double )svnk;//dish to have at table
        }
        //save dish index for table; table index should be the same in all processes
        dnxmap[newt+offset] = newk;
    }
}





// =========================================================================
// MasterProcessRestTables: process all tables in given restaurants
//
bool HDPsampler::MasterProcessRestTables( int n, int nn )
{
    if( !mMPIIAmMaster())
        return false;

    if( mMPIMyRank() < 0 || mMPIRingSize() < 1 )
        throw myruntime_error( "MasterProcessRestVectors: Invalid MPI rank or ring size." );

    if( GetBasin() == NULL || GetMenu() == NULL || GetChain() == NULL )
        return false;

    mystring    preamb = "MasterProcessRestTables: ";
    char        strbuf[BUF_MAX];
    const int   nopus = mMPIRingSize() - 1;//number of processing units
    int         norests = GetChain()->GetActualSize();
    double*     values = GetResMBufferValues();
    mystring    merror;
    int     len, err;
    int     p, pp, cc = 0;

    if( n < 0 || norests <= n || nn < 0 || norests <= nn )
        throw myruntime_error( preamb + "Invalid restaurant indices.");
    if( values == NULL )
        throw myruntime_error( preamb + "Null address of buffer.");

    SetMBufferHead( THeadMsg );
    SetMBufferCmd( TComInitMtx );
    SetMBufferNoVals( 2 );
    SetMBufferValueAt( 0, ( double )n );
    SetMBufferValueAt( 1, ( double )nn );
    len = SetMBufferCRC();
    if( GetMaxSizeOfResMBuffer() <= len || GetMaxSizeOfDataMessage() <= len ) {
        error("MasterProcessRestTables: Short of size in buffers. Increase it and recompile.");
        return false;
    }
//     if(( err = BcastMPIMessage( GetMBuffer(), /*len*/GetMaxSizeOfDataMessage(), 
//             false/*throw*/)) != MOK )
//         return false;
    if(( err = BcastlenMPIMessage( GetMBuffer(), len, GetMaxSizeOfDataMessage(), 
            false/*throw*/)) != MOK )
        return false;

    try {
        //iterate until slaves send empty messages
        for( ;; ) {
            //wait for data (table and dish indices)
            for( p = 0, pp = 0; p < nopus; p++ ) {
                //all processes should reply
                err = ReceiveMPIMessage( GetMBuffer(), GetMaxSizeOfDataMessage(), 
                                          false/*throw*/);
                if( err != MOK ) {
                    sprintf( strbuf, "MasterProcessRestTables: Receive failed: Error code %d.", err );
                    merror = strbuf;
                    continue;
                }
                if( !merror.empty())
                    continue;
                try {
                    MasterProcessCmdValsMtx( n, nn, values, &pp );
                } catch( myexception const& ex ) {
                    merror = ex.what();
                    continue;
                }
            }
            if( !merror.empty()) {
                error( merror.c_str());
                break;
            }
            hdmessage( 0 );
            if( pp <= 0 ) {
                //move on to next restaurants
                hdmessage( -1 );
                message( NULL );
                break;
            }

            //broadcast message of the vectors have migrated
            SetResMBufferHead( THeadDat );
            SetResMBufferCmd( TComMigrTbl );
            SetResMBufferNoVals( pp );
            len = SetResMBufferCRC();
            if( GetMaxSizeOfResMBuffer() <= len || GetMaxSizeOfDataMessage() <= len )
                throw myruntime_error( preamb + "Short of size in buffers. Increase it and recompile.");

//             if(( err = BcastMPIMessage( GetResMBuffer(), /*len*/GetMaxSizeOfDataMessage(), 
//                           false/*throw*/)) != MOK )
//                 throw myruntime_error( preamb + "Broadcast failed.");
            if(( err = BcastlenMPIMessage( GetResMBuffer(), len, GetMaxSizeOfDataMessage(), 
                          false/*throw*/)) != MOK )
                throw myruntime_error( preamb + "Broadcast failed.");

            if( !merror.empty()) {
                error( merror.c_str());
                break;
            }
        }
    } catch( myexception const& ex ) {
        merror = ex.what();
    }

    if( !merror.empty()) {
        error( merror.c_str());
        return false;
    }

    return true;
}





// =========================================================================
// MasterProcessCmdValsMtx: process command of sampling from matrix
//  probabilities
//
void HDPsampler::MasterProcessCmdValsMtx( int rb, int re, double* values, int* pp )
{
    if( !mMPIIAmMaster())
        return;

    if( mMPIMyRank() < 0 || mMPIRingSize() < 1 )
        throw myruntime_error( "MasterProcessCmdValsMtx: Invalid MPI rank or ring size." );
    if( values == NULL )
        throw myruntime_error("MasterProcessCmdValsMtx: Null address of buffer.");
    if( pp == NULL || *pp < 0  )
        throw myruntime_error( "MasterProcessCmdValsMtx: Memory access error." );

    if( GetBasin() == NULL || GetMenu() == NULL || GetChain() == NULL )
        return;

    int*        dnxmap = GetLocIntBuf();//map of dish indices
    Restaurant* rest = NULL;
    Table*      tbl = NULL;
    Dish*       newdsh = NULL;
    int head, cmd, novals;
    int r, t, k;
    int newk, svnk, updk;
    int mm;

    if( dnxmap == NULL )
        throw myruntime_error("MasterProcessCmdValsMtx: Null local buffer of integers.");
    if( GetSizeLocIntBuf() < GetBasin()->GetSize())
        throw myruntime_error("MasterProcessCmdValsMtx: Short of size of local buffer of integers.");
    ResetLocIntBuf( -1 );

    head = GetMBufferHead();
    cmd = GetMBufferCmd();
    novals = GetMBufferNoVals();
    if( !CheckMBufferCRC())
        throw myruntime_error("MasterProcessCmdValsMtx: Invalid CRC.");

    if( head != THeadDat )
        throw myruntime_error("MasterProcessCmdValsMtx: Invalid header.");
    if( cmd != TComValsMtx )
        throw myruntime_error("MasterProcessCmdValsMtx: Invalid command.");

    if( novals <= 0 )
        return;

    for( mm = 0; mm < novals; ) {
        tbl = NULL;
        newdsh = NULL;
        r = ( int )GetMBufferValueAt( mm++ );//restaurant index
        t = ( int )GetMBufferValueAt( mm++ );//table index
        svnk = newk = ( int )GetMBufferValueAt( mm++ );//new dish of table
        updk = ( int )GetMBufferValueAt( mm++ );//dish index after update by slave

        if( novals < mm )
            throw myruntime_error("MasterProcessCmdValsMtx: Unexpected number of values.");

        //chain stays static in process
        if( r < rb || re < r || GetChain()->GetActualSize() <= r ||
          ( rest = GetChain()->GetRestaurantAt( r )) == NULL )
            throw myruntime_error("MasterProcessCmdValsMtx: Invalid restaurant index received.");

        if( GetSizeLocIntBuf() <= GetMenu()->GetSize())
            throw myruntime_error("MasterProcessCmdValsMtx: Short of size of local buffer of integers.");

        if( t < 0 || rest->GetSize() <= t ||
        ( tbl = rest->GetTableAt( t )) == NULL )
            throw myruntime_error("MasterProcessCmdValsMtx: Invalid table index received.");

        if(( k = tbl->GetDishIndex()) < 0 || GetMenu()->GetSize() <= k )
            throw myruntime_error("MasterProcessCmdValsMtx: Invalid table's dish index.");

        if( 0 <= newk && updk != newk )
            throw myruntime_error("MasterProcessCmdValsMtx: Invalid dish index updated by slave.");
        if( updk < 0 || GetBasin()->GetSize() <= updk )
            throw myruntime_error("MasterProcessCmdValsMtx: Invalid dish index updated by slave.");

        //slave and master dish indices may differ when creating new table with new dish;
        //dishes are global across restaurants, 
        //no need for reinitialization of dnxmap -- processing information from single slave process
        if( 0 <= dnxmap[updk])
            svnk = newk = dnxmap[updk];

        if( 0 <= newk && ( GetMenu()->GetSize() <= newk ||
          ( newdsh = GetMenu()->GetDishAt( newk )) == NULL )) {
              MasterMessage( "New table's dish instead of vanished one." );
              newdsh = NULL;
              svnk = newk = -1;
        }

        if( newk != k ) {
            //value of newk will be changed if it's <0
            if( !MoveTable( rest, t, k, tbl, newk ))
                throw myruntime_error("MasterProcessCmdValsMtx: Moving of table failed.");

            values[(*pp)++] = ( double )rb;//begin index of restaurants
            values[(*pp)++] = ( double )re;//end index of restaurants
            values[(*pp)++] = ( double )r;//restaurant index
            values[(*pp)++] = ( double )t;//table index
            values[(*pp)++] = ( double )svnk;//new (or existing) dish for table
        }
        //map dish indices between slave and master after update
        dnxmap[updk] = newk;
    }
}





// =========================================================================
// MasterprocessBcastKappa0: broadcast hyperparameter kappa0 value
//
bool HDPsampler::MasterprocessBcastKappa0()
{
    if( !mMPIIAmMaster())
        return false;

    if( GetMenu() == NULL )
        throw myruntime_error( "MasterprocessBcastKappa0: Null Menu." );

    if( mMPIMyRank() < 0 || mMPIRingSize() < 1 )
        throw myruntime_error( "MasterprocessBcastKappa0: Invalid MPI rank or ring size." );

    double* values = GetResMBufferValues();
    int     len, err;

    if( values == NULL )
        throw myruntime_error("MasterprocessBcastKappa0: Null address of buffer.");

    SetMBufferHead( THeadDat );
    SetMBufferCmd( TComValsKappa0 );
    SetMBufferNoVals( 1 );
    SetMBufferValueAt( 0, GetMenu()->GetKappa0());
    len = SetMBufferCRC();
    if( GetMaxSizeOfResMBuffer() <= len || GetMaxSizeOfDataMessage() <= len ) {
        error("MasterprocessBcastKappa0: Short of size in buffers.");
        return false;
    }
//     if(( err = BcastMPIMessage( GetMBuffer(), /*len*/GetMaxSizeOfDataMessage(), 
//             false/*throw*/)) != MOK )
//         return false;
    if(( err = BcastlenMPIMessage( GetMBuffer(), len, GetMaxSizeOfDataMessage(), 
            false/*throw*/)) != MOK )
        return false;

    return true;
}





// =========================================================================
// MasterprocessBcastNu0: broadcast hyperparameter nu0 value
//
bool HDPsampler::MasterprocessBcastNu0()
{
    if( !mMPIIAmMaster())
        return false;

    if( GetMenu() == NULL )
        throw myruntime_error( "MasterprocessBcastNu0: Null Menu." );

    if( mMPIMyRank() < 0 || mMPIRingSize() < 1 )
        throw myruntime_error( "MasterprocessBcastNu0: Invalid MPI rank or ring size." );

    double* values = GetResMBufferValues();
    int     len, err;

    if( values == NULL )
        throw myruntime_error("MasterprocessBcastNu0: Null address of buffer.");

    SetMBufferHead( THeadDat );
    SetMBufferCmd( TComValsNu0 );
    SetMBufferNoVals( 1 );
    SetMBufferValueAt( 0, GetMenu()->GetNu0());
    len = SetMBufferCRC();
    if( GetMaxSizeOfResMBuffer() <= len || GetMaxSizeOfDataMessage() <= len ) {
        error("MasterprocessBcastNu0: Short of size in buffers.");
        return false;
    }
//     if(( err = BcastMPIMessage( GetMBuffer(), /*len*/GetMaxSizeOfDataMessage(), 
//             false/*throw*/)) != MOK )
//         return false;
    if(( err = BcastlenMPIMessage( GetMBuffer(), len, GetMaxSizeOfDataMessage(), 
            false/*throw*/)) != MOK )
        return false;

    return true;
}





// =========================================================================
// MasterprocessBcastConcParams: broadcast values of concentration 
//  parameters
//
bool HDPsampler::MasterprocessBcastConcParams()
{
    if( !mMPIIAmMaster())
        return false;

    if( mMPIMyRank() < 0 || mMPIRingSize() < 1 )
        throw myruntime_error( "MasterprocessBcastConcParams: Invalid MPI rank or ring size." );

    double* values = GetResMBufferValues();
    int     len, err;

    if( values == NULL )
        throw myruntime_error("MasterprocessBcastConcParams: Null address of buffer.");

    SetMBufferHead( THeadDat );
    SetMBufferCmd( TComValsConc );
    SetMBufferNoVals( 2 );
    SetMBufferValueAt( 0, GetDPMTau());
    SetMBufferValueAt( 1, GetDPMGamma());
    len = SetMBufferCRC();
    if( GetMaxSizeOfResMBuffer() <= len || GetMaxSizeOfDataMessage() <= len ) {
        error("MasterprocessBcastConcParams: Short of size in buffers.");
        return false;
    }
//     if(( err = BcastMPIMessage( GetMBuffer(), /*len*/GetMaxSizeOfDataMessage(), 
//             false/*throw*/)) != MOK )
//         return false;
    if(( err = BcastlenMPIMessage( GetMBuffer(), len, GetMaxSizeOfDataMessage(), 
            false/*throw*/)) != MOK )
        return false;

    return true;
}





// =========================================================================
// ss_hdphyperkappa0: function to calculate 
//  log p(kappa0|observations,state vector,other hyperparameters)
//
int ss_hdphyperkappa0( double x, double* lfx, void* params )
{
    mystring preamb = "ss_hdphyperkappa0: ";
    if( lfx == NULL )
        throw myruntime_error( preamb + "Memory access error.");
    if( params == NULL )
        throw myruntime_error( preamb + "Null parameters.");

    void**      ppars = ( void** )params;
    HDPsampler* pthis = ( HDPsampler* )( ppars[0]);
    SPDmatrix** devmtcs = ( SPDmatrix** )( ppars[1]);

    if( pthis == NULL || devmtcs == NULL )
        throw myruntime_error( preamb + "Invalid parameters.");

    const Menu* menu = pthis->GetMenu();
    const Basin* basin = pthis->GetBasin();

    if( menu == NULL || basin == NULL )
        throw myruntime_error( preamb + "Null menu structures.");

    const int       dim = menu->GetDim();
    const double    do2 = ( double )dim * 0.5;
    const double    nu0 = menu->GetNu0();
    const double    nu0o2 = nu0 * 0.5;
    const double    no2 = nu0o2 + do2*( 1.0+pthis->GetDegFAdjustment());
    const double    k0 = menu->GetKappa0();
    const double    a_k0 = menu->GetKappa0_pp_a();
    const double    b_k0 = menu->GetKappa0_pp_b();
    const double    nodshs = ( double )menu->GetActualSize();//number of dishes

    if( dim < 1 )
        throw myruntime_error( preamb + "Invalid dimensions.");

    if( a_k0 <= 0.0 || b_k0 <= 0.0 || k0 <= 0.0 || nu0 <= 0.0 )
        throw myruntime_error( preamb + "Invalid hyperparameters.");

    if( nodshs < 1.0 )
        throw myruntime_error( preamb + "No dishes.");

    if( x <= 0.0 )
        throw myruntime_error( preamb + "Invalid argument.");

    const Dish*         dish;
    const SPDmatrix*    L_k;//dish scale matrix
    double              n_k;//number of samples in dish
    SPDmatrix           Last( dim );//adjusted scale matrix
    double  fact, ldet;
    double  term;
    int     k, err;

    *lfx = 0.0;
    term = 0.0;

    for( k = 0; k < menu->GetSize(); k++ ) {
        dish = menu->GetDishAt( k );
        if( dish == NULL )
            continue;
        n_k = ( double )dish->GetActualSize();
        L_k = menu->GetSMatrixAt( k );
        ldet = menu->GetLDetSMAt( k );
        if( n_k < 1 )
            throw myruntime_error( preamb + "Invalid no. samples in dish.");
        if( L_k == NULL )
            throw myruntime_error( preamb + "Null dish scale matrix.");
        if( devmtcs[k] == NULL )
            throw myruntime_error( preamb + "Null deviation matrix.");
        //calculate adjusted scale matrix for dish
        Last = *devmtcs[k];
        fact = x*n_k/(x+n_k) - k0*n_k/(k0+n_k);
        if( fact ) {
            if(( err = Last.Scale( fact )) != 0 )
                throw myruntime_error( preamb + TranslatePSLError( err ));
            if(( err = Last.Add( *L_k )) != 0 )
                throw myruntime_error( preamb + TranslatePSLError( err ));
            //calculate determinant of adjusted scale matrix
            if(( err = Last.CholeskyDecompose()) != 0 )
                throw myruntime_error( preamb + TranslatePSLError( err ));
            if(( err = Last.CDedLogDet( &ldet )) != 0 )
                throw myruntime_error( preamb + TranslatePSLError( err ));
        }
        //add to log probability
        *lfx -= ( no2 + n_k*0.5 )* ldet;
        term += log( x + n_k );
    }

    *lfx -= do2 * term;
    *lfx -= x / b_k0;
    *lfx += ( a_k0-1.0 + nodshs*do2 ) * log(x);

    if( !isfinite( *lfx ))
        throw myruntime_error( preamb + "Invalid final log probability value.");

    return 0;
}

// -------------------------------------------------------------------------
// MasterprocessSampleKappa0: sample hyperparameter k_0
//
bool HDPsampler::MasterprocessSampleKappa0()
{
    if( GetMenu() == NULL || GetBasin() == NULL )
        throw myruntime_error("MasterprocessSampleKappa0: Null menu structures.");

    SliceSmpl   sls;
    const int   slsmpar = 2;
    const int   burnin = 10;
    void*       params[2] = {( void* )this,NULL};

    mystring    preamb = "MasterprocessSampleKappa0: ";
    mystring    strerr;
    const int   dim = GetMenu()->GetDim();
    int         totvs = GetBasin()->GetActualSize();//total number of observations
    int         nodshs = GetMenu()->GetActualSize();//number of dishes
    int         szmenu = GetMenu()->GetSize();
    double      k0 = GetMenu()->GetKappa0();
    double      nk;//number of samples in dish
    SPDmatrix**         devmtcs = NULL;
    const Pslvector*    mu0 = GetMenu()->GetMu0();
    const Pslvector*    pnts;
    const Dish*         dish;
    const Pslvector*    mu;
    Pslvector           eta_k;//mean of sample vectors for dish k
    double      left, right;
    double      smpl;
    int k, err;

    if( mu0 == NULL )
        throw myruntime_error( preamb + "Null prior mean vector.");

    if( dim < 1 || nodshs < 1 )
        throw myruntime_error( preamb + "No dishes.");

    devmtcs = ( SPDmatrix** )malloc( szmenu * sizeof( void* ));
    if( devmtcs == NULL )
        throw myruntime_error( preamb + "Not enough memory.");
    params[1] = ( void* )devmtcs;

    //calculate deviation matrix for scale matrix adjustment for each dish
    for( k = 0; k < szmenu; k++ ) {
        devmtcs[k] = NULL;
        dish = GetMenu()->GetDishAt( k );
        if( dish == NULL )
            continue;
        devmtcs[k] = new SPDmatrix( dim );
        if( devmtcs[k] == NULL )
            throw myruntime_error( preamb + "Not enough memory.");
        //calculate mean of sample vectors for dish d
        nk = ( double )dish->GetActualSize();
        mu = GetMenu()->GetMuVectorAt( k );
        if( mu == NULL || nk < 1 )
            throw myruntime_error( preamb + "Null mean vector of dish.");
        eta_k = *mu;
        if(( err = eta_k.Superposition( -1.0, *mu0 )) != 0 )
            throw myruntime_error( preamb + TranslatePSLError( err ));
        if(( err = eta_k.MultiplyBy( (k0+nk)/nk )) != 0 )
            throw myruntime_error( preamb + TranslatePSLError( err ));
        //NOTE:eta_k is now mean of samples - mu0
        if(( err = devmtcs[k]->Mul( eta_k, eta_k )) != 0 )
            throw myruntime_error( preamb + TranslatePSLError( err ));
    }

    right = k0;
    left = 1.e-5;//left bound for sampler

    sls.SetParams(( void* )params );
    sls.SetFEval( ss_hdphyperkappa0 );
    sls.SetLeftX( left );
    sls.SetX0( right );
    sls.SetBurnin( burnin );
    sls.SetW( SLC_MAX( 1.0, TIMES2( right )));//set width dependent on prev. value
    sls.SetM( slsmpar );

    try {
        if(( err = sls.Run()) != 0 )
            throw myruntime_error( preamb + TranslateSliceSmplError( err ));
    } catch( myexception const& ex ) {
        strerr = ex.what();
    }

    if( devmtcs ) {
        for( k = 0; k < szmenu; k++ )
            if( devmtcs[k])
                delete devmtcs[k];
        free( devmtcs );
        devmtcs = NULL;
    }

    if( !strerr.empty())
        throw myruntime_error( preamb + strerr );

    pnts = const_cast<const SliceSmpl&>( sls ).GetLastSamples();
    if( pnts == NULL || pnts->GetSize() < 1 )
        throw myruntime_error( preamb + "Null sampled points.");
    smpl = pnts->GetValueAt( pnts->GetSize()-1 );
    if( smpl <= 0.0 || ( double )totvs * 0.5 < smpl )
        throw myruntime_error( preamb + "Too large sampled point.");

    //NOTE: do not set new value of kappa0 here; 
    //current value is required for adjustments
    if( !ProcessCmdNewKappa0( smpl ))
        return false;
    return true;
}





// =========================================================================
// ss_hdphypernu0: function to calculate 
//  log p(nu0|observations,state vector,other hyperparameters)
//
int ss_hdphypernu0( double x, double* lfx, void* params )
{
    mystring preamb = "ss_hdphypernu0: ";
    if( lfx == NULL )
        throw myruntime_error( preamb + "Memory access error.");
    if( params == NULL )
        throw myruntime_error( preamb + "Null parameters.");

    void**      ppars = ( void** )params;
    HDPsampler* pthis = ( HDPsampler* )( ppars[0]);

    if( pthis == NULL )
        throw myruntime_error( preamb + "Invalid parameters.");

    const Menu* menu = pthis->GetMenu();
    const Basin* basin = pthis->GetBasin();

    if( menu == NULL || basin == NULL )
        throw myruntime_error( preamb + "Null menu structures.");

    const double    dim = menu->GetDim();
    const double    do2 = dim * 0.5;
    const double    nu0 = x;
    const double    nu0o2 = nu0 * 0.5;
    const double    no2 = nu0o2 + do2*( 1.0+pthis->GetDegFAdjustment());
    const double    np1o2 = no2 + 0.5;//(nu0+1)/2
    const double    nodshs = ( double )menu->GetActualSize();//number of dishes
    const double    L0ldet = menu->GetLDetS0();

    if( dim < 1 )
        throw myruntime_error( preamb + "Invalid dimensions.");

    if( no2 < do2 )
        throw myruntime_error( preamb + "Invalid hyperparameters.");

    if( nodshs < 1.0 )
        throw myruntime_error( preamb + "No dishes.");

    if( x <= 0.0 )
        throw myruntime_error( preamb + "Invalid argument.");

    const Dish*         dish;
    double              ldet;//log det of dish's scale matrix
    double              nko2;//number of samples in dish over 2
    double  lgx, pg1x, lgxn, pg1xn;
    double  lgxsum, pg1xsum, lgxnsum, pg1xnsum;
    double  j, ldtsum, operr;
    int     k, err;

    *lfx = 0.0;
    ldtsum = 0.0;
    lgxsum = pg1xsum = 0.0;
    lgxnsum = pg1xnsum = 0.0;

    //calculate sums of log gamma and polygamma functions
    for( j = 0.5; j <= do2; j += 0.5 ) {
        if(( err = psl_lngamma_e( np1o2 - j, &lgx, &operr )) != 0 )
            throw myruntime_error( preamb + TranslatePSLError( err ));
        if(( err = psl_psi_n_e( 1, np1o2 - j, &pg1x, &operr )) != PSL_SUCCESS )
            throw myruntime_error( preamb + TranslatePSLError( err ));
        lgxsum += lgx;
        pg1xsum += pg1x;
    }

    for( k = 0; k < menu->GetSize(); k++ ) {
        dish = menu->GetDishAt( k );
        if( dish == NULL )
            continue;
        ldet = menu->GetLDetSMAt( k );
        nko2 = ( double )dish->GetActualSize();
        if( nko2 < 1 )
            throw myruntime_error( preamb + "Invalid no. samples in dish.");
        nko2 *= 0.5;
        //add factorized log det
        ldtsum += ldet * ( no2+nko2 );
        //calculate sums of log gamma and polygamma functions
        for( j = 0.5; j <= do2; j += 0.5 )
        {
            if(( err = psl_lngamma_e( np1o2 + nko2 - j, &lgxn, &operr )) != 0 )
                throw myruntime_error( preamb + TranslatePSLError( err ));
            if(( err = psl_psi_n_e( 1, np1o2 + nko2 - j, &pg1xn, &operr )) != PSL_SUCCESS )
                throw myruntime_error( preamb + TranslatePSLError( err ));
            lgxnsum += lgxn;
            pg1xnsum += pg1xn;
        }
    }

    //add to final log probability
    *lfx += 0.5 * log(0.25*(nodshs*pg1xsum - pg1xnsum));
    *lfx += lgxnsum - nodshs*lgxsum;
    *lfx += nodshs*no2*L0ldet - ldtsum;

    if( !isfinite( *lfx ))
        throw myruntime_error( preamb + "Invalid final log probability value.");

    return 0;
}

// -------------------------------------------------------------------------
// MasterprocessSampleNu0: sample hyperparameter nu_0
//
bool HDPsampler::MasterprocessSampleNu0()
{
    if( GetMenu() == NULL || GetBasin() == NULL )
        throw myruntime_error("MasterprocessSampleNu0: Null menu structures.");

    SliceSmpl   sls;
    const int   slsmpar = 2;
    const int   burnin = 10;
    void*       params[2] = {( void* )this,NULL};

    mystring    preamb = "MasterprocessSampleNu0: ";
    mystring    strerr;
    const int   dim = GetMenu()->GetDim();
    int         totvs = GetBasin()->GetActualSize();//total number of observations
    int         nodshs = GetMenu()->GetActualSize();//number of dishes
    double      nu0 = GetMenu()->GetNu0();//number of samples in dish
    const Pslvector*    pnts;
    double      left, right;
    double      smpl;
    int err;

    if( dim < 1 || nodshs < 1 )
        throw myruntime_error( preamb + "No dishes.");

    right = nu0;
    left = 1.e-5;//left bound for sampler
    if( GetDegFAdjustment() < 0 ) {
        left += ( double )dim;
        if( right < ( double )dim )
            right += ( double )dim;
    }

    sls.SetParams(( void* )params );
    sls.SetFEval( ss_hdphypernu0 );
    sls.SetLeftX( left );
    sls.SetX0( right );
    sls.SetBurnin( burnin );
    sls.SetW( SLC_MAX( 1.0, TIMES2( right )));//set width dependent on prev. value
    sls.SetM( slsmpar );

    try {
        if(( err = sls.Run()) != 0 )
            throw myruntime_error( preamb + TranslateSliceSmplError( err ));
    } catch( myexception const& ex ) {
        strerr = ex.what();
    }

    if( !strerr.empty())
        throw myruntime_error( preamb + strerr );

    pnts = const_cast<const SliceSmpl&>( sls ).GetLastSamples();
    if( pnts == NULL || pnts->GetSize() < 1 )
        throw myruntime_error( preamb + "Null sampled points.");
    smpl = pnts->GetValueAt( pnts->GetSize()-1 );
    if( smpl < left || ( double )totvs * 0.5 < smpl )
        throw myruntime_error( preamb + "Too large sampled point.");

    //NOTE: do not set new value of nu0 here
    if( !ProcessCmdNewNu0( smpl ))
        return false;
    return true;
}





// =========================================================================
// ss_hdpconctau: function to calculate log p(tau|m_1,...,m_J;n_1,...,n_J)
//
int ss_hdpconctau( double x, double* lfx, void* params )
{
    mystring preamb = "ss_hdpconctau: ";
    if( lfx == NULL )
        throw myruntime_error( preamb + "Memory access error.");
    if( params == NULL )
        throw myruntime_error( preamb + "Null parameters.");

    void**      ppars = ( void** )params;
    HDPsampler* pthis = ( HDPsampler* )( ppars[0] );

    if( pthis == NULL )
        throw myruntime_error( preamb + "Invalid parameters.");

    const Menu* menu = pthis->GetMenu();
    const RntChain* chain = pthis->GetChain();
    const Restaurant* rest;

    if( menu == NULL || chain == NULL )
        throw myruntime_error( preamb + "Null restaurant structures.");

    int     norsts = chain->GetSize();
    double  J = norsts;

    if( norsts != chain->GetActualSize() || norsts < 1 )
        throw myruntime_error( preamb + "Inconsistent size of restaurant chain.");

    if( x <= 0.0 )
        throw myruntime_error( preamb + "Invalid argument.");

    double  cerr;
    double  klogk;
    double  lgk, lgkn, lx;
    double  lgx, lgxn, psix, psixn;
    double  pg1x, pg1xn;
    double  lgxnsum, psixnsum, pg1xnsum;
    int code;
    int r, n, k;
    int nsum = 0, ksum = 0; 

    if(( code = psl_lngamma_e( x, &lgx, &cerr )) != PSL_SUCCESS )
        throw myruntime_error( preamb + TranslatePSLError( code ));
    if(( code = psl_psi_e( x, &psix, &cerr )) != PSL_SUCCESS )
        throw myruntime_error( preamb + TranslatePSLError( code ));

    if(( code = psl_psi_n_e( 1, x, &pg1x, &cerr )) != PSL_SUCCESS )
        throw myruntime_error( preamb + TranslatePSLError( code ));

    lgxnsum = psixnsum = pg1xnsum = 0.0;

    for( r = 0; r < norsts; r++ )
    {
        rest = chain->GetRestaurantAt( r );
        if( rest == NULL )
            throw myruntime_error( preamb + "Null restaurant.");

        n = rest->GetNoVectors();//number of vectors in restaurant
        k = rest->GetActualSize();//number of tables in restaurant

        if( n < 1 ) throw myruntime_error( preamb + "No vectors in restaurant.");
        if( k < 1 ) throw myruntime_error( preamb + "No tables in restaurant.");

        nsum += n;
        ksum += k;

        if(( code = psl_lngamma_e( x + n, &lgxn, &cerr )) != PSL_SUCCESS )
            throw myruntime_error( preamb + TranslatePSLError( code ));
        if(( code = psl_psi_e( x + n, &psixn, &cerr )) != PSL_SUCCESS )
            throw myruntime_error( preamb + TranslatePSLError( code ));
        lgxnsum += lgxn;
        psixnsum += psixn;

        if(( code = psl_psi_n_e( 1, x + n, &pg1xn, &cerr )) != PSL_SUCCESS )
            throw myruntime_error( preamb + TranslatePSLError( code ));
        pg1xnsum += pg1xn;
    }

//     if( ksum < 1 )
//         throw myruntime_error( preamb + "Invalid overall number of tables.");
//     klogk = ( double )ksum * log(( double )ksum / J );
//     if(( code = psl_lngamma_e( ksum, &lgk, &cerr )) != PSL_SUCCESS )
//         throw myruntime_error( preamb + TranslatePSLError( code ));
//     if(( code = psl_lngamma_e( ksum + nsum, &lgkn, &cerr )) != PSL_SUCCESS )
//         throw myruntime_error( preamb + TranslatePSLError( code ));

    *lfx = 0.0;
//     *lfx = lgkn-lgk-klogk;//to compesnate for underflow
    lx = pg1xnsum-J*pg1x+(psixnsum-J*psix)/x;
    if( 0.0 < lx )
        *lfx += 0.5 * log(lx);
    else
        *lfx -= 0.5 * SLC_LOG_DBL_MAX;
    *lfx += (( double )ksum)*log(x) + J*lgx - lgxnsum;

    return 0;
}

// =========================================================================
// MasterprocessSampleTau: sample concentration parameter tau (inner 
//      hierarchical level)
//
bool HDPsampler::MasterprocessSampleTau()
{
    if( GetMenu() == NULL || GetChain() == NULL )
        throw myruntime_error("MasterprocessSampleTau: Null menu.");

    SliceSmpl   sls;
    const int   slsmpar = 2;
    const int   burnin = 10;
    void*       params[] = {( void* )this};

    mystring    preamb = "MasterprocessSampleTau: ";
    double      notbls = GetMenu()->GetNoTables();//overall number of tables
    int         norsts = GetChain()->GetActualSize();//number of restaurants
    const Pslvector* pnts;
    double      left, right;
    double      smpl;
    int err;

    if( norsts < 1 )
        throw myruntime_error( preamb + "No restaurants.");

    right = notbls /( double )norsts;
    if( right < 1.0 ) {
        warning("Average no. tables per restaurant <1.", false );
        right = 1.0;
    }
    left = 0.0001;//left bound for sampler

    sls.SetParams(( void* )params );
    sls.SetFEval( ss_hdpconctau );
    sls.SetLeftX( left );
    sls.SetX0( right );
    sls.SetBurnin( burnin );
    sls.SetW( right );//set width dependent on avg. no. tables over all restautants
    sls.SetM( slsmpar );

    if(( err = sls.Run()) != 0 )
        throw myruntime_error( preamb + TranslateSliceSmplError( err ));

    pnts = const_cast<const SliceSmpl&>( sls ).GetLastSamples();
    if( pnts == NULL || pnts->GetSize() < 1 )
        throw myruntime_error( preamb + "Null sampled points.");
    smpl = pnts->GetValueAt( pnts->GetSize()-1 );
    if( smpl <= 0.0 )
        throw myruntime_error( preamb + "Invalid sampled point.");

    SetDPMTau( smpl );
    return true;
}





// =========================================================================
// ss_hdpconcgamma: function to calculate log p(gamma|k,n)
//
int ss_hdpconcgamma( double x, double* lfx, void* params )
{
    mystring preamb = "ss_hdpconcgamma: ";
    if( lfx == NULL )
        throw myruntime_error( preamb + "Memory access error.");
    if( params == NULL )
        throw myruntime_error( preamb + "Null parameters.");

    void**      ppars = ( void** )params;
    HDPsampler* pthis = ( HDPsampler* )( ppars[0] );

    if( pthis == NULL )
        throw myruntime_error( preamb + "Invalid parameters.");

    Menu*       menu = pthis->GetMenu();

    if( menu == NULL )
        throw myruntime_error( preamb + "Null menu.");

    double  n = menu->GetNoTables();//overall number of tables
    double  k = menu->GetActualSize();//number of dishes

    if( x <= 0.0 ) throw myruntime_error( preamb + "Invalid argument.");
    if( n <= 0 ) throw myruntime_error( preamb + "No tables.");
    if( k <= 0 ) throw myruntime_error( preamb + "No dishes.");

    double  cerr;
//     double  klogk = k * log( k );
    double  lgk, lgkn, lx;
    double  lgx, lgxn, psix, psixn;
    double  pg1x, pg1xn;
    int code;

    if(( code = psl_lngamma_e( x, &lgx, &cerr )) != PSL_SUCCESS )
        throw myruntime_error( preamb + TranslatePSLError( code ));
    if(( code = psl_lngamma_e( x + n, &lgxn, &cerr )) != PSL_SUCCESS )
        throw myruntime_error( preamb + TranslatePSLError( code ));

//     if(( code = psl_lngamma_e( k, &lgk, &cerr )) != PSL_SUCCESS )
//         throw myruntime_error( preamb + TranslatePSLError( code ));
//     if(( code = psl_lngamma_e( k + n, &lgkn, &cerr )) != PSL_SUCCESS )
//         throw myruntime_error( preamb + TranslatePSLError( code ));

    if(( code = psl_psi_e( x, &psix, &cerr )) != PSL_SUCCESS )
        throw myruntime_error( preamb + TranslatePSLError( code ));
    if(( code = psl_psi_e( x + n, &psixn, &cerr )) != PSL_SUCCESS )
        throw myruntime_error( preamb + TranslatePSLError( code ));

    if(( code = psl_psi_n_e( 1, x, &pg1x, &cerr )) != PSL_SUCCESS )
        throw myruntime_error( preamb + TranslatePSLError( code ));
    if(( code = psl_psi_n_e( 1, x + n, &pg1xn, &cerr )) != PSL_SUCCESS )
        throw myruntime_error( preamb + TranslatePSLError( code ));

    *lfx = 0.0;
//     *lfx += (lgkn-lgk-klogk);//to compesnate for underflow
    lx = pg1xn-pg1x+(psixn-psix)/x;
    if( 0.0 < lx )
        *lfx += 0.5 * log(lx);
    else
        *lfx -= 0.5 * SLC_LOG_DBL_MAX;
    *lfx += k*log(x) + lgx - lgxn;

    return 0;
}

// -------------------------------------------------------------------------
// MasterprocessSampleGamma: sample concentration parameter gamma
//
bool HDPsampler::MasterprocessSampleGamma()
{
    if( GetMenu() == NULL )
        throw myruntime_error("MasterprocessSampleGamma: Null menu.");

    SliceSmpl   sls;
    const int   slsmpar = 2;
    const int   burnin = 10;
    void*       params[] = {( void* )this};

    mystring    preamb = "MasterprocessSampleGamma: ";
    int         notbls = GetMenu()->GetNoTables();//overall number of tables
    int         nodshs = GetMenu()->GetActualSize();//number of dishes
    const Pslvector* pnts;
    double      left, right;
    double      smpl;
    int err;

    if( nodshs < 1 )
        throw myruntime_error( preamb + "No dishes.");

    right = ( double )nodshs;
    left = 0.0001;//left bound for sampler

    sls.SetParams(( void* )params );
    sls.SetFEval( ss_hdpconcgamma );
    sls.SetLeftX( left );
    sls.SetX0( right );
    sls.SetBurnin( burnin );
    sls.SetW( right );//set width dependent on no. dishes
    sls.SetM( slsmpar );

    if(( err = sls.Run()) != 0 )
        throw myruntime_error( preamb + TranslateSliceSmplError( err ));

    pnts = const_cast<const SliceSmpl&>( sls ).GetLastSamples();
    if( pnts == NULL || pnts->GetSize() < 1 )
        throw myruntime_error( preamb + "Null sampled points.");
    smpl = pnts->GetValueAt( pnts->GetSize()-1 );
    if( smpl <= 0.0 )
        throw myruntime_error( preamb + "Invalid sampled point.");

    SetDPMGamma( smpl );
    return true;
}





// ================================OBS======================================
// ================================Tau======================================
// ars_hdpconctau_helper: function to calculate approx. of
//      d/dtau log p(tau|m_1,...,m_J;n_1,...,n_J) and 
//      d^2/dtau^2 log p(tau|m_1,...,m_J;n_1,...,n_J)
//
int ars_hdpconctau_helper( double x, double* hx, double* hdx, void* params )
{
    if( hx == NULL || hdx == NULL )
        throw myruntime_error("ars_hdpconctau_helper: Memory access error.");
    if( params == NULL )
        throw myruntime_error("ars_hdpconctau_helper: Null parameters.");

    void**      ppars = ( void** )params;
    HDPsampler* pthis = ( HDPsampler* )( ppars[0] );
    int         approx = ( int )( size_t )( ppars[1] );//approximate calculation

    if( pthis == NULL )
        throw myruntime_error("ars_hdpconctau_helper: Invalid parameters.");

    const Menu* menu = pthis->GetMenu();
    const RntChain* chain = pthis->GetChain();
    const Restaurant* rest;

    if( menu == NULL || chain == NULL )
        throw myruntime_error("ars_hdpconctau_helper: Null restaurant structures.");

    int     norsts = chain->GetSize();
    double  J = norsts;

    if( norsts != chain->GetActualSize() || norsts < 1 )
        throw myruntime_error("ars_hdpconctau_helper: Not consistent size of restaurant chain.");

    if( x <= 0.0 )
        throw myruntime_error("ars_hdpconctau_helper: Invalid argument.");

    double  err;
    double  x2 = x * x;
    double  psix, psixn;
    double  pg1x, pg1xn;
    double  psixnsum, pg1xnsum;
    int code;
    int r, n, k;
    int nsum = 0, ksum = 0; 

    if(( code = psl_psi_e( x, &psix, &err )) != PSL_SUCCESS )
        throw myruntime_error( TranslatePSLError( code ));
    if(( code = psl_psi_n_e( 1, x, &pg1x, &err )) != PSL_SUCCESS )
        throw myruntime_error( TranslatePSLError( code ));

    psixnsum = pg1xnsum = 0.0;

    for( r = 0; r < norsts; r++ )
    {
        rest = chain->GetRestaurantAt( r );
        if( rest == NULL )
            throw myruntime_error("ars_hdpconctau_helper: Null restaurant.");

        n = rest->GetNoVectors();//number of vectors in restaurant
        k = rest->GetActualSize();//number of tables in restaurant

        if( n < 1 ) throw myruntime_error("ars_hdpconctau_helper: No. vectors in restaurant is 0.");
        if( k < 1 ) throw myruntime_error("ars_hdpconctau_helper: No. tables in restaurant is 0.");

        nsum += n;
        ksum += k;

        if(( code = psl_psi_e( x + n, &psixn, &err )) != PSL_SUCCESS )
            throw myruntime_error( TranslatePSLError( code ));
        if(( code = psl_psi_n_e( 1, x + n, &pg1xn, &err )) != PSL_SUCCESS )
            throw myruntime_error( TranslatePSLError( code ));
        psixnsum += psixn;
        pg1xnsum += pg1xn;
    }
    //unable to optimize when k=1: non (log-)cocave
    if( ksum < 2 )
        ksum = 2;

    *hx = ( double )ksum/x + J*psix - psixnsum;
    *hdx = -( double )ksum/x2 + J*pg1x - pg1xnsum;

    return 0;
}

// -------------------------------------------------------------------------
// ars_hdpconctau: function to calculate log p(tau|m_1,...,m_J;n_1,...,n_J) and 
//      d/dtau log p(tau|m_1,...,m_J;n_1,...,n_J)
//
int ars_hdpconctau( double x, double* hx, double* hdx, void* params )
{
    if( hx == NULL || hdx == NULL )
        throw myruntime_error("ars_hdpconctau: Memory access error.");
    if( params == NULL )
        throw myruntime_error("ars_hdpconctau: Null parameters.");

    void**      ppars = ( void** )params;
    HDPsampler* pthis = ( HDPsampler* )( ppars[0] );
    int         approx = ( int )( size_t )( ppars[1] );//approximate calculation

    if( pthis == NULL )
        throw myruntime_error("ars_hdpconctau: Invalid parameters.");

    const Menu* menu = pthis->GetMenu();
    const RntChain* chain = pthis->GetChain();
    const Restaurant* rest;

    if( menu == NULL || chain == NULL )
        throw myruntime_error("ars_hdpconctau: Null restaurant structures.");

    int     norsts = chain->GetSize();
    double  J = norsts;

    if( norsts != chain->GetActualSize() || norsts < 1 )
        throw myruntime_error("ars_hdpconctau: Not consistent size of restaurant chain.");

    if( x <= 0.0 )
        throw myruntime_error("ars_hdpconctau: Invalid argument.");

    double  err;
    double  x2 = x * x;
    double  klogk;
    double  lgk, lgkn;
    double  lgx, lgxn, psix, psixn;
    double  pg1x, pg1xn, pg2x, pg2xn;
    double  lgxnsum, psixnsum, pg1xnsum, pg2xnsum;
    int code;
    int r, n, k;
    int nsum = 0, ksum = 0; 

    if(( code = psl_lngamma_e( x, &lgx, &err )) != PSL_SUCCESS )
        throw myruntime_error( TranslatePSLError( code ));
    if(( code = psl_psi_e( x, &psix, &err )) != PSL_SUCCESS )
        throw myruntime_error( TranslatePSLError( code ));

    if( !approx ) {
        if(( code = psl_psi_n_e( 1, x, &pg1x, &err )) != PSL_SUCCESS )
            throw myruntime_error( TranslatePSLError( code ));
        if(( code = psl_psi_n_e( 2, x, &pg2x, &err )) != PSL_SUCCESS )
            throw myruntime_error( TranslatePSLError( code ));
    }

    lgxnsum = psixnsum = pg1xnsum = pg2xnsum = 0.0;

    for( r = 0; r < norsts; r++ )
    {
        rest = chain->GetRestaurantAt( r );
        if( rest == NULL )
            throw myruntime_error("ars_hdpconctau: Null restaurant.");

        n = rest->GetNoVectors();//number of vectors in restaurant
        k = rest->GetActualSize();//number of tables in restaurant

        if( n < 1 ) throw myruntime_error("ars_hdpconctau: No. vectors in restaurant is 0.");
        if( k < 1 ) throw myruntime_error("ars_hdpconctau: No. tables in restaurant is 0.");

        nsum += n;
        ksum += k;

        if(( code = psl_lngamma_e( x + n, &lgxn, &err )) != PSL_SUCCESS )
            throw myruntime_error( TranslatePSLError( code ));
        if(( code = psl_psi_e( x + n, &psixn, &err )) != PSL_SUCCESS )
            throw myruntime_error( TranslatePSLError( code ));
        lgxnsum += lgxn;
        psixnsum += psixn;

        if( !approx ) {
            if(( code = psl_psi_n_e( 1, x + n, &pg1xn, &err )) != PSL_SUCCESS )
                throw myruntime_error( TranslatePSLError( code ));
            if(( code = psl_psi_n_e( 2, x + n, &pg2xn, &err )) != PSL_SUCCESS )
                throw myruntime_error( TranslatePSLError( code ));
            pg1xnsum += pg1xn;
            pg2xnsum += pg2xn;
        }
    }
    //unable to sample and optimize when k=1: non (log-)cocave
    if( ksum < 2 )
        ksum = 2;

    klogk = ( double )ksum * log(( double )ksum / J );

    if(( code = psl_lngamma_e( ksum, &lgk, &err )) != PSL_SUCCESS )
        throw myruntime_error( TranslatePSLError( code ));
    if(( code = psl_lngamma_e( ksum + nsum, &lgkn, &err )) != PSL_SUCCESS )
        throw myruntime_error( TranslatePSLError( code ));

    //(lgkn-lgk-klogk): to compesnate for underflow
    *hx = (lgkn-lgk-klogk) + (( double )ksum)*log(x) + J*lgx - lgxnsum;
    if( !approx )
        *hx += 0.5*log(pg1xnsum-J*pg1x+(psixnsum-J*psix)/x);

    *hdx = ( double )ksum/x + J*psix - psixnsum;
    if( !approx )
        *hdx += 0.5*(x2*(pg2xnsum-J*pg2x)+x*(pg1xnsum-J*pg1x)-(psixnsum-J*psix))/
                    (x2*(pg1xnsum-J*pg1x)+x*(psixnsum-J*psix));

    return 0;
}

// -------------------------------------------------------------------------
// MasterprocessSampleTau: sample concentration parameter tau
//
bool HDPsampler::MasterprocessSampleTauObs()
{
    CARS    ars;

    if( GetMenu() == NULL || GetChain() == NULL )
        throw myruntime_error("MasterprocessSampleTau: Null menu.");

    void*   params[] = {( void* )this,( void* )0 };
    void*   parapp[] = {( void* )this,( void* )1 };

    const double    rtaccuracy = 0.01;//required accuracy of root
    const int       rtmaxits = 25;//max no. iterations for findiing root

    double  notbls = GetMenu()->GetNoTables();//number of tables over all dishes
    int     norsts = GetChain()->GetActualSize();//number of restaurants
    const Pslvector* points;
    double  left, right;
    double  smpl;
    int     code;
    int     burnin = 10;
    int     maxits = 10;
    int     its = ( int )( RNGc.GetDouble() * maxits ) + 1;

    if( norsts < 1 )
        throw myruntime_error("MasterprocessSampleTau: No restaurants.");

    right = notbls /( double )norsts;
    if( right < 1.0 ) {
        warning("MasterprocessSampleTau: Average no. tables per restaurant <1.", false );
        right = 1.0;
    }
    left = right / 100.0;

    ars.SetParams(( void* )params );
    ars.SetHEval( ars_hdpconctau );
    ars.SetLowerX( 0.0 );
    ars.SetNoIterats( burnin + its );

    try { ars.PushInitialX( left );
          ars.PushInitialX( right );
          ars.Run();
    } catch( myexception const& ex ) {
        warning( ex.what(), false );
        message( "Solving for suboptimal value...", false );
        for( ;; ) {
            code = rootByNRandBisection(
                ars_hdpconctau_helper,
                left,
                right,
                &right,/*root*/
                rtaccuracy/*accuracy*/,
                rtmaxits/*max no. iterations*/,
                ( void* )parapp 
            );
            if( code != PRT_ERR_DOMAIN || notbls <= right )
                break;
            left = right;
            right = notbls;
        }
        if( code != PRT_OK && code != PRT_MAXITERATS ) {
            warning( TranslatePRTError( code ), false );
            message( "Resetting to fallback value...", false );
            SetDPMTau( GetDefDPMTau());
            return true;
        }
        if( code != PRT_OK )
            warning( TranslatePRTError( code ), false );
        message( "Sampling again...", false );

        try { ars.PushInitialX( left );
              ars.PushInitialX( right );
              ars.Run();
        } catch( myexception const& ex2 ) {
            warning( ex2.what(), false );
            message( "Setting suboptimal value...", false );
            SetDPMTau( right );
            return true;
        }
    }

    points = const_cast<const CARS&>( ars ).GetResPoints();
    if( points == NULL )
        throw myruntime_error("MasterprocessSampleTau: Null sampled points.");
    if( points->GetSize() < burnin || points->GetSize() < 1 )
        throw myruntime_error("MasterprocessSampleTau: Invalid size of points.");
    smpl = points->GetValueAt( points->GetSize()-1 );
    if( smpl <= 0.0 )
        throw myruntime_error("MasterprocessSampleTau: Invalid sampled point.");

    SetDPMTau( smpl );
    return true;
}





// ===============================Gamma=====================================
// ars_hdpconcgamma_helper: function to calculate approx. of 
//      d/dgamma log p(gamma|k,n) and 
//      d^2/dgamma^2 log p(gamma|k,n)
//
int ars_hdpconcgamma_helper( double x, double* hx, double* hdx, void* params )
{
    if( params == NULL )
        throw myruntime_error("ars_hdpconcgamma_helper: Null parameters.");

    void**      ppars = ( void** )params;
    HDPsampler* pthis = ( HDPsampler* )( ppars[0] );
    int         approx = ( int )( size_t )( ppars[1] );//approximate calculation

    if( pthis == NULL )
        throw myruntime_error("ars_hdpconcgamma_helper: Invalid parameters.");

    Menu*       menu = pthis->GetMenu();

    if( menu == NULL )
        throw myruntime_error("ars_hdpconcgamma_helper: Null menu.");

    double  n = menu->GetNoTables();//number of tables over all dishes
    double  k = menu->GetActualSize();//number of dishes

    if( x <= 0.0 ) throw myruntime_error("ars_hdpconcgamma_helper: Invalid argument.");
    if( n <= 0 ) throw myruntime_error("ars_hdpconcgamma_helper: No. tables is 0.");
    if( k <= 0 ) throw myruntime_error("ars_hdpconcgamma_helper: No. dishes is 0.");

    //unable to sample and optimize when k=1: non (log-)cocave
    if( k < 2.0 )
        k = 2.0;

    double  err;
    double  x2 = x * x;
    double  psix, psixn;
    double  pg1x, pg1xn;
    int code;

    if(( code = psl_psi_e( x, &psix, &err )) != PSL_SUCCESS )
        throw myruntime_error( TranslatePSLError( code ));
    if(( code = psl_psi_e( x + n, &psixn, &err )) != PSL_SUCCESS )
        throw myruntime_error( TranslatePSLError( code ));

    if(( code = psl_psi_n_e( 1, x, &pg1x, &err )) != PSL_SUCCESS )
        throw myruntime_error( TranslatePSLError( code ));
    if(( code = psl_psi_n_e( 1, x + n, &pg1xn, &err )) != PSL_SUCCESS )
        throw myruntime_error( TranslatePSLError( code ));

    if( hx ) *hx = k/x + psix - psixn;
    if( hdx ) *hdx = -k/x2 + pg1x - pg1xn;

    return 0;
}

// -------------------------------------------------------------------------
// ars_hdpconcgamma: function to calculate log p(gamma|k,n) and 
//      d/dgamma log p(gamma|k,n)
//
int ars_hdpconcgamma( double x, double* hx, double* hdx, void* params )
{
    if( params == NULL )
        throw myruntime_error("ars_hdpconcgamma: Null parameters.");

    void**      ppars = ( void** )params;
    HDPsampler* pthis = ( HDPsampler* )( ppars[0] );
    int         approx = ( int )( size_t )( ppars[1] );//approximate calculation

    if( pthis == NULL )
        throw myruntime_error("ars_hdpconcgamma: Invalid parameters.");

    Menu*       menu = pthis->GetMenu();

    if( menu == NULL )
        throw myruntime_error("ars_hdpconcgamma: Null menu.");

    double  n = menu->GetNoTables();//number of tables over all dishes
    double  k = menu->GetActualSize();//number of dishes

    if( x <= 0.0 ) throw myruntime_error("ars_hdpconcgamma: Invalid argument.");
    if( n <= 0 ) throw myruntime_error("ars_hdpconcgamma: No. tables is 0.");
    if( k <= 0 ) throw myruntime_error("ars_hdpconcgamma: No. dishes is 0.");

    //unable to sample and optimize when k=1: non (log-)concave
    if( k < 2.0 )
        k = 2.0;

    double  err;
    double  x2 = x * x;
    double  klogk = k * log( k );
    double  lgk, lgkn;
    double  lgx, lgxn, psix, psixn;
    double  pg1x, pg1xn, pg2x, pg2xn;
    int code;

    if(( code = psl_lngamma_e( x, &lgx, &err )) != PSL_SUCCESS )
        throw myruntime_error( TranslatePSLError( code ));
    if(( code = psl_lngamma_e( k, &lgk, &err )) != PSL_SUCCESS )
        throw myruntime_error( TranslatePSLError( code ));
    if(( code = psl_lngamma_e( k + n, &lgkn, &err )) != PSL_SUCCESS )
        throw myruntime_error( TranslatePSLError( code ));
    if(( code = psl_lngamma_e( x + n, &lgxn, &err )) != PSL_SUCCESS )
        throw myruntime_error( TranslatePSLError( code ));

    if(( code = psl_psi_e( x, &psix, &err )) != PSL_SUCCESS )
        throw myruntime_error( TranslatePSLError( code ));
    if(( code = psl_psi_e( x + n, &psixn, &err )) != PSL_SUCCESS )
        throw myruntime_error( TranslatePSLError( code ));

    if( !approx ) {
        if(( code = psl_psi_n_e( 1, x, &pg1x, &err )) != PSL_SUCCESS )
            throw myruntime_error( TranslatePSLError( code ));
        if(( code = psl_psi_n_e( 1, x + n, &pg1xn, &err )) != PSL_SUCCESS )
            throw myruntime_error( TranslatePSLError( code ));

        if(( code = psl_psi_n_e( 2, x, &pg2x, &err )) != PSL_SUCCESS )
            throw myruntime_error( TranslatePSLError( code ));
        if(( code = psl_psi_n_e( 2, x + n, &pg2xn, &err )) != PSL_SUCCESS )
            throw myruntime_error( TranslatePSLError( code ));
    }

    if( hx ) {
        //(lgkn-lgk-klogk): to compesnate for underflow
        *hx = (lgkn-lgk-klogk) + k*log(x) + lgx - lgxn;
        if( !approx )
            *hx += 0.5*log(pg1xn-pg1x+(psixn-psix)/x);
    }
    if( hdx ) {
        *hdx = k/x + psix - psixn;
        if( !approx )
            *hdx += 0.5*(x2*(pg2xn-pg2x)+x*(pg1xn-pg1x)-(psixn-psix))/
                        (x2*(pg1xn-pg1x)+x*(psixn-psix));
    }
    return 0;
}

// -------------------------------------------------------------------------
// MasterprocessSampleGamma: sample concentration parameter gamma
//
bool HDPsampler::MasterprocessSampleGammaObs()
{
    CARS    ars;

    if( GetMenu() == NULL )
        throw myruntime_error("MasterprocessSampleGamma: Null menu.");

    void*   params[] = {( void* )this,( void* )0 };
    void*   parapp[] = {( void* )this,( void* )1 };

    const double    rtaccuracy = 0.01;//required accuracy of root
    const int       rtmaxits = 25;//max no. iterations for findiing root

    int     notbls = GetMenu()->GetNoTables();//number of tables over all dishes
    int     nodshs = GetMenu()->GetActualSize();//number of dishes
    const Pslvector*    points;
    double  left, right;
    double  smpl;
    int     code;
    int     burnin = 10;
    int     maxits = 10;
    int     its = ( int )( RNGc.GetDouble() * maxits ) + 1;

    if( nodshs < 2 ) {
        warning("MasterprocessSampleGamma: No. dishes <2.", false );
        nodshs = 2;
    }

    right = nodshs;
    left = right / 100.0;

    ars.SetParams(( void* )params );
    ars.SetHEval( ars_hdpconcgamma );
    ars.SetLowerX( 0.0 );
    ars.SetNoIterats( burnin + its );

    try { ars.PushInitialX( left );
          ars.PushInitialX( right );
          ars.Run();
    } catch( myexception const& ex ) {
        warning( ex.what(), false );
        message( "Solving for suboptimal value...", false );
        code = rootByNRandBisection(
                ars_hdpconcgamma_helper,
                left,
                right,
                &right,/*root*/
                rtaccuracy/*accuracy*/,
                rtmaxits/*max no. iterations*/,
                ( void* )parapp 
        );
        if( code != PRT_OK && code != PRT_MAXITERATS ) {
            warning( TranslatePRTError( code ), false );
            message( "Resetting to fallback value...", false );
            SetDPMGamma( GetDefDPMGamma());
            return true;
        }
        if( code != PRT_OK )
            warning( TranslatePRTError( code ), false );
        message( "Sampling again...", false );

        try { ars.PushInitialX( left );
              ars.PushInitialX( right );
              ars.Run();
        } catch( myexception const& ex2 ) {
            warning( ex2.what(), false );
            message( "Setting suboptimal value...", false );
            SetDPMGamma( right );
            return true;
        }
    }

    points = const_cast<const CARS&>( ars ).GetResPoints();
    if( points == NULL )
        throw myruntime_error("MasterprocessSampleGamma: Null sampled points.");
    if( points->GetSize() < burnin || points->GetSize() < 1 )
        throw myruntime_error("MasterprocessSampleGamma: Invalid size of points.");
    smpl = points->GetValueAt( points->GetSize()-1 );
    if( smpl <= 0.0 )
        throw myruntime_error("MasterprocessSampleGamma: Invalid sampled point.");

    SetDPMGamma( smpl );
    return true;
}





// =========================================================================
// MasterCmdInitTest: initialize parallel testing by sending appropriate 
//  command to slaves
//
bool HDPsampler::MasterCmdInitTest()
{
    if( !mMPIIAmMaster())
        return false;
    if( GetBasin() == NULL || GetMenu() == NULL )
        return false;

    const int   nopus = mMPIRingSize() - 1;//number of processing units
    double*     data = GetMBufferValues();
    int len, code;
    mystring preamb = "MasterCmdInitTest: ";
    mystring merror;

    if( data == NULL )
        throw myruntime_error( preamb + "Null data buffer.");
    if( nopus < 1 )
        throw myruntime_error( preamb + "Invalid MPI ring size.");

    SetMBufferHead( THeadMsg );
    SetMBufferCmd( TComInitTest );
    SetMBufferNoVals( 0 );
    len = SetMBufferCRC();
    if( GetMaxSizeOfResMBuffer() <= len || GetMaxSizeOfDataMessage() <= len ) {
        error("MasterCmdInitTest: Short of size in buffers. Increase it and recompile.");
        return false;
    }
//     if(( code = BcastMPIMessage( GetMBuffer(), /*len*/GetMaxSizeOfDataMessage(), 
//                   false/*throw*/)) != MOK )
//         return false;
    if(( code = BcastlenMPIMessage( GetMBuffer(), len, GetMaxSizeOfDataMessage(), 
                  false/*throw*/)) != MOK )
        return false;
 
    return true;
}

// -------------------------------------------------------------------------
// MasterCmdValsTest: receive data from slaves and compare it with that of 
//  master
//
bool HDPsampler::MasterCmdValsTest()
{
    if( !mMPIIAmMaster())
        return false;

    char        strbuf[BUF_MAX];
    const int   nopus = mMPIRingSize() - 1;//number of processing units
    double*     data = GetMBufferValues();
    double*     values = GetResMBufferValues();
    mystring    preamb = "MasterCmdValsTest: ";
    mystring    merror;
    Restaurant* rest;
    Table*      tbl;
    Pslvector*  vec;
    Dish*       dish;
    int head, cmd, novals;
    int len, err;
    int veck, veckk; //vector's dish
    int r, t, v, b, bb;
    int p, n;

    if( mMPIMyRank() < 0 )
        throw myruntime_error( preamb + "Invalid MPI rank.");
    if( nopus < 1 )
        throw myruntime_error( preamb + "Invalid ring size.");
    if( data == NULL || values == NULL )
        throw myruntime_error( preamb + "Null buffers.");

    //wait for data (dish and vector indices)
    for( p = 0; p < nopus; p++ ) {
        //all processes should reply
        err = ReceiveMPIMessage( GetResMBuffer(), GetMaxSizeOfResMBuffer(), 
                                  false/*throw*/);
        if( err != MOK ) {
            sprintf( strbuf, "MasterCmdValsTest: Receive failed: Error code %d.", err );
            merror = strbuf;
            continue;
        }
        if( !merror.empty())
            continue;
        try {
            head = GetResMBufferHead();
            cmd = GetResMBufferCmd();
            novals = GetResMBufferNoVals();
            if( !CheckResMBufferCRC())
                throw myruntime_error("Invalid CRC.");

            if( head != THeadDat )
                throw myruntime_error("Invalid header.");
            if( cmd != TComValsTest )
                throw myruntime_error("Invalid command.");

            for( n = 0; n < novals; ) {
                r = ( int )GetResMBufferValueAt( n++ );//restaurant index
                t = ( int )GetResMBufferValueAt( n++ );//table index
                v = ( int )GetResMBufferValueAt( n++ );//vector index
                b = ( int )GetResMBufferValueAt( n++ );//vector index in basin
                veck = ( int )GetResMBufferValueAt( n++ );//vector's dish index

                rest = GetChain()->GetRestaurantAt( r );
                if( rest == NULL )
                    throw myruntime_error("Invalid restaurant index.");
                tbl = rest->GetTableAt( t );
                if( tbl == NULL )
                    throw myruntime_error("Invalid table index.");
                bb = tbl->GetVectorNIndAt( v );
                if( bb < 0 || bb != b )
                    throw myruntime_error("Invalid vector index in table.");

                vec = tbl->GetVectorNAt( v );
                if( vec == NULL )
                    throw myruntime_error("Invalid vector index.");
                veckk = tbl->GetDishIndex();
                if( veckk < 0 || veckk != veck )
                    throw myruntime_error("Invalid vector's dish index.");
                dish = GetMenu()->GetDishAt( veckk );
                if( dish == NULL )
                    throw myruntime_error("Null vector's dish.");
                ;
            }

        } catch( myexception const& ex ) {
            merror = preamb + ex.what();
            continue;
        }
    }

    if( !merror.empty()) {
        error( merror.c_str());
        return false;
    }

    return true;
}
