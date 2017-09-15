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
#include "ext/rv/rvdisc.h"
#include "HDPsampler.h"

//random number generators
MTRng   RNGsp;
MTRng   RNGsm;
MTRng   RNGsa;

// =========================================================================
// SlaveProcess: slave process
//
bool HDPsampler::SlaveProcess()
{
    if( mMPIIAmMaster())
        return false;

    if( mMPIMyRank() < 0 || mMPIRingSize() < 1 )
        throw myruntime_error( "MasterProcess: Invalid MPI rank or ring size." );

    if( GetBasin() == NULL || GetMenu() == NULL || GetChain() == NULL )
        throw myruntime_error("HDPsampler: SlaveProcess: Memory access error.");

    const int rank = mMPIMyRank();//rank of process
    const int nopus = mMPIRingSize() - 1;//number of processing units
    int norests = GetChain()->GetActualSize();
    int head, cmd, novals;
    int len, err;
    int n, nn;

    while( 1 ) {
//         if(( err = BcastMPIMessage( GetMBuffer(), GetMaxSizeOfDataMessage(), 
//                 true/*throw*/)) != MOK )
//             return false;
        if(( err = BcastlenMPIMessage( GetMBuffer(), len, GetMaxSizeOfDataMessage(), 
                true/*throw*/)) != MOK )
            return false;

        head = GetMBufferHead();
        cmd = GetMBufferCmd();
        novals = GetMBufferNoVals();
        if( !CheckMBufferCRC())
            throw myruntime_error("HDPsampler: SlaveProcess: Invalid CRC.");

        if( head == THeadDat ) {
            if( cmd == TComInitMCMC ) {
                //sends via calling 
                if( !SlaveCmdInitMCMC())
                    return false;
            }
            else if( cmd == TComMCMCcomplete ) {
                if( !ProcessCmdMCMCcomplete())
                    return false;
                //blocks untill all processes reaches this routine!
//                 if( !BlockUntilAllReady( true ))
//                     return false;
            }
            else if( cmd == TComMigrVec ) {
                //sends via calling SlaveProcessCmdSampleVec
                if( !SlaveProcessCmdMigrVec())
                    return false;
            }
            else if( cmd == TComMigrTbl ) {
                //sends via calling SlaveProcessCmdSampleMtx
                if( !SlaveProcessCmdMigrTbl())
                    return false;
            }
            else if( cmd == TComValsKappa0 ) {
                //does not send
                if( !SlaveProcessCmdValsKappa0())
                    return false;
            }
            else if( cmd == TComValsNu0 ) {
                //does not send
                if( !SlaveProcessCmdValsNu0())
                    return false;
            }
            else if( cmd == TComValsConc ) {
                //does not send
                if( !SlaveProcessCmdValsConc())
                    return false;
#ifdef HDPTESTINL
                if( !TestMenu())
                    throw myruntime_error("HDPsampler: SlaveProcess: Test failed.");
#endif
            }
            else
                throw myruntime_error("HDPsampler: SlaveProcess: Invalid command.");
        }
        else if( head == THeadMsg ) {
            if( cmd == TComTerminate )
                break;
            else if( cmd == TComInitVec ) {
                //sends data
                if( !SlaveProcessCmdSampleVec())
                    return false;
            }
            else if( cmd == TComInitMtx ) {
                //sends data
                if( !SlaveProcessCmdSampleMtx())
                    return false;
            }
            else if( cmd == TComInitTest ) {
                //sends data
                if( !SlaveCmdInitTest())
                    return false;
            }
            else
                throw myruntime_error("HDPsampler: SlaveProcess: Invalid command.");
        }
        else
            throw myruntime_error("HDPsampler: SlaveProcess: Invalid header.");
    }
    return true;
}





// =========================================================================
// SlaveCmdInitMCMC: process command start MCMC procedure
//
bool HDPsampler::SlaveCmdInitMCMC()
{
    if( GetBasin() == NULL || GetMenu() == NULL )
        throw myruntime_error("SlaveCmdInitMCMC: Null basin and menu.");

    const int   rank = mMPIMyRank();//rank of process
    const int   nopus = mMPIRingSize() - 1;//number of processing units
    int         novals = GetMBufferNoVals();
    double*     values = GetResMBufferValues();
    Dish*       dsh = NULL;
    Pslvector*  vec = NULL;
    bool    acpt;
    int mm, rasgn;
    int head, cmd;
    int len, err;
    int d, v, dd, vv;

    if( nopus < 1 )
        throw myruntime_error("SlaveCmdInitMCMC: Invalid ring size.");
    if( rank < 1 )
        throw myruntime_error("SlaveCmdInitMCMC: Invalid rank.");
    if( values == NULL )
        throw myruntime_error("SlaveCmdInitMCMC: Null values.");

    SetResMBufferNoVals( 0 );

    for( mm = 0; mm < novals; ) {
        d = ( int )GetMBufferValueAt( mm++ );//dish index
        v = ( int )GetMBufferValueAt( mm++ );//vector index in dish d
        dd = ( int )GetMBufferValueAt( mm++ );//index of the 2nd dish 
        vv = ( int )GetMBufferValueAt( mm++ );//vector index in dish dd

        if( novals < mm )
            throw myruntime_error("SlaveCmdInitMCMC: Invalid no. values in message from master.");

        rasgn = 1 + (( mm >> 2 ) - 1 )% nopus;
        if( rasgn != rank )
            //this pair of indices does not belong to this rank
            continue;

        if( d < 0 || GetMenu()->GetSize() <= d ||
            dd < 0 || GetMenu()->GetSize() <= dd )
            throw myruntime_error("SlaveCmdInitMCMC: Invalid dish indices received.");

        if( d == dd && v == vv )
            throw myruntime_error("SlaveCmdInitMCMC: Vector indices coincide.");

        //if accepted, values will be composed in-line
        if( GetModJNMCMC())
              acpt = MHupdate( d, v, dd, vv );
        else  acpt = MHupdateS( d, v, dd, vv );
        if( !acpt ) {
            //not accepted
            SetResMBufferNoVals( 2 );
            SetResMBufferValueAt( 0, ( double )d );
            SetResMBufferValueAt( 1, ( double )dd );
        }

        //quit once processed at most one pair
        break;
    }

    //compose message
    SetResMBufferHead( THeadDat );
    SetResMBufferCmd( TComMCMCaccpd );
    len = SetResMBufferCRC();
    if( GetMaxSizeOfResMBuffer() <= len || GetMaxSizeOfDataMessage() <= len )
        throw myruntime_error(
            "SlaveCmdInitMCMC: Short of size in buffers. Increase it and recompile.");
    if(( err = SendMPIMessage( GetResMBuffer(), len/*GetMaxSizeOfResMBuffer()*/, 
                  true/*throw*/)) != MOK )
        return false;

    return true;
}





// =========================================================================
// SlaveProcessCmdSampleVec: calculate multivariate probabilities
//  ttt -- table index in restaurant to start with
//  vvv -- vector index in table to start with
//
bool HDPsampler::SlaveProcessCmdSampleVec( int ttt, int vvv )
{
    mystring    preamb = "SlaveProcessCmdSampleVec: ";
    mystring    merror;
    const int rank = mMPIMyRank();//rank of process
    const int nopus = mMPIRingSize() - 1;//number of processing units
    int novals = GetMBufferNoVals();
    int resv, err;  //reserved, error code
    int r, n, nn;

    SetMBufferHead( THeadDat );
    SetMBufferCmd( TComValsVec );

    if( rank < 1 )
        throw myruntime_error( preamb + "Invalid rank.");
    if( nopus < 1 )
        throw myruntime_error( preamb + "Invalid ring size.");

    if( novals == 2 ) {
        n = ( int )GetMBufferValueAt( 0 );
        nn = ( int )GetMBufferValueAt( 1 );
        n += rank - 1;//0th rank is master's
    }
    if( nn < n || novals == 0 ) {
        //anyway, send a message
        SetMBufferNoVals( 0 );
        resv = SetMBufferCRC();
        if( GetMaxSizeOfResMBuffer() <= resv || GetMaxSizeOfDataMessage() <= resv )
            throw myruntime_error(
                preamb + "Short of size in buffers. Increase it and recompile.");
        if(( err = SendMPIMessage( GetMBuffer(), resv/*GetMaxSizeOfDataMessage()*/, 
                      true/*throw*/)) != MOK )
            return false;
        return true;
    }

    if( novals != 2 )
        throw myruntime_error( preamb + "Wrong number of values.");
    if( ttt < 0 || vvv < 0 )
        throw myruntime_error( preamb + "Invalid table and vector indices.");
    if( GetBasin() == NULL || GetMenu() == NULL || GetChain() == NULL )
        throw myruntime_error( preamb + "Memory access error.");

    int norests = GetChain()->GetActualSize();//restaurants do not alter in time
    if( norests <= n )
        throw myruntime_error( preamb + "Invalid restaurant identifier.");

    const int   lastgrp = ( GetParallelProc() == 3 )? nn: n;
    Restaurant* rest, *chrest;
    Table*      tbl, *tblaux, *newtbl;
    Pslvector*  vec;
    bool        processed = false;
    Dish*       dish;
    double      finorm, pinorm;
    double      prob, prprob, newtblprob;
    double*     data = GetMBufferValues();
    double*     values = GetResMBufferValues();
    const int   prd = 5;
    int nods;//no.dishes
    int veck, newk, newt; //dish the vector belongs to, new dish and table indices for vector
    int tbl_asz;//table's actual size
    int taux, t, v, b;
    int d, mm, pp, i;
    int code;

    if( data == NULL || values == NULL )
        throw myruntime_error( preamb + "Null buffer values.");

    //NOTE: STORE menu configuration here (make copy )
    Menu*       orgmenu = GetMenu();
    Menu*       menu = new Menu( *GetMenu());
    if( menu == NULL )
        throw myruntime_error( preamb + "Not enough memory.");
    ReplaceMenu( menu );

    try {
        for( d = 0; d < GetMenu()->GetSize(); d++ ) {
            dish = GetMenu()->GetDishAt( d );
            if( dish == NULL )
                continue;
            dish->SetTmpValue( 0.0 );
        }

        mm = 0;
        for( r = n; r <= lastgrp; r += nopus, ttt = 0 ) {
            chrest = GetChain()->GetRestaurantAt( r );
            if( chrest == NULL )
                throw myruntime_error( preamb + "Null restaurant.");
            //
            //NOTE: make COPY of Restaurant configuration here
            rest = new Restaurant( *chrest );
            if( rest == NULL )
                throw myruntime_error( preamb + "Not enough memory.");
            //
            //set all tables in this restaurant unprocessed
            if( ttt == 0 && vvv == 0 )
                //to start processing
                for( t = 0; t < rest->GetSize(); t++ ) {
                    tbl = rest->GetTableAt( t );
                    if( tbl == NULL )
                        continue;
                    for( v = 0; v < tbl->GetSize(); v++ )
                        tbl->SetProcessedAt( v, false );
                }
            //
            for( t = ttt; t < rest->GetSize() && !processed; t++, vvv = 0 ) {
                tbl = rest->GetTableAt( t );
                if( tbl == NULL )
                    continue;
                tbl->SetMenu( GetMenu());//set temporary menu in temporary restaurant
                for( v = vvv; v < tbl->GetSize() && !processed; v++ ) {
                    //table may have chahge its size
                    tbl_asz = tbl->GetActualSize();
                    b = tbl->GetVectorNIndAt( v );
                    if( b < 0 )//vacant
                        continue;
                    vec = tbl->GetVectorNAt( v );
                    if( vec == NULL )
                        throw myruntime_error( preamb + "Null vector.");
                    if( tbl->GetProcessedAt( v ))
                        continue;//this vector has been processed
                    veck = tbl->GetDishIndex();

                    data[mm++] = ( double )r;//restaurant index
                    data[mm++] = ( double )t;//table index
                    data[mm++] = ( double )v;//vector index
                    data[mm++] = ( double )b;//vector index in basin

                    nods = GetMenu()->GetActualSize();

                    //calculate probability of new table
                    for( d = 0, pp = 0; d < GetMenu()->GetSize(); d++ ) {
                        dish = GetMenu()->GetDishAt( d );
                        if( dish == NULL )
                            continue;
                        if( d == veck )
                            ProbVecOfDishExc( vec, d, &prob );
                        else
                            ProbVecOfDish( vec, d, &prob );
                        dish->SetTmpValue( prob );
                        values[pp++] = -1.0;//new table
                        values[pp++] = ( double )d;//dish index
                        values[pp++] = ( double )( dish->GetNoTables()-(( d == veck && tbl_asz < 2 )? 1: 0 ));
                        values[pp++] = prob;
                        values[pp++] = prob;//for normalized prob. value
                    }
                    if( pp != TIMES5( nods ))
                        throw myruntime_error( preamb + "Inconsistent menu size.");

                    PriorProbVec( vec, &prprob );
                    values[pp++] = -1.0;//new table
                    values[pp++] = -1.0;//new dish
                    values[pp++] = GetDPMGamma();//factor for new table
                    values[pp++] = prprob;
                    values[pp++] = prprob;//for original prob. value

                    //calculate normalized weighted probabilities of dishes (don't sample)
                    if( !SlaveProcessSampleFromProbs( values, pp, &newt, &newk, false/*sample*/ ))
                        return false;

                    newtblprob = 0.0;
                    for( i = 3; i < pp; i += prd ) {
                        finorm = values[i-1];//normalized factor
                        pinorm = values[i+1];//normalized probability
                        if( finorm < 0.0 || 1.0 < finorm || pinorm < 0.0 || 1.0 < pinorm )
                            throw myruntime_error( preamb + "Invalid probabilities in expectation of new table.");
                        newtblprob += finorm * pinorm;//calculate probability
                        d = ( int )values[i-2];
                        if( i + prd < pp ) {
                            if( d < 0 || GetMenu()->GetSize() <= d ||
                              ( dish = GetMenu()->GetDishAt( d )) == NULL )
                                throw myruntime_error( preamb + "Invalid dish index in expectation of new table.");
                            dish->SetTmpValue( pinorm );//save normalized dish probabilities
                        }
                        else
                            prprob = pinorm;//normalized prior probability
                    }

                    if( newtblprob < 0.0 || 1.0 < newtblprob )
                        throw myruntime_error( preamb + "Invalid expectation of new table.");


                    //probabilities for existing tables in restaurant
                    for( taux = 0, pp = 0; taux < rest->GetSize(); taux++ ) {
                        tblaux = rest->GetTableAt( taux );
                        if( tblaux == NULL )
                            continue;
                        dish = tblaux->GetDish();
                        if( dish == NULL )
                            throw myruntime_error( preamb + "Null dish in table.");
                        d = tblaux->GetDishIndex();
                        prob = dish->GetTmpValue();
                        if( prob < 0.0 || 1.0 < prob )
                            throw myruntime_error( preamb + "Invalid normalized probabilities.");
                        ;
                        values[pp++] = ( double )taux;
                        values[pp++] = ( double )d;//table's dish index
                        values[pp++] = ( double )( tblaux->GetActualSize()-(( tbl == tblaux )? 1: 0 ));
                        values[pp++] = prob;
                        values[pp++] = prob;//for normalized prob. value
                    }

                    //NOTE:probability to create new table
                    values[pp++] = -1.0;//new t (table)
                    values[pp++] = -1.0;//new dish for a while (sampling next)
                    values[pp++] = GetDPMTau();//factor for new table
                    values[pp++] = newtblprob;
                    values[pp++] = newtblprob;//for original prob. value

                    if( !SlaveProcessSample( values, pp, &newt, &newk, false/*logs*/ ))
                        return false;

                    if( newt < 0 ) {
                        //sample dish for new table;
                        for( d = 0, pp = 0; d < GetMenu()->GetSize(); d++ ) {
                            dish = GetMenu()->GetDishAt( d );
                            if( dish == NULL )
                                continue;
                            prob = dish->GetTmpValue();
                            if( prob < 0.0 || 1.0 < prob )
                                throw myruntime_error( preamb + "Invalid normalized probabilities.");
                            ;
                            values[pp++] = -1.0;//new table
                            values[pp++] = ( double )d;//dish index
                            values[pp++] = ( double )( dish->GetNoTables()-(( d == veck && tbl_asz < 2 )? 1: 0 ));
                            values[pp++] = prob;
                            values[pp++] = prob;//for normalized prob. value
                        }
                        if( pp != TIMES5( nods ))
                            throw myruntime_error( preamb + "Inconsistent menu size.");

                        //prior for new dish
                        values[pp++] = -1.0;//new t (table)
                        values[pp++] = -1.0;//new k (dish)
                        values[pp++] = GetDPMGamma();//factor of new dish for new table 
                        values[pp++] = prprob;
                        values[pp++] = prprob;//for normalized prob. value

                        if( !SlaveProcessSample( values, pp, &newt, &newk, false/*logs*/ ))
                            return false;
                    }

                    data[mm++] = ( double )newt;//sampled table index for vector
                    data[mm++] = ( double )newk;//sampled dish index for vector

                    //move vector here to ensure valid MCMC chain when GetParallelProc()>1
                    newtbl = NULL;
                    if( 0 <= newt && ( newtbl = rest->GetTableAt( newt )) == NULL )
                        throw myruntime_error( preamb + "Invalid new table index.");
                    if( !MoveVector( rest, t, veck, tbl, v, vec, newt, newk, newtbl ))
                        throw myruntime_error( preamb + "Moving of vector failed.");

                    if( GetParallelProc() == 1 )
                        processed = true;

                    if( tbl == NULL )
                        //this table has been removed (had single vector)
                        break;
                    //tbl->SetProcessedAt( v, true );//has been set in MoveVector
                }
            }
            //
            if( !processed )
                //check whether all vectors in restaurant have been processed
                for( t = 0; t < rest->GetSize(); t++ ) {
                    tbl = rest->GetTableAt( t );
                    if( tbl == NULL )
                        continue;
                    for( v = 0; v < tbl->GetSize(); v++ ) {
                        if( tbl->GetVectorNIndAt( v ) < 0 )
                            continue;
                        if( !tbl->GetProcessedAt( v ))
                          throw myruntime_error( preamb + "Not all vectors processed.");
                    }
                }
            //NOTE: DELETE copy of restaurant
            if( rest ) { delete rest; rest = NULL; }
        }
    } catch( myexception const& ex ) {
        merror = ex.what();
    }

    //NOTE: RESTORE menu configuration
    SetMenu( orgmenu );
    if( rest ) { delete rest; rest = NULL; }

    if( !merror.empty()) {
        throw myruntime_error( preamb + merror );
        return false;
    }

    SetMBufferNoVals( mm );
    resv = SetMBufferCRC();
    if( GetMaxSizeOfResMBuffer() <= resv || GetMaxSizeOfDataMessage() <= resv )
        throw myruntime_error(
            preamb + "Short of size in buffers. Increase it and recompile.");
    if(( err = SendMPIMessage( GetMBuffer(), resv/*GetMaxSizeOfDataMessage()*/, 
                  true/*throw*/)) != MOK )
        return false;
    return true;
}

// =========================================================================
// SlaveProcessCmdSampleVec: calculate multivariate probabilities
//  ttt -- table index in restaurant to start with
//  vvv -- vector index in table to start with
//
bool HDPsampler::SlaveProcessCmdSampleVecObs( int ttt, int vvv )
{
    const int rank = mMPIMyRank();//rank of process
    const int nopus = mMPIRingSize() - 1;//number of processing units
    int novals = GetMBufferNoVals();
    int resv, err;  //reserved, error code
    int n, nn;

    SetMBufferHead( THeadDat );
    SetMBufferCmd( TComValsVec );

    if( rank <= 0 )
        throw myruntime_error("SlaveProcessCmdSampleVec: Invalid rank.");

    if( novals == 2 ) {
        n = ( int )GetMBufferValueAt( 0 );
        nn = ( int )GetMBufferValueAt( 1 );
        n += rank - 1;//0th rank is master's
    }
    if( nn < n || novals == 0 ) {
        //anyway, send a message
        SetMBufferNoVals( 0 );
        resv = SetMBufferCRC();
        if( GetMaxSizeOfResMBuffer() <= resv || GetMaxSizeOfDataMessage() <= resv )
            throw myruntime_error(
                "SlaveProcessCmdSampleVec: Short of size in buffers. Increase it and recompile.");
        if(( err = SendMPIMessage( GetMBuffer(), resv/*GetMaxSizeOfDataMessage()*/, 
                      true/*throw*/)) != MOK )
            return false;
        return true;
    }

    if( novals != 2 )
        throw myruntime_error("SlaveProcessCmdSampleVec: Wrong number of values.");
    if( ttt < 0 || vvv < 0 )
        throw myruntime_error("SlaveProcessCmdSampleVec: Invalid table and vector indices.");
    if( GetBasin() == NULL || GetMenu() == NULL || GetChain() == NULL )
        throw myruntime_error("SlaveProcessCmdSampleVec: Memory access error.");

    //calculate total number of tables over all dishes
    GetMenu()->CalcNoTables();

    int norests = GetChain()->GetActualSize();//restaurants do not alter in time
    if( norests <= n )
        throw myruntime_error("SlaveProcessCmdSampleVec: Invalid restaurant identifier.");

    Restaurant* rest = GetChain()->GetRestaurantAt( n );
    Table*      tbl, *tblaux;
    Pslvector*  vec;
    bool        processed = false;
    Dish*       dish;
    double*     data = GetMBufferValues();
    double*     values = GetResMBufferValues();
    int nods = GetMenu()->GetActualSize();//no.dishes
    int nots = GetMenu()->GetNoTables();//no.tables over all dishes
    int veck, newk, newt; //dish the vector belongs to, new dish and table indices for vector
    int wels = 6;   //elements to write
    int taux, t, v, b;
    int d, mm, pp;

    double  fact = ( double )nots + GetDPMGamma();//factor for size of a table
    double  prob;

    if( rest == NULL )
        throw myruntime_error("SlaveProcessCmdSampleVec: Null restaurant.");
    if( data == NULL || values == NULL )
        throw myruntime_error("SlaveProcessCmdSampleVec: Null buffer values.");

    if( ttt == 0 && vvv == 0 ) {
        //this restaurant is about to start being processed
        for( t = 0; t < rest->GetSize(); t++ ) {
            tbl = rest->GetTableAt( t );
            if( tbl == NULL )
                continue;
            for( v = 0; v < tbl->GetSize(); v++ )
                tbl->SetProcessedAt( v, false );
        }
    }

    mm = 0;
    for( t = ttt; t < rest->GetSize() && !processed; t++, vvv = 0 ) {
        tbl = rest->GetTableAt( t );
        if( tbl == NULL )
            continue;
        for( v = vvv; v < tbl->GetSize() && !processed; v++ ) {
            b = tbl->GetVectorNIndAt( v );
            if( b < 0 )//vacant
                continue;
            vec = tbl->GetVectorNAt( v );
            if( vec == NULL )
                throw myruntime_error("SlaveProcessCmdSampleVec: Null vector.");
            if( tbl->GetProcessedAt( v ))
                continue;//this vector has been processed
            veck = tbl->GetDishIndex();

            data[mm++] = ( double )n;//restaurant index
            data[mm++] = ( double )t;//table index
            data[mm++] = ( double )v;//vector index
            data[mm++] = ( double )b;//vector index in basin

            //calculate probabilities for each dish
            for( d = 0, pp = 0; d < GetMenu()->GetSize(); d++ ) {
                dish = GetMenu()->GetDishAt( d );
                if( dish == NULL )
                    continue;
                if( d == veck )
                    ProbVecOfDishExc( vec, d, &prob );
                else
                    ProbVecOfDish( vec, d, &prob );
                dish->SetTmpValue( prob );

                values[pp++] = -1.0;//new table
                values[pp++] = ( double )d;//dish index
                values[pp++] = ( double )dish->GetNoTables() * GetDPMTau();//factor of new table with this dish
                values[pp++] = prob;
                values[pp++] = prob;//for original prob. value
            }
            if( pp != TIMES5( nods ))
                throw myruntime_error("SlaveProcessCmdSampleVec: Inconsistent menu size.");

            //prior for new table new dish
            PriorProbVec( vec, &prob );
            values[pp++] = -1.0;//new t (table)
            values[pp++] = -1.0;//new k (dish)
            values[pp++] = GetDPMTau() * GetDPMGamma();//factor of new table new dish
            values[pp++] = prob;
            values[pp++] = prob;//for original prob. value

            //probabilities for existing tables in restaurant
            for( taux = 0; taux < rest->GetSize(); taux++ ) {
                tblaux = rest->GetTableAt( taux );
                if( tblaux == NULL )
                    continue;
                dish = tblaux->GetDish();
                if( dish == NULL )
                    throw myruntime_error("SlaveProcessCmdSampleVec: Null dish in table.");
                values[pp++] = ( double )taux;
                values[pp++] = ( double )tblaux->GetDishIndex();
                values[pp++] =(( double )tblaux->GetActualSize()-(( tbl == tblaux )? 1: 0)) * fact;
                values[pp++] = dish->GetTmpValue();
                values[pp++] = dish->GetTmpValue();//for original prob. value
            }

            if( !SlaveProcessSample( values, pp, &newt, &newk ))
                return false;

            data[mm++] = ( double )newt;//sampled table index for vector
            data[mm++] = ( double )newk;//sampled dish index for vector

            tbl->SetProcessedAt( v, true );
            if( GetParallelProc() == 1 )
                processed = true;
        }
    }

    if( !processed ) {
        //check whether all vectors in restaurant have been processed
        for( t = 0; t < rest->GetSize(); t++ ) {
            tbl = rest->GetTableAt( t );
            if( tbl == NULL )
                continue;
            for( v = 0; v < tbl->GetSize(); v++ ) {
                if( tbl->GetVectorNIndAt( v ) < 0 )
                    continue;
                if( !tbl->GetProcessedAt( v ))
                  throw myruntime_error("SlaveProcessCmdSampleVec: Not all vectors processed.");
            }
        }
    }

    SetMBufferNoVals( mm );
    resv = SetMBufferCRC();
    if( GetMaxSizeOfResMBuffer() <= resv || GetMaxSizeOfDataMessage() <= resv )
        throw myruntime_error(
            "SlaveProcessCmdSampleVec: Short of size in buffers. Increase it and recompile.");
    if(( err = SendMPIMessage( GetMBuffer(), resv/*GetMaxSizeOfDataMessage()*/, 
                  true/*throw*/)) != MOK )
        return false;
    return true;
}





// -------------------------------------------------------------------------
// SlaveProcessSample: sample from multivariate probabilities
//  values -- log probability values
//  mp -- number of values
//  newt -- index of sampled table
//  newk -- index of sampled dish
//  logs -- `values' contain log probability values
//
bool HDPsampler::SlaveProcessSample( 
    double* values, int mp, int* newt, int* newk, bool logs )
{
    return SlaveProcessSampleFromProbs( values, mp, newt, newk, true, NULL, logs );
}

// -------------------------------------------------------------------------
// SlaveProcessSampleFromProbs: sample from multivariate discrete distribution
//  values -- log probability values
//  mp -- number of values
//  newt -- index of sampled table
//  newk -- index of sampled dish
//  sample -- flag of sampling (whether to sample)
//  ndx -- address to save index at
//  logs -- `values' contain log probability values
//
bool HDPsampler::SlaveProcessSampleFromProbs( 
    double* values, int mp, int* newt, int* newk, bool sample, int* ndx, bool logs )
{
    mystring    preamb = "SlaveProcessSampleFromProbs: ";
    const int   prd = 5;
    double  lprob;
    double  minlp = 0.0;
    double  maxlp = -LOC_DBL_MAX;
    double  normc = 0.0, fnorm = 0.0, pnorm = 0.0;
    double  F = 0.0, pf;//distribution function value
    double  ru = RNGsp.GetDouble1();//r.n. excluding 1
    double  fact;
    int n;

    if( values == NULL || newt == NULL || newk == NULL )
        throw myruntime_error( preamb + "Memory access error.");
    if( ru < 0.0 || 1.0 < ru )
        throw myruntime_error( preamb + "Invalid RNG value.");

    for( n = 3; n < mp; n += prd ) {
        lprob = values[n];
        //density values can be >1
        if( !isfinite( lprob ))
            throw myruntime_error( preamb + "Invalid probability value.");
        if( lprob < minlp )
            minlp = lprob;
        if( maxlp < lprob )
            maxlp = lprob;
    }
    for( n = 3; n < mp; n += prd ) {
        fact = values[n-1];
        if( fact < 0 )
            throw myruntime_error( preamb + "Invalid factor.");
        fnorm += fact;
        if( logs ) {
            values[n] -= maxlp;
            if( values[n] < cdp_LOG_DBL_MIN || !isfinite(values[n]))
                values[n] = 0.0;
            else
                values[n] = exp( values[n]);
        }
        pnorm += values[n];
    }

    if( fnorm <= 0.0 || pnorm <= 0.0 )
        throw myruntime_error( preamb + "Failed to normalize probabilities.");
    for( n = 3; n < mp; n += prd ) {
        values[n-1] /= fnorm;//save normalized factors
        values[n+1] = values[n] / pnorm;//save normalized probabilities
        values[n] = values[n-1] * values[n+1];
        normc += values[n];
    }
    if( normc <= 0.0 )
        throw myruntime_error( preamb + "Failed to normalize probabilities.");

    for( n = 3; n < mp; n += prd ) {
        pf = F;
        if( values[n])
            F += ( values[n] /= normc );
        if( sample && pf <= ru && ru < F ) {
            *newt = ( int )values[n-3];//table index
            *newk = ( int )values[n-2];//dish index
            if( ndx )
                *ndx = n - 3;
            return true;
        }
    }
    if( !sample )
        return true;
    throw myruntime_error( preamb + "Sampling failed.");
    return false;
}





// =========================================================================
// SlaveProcessCmdMigrVec: process command of migration of vectors
//
bool HDPsampler::SlaveProcessCmdMigrVec()
{
    if( GetBasin() == NULL || GetMenu() == NULL || GetChain() == NULL )
        throw myruntime_error("SlaveProcessCmdMigrVec: Memory access error.");

    mystring    preamb = "SlaveProcessCmdMigrVec: ";
    const int   rank = mMPIMyRank();//rank of process
    const int   nopus = mMPIRingSize() - 1;//number of processing units
    int         novals = GetMBufferNoVals();
    int         norests = GetChain()->GetActualSize();
    double*     values = GetResMBufferValues();
    const int   slaveid = rank - 1;
    Restaurant* rest = NULL;
    Table*      tbl = NULL;
    Table*      newtbl = NULL;
    Dish*       newdsh = NULL;
    Pslvector*  vec = NULL;
    int head, cmd;
    int len, err;
    int rb, re, r, t, v, b, k;
    int newt, newk;
    int d, mm, mp;

    if( rank < 1 )
        throw myruntime_error( preamb + "Invalid rank.");
    if( nopus < 1 )
        throw myruntime_error( preamb + "Invalid ring size.");

    values[0] = values[1] = values[2] = values[3] = -1;

    for( mm = 0; mm < novals; ) {
        rest = NULL;
        tbl = newtbl = NULL;
        newdsh = NULL;

        rb = ( int )GetMBufferValueAt( mm++ );//begin index of restaurants
        re = ( int )GetMBufferValueAt( mm++ );//end index of restaurants
        r = ( int )GetMBufferValueAt( mm++ );//restaurant index
        t = ( int )GetMBufferValueAt( mm++ );//table index
        v = ( int )GetMBufferValueAt( mm++ );//vector index
        b = ( int )GetMBufferValueAt( mm++ );//vector index in basin
        newt = ( int )GetMBufferValueAt( mm++ );//table for vector to migrate to
        newk = ( int )GetMBufferValueAt( mm++ );//dish for vector to have at table

        if( novals < mm )
            throw myruntime_error( preamb + "Invalid no. values in message from master.");

        if( r < 0 || GetChain()->GetSize() <= r ||
          ( rest = GetChain()->GetRestaurantAt( r )) == NULL )
            throw myruntime_error( preamb + "Invalid restaurant index received.");

        if( t < 0 || rest->GetSize() <= t ||
          ( tbl = rest->GetTableAt( t )) == NULL )
            throw myruntime_error( preamb + "Invalid table index.");

        if(( k = tbl->GetDishIndex()) < 0 || GetMenu()->GetSize() <= k )
            throw myruntime_error( preamb + "Invalid table's dish index.");

        if( v < 0 || tbl->GetSize() <= v ||
          ( vec = tbl->GetVectorNAt( v )) == NULL )
            throw myruntime_error( preamb + "Invalid vector index received.");

        if( b < 0 || GetBasin()->GetSize() <= b ||
            vec != GetBasin()->GetValueAt( b ))
            throw myruntime_error( preamb + "Inconsistent vector index.");

        if( 0 <= newt && ( rest->GetSize() <= newt ||
          ( newtbl = rest->GetTableAt( newt )) == NULL ))
            throw myruntime_error( preamb + "Invalid table index received.");

        if( 0 <= newk && ( GetMenu()->GetSize() <= newk ||
          ( newdsh = GetMenu()->GetDishAt( newk )) == NULL ))
            throw myruntime_error( preamb + "Invalid dish index received.");

        if( newt == t && newk == k )
            //proceed further to on the last index processed
            ;//continue;

        if( !MoveVector( rest, t, k, tbl, v, vec, newt, newk, newtbl ))
            throw myruntime_error( preamb + "Moving of vector failed.");

        if( r == rb + rank - 1 ) {
//             if( 0 <= values[0])
//                 throw myruntime_error(
//                     "SlaveProcessCmdMigrVec: Two restaurant indices at a time.");
            values[0] = ( double )rb;
            values[1] = ( double )re;
            //other values serve as params
            values[2] = ( double )t;
            values[3] = ( double )( v + 1 );//start with next vector in the table
        }
    }

    if( novals != mm )
        throw myruntime_error( preamb + "Unexpected number of values.");

    //pretend true message
    SetMBufferHead( THeadMsg );
    SetMBufferCmd( TComInitVec );

    if( values[0] < 0 || GetParallelProc() != 1 )
        //processed either one or all groups (restaurants)
        SetMBufferNoVals( 0 );
    else {
        SetMBufferNoVals( 2 );
        SetMBufferValueAt( 0, values[0] );
        SetMBufferValueAt( 1, values[1] );
    }
    //don't set CRC intentionally

    if( !SlaveProcessCmdSampleVec(( int )values[2], ( int )values[3]))
        return false;

    return true;
}





// =========================================================================
// SlaveProcessCmdSampleMtx: calculate matrix probabilities
//  ttt -- table index in restaurant to start with
//
bool HDPsampler::SlaveProcessCmdSampleMtx( int ttt )
{
    mystring    preamb = "SlaveProcessCmdSampleMtx: ";
    mystring    merror;
    const int rank = mMPIMyRank();//rank of process
    const int nopus = mMPIRingSize() - 1;//number of processing units
    int novals = GetMBufferNoVals();
    int resv, err;  //reserved, error code
    int r, n, nn;

    SetMBufferHead( THeadDat );
    SetMBufferCmd( TComValsMtx );

    if( rank < 1 )
        throw myruntime_error( preamb + "Invalid rank.");
    if( nopus < 1 )
        throw myruntime_error( preamb + "Invalid ring size.");

    if( novals == 2 ) {
        n = ( int )GetMBufferValueAt( 0 );
        nn = ( int )GetMBufferValueAt( 1 );
        n += rank - 1;//0th rank is master's
    }
    if( nn < n || novals == 0 ) {
        //anyway, send a message
        SetMBufferNoVals( 0 );
        resv = SetMBufferCRC();
        if( GetMaxSizeOfResMBuffer() <= resv || GetMaxSizeOfDataMessage() <= resv )
            throw myruntime_error(
                preamb + "Short of size in buffers. Increase it and recompile.");
        if(( err = SendMPIMessage( GetMBuffer(), resv/*GetMaxSizeOfDataMessage()*/, 
                      true/*throw*/)) != MOK )
            return false;
        return true;
    }

    if( novals != 2 )
        throw myruntime_error( preamb + "Wrong number of values.");
    if( ttt < 0 )
        throw myruntime_error( preamb + "Invalid table index.");
    if( GetBasin() == NULL || GetMenu() == NULL || GetChain() == NULL )
        throw myruntime_error( preamb + "Memory access error.");

    int norests = GetChain()->GetActualSize();//restaurants do not alter in time
    if( norests <= n )
        throw myruntime_error( preamb + "Invalid restaurant identifier.");

    const int   lastgrp = ( GetParallelProc() == 3 )? nn: n;
    Restaurant* rest, *chrest;
    Table*      tbl, *tblaux;
//     Pslvector*  vec;
    bool        processed = false;
    Dish*       dish;
    double      prob;
    double*     data = GetMBufferValues();
    double*     values = GetResMBufferValues();
    int nods;//no.dishes
    int tblk, svnk, newk; //table dish, new dish index for table
    int t, newt;
    int d, mm, pp;

    if( data == NULL || values == NULL )
        throw myruntime_error( preamb + "Null buffer values.");

    //NOTE: STORE menu configuration here (make copy )
    Menu*       orgmenu = GetMenu();
    Menu*       menu = new Menu( *GetMenu());
    if( menu == NULL )
        throw myruntime_error( preamb + "Not enough memory.");
    ReplaceMenu( menu );

    try {
        for( d = 0; d < GetMenu()->GetSize(); d++ ) {
            dish = GetMenu()->GetDishAt( d );
            if( dish == NULL )
                continue;
            dish->SetTmpValue( 0.0 );
        }

        mm = 0;
        for( r = n; r <= lastgrp; r += nopus, ttt = 0 ) {
            chrest = GetChain()->GetRestaurantAt( r );
            if( chrest == NULL )
                throw myruntime_error( preamb + "Null restaurant.");
            //
            //NOTE: make COPY of Restaurant configuration here
            rest = new Restaurant( *chrest );
            if( rest == NULL )
                throw myruntime_error( preamb + "Not enough memory.");
            //
            //set all tables in this restaurant unprocessed
            if( ttt == 0 )
                //to start processing
                for( t = 0; t < rest->GetSize(); t++ ) {
                    tbl = rest->GetTableAt( t );
                    if( tbl == NULL )
                        continue;
                    tbl->SetProcessed( false );
                }
            //
            for( t = ttt; t < rest->GetSize() && !processed; t++ ) {
                tbl = rest->GetTableAt( t );
                if( tbl == NULL )
                    continue;

                tbl->SetMenu( GetMenu());//set temporary menu in temporary restaurant
                if( tbl->GetProcessed())
                    continue;//this table has been processed
                tblk = tbl->GetDishIndex();

                data[mm++] = ( double )r;//restaurant index
                data[mm++] = ( double )t;//table index

                nods = GetMenu()->GetActualSize();

                //calculate probabilities for each dish
                for( d = 0, pp = 0; d < GetMenu()->GetSize(); d++ ) {
                    dish = GetMenu()->GetDishAt( d );
                    if( dish == NULL )
                        continue;
                    if( d == tblk )
                        ProbMtxOfDishExc( tbl, d, &prob );
                    else
                        ProbMtxOfDish( tbl, d, &prob );
                    dish->SetTmpValue( prob );

                    if( dish->GetNoTables() <= 0 )
                        throw myruntime_error( preamb + "Dish not served at any of tables.");
                    values[pp++] = -1.0;//fake table index (for comp.)
                    values[pp++] = ( double )d;//dish index
                    values[pp++] = ( double )( dish->GetNoTables() - (( d == tblk )? 1: 0 ));
                    values[pp++] = prob;
                    values[pp++] = prob;//for normalized prob. value
                }
                if( pp != TIMES5( nods ))
                    throw myruntime_error( preamb + "Inconsistent menu size.");

                //prior for new dish of table
                PriorProbMtx( tbl, &prob );
                values[pp++] = -1.0;//fake table index (for comp.)
                values[pp++] = -1.0;//new k (dish)
                values[pp++] = GetDPMGamma();//factor of new dish
                values[pp++] = prob;
                values[pp++] = prob;//for original prob. value

                if( !SlaveProcessSample( values, pp, &newt/*comp.*/, &newk ))
                    return false;

                data[mm++] = ( double )newk;//sampled dish index for table

                //move table here to ensure valid MCMC chain when GetParallelProc()>1
                if( !MoveTable( rest, t, tblk, tbl, newk ))
                    throw myruntime_error( preamb + "Moving of table failed.");

                //newk changed if it's been <0; write assigned value for new dish
                data[mm++] = ( double )newk;

                if( GetParallelProc() == 1 )
                    processed = true;
                //tbl->SetProcessed( true );//has been set in MoveTable
            }

            if( !processed )
                //check whether all tables in restaurant have been processed
                for( t = 0; t < rest->GetSize(); t++ ) {
                    tbl = rest->GetTableAt( t );
                    if( tbl == NULL )
                        continue;
                    if( !tbl->GetProcessed())
                      throw myruntime_error( preamb + "Not all tables processed.");
                }
            //NOTE: DELETE copy of restaurant
            if( rest ) { delete rest; rest = NULL; }
        }
    } catch( myexception const& ex ) {
        merror = ex.what();
    }

    //NOTE: RESTORE menu configuration
    SetMenu( orgmenu );
    if( rest ) { delete rest; rest = NULL; }

    if( !merror.empty()) {
        throw myruntime_error( preamb + merror );
        return false;
    }

    SetMBufferNoVals( mm );
    resv = SetMBufferCRC();
    if( GetMaxSizeOfResMBuffer() <= resv || GetMaxSizeOfDataMessage() <= resv )
        throw myruntime_error(
            preamb + "Short of size in buffers. Increase it and recompile.");
    if(( err = SendMPIMessage( GetMBuffer(), resv/*GetMaxSizeOfDataMessage()*/, 
                  true/*throw*/)) != MOK )
        return false;
    return true;
}





// =========================================================================
// SlaveProcessCmdMigrTbl: process command of migration of table
//
bool HDPsampler::SlaveProcessCmdMigrTbl()
{
    if( GetBasin() == NULL || GetMenu() == NULL || GetChain() == NULL )
        throw myruntime_error("SlaveProcessCmdMigrTbl: Memory access error.");

    mystring    preamb = "SlaveProcessCmdMigrTbl: ";
    const int   rank = mMPIMyRank();//rank of process
    const int   nopus = mMPIRingSize() - 1;//number of processing units
    int         novals = GetMBufferNoVals();
    int         norests = GetChain()->GetActualSize();
    double*     values = GetResMBufferValues();
    Restaurant* rest = NULL;
    Table*      tbl = NULL;
    Table*      newtbl = NULL;
    Dish*       newdsh = NULL;
    Pslvector*  vec = NULL;
    int head, cmd;
    int len, err;
    int rb, re, r, t, v, b, k;
    int newt, newk;
    int d, mm, mp;

    if( rank < 1 )
        throw myruntime_error( preamb + "Invalid rank.");
    if( nopus < 1 )
        throw myruntime_error( preamb + "Invalid ring size.");

    values[0] = values[1] = values[2] = -1;

    for( mm = 0; mm < novals; ) {
        rest = NULL;
        tbl = newtbl = NULL;
        newdsh = NULL;

        rb = ( int )GetMBufferValueAt( mm++ );//begin index of restaurants
        re = ( int )GetMBufferValueAt( mm++ );//end index of restaurants
        r = ( int )GetMBufferValueAt( mm++ );//restaurant index
        t = ( int )GetMBufferValueAt( mm++ );//table index
        newk = ( int )GetMBufferValueAt( mm++ );//dish to have at table

        if( novals < mm )
            throw myruntime_error( preamb + "Invalid no. values in message from master.");

        if( r < 0 || GetChain()->GetSize() <= r ||
          ( rest = GetChain()->GetRestaurantAt( r )) == NULL )
            throw myruntime_error( preamb + "Invalid restaurant index received.");

        if( t < 0 || rest->GetSize() <= t ||
          ( tbl = rest->GetTableAt( t )) == NULL )
            throw myruntime_error( preamb + "Invalid table index.");

        if(( k = tbl->GetDishIndex()) < 0 || GetMenu()->GetSize() <= k )
            throw myruntime_error( preamb + "Invalid table's dish index.");

        if( 0 <= newk && ( GetMenu()->GetSize() <= newk ||
          ( newdsh = GetMenu()->GetDishAt( newk )) == NULL ))
            throw myruntime_error( preamb + "Invalid dish index received.");

        if( newk == k )
            //proceed further to on the last index processed
            ;//continue;

        if( !MoveTable( rest, t, k, tbl, newk ))
            throw myruntime_error( preamb + "Moving of table failed.");

        if( r == rb + rank - 1 ) {
//             if( 0 <= values[0])
//                 throw myruntime_error(
//                     "SlaveProcessCmdMigrTbl: Two restaurant indices at a time.");
            values[0] = ( double )rb;
            values[1] = ( double )re;
            //other values serve as params
            values[2] = ( double )( t + 1 );
        }
    }

    if( novals != mm )
        throw myruntime_error( preamb + "Unexpected number of values.");

    //pretend true message
    SetMBufferHead( THeadMsg );
    SetMBufferCmd( TComInitMtx );
    if( values[0] < 0 || GetParallelProc() != 1 )
        SetMBufferNoVals( 0 );
    else {
        SetMBufferNoVals( 2 );
        SetMBufferValueAt( 0, values[0] );
        SetMBufferValueAt( 1, values[1] );
    }
    //don't set CRC intentionally

    if( !SlaveProcessCmdSampleMtx(( int )values[2]))
        return false;

    return true;
}





// =========================================================================
// SlaveProcessCmdValsKappa0: process command of receiving new value of 
//  hyperparameter kappa0
//
bool HDPsampler::SlaveProcessCmdValsKappa0()
{
    const int   rank = mMPIMyRank();//rank of process
    const int   nopus = mMPIRingSize() - 1;//number of processing units
    double      newk0;
    int         novals = GetMBufferNoVals();

    if( rank <= 0 )
        throw myruntime_error("SlaveProcessCmdValsKappa0: Invalid rank.");
    if( novals != 1 )
        throw myruntime_error("SlaveProcessCmdValsKappa0: Invalid number of values.");

    newk0 = GetMBufferValueAt( 0 );

    if( newk0 <= 0.0 )
        throw myruntime_error("SlaveProcessCmdValsKappa0: Invalid new value of kappa0.");

    //NOTE: do not set new value of kappa0 here; 
    //current value is required for adjustments
    if( !ProcessCmdNewKappa0( newk0 ))
        return false;
    return true;
}





// =========================================================================
// SlaveProcessCmdValsNu0: process command of receiving new value of 
//  hyperparameter nu0
//
bool HDPsampler::SlaveProcessCmdValsNu0()
{
    const int   rank = mMPIMyRank();//rank of process
    const int   nopus = mMPIRingSize() - 1;//number of processing units
    double      newnu0;
    int         novals = GetMBufferNoVals();

    if( rank <= 0 )
        throw myruntime_error("SlaveProcessCmdValsNu0: Invalid rank.");
    if( novals != 1 )
        throw myruntime_error("SlaveProcessCmdValsNu0: Invalid number of values.");

    newnu0 = GetMBufferValueAt( 0 );

    if( newnu0 <= 0.0 )
        throw myruntime_error("SlaveProcessCmdValsNu0: Invalid new value of nu0.");

    //NOTE: do not set new value of nu0 here
    if( !ProcessCmdNewNu0( newnu0 ))
        return false;
    return true;
}





// =========================================================================
// SlaveProcessCmdValsConc: process command of receiving values of 
//  concentration parameters
//
bool HDPsampler::SlaveProcessCmdValsConc()
{
    const int   rank = mMPIMyRank();//rank of process
    const int   nopus = mMPIRingSize() - 1;//number of processing units
    double      tau, gamma;
    int         novals = GetMBufferNoVals();

    if( rank <= 0 )
        throw myruntime_error("SlaveProcessCmdValsConc: Invalid rank.");

    tau = GetMBufferValueAt( 0 );
    gamma = GetMBufferValueAt( 1 );

    if( tau <= 0.0 )
        throw myruntime_error("SlaveProcessCmdValsConc: Invalid parameter tau.");
    if( gamma <= 0.0 )
        throw myruntime_error("SlaveProcessCmdValsConc: Invalid parameter gamma.");

    SetDPMTau( tau );
    SetDPMGamma( gamma );

    return true;
}





// =========================================================================
// SlaveCmdInitTest: sollect HDP configuration data and send it to master
//
bool HDPsampler::SlaveCmdInitTest()
{
    const int rank = mMPIMyRank();//rank of process
    const int nopus = mMPIRingSize() - 1;//number of processing units
    int novals = GetMBufferNoVals();
    int resv, err;  //reserved, error code

    SetMBufferHead( THeadDat );
    SetMBufferCmd( TComValsTest );

    if( rank <= 0 )
        throw myruntime_error("SlaveCmdInitTest: Invalid rank.");

    if( novals != 0 )
        throw myruntime_error("SlaveCmdInitTest: Invalid no. values.");

    if( GetBasin() == NULL || GetMenu() == NULL || GetChain() == NULL )
        throw myruntime_error("SlaveCmdInitTest: Null HDP structures.");

    const int   maxdatsize = 400000;//maximum data size in bytes
    int wels = 5;//elements to write
    int period = 1;
    static int  s_pos = 0;

    Restaurant* rest;
    Table*      tbl;
    Pslvector*  vec;
    Dish*       dish;
    double*     data = GetMBufferValues();
    double*     values = GetResMBufferValues();
    int totv = GetBasin()->GetActualSize();
    int nods = GetMenu()->GetActualSize();//no.dishes
    int nots = GetMenu()->GetNoTables();//no.tables over all dishes
    int veck; //vector's dish
    int r, t, v, b;
    int mm;

    if( data == NULL || values == NULL )
        throw myruntime_error("SlaveCmdInitTest: Null buffer values.");

    //set step (period) to get vectors at periodically
    if( maxdatsize < ( mm = totv * wels * sizeof( double ))) {
        period = SLC_MAX( 1, mm / maxdatsize );
        if( period <= s_pos )
            s_pos = 0;
    }

    mm = 0;
    for( r = 0; r < GetChain()->GetSize(); r++ ) {
        rest = GetChain()->GetRestaurantAt( r );
        if( rest == NULL )
            throw myruntime_error("SlaveCmdInitTest: Null restaurant.");
        for( t = 0; t < rest->GetSize(); t++ ) {
            tbl = rest->GetTableAt( t );
            if( tbl == NULL )
                continue;
            for( v = 0; v < tbl->GetSize(); v++ ) {
                b = tbl->GetVectorNIndAt( v );
                if( b < 0 )//vacant
                    continue;
                vec = tbl->GetVectorNAt( v );
                if( vec == NULL )
                    throw myruntime_error("SlaveCmdInitTest: Null vector.");
                veck = tbl->GetDishIndex();
                if( veck < 0 )
                    throw myruntime_error("SlaveCmdInitTest: Unassoc. vector's dish index.");
                dish = GetMenu()->GetDishAt( veck );
                if( dish == NULL )
                    throw myruntime_error("SlaveCmdInitTest: Null vector's dish.");

                if( GetMaxSizeOfDataMessage() <= ( mm + 5 )* sizeof( double ))
                    throw myruntime_error(
                        "SlaveCmdInitTest: Short of size in buffers. Increase it and recompile.");

                if( b % period - s_pos == 0 ) {
                    //pass only fraction of vectors to lower data transfer load
                    data[mm++] = ( double )r;//restaurant index
                    data[mm++] = ( double )t;//table index
                    data[mm++] = ( double )v;//vector index
                    data[mm++] = ( double )b;//vector index in basin
                    data[mm++] = ( double )veck;//vector's dish index
                }
            }
        }
    }

    //increase phase (position) to take next fraction later
    s_pos++;

    SetMBufferNoVals( mm );
    resv = SetMBufferCRC();
    if( GetMaxSizeOfResMBuffer() <= resv || GetMaxSizeOfDataMessage() <= resv )
        throw myruntime_error(
            "SlaveCmdInitTest: Short of size in buffers. Increase it and recompile.");
    if(( err = SendMPIMessage( GetMBuffer(), resv/*GetMaxSizeOfDataMessage()*/, 
                  true/*throw*/)) != MOK )
        return false;
    return true;
}
