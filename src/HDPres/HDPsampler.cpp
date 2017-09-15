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
#include "HDPsampler.h"

const double    cdp_LOG_DBL_MIN = log( LOC_DBL_MIN );
const double    cdp_LOG_DBL_MAX = log( LOC_DBL_MAX );

// -------------------------------------------------------------------------
// constructor: initialization
//
HDPsampler::HDPsampler()
:   filename_( NULL ),
    outputfile_( NULL ),
    parsread_( false ),
    parproc_( 1 ),
    noiters_( 0 ),
    noMHRGSscans_( 0 ),
    noMHups_( 0 ),
    smprops_( 0 ),
    jnmodsmpl_( false ),
    unfdishMHups_( false ),
    unfvectMHups_( false ),
    notblsmpl_( false ),
    taupreset_( -1.0 ),
    gammapreset_( -1.0 ),
    kappa0preset_( -1.0 ),
    nu0preset_( -1.0 ),
    mbuffer_( NULL ),
    resmbuffer_( NULL ),
    locintbuf_( NULL ),
    szlocintbuf_( 0 )
{
    SetFMinMsgSize( &HDPsampler::GetMinSizeOfMessage );
    SetFIsMsgValid( &HDPsampler::AreDataValid );
}

// -------------------------------------------------------------------------
// destructor:
//
HDPsampler::~HDPsampler()
{
    DestroyLocIntBuf();
    DestroyMBuffer();
    DestroyResMBuffer();
}

// -------------------------------------------------------------------------
// TellSlavesToTerminate: broadcast message to slave processes to tell of
//     termination
//
void HDPsampler::TellSlavesToTerminate( bool throw_on_error )
{
    int code, len;
    if( !mMPIIAmMaster())
        return;
    len = ( int )FormatTermMessage();
    code = BcastlenMPIMessage(
                GetMBuffer(),
                len,
                GetMaxSizeOfDataMessage(),//lengths of master and slave messages must match
                throw_on_error
    );
    if( code != MOK )
        if( throw_on_error )
            throw myruntime_error( "HDPsampler: Terminating MPI broadcast failed." );
}

// -------------------------------------------------------------------------
// Run: begin sampling procedure
//
long HDPsampler::Run()
{
    PrivateMPIInit();

    char        strbuf[BUF_MAX];
    bool        needsending = false;            //no need for sending of message on error
    static char s_errmsg[MSG_MAX];              //string buffer for error message
    int         s_szemsg = 0;                   //size of message
    int         szlimit = MSG_MAX;
    mystring    throwstr;

    const double    tol = 0.1;//error tolerance for target probabilities
    int         ctx = 1;
    int         dim = 1;
    Pslvector   mvec( dim );
    time_t      tm;
    int a, c;

    time( &tm );

    sprintf( s_errmsg, "HDPsampler(%d): ", mMPIMyRank());
    szlimit -= strlen( s_errmsg );

    if( GetDegFAdjustment() < -1.0 )
        throw myruntime_error("Invalid value of adjustment to deg. of freedom.");

    try {
        MasterMessage("Read...");
        SetParsRead( false );
        Read();

        //initialize buffers after read to have definite sizes
        InitMBuffer();
        InitResMBuffer();
        if( GetBasin())
            InitLocIntBuf( TIMES2( GetBasin()->GetSize()));
        RNGi.Set(( unsigned long )( size_t )this +( unsigned long )( size_t )( &RNGi ) +( unsigned long )tm);
        RNGc.Set(( unsigned long )( size_t )this +( unsigned long )( size_t )( &RNGc ) +( unsigned long )tm);
        RNGsp.Set(( unsigned long )( size_t )this +( unsigned long )( size_t )( &RNGsp ) +( unsigned long )tm);
        RNGsm.Set(( unsigned long )( size_t )this +( unsigned long )( size_t )( &RNGsm ) +( unsigned long )tm);
        RNGsa.Set(( unsigned long )( size_t )this +( unsigned long )( size_t )( &RNGsa ) +( unsigned long )tm);
        RNGjl.Set(( unsigned long )( size_t )this +( unsigned long )( size_t )( &RNGjl ) +( unsigned long )tm);
        RNGja.Set(( unsigned long )( size_t )this +( unsigned long )( size_t )( &RNGja ) +( unsigned long )tm);

        if( !GetParsRead()) {
            //NOTE: get dimensions and context size after read!
            ctx = GetCtxtSize();
            if( GetMenu()) {
                dim = GetMenu()->GetDim();
            }
            if( dim < 1 || ctx < 1 )
                throw myruntime_error("Invalid dimensions and/or context size.");

            if( GetUninfPrior()) {
                MasterMessage("Setting uninformative prior parameters...");
                Pslvector   mean( dim );
                //set mean vector to 0 meaning that vectors in 
                //frequency space at all context positions have all equal terms
                mean.SetAllToValue( 0.0 );
                SetUninfPriorParamsS0I( dim, ctx, mean, GetKappa0Preset(), GetNu0Preset());
            }
            else {
                MasterMessage("Calculating prior parameters...");
                CalcPriorParams();
            }

            MasterMessage("Calculating parameters of initial clusters...");
            RecalcMenuParams();
        }

        MasterMessage("Sampling...");
        MSSampling();

    } catch( myexception const& ex )
    {
        if( mMPIIAmMaster()) {
            //the throw will terminate slave processes
            error( ex.what());
            throwstr = ex.what();
        } else {
            //inform the master on error if needed and silently exit
            if( needsending ) {
                s_szemsg = strlen( ex.what());
                if( szlimit <= s_szemsg )
                    s_szemsg = szlimit - 1;

                memcpy( s_errmsg + strlen( s_errmsg ), ex.what(), s_szemsg );
                s_errmsg[s_szemsg] = 0;
                if( SendMPIMessage( s_errmsg, s_szemsg, false/*throw_on_error*/ ) != MOK )
                    throwstr = ex.what();
            }
            else
                throwstr = ex.what();
        }
    }

    if( !throwstr.empty())
        throw myruntime_error( throwstr );

    PrivateMPITerminate();
    return 0;
}

// -------------------------------------------------------------------------
// MSSampling: master-slave sampling
//
bool HDPsampler::MSSampling()
{
    bool code = true;
    if( mMPIIAmMaster()) {
        code = MasterProcess();
        if( code )
            MasterMessage("Finished.");
        TellSlavesToTerminate( true );
    }
    else {
        code = SlaveProcess();
    }
    return code;
}





// =========================================================================
// MoveVector: move vector from one table to another;
//  n,nn -- table indices in restaurant `rest';
//  v -- vector index in table `from';
//  d,dd -- respective table dishes;
//  adjust -- adjust dish parameters
//
bool HDPsampler::MoveVector( 
    Restaurant* rest, 
    int n, int& d, Table*& from, const int v, const Pslvector* vec, 
    int& nn, int& dd, Table*& to, bool adjust )
{
    int     ndx, ddx;
    int     mloc, rloc, tloc, dloc;
    Dish*   dish = NULL;

    if( GetMenu() == NULL )
        throw myruntime_error( "MoveVector: Null menu." );
    if( rest == NULL || from == NULL || vec == NULL ||
        ( 0 <= nn && to == NULL ) || ( nn < 0 && to != NULL ))
        throw myruntime_error("MoveVector: Memory access error.");
    if( n < 0 || rest->GetSize() <= n || rest->GetSize() <= nn )
        throw myruntime_error("MoveVector: Invalid table indices.");
    if( rest->GetTableAt( n ) != from || ( 0 <= nn && rest->GetTableAt( nn ) != to ))
        throw myruntime_error("MoveVector: Invalid table adresses.");
    if( from->GetDishIndex() != d || ( to && to->GetDishIndex() != dd ))
        throw myruntime_error("MoveVector: Invalid table dish ids.");
    if( v < 0 || from->GetSize() <= v )
        throw myruntime_error("MoveVector: Invalid vector index.");
    if( from->GetVectorNAt( v ) != vec )
        throw myruntime_error("MoveVector: Invalid vector address.");

    if( from == to ) {
        from->SetProcessedAt( v, true );
        return true;
    }

    ndx = from->GetVectorNIndAt( v );
    ddx = from->GetVectorNDishIndAt( v );
    if( 1 < from->GetDish()->GetActualSize()) {
        //update if dish has more than 1 member, otherwise it'll be removed
        //IMPORTANT: remove vector from 1.table and 2.dish after updating dish params
        if( adjust )
            AdjustDishParams( d, vec, false/*subtract*/);
    }
    else {
        if( from->GetActualSize() != 1 )
            throw myruntime_error("MoveVector: Table is expected to have a single vector.");
        //table and dish have one vector each;
        //if that single vector is to be moved to the same dish, do nothing
        if( dd == d ) {
            from->SetProcessedAt( v, true );
            return true;
        }
    }
    from->RemValueAt( v, vec );
    from->GetDish()->RemValueAt( ddx, vec );
    //if there're no vectors left in table
    if( from->GetActualSize() <= 0 ) {
        rest->RemTableAt( n, from );
        //if this table is the only one having dish k
        if( from->GetDish()->GetActualSize() <= 0 ) {
            if( 0 < from->GetDish()->GetNoTables())
                throw myruntime_error("MoveVector: Dish's not expected to be served on tables.");
            GetMenu()->RemDishAt( from->GetDishIndex(), from->GetDish());
//             delete from->GetDish();//RemDishAt will delete dish
            if( dd == d )
                dd = -1;
            d = -1;
        }
        delete from;
        from = NULL;
        n = -1;
    }
    else if( from->GetDish()->GetActualSize() <= 0 )
            throw myruntime_error("MoveVector: Null dish while table is not.");
    rloc = nn;
    if( to == NULL ) {
        //create new table
        to = new Table( GetDefTableSize());
        if( to == NULL )
            throw myruntime_error("MoveVector: Not enough memory for new table.");
        to->SetBasin( GetBasin());
        to->SetMenu( GetMenu());
        mloc = dd;
        if( dd < 0 ) {
            //create new dish
            dish = new Dish( GetDefDishSize());
            if( dish == NULL )
                throw myruntime_error( "MoveVector: Not enough memory for new dish." );
            dish->SetBasin( GetBasin());
            dd = mloc = GetMenu()->NewDish( dish );
        }
        to->SetDishIndex( dd );
        nn = rloc = rest->NewTable( to );
    }
    //IMPORTANT: add vector to 1.dish and 2.table after updating dish params
    if( adjust )
        AdjustDishParams( dd, vec, true/*add*/);
    dloc = to->GetDish()->NewVectorNInd( ndx );
    tloc = to->NewVectorNInd( ndx, dloc );
    to->SetProcessedAt( tloc, true );
    return true;
}

// -------------------------------------------------------------------------
// AdjustDishParams: adjust dish's k parameters assuming vector `vec' has 
//  been added (`add') or subtracted from the set of vectors assigned to 
//  that dish k
//
void HDPsampler::AdjustDishParams( int k, const Pslvector* vec, bool add )
{
    if( GetMenu() == NULL || vec == NULL )
        throw myruntime_error( "AdjustDishParams: Memory access error." );

    const int           dim = GetMenu()->GetDim();
    const double        kp0 = GetMenu()->GetKappa0();
    const Pslvector*    mu0 = GetMenu()->GetMu0();
    const SPDmatrix*    L0 = GetMenu()->GetS0();
    const SPDmatrix*    L0inv = GetMenu()->GetInvS0();
    const double        ldet0 = GetMenu()->GetLDetS0();
    Dish*               dish = NULL;

    if( mu0 == NULL || L0 == NULL )
        throw myruntime_error( "AdjustDishParams: Null prior parameters." );
    if( dim <= 0 || kp0 <= 0/*dim*/ )
        throw myruntime_error( "AdjustDishParams: Invalid dimensionality." );

    if( k < 0 || GetMenu()->GetSize() <= k || 
      ( dish = GetMenu()->GetDishAt( k )) == NULL )
        throw myruntime_error( "AdjustDishParams: Memory access error." );

    Pslvector*          mu = GetMenu()->GetMuVectorAt( k );
    SPDmatrix*          L = GetMenu()->GetSMatrixAt( k );
    SPDmatrix           S( dim );
    SPDmatrix*          Linv = &S;
    Pslvector           dev;
    int                 nk = dish->GetActualSize();
    double              kpnk = ( double )nk + kp0;
    double              kpnk1 = kpnk + ( add? 1.0: -1.0 );
    double              kpnkokpnk1 = kpnk / kpnk1;
    mystring    preamb = "AdjustDishParams: ";
    mystring    errstr;
    double      ldet = 0.0;
    int         err;

    if( nk < 0 )
        throw myruntime_error( "AdjustDishParams: Invalid number of vectors." );

#ifdef PROBMTXINVSM
    Linv = GetMenu()->GetInvSMatrixAt( k );
#endif
    if( nk == 0 ) {
        //there're no parameters set for this dish yet
        GetMenu()->SetMuVectorAt( k, mu = new Pslvector( dim ));
        GetMenu()->SetSMatrixAt( k, L = new SPDmatrix( dim ));
#ifdef PROBMTXINVSM
        GetMenu()->SetInvSMatrixAt( k, Linv = new SPDmatrix( dim ));
#endif
        if( mu == NULL || L == NULL || Linv == NULL )
            throw myruntime_error( "AdjustDishParams: Not enough memory." );
        *mu = *mu0;
        *L = *L0;
    }

    if( mu == NULL || L == NULL || Linv == NULL )
        throw myruntime_error( "AdjustDishParams: Null dish parameters." );

    try {
        if( add ) {
            //ADD vector to dish
            //first, update scale matrix since we need not updated mean vector `mu'
            dev = *mu;
            if(( err = dev.Superposition( -1.0, *vec )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            if(( err = S.Mul( dev, dev )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            if(( err = S.Scale( kpnkokpnk1 )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            if(( err = L->Add( S )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            //update mean vector `mu'
            if(( err = mu->MultiplyBy( kpnkokpnk1 )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
            if(( err = mu->Superposition( 1.0 / kpnk1, *vec )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
        }
        else {
            //SUBTRACT vector from dish
            if( nk == 0 )
                throw myruntime_error( "Subtraction of vector from null dish." );
            if( nk == 1 ) {
                //IMPORTANT: assuming vec is equal to the only one in dish
                *mu = *mu0;
                *L = *L0;
                GetMenu()->SetLDetSMAt( k, ldet0 );
#ifdef PROBMTXINVSM
                *Linv = *L0inv;
#endif
                return;
            }

            //first, update mean vector `mu' to mu_{n_k-1}
            if(( err = mu->MultiplyBy( kpnkokpnk1 )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
            if(( err = mu->Superposition( -1.0 / kpnk1, *vec )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            //update scale matrix using updated mean vector mu_{n_k-1}
            dev = *mu;
            if(( err = dev.Superposition( -1.0, *vec )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            if(( err = S.Mul( dev, dev )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            if(( err = S.Scale( kpnk1 / kpnk )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            if(( err = L->Sub( S )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
        }

        //solve for inverse matrix and calculate determinant
        *Linv = *L;
        if(( err = Linv->CholeskyDecompose()) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = Linv->CDedLogDet( &ldet )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));
#ifdef PROBMTXINVSM
        if(( err = Linv->CDedInvert()) != 0 )
            throw myruntime_error( TranslatePSLError( err ));
#endif
        //SAVE
        //save log-determinant only;
        //mean vector, scale matrix, and its inverse have already been updated
        GetMenu()->SetLDetSMAt( k, ldet );

    } catch( myexception const& ex ) {
        errstr = preamb + ex.what();
    }
    if( !errstr.empty())
        throw myruntime_error( errstr );
}

// -------------------------------------------------------------------------
// AdjustDishParams: adjust dish's k parameters assuming table `tbl' has 
//  been added (`add') or subtracted from the set of vectors assigned to 
//  that dish k
//
void HDPsampler::AdjustDishParams( int k, const Table* tbl, bool add )
{
    if( GetMenu() == NULL || tbl == NULL )
        throw myruntime_error( "AdjustDishParams: Memory access error." );

    const int           dim = GetMenu()->GetDim();
    const double        kp0 = GetMenu()->GetKappa0();
    const Pslvector*    mu0 = GetMenu()->GetMu0();
    const SPDmatrix*    L0 = GetMenu()->GetS0();
    const SPDmatrix*    L0inv = GetMenu()->GetInvS0();
    const double        ldet0 = GetMenu()->GetLDetS0();
    Dish*               dish = NULL;

    if( mu0 == NULL || L0 == NULL )
        throw myruntime_error( "AdjustDishParams: Null prior parameters." );
    if( dim <= 0 || kp0 <= 0/*dim*/ )
        throw myruntime_error( "AdjustDishParams: Invalid dimensionality." );

    if( k < 0 || GetMenu()->GetSize() <= k || 
      ( dish = GetMenu()->GetDishAt( k )) == NULL )
        throw myruntime_error( "AdjustDishParams: Memory access error." );

    Pslvector*          mu = GetMenu()->GetMuVectorAt( k );
    SPDmatrix*          L = GetMenu()->GetSMatrixAt( k );
    SPDmatrix           S( dim );
    SPDmatrix*          Linv = &S;
    Pslvector*          vec;
    Pslvector           dev( dim ), mnew( dim );
    const int           tsz = tbl->GetSize();
    const int           nov = tbl->GetActualSize();//number of vectors in table
    const int           nk = dish->GetActualSize();
    const double        kp = ( double )nk + kp0;
    const double        kppv = kp + ( double )( add? nov: -nov );
    const double        kpokppv = kp / ( kppv? kppv: 1.0 );
    mystring    preamb = "AdjustDishParams: ";
    mystring    errstr;
    double      ldet = 0.0;
    int         v, err;

    if( nk < 0 )
        throw myruntime_error( "AdjustDishParams: Invalid number of vectors in dish." );
    if( nov <= 0 )
        throw myruntime_error( "AdjustDishParams: Invalid number of vectors in table." );
    if( nov == 1 ) {
        //if table has only one memeber
        for( v = 0; v < tsz; v++ ) {
            if( tbl->GetVectorNIndAt( v ) < 0 )
                continue;
            vec = tbl->GetVectorNAt( v );
            if( vec == NULL )
                continue;
            AdjustDishParams( k, vec, add );
            break;
        }
        return;
    }


#ifdef PROBMTXINVSM
    Linv = GetMenu()->GetInvSMatrixAt( k );
#endif
    if( nk == 0 ) {
        //there're no parameters set for this dish yet
        GetMenu()->SetMuVectorAt( k, mu = new Pslvector( dim ));
        GetMenu()->SetSMatrixAt( k, L = new SPDmatrix( dim ));
#ifdef PROBMTXINVSM
        GetMenu()->SetInvSMatrixAt( k, Linv = new SPDmatrix( dim ));
#endif
        if( mu == NULL || L == NULL || Linv == NULL )
            throw myruntime_error( "AdjustDishParams: Not enough memory." );
        *mu = *mu0;
        *L = *L0;
    }

    if( mu == NULL || L == NULL || Linv == NULL )
        throw myruntime_error( "AdjustDishParams: Null dish parameters." );

    if( add == false ) {
        //if subtract
        if( nk < nov )
            throw myruntime_error( "AdjustDishParams: Subtraction of table from smaller dish." );
        if( nk == 0 )
            throw myruntime_error( "Subtraction of table from null dish." );
        if( nk == nov ) {
            //IMPORTANT: assuming `tbl' is equal to the only table a dish has
            *mu = *mu0;
            *L = *L0;
            GetMenu()->SetLDetSMAt( k, ldet0 );
#ifdef PROBMTXINVSM
            *Linv = *L0inv;
#endif
            return;
        }
    }

    try {
        //precalculate mean vector of table
        for( v = 0; v < tsz; v++ )
        {
            if( tbl->GetVectorNIndAt( v ) < 0 )
                continue;
            vec = tbl->GetVectorNAt( v );
            if( vec == NULL )
                continue;
            if(( err = mnew.Superposition( 1.0, *vec )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
        }

        if( add ) {
            //ADD table to dish
            //first, update scale matrix since we need not updated mean vector `mu'
            for( v = 0; v < tsz; v++ )
            {
                if( tbl->GetVectorNIndAt( v ) < 0 )
                    continue;
                vec = tbl->GetVectorNAt( v );
                if( vec == NULL )
                    continue;
                dev = *vec;
                if(( err = dev.Superposition( -1.0, *mu )) != 0 )
                    throw myruntime_error( TranslatePSLError( err ));

                if(( err = S.Mul( dev, dev )) != 0 )
                    throw myruntime_error( TranslatePSLError( err ));

                if(( err = L->Add( S )) != 0 )
                    throw myruntime_error( TranslatePSLError( err ));
            }

            dev = mnew;
            if(( err = dev.Superposition(( double )-nov, *mu )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            if(( err = S.Mul( dev, dev )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            if(( err = S.Scale( -1.0 / kppv )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            if(( err = L->Add( S )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            //update mean vector `mu'
            if(( err = mu->MultiplyBy( kpokppv )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
            if(( err = mu->Superposition( 1.0 / kppv, mnew )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
        }
        else {
            //SUBTRACT table from dish
            //first, update mean vector `mu' to mu_{n_k-L}
            if(( err = mu->MultiplyBy( kpokppv )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
            if(( err = mu->Superposition( -1.0 / kppv, mnew )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            //update scale matrix using updated mean vector `mu'
            for( v = 0; v < tsz; v++ )
            {
                if( tbl->GetVectorNIndAt( v ) < 0 )
                    continue;
                vec = tbl->GetVectorNAt( v );
                if( vec == NULL )
                    continue;
                dev = *vec;
                if(( err = dev.Superposition( -1.0, *mu )) != 0 )
                    throw myruntime_error( TranslatePSLError( err ));

                if(( err = S.Mul( dev, dev )) != 0 )
                    throw myruntime_error( TranslatePSLError( err ));

                if(( err = L->Sub( S )) != 0 )
                    throw myruntime_error( TranslatePSLError( err ));
            }

            dev = mnew;
            if(( err = dev.Superposition(( double )-nov, *mu )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            if(( err = S.Mul( dev, dev )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            if(( err = S.Scale( 1.0 / kp )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));

            if(( err = L->Add( S )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
        }

        //solve for inverse matrix and calculate determinant
        *Linv = *L;
        if(( err = Linv->CholeskyDecompose()) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = Linv->CDedLogDet( &ldet )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));
#ifdef PROBMTXINVSM
        if(( err = Linv->CDedInvert()) != 0 )
            throw myruntime_error( TranslatePSLError( err ));
#endif
        //SAVE
        //save log-determinant only;
        //mean vector, scale matrix, and its inverse have already been updated
        GetMenu()->SetLDetSMAt( k, ldet );

    } catch( myexception const& ex ) {
        errstr = preamb + ex.what();
    }
    if( !errstr.empty())
        throw myruntime_error( errstr );
}





// =========================================================================
// MoveTable: assign dish for table;
//  n -- table index in restaurant `rest';
//  d,dd -- table's current and new dishes;
//
bool HDPsampler::MoveTable( 
    Restaurant* rest, int n, int d, Table* from, int& dd, bool adjust )
{
    int     ndx, ddx, v;
    int     mloc, dloc;
    int     notbls;
    Dish*   dish = NULL;
    Pslvector* vec = NULL;

    if( GetMenu() == NULL )
        throw myruntime_error( "MoveTable: Null menu." );
    if( rest == NULL || from == NULL )
        throw myruntime_error("MoveTable: Memory access error.");
    if( n < 0 || rest->GetSize() <= n )
        throw myruntime_error("MoveTable: Invalid table index.");
    if( rest->GetTableAt( n ) != from )
        throw myruntime_error("MoveTable: Invalid table adress.");
    if( d < 0 || GetMenu()->GetSize() <= d || GetMenu()->GetSize() <= dd )
        throw myruntime_error("MoveTable: Invalid dish indices.");
    if( from->GetDishIndex() != d )
        throw myruntime_error("MoveTable: Invalid table's dish id.");
    if( 0 <= dd && ( dish = GetMenu()->GetDishAt( dd )) == NULL )
        throw myruntime_error("MoveTable: Null new table's dish.");

    if( d == dd ) {
        from->SetProcessed( true );
        return true;
    }
    notbls = from->GetDish()->GetNoTables();
    if( notbls < 1 )
        throw myruntime_error("MoveTable: Invalid no. tables with the given dish.");
    if( from->GetDish()->GetActualSize() < from->GetActualSize())
        throw myruntime_error("MoveTable: Table's size greater than that of dish.");
    if( notbls == 1 ) {
        //there's only one table (`from') the dish is served on
        if( dd < 0 ) {
            //the only table of the dish has to be served with new dish;
            //this is the same as to leave it unchanged
            from->SetProcessed( true );
            dd = d;
            return true;
        }
        if( from->GetDish()->GetActualSize() != from->GetActualSize())
            throw myruntime_error("MoveTable: Inconsistency of sizes of table and dish.");
        GetMenu()->RemDishAt( from->GetDishIndex(), from->GetDish());
//         delete from->GetDish();//RemDishAt will delete dish
        if( dd == d )
            dd = -1;//the dish has been removed
    }
    else {
        //adjust dish params only if no. tables >1; 
        //otherwise dish is removed
        if( adjust )
            AdjustDishParams( d, from, false/*subtract*/);
        //remove all vectors at the table from the dish
        for( v = 0; v < from->GetSize(); v++ ) {
            ndx = from->GetVectorNIndAt( v );
            if( ndx < 0 )
                continue;
            vec = from->GetVectorNAt( v );
            if( vec == NULL )
                throw myruntime_error("MoveTable: Null vector at table.");
            ddx = from->GetVectorNDishIndAt( v );
            from->GetDish()->RemValueAt( ddx, vec );
        }
        from->GetDish()->DecNoTables();
        if( from->GetDish()->GetNoTables() <= 0 )
            throw myruntime_error("MoveTable: Invalid no. tables of dish obtained.");
    }

    if( dd < 0 ) {
        //create new dish
        dish = new Dish( GetDefDishSize());
        if( dish == NULL )
            throw myruntime_error( "MoveTable: Not enough memory for new dish." );
        dish->SetBasin( GetBasin());
        dd = mloc = GetMenu()->NewDish( dish );
    }
    from->SetDishIndex( dd );//FIRST, change dish index (move table to new dish)
    if( adjust )
        AdjustDishParams( dd, from, true/*add*/);
    for( v = 0; v < from->GetSize(); v++ ) {
        ndx = from->GetVectorNIndAt( v );
        if( ndx < 0 )
            continue;
        vec = from->GetVectorNAt( v );
        if( vec == NULL )
            throw myruntime_error("MoveTable: Null vector at table.");
        dloc = from->GetDish()->NewVectorNInd( ndx );//add vector to dish
        from->SetVectorNDishIndAt( v, dloc );//change dish index for vector
    }
    from->GetDish()->IncNoTables();
    from->SetProcessed( true );
    return true;
}





// =========================================================================
// ProcessCmdMCMCcomplete: perform split/merge of dish(es) determined by 
//  MCMC procedure
//
bool HDPsampler::ProcessCmdMCMCcomplete( int* noups )
{
    double* data = GetMBufferValues();
    int     totvals = GetMBufferNoVals();
    int novals;
    int n;

    if( data == NULL )
        throw myruntime_error("ProcessCmdMCMCcomplete: Null data buffer.");
    if( noups )
        *noups = 0;
    for( n = 0; n < totvals; ) {
        novals = ( int )data[n++];
        if( !SplitMergeDish( data + n, novals ))
            return false;
        n += novals;
        if( noups && 2 < novals )
            (*noups)++;
    }

    return true;
}

// -------------------------------------------------------------------------
// SplitMergeDish: perform explicit split/merge of dish(es); required 
//  indices are given in data buffer
//
bool HDPsampler::SplitMergeDish( double* values, int nn )
{
    bool    adjinturn = true;//adjust dish params in turn of each vector
    mystring errstr;
    Restaurant* rest;
    Table*  tbl, *tblto;
    Dish*   dfrom, *ddto = NULL;
    const Pslvector* vec;
    bool    split;
    int*    varray = NULL;
    int     d, dd;//dish from and dish to
    int     bn, vtn;
    int     r, t, tt;
    int n = 0;

    if( GetBasin() == NULL || GetMenu() == NULL || GetChain() == NULL )
        throw myruntime_error("SplitMergeDish: Null restaurant structures.");
    if( values == NULL )
        throw myruntime_error("SplitMergeDish: Null values.");
    if( nn <= 0 )
        return true;
    if( nn < 2 )
        throw myruntime_error("SplitMergeDish: Invalid size of values.");

    d = ( int )values[n++];
    dd = ( int )values[n++];

    split = ( d != dd && dd < 0 );

    if( d < 0 || GetMenu()->GetSize() <= d ||
        GetMenu()->GetSize() <= dd )
        throw myruntime_error("SplitMergeDish: Invalid dish indices in values.");

    if(( dfrom = GetMenu()->GetDishAt( d )) == NULL ||
       ( 0 <= dd && ( ddto = GetMenu()->GetDishAt( dd )) == NULL ))
        throw myruntime_error("SplitMergeDish: Null dishes at index locations.");

    if( nn <= 2 )
        return true;

    //{{dicide how to update dish parameters (to get favourable complexity)
    //(1 call of `RecalcDishParams' compares favourably to 2N calls of `AdjustDishParams')
    int novmv = nn - 2;//number of vectors to move
    int szdfrom = dfrom->GetActualSize();
    int szddto = ( ddto? ddto->GetActualSize(): 0 );
    if( szdfrom + szddto <= novmv + novmv )
        //if twice the no. vectors (2 calls of `AdjustDishParams' per vector) exceeds 
        //overall no. vectors in dishes, then
        //turn off adjustment of params in turn of each vector;
        //when MERGE state is proposed, one of dishes will be merged so that
        //just 1 call of `RecalcDishParams' suffices 
        //(complexity of `AdjustDishParams(add/sub--table)' and `RecalcDishParams' 
        // matches only if table size is comparable to that of dish)
        adjinturn = false;
    //}}

    varray = ( int* )malloc( GetBasin()->GetSize() * sizeof( int ));
    if( varray == NULL )
        throw myruntime_error( "SplitMergeDish: Not enough memory.");
    memset( varray, 0, GetBasin()->GetSize() * sizeof( int ));

    try {
        for( ; n < nn; n++ ) {
            bn = ( int )values[n];
            if( bn < 0 || GetBasin()->GetSize() <= bn )
                throw myruntime_error("SplitMergeDish: Invalid vector index in values.");
            if(( vec = GetBasin()->GetValueAt( bn )) == NULL )
                throw myruntime_error("SplitMergeDish: Null vector at index locations.");
            varray[bn] = 1;
        }
        //iterate over all restaurants in chain
        for( r = 0; r < GetChain()->GetSize(); r++ ) {
            rest = GetChain()->GetRestaurantAt( r );
            if( rest == NULL )
                continue;
//             //NOTE: obsolete scheme
//             //one table per restaurant is assigned for split/merge vectors
//             // (even in case of merge state)
//             tblto = NULL;
//             tt = -1;
            for( t = 0; t < rest->GetSize(); t++ ) {
                //NOTE:table is created for each existing table processed when
                //splitting to new dish
                tblto = NULL;
                tt = -1;
                tbl = rest->GetTableAt( t );
                if( tbl == NULL || tbl == tblto )
                    continue;
                if( d != tbl->GetDishIndex())
                    continue;
                for( vtn = 0; vtn < tbl->GetSize(); vtn++ ) {
                    bn = tbl->GetVectorNIndAt( vtn );
                    if( bn < 0 )
                        continue;
                    vec = tbl->GetVectorNAt( vtn );
                    if( GetBasin()->GetSize() <= bn )
                        throw myruntime_error("SplitMergeDish: Invalid vector index at table.");
                    ;
                    if( !split ) {
                        //in merge state, do nothing with vectors
                        if( !varray[bn])
                            throw myruntime_error("SplitMergeDish: Vector unexpectedly not to be moved.");
                        varray[bn] = 0;
                        continue;
                    }
                    if( !varray[bn])
                        continue;
                    ;
                    //dish index `d' will be set to -1 once it is removed
                    //dish index `dd' will be updated to new dish index if it's <0
                    if( !MoveVector( rest, t, d, tbl, vtn, vec,  tt, dd, tblto, adjinturn ))
                        throw myruntime_error("SplitMergeDish: Failed to move vector.");
                    varray[bn] = 0;
                    nn--;
                    if( tbl == NULL )
                        //the table had the only vector and has been removed
                        break;
                }
                if( !split ) {
                    //move whole table (relabel table to belong to other dish)
                    if( !MoveTable( rest, t, d, tbl, dd, adjinturn ))
                        throw myruntime_error("SplitMergeDish: Failed to relabel tabel.");
                }
            }
        }
        if( split ) {
            if( nn != 2 )
                throw myruntime_error("SplitMergeDish: Not all vectors processed.");
        }
        else {
            if( GetMenu()->GetDishAt( d ) != NULL )
                //if dish d still has members
                throw myruntime_error("SplitMergeDish: Dish after merging expected to have no members.");
            d = -1;
        }
        ;
        if( !adjinturn ) {
            //adjust dish parameters
            if( 0 <= d )
                //if this dish hasn't been removed
                RecalcDishParams( d );//a dish vectors moved from
            RecalcDishParams( dd );//a dish vectors moved to
        }
    } catch( myexception const& ex ) {
        errstr = ex.what();
    }

    if( varray ) {
        free( varray );
        varray = NULL;
    }
    if( !errstr.empty())
        throw myruntime_error( errstr );

    return true;
}





// =========================================================================
// ProcessCmdNewKappa0: adjust dish parameters according to new value of 
//  kappa0 and save this new value for kappa0
//
bool HDPsampler::ProcessCmdNewKappa0( double newk0 )
{
    if( GetMenu() == NULL || GetBasin() == NULL )
        throw myruntime_error("ProcessCmdNewKappa0: Null menu structures.");

    mystring    preamb = "ProcessCmdNewKappa0: ";
    const int   dim = GetMenu()->GetDim();
    int         nodshs = GetMenu()->GetActualSize();//number of dishes
    int         szmenu = GetMenu()->GetSize();
    double      k0 = GetMenu()->GetKappa0();
    double      nk;//number of samples in dish
    double      fact, ldet;
    const Pslvector*    mu0 = GetMenu()->GetMu0();
    Pslvector*          mu;//dish mean vector
    SPDmatrix*          L, *Linv;//dish scale matrix and its inverse
    const Dish*         dish;
    int k, err;

    if( mu0 == NULL )
        throw myruntime_error( preamb + "Null prior mean vector.");

    if( dim < 1 || nodshs < 1 )
        throw myruntime_error( preamb + "No dishes.");

    if( newk0 == k0 )
        //nothing to do
        return true;

    SPDmatrix       devmtx( dim );
    Pslvector       eta_k;//mean of sample vectors for dish d

    for( k = 0; k < szmenu; k++ ) {
        dish = GetMenu()->GetDishAt( k );
        if( dish == NULL )
            continue;
        //calculate mean of sample vectors for dish d
        nk = ( double )dish->GetActualSize();
        mu = GetMenu()->GetMuVectorAt( k );
        L = GetMenu()->GetSMatrixAt( k );
        Linv = &devmtx;
#ifdef PROBMTXINVSM
        Linv = GetMenu()->GetInvSMatrixAt( k );
        if( Linv == NULL )
            throw myruntime_error( preamb + "Null inverse of scale matrix for dish.");
#endif
        if( mu == NULL || nk < 1 )
            throw myruntime_error( preamb + "Null mean vector of dish.");
        if( L == NULL )
            throw myruntime_error( preamb + "Null scale matrix of dish.");
        eta_k = *mu;
        if(( err = eta_k.Superposition( -1.0, *mu0 )) != 0 )
            throw myruntime_error( preamb + TranslatePSLError( err ));
        if(( err = eta_k.MultiplyBy( (k0+nk)/nk )) != 0 )
            throw myruntime_error( preamb + TranslatePSLError( err ));
        //NOTE:eta_k is now mean of samples - mu0
        //set adjustment matrix
        if(( err = devmtx.Mul( eta_k, eta_k )) != 0 )
            throw myruntime_error( preamb + TranslatePSLError( err ));
        fact = newk0*nk/(newk0+nk) - k0*nk/(k0+nk);
        if(( err = devmtx.Scale( fact )) != 0 )
            throw myruntime_error( preamb + TranslatePSLError( err ));
        //adjust scale matrix and calculate its inverse and determinant
        if(( err = L->Add( devmtx )) != 0 )
            throw myruntime_error( preamb + TranslatePSLError( err ));
        //solve for inverse matrix and calculate determinant
        *Linv = *L;
        if(( err = Linv->CholeskyDecompose()) != 0 )
            throw myruntime_error( preamb + TranslatePSLError( err ));
        if(( err = Linv->CDedLogDet( &ldet )) != 0 )
            throw myruntime_error( preamb + TranslatePSLError( err ));
#ifdef PROBMTXINVSM
        if(( err = Linv->CDedInvert()) != 0 )
            throw myruntime_error( preamb + TranslatePSLError( err ));
#endif
        //save log-determinant; matrices adjusted above
        GetMenu()->SetLDetSMAt( k, ldet );
        //adjust mean vector for dish
        //NOTE:eta_k is mean of samples - mu0
        *mu = eta_k;
        if(( err = mu->MultiplyBy( nk/(newk0+nk) )) != 0 )
            throw myruntime_error( preamb + TranslatePSLError( err ));
        if(( err = mu->Superposition( 1.0, *mu0 )) != 0 )
            throw myruntime_error( preamb + TranslatePSLError( err ));
    }

    //save new value of kappa0
    GetMenu()->SetKappa0( newk0 );
    //recalculate factor for prior probabilities of vectors
    CalcPriorProbFact();
    return true;
}





// =========================================================================
// ProcessCmdNewNu0: adjust dependent parameters according to new value of 
//  nu0 and save this new value for nu0
//
bool HDPsampler::ProcessCmdNewNu0( double newnu0 )
{
    if( GetMenu() == NULL || GetBasin() == NULL )
        throw myruntime_error("ProcessCmdNewNu0: Null menu structures.");

    mystring    preamb = "ProcessCmdNewNu0: ";
    const int   dim = GetMenu()->GetDim();
    double      nu0 = GetMenu()->GetNu0();

    if( dim < 1 )
        throw myruntime_error( preamb + "Invalid dimensions.");

    if( newnu0 == nu0 )
        //nothing to do
        return true;

    //save new value of nu0
    GetMenu()->SetNu0( newnu0 );
    //recalculate factor for prior probabilities of vectors
    CalcPriorProbFact();
    return true;
}





// =========================================================================
// FORMATING of MESSAGES ROUTINES
//
// FormatMessage: format message to be sent;
//     returns number of bytes written to the string stream
//
// NOTE: it is assumed that space enough to keep formatted data is
//     preallocated
// -------------------------------------------------------------------------

size_t HDPsampler::FormatMessage( 
    char* strstream, int header, int cmd, 
    int novals, double* values )
{
    size_t  written;
    int     n;

    if( !strstream )
        throw myruntime_error( "FormatMessage: Unable to format message." );

    if( novals < 0 || values == NULL )
        novals = 0;

    SetBufferHead( strstream, header );
    SetBufferCmd( strstream, cmd );
    SetBufferNoVals( strstream, novals );

    for( n = 0; n < novals; n++ )
        SetBufferValueAt( strstream, n, values[n]);

    written = SetBufferCRC( strstream );
    return written;
}
