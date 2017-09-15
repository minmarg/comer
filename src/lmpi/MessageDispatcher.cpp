/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <stdlib.h>
#include <string.h>

#include "mpiloc.h"
#include "MessageDispatcher.h"

#include "mystring.h"
#include "myexcept.h"



int MessageDispatcher::mpi_rank_ = -1;       //rank of this process in the MPI ring
int MessageDispatcher::mpi_size_ = -1;       //number of processes running in the MPI ring

MessageDispatcher::TMinMsgSizeDet   MessageDispatcher::fGetMinMsgSize_ = NULL;
MessageDispatcher::TMsgValidator    MessageDispatcher::fIsMsgValid_ = NULL;


////////////////////////////////////////////////////////////////////////////
// CLASS MessageDispatcher
//
// Constructor
//
MessageDispatcher::MessageDispatcher()
{
}

// Destructor
//
MessageDispatcher::~MessageDispatcher()
{
}

// -------------------------------------------------------------------------
// PrivateMPIInit: initializes MPI envorinment
// -------------------------------------------------------------------------

void MessageDispatcher::PrivateMPIInit()
{
    if( mpi_rank_ != -1 || mpi_size_ != -1 )
            throw myruntime_error( mystring( "MessageDispatcher: Unable to initialize MPI environment." ));

    if( MPI_Init( __PARGC__, __PARGV__ ) != MPI_SUCCESS )
            throw myruntime_error( mystring( "MessageDispatcher: Failed to initialize MPI environment." ));

    if( MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank_ ) != MPI_SUCCESS )
            throw myruntime_error( mystring( "MessageDispatcher: Failed to determine rank in the MPI ring." ));

    if( MPI_Comm_size( MPI_COMM_WORLD, &mpi_size_ ) != MPI_SUCCESS )
            throw myruntime_error( mystring( "MessageDispatcher: Failed to determine size of the MPI ring." ));

    if( mpi_size_ < 2 )
            throw myruntime_error( mystring( "MessageDispatcher: MPI ring size is too small for computations." ));
}

// -------------------------------------------------------------------------
// PrivateMPITerminate: finalizes MPI envorinment. Must be called for each
//     PrivateMPIInit call, respectively!
// -------------------------------------------------------------------------

void MessageDispatcher::PrivateMPITerminate()
{
    if( MPI_Finalize() != MPI_SUCCESS )
        error( "MessageDispatcher: Failed to finalize MPI environment." );

    mpi_rank_ = -1;
    mpi_size_ = -1;
}

// -------------------------------------------------------------------------
// BcastMPIMessage: broadcasts a message from the master process to the
//     slave processes; the slaves block until receive the message;
//     on successful return code is MPI_SUCCESS
// -------------------------------------------------------------------------

int MessageDispatcher::BcastMPIMessage( char* msg, size_t length, bool throw_on_error )
{
    int         code = MPI_SUCCESS;
    static char errbuf[KBYTE];

    code = MPI_Bcast( msg, length, MPI_CHAR, mMPIMasterRank(), MPI_COMM_WORLD );

    if( code == MPI_SUCCESS )
        return MOK;

    sprintf( errbuf, "MessageDispatcher(%d): Erroneous MPI broadcast - %d.", mMPIMyRank(), code );

    if( throw_on_error )
        throw myruntime_error( mystring( errbuf ));
    else if( mMPIIAmMaster())
            error( errbuf );

    return code;
}

// -------------------------------------------------------------------------
// BcastlenMPIMessage: broadcast a message with its length at first from the 
//      master process to the slave processes; the slaves block until 
//      receive the message;
//      on successful return code is MPI_SUCCESS
// -------------------------------------------------------------------------

int MessageDispatcher::BcastlenMPIMessage( char* msg, size_t length, size_t max_len, bool throw_on_error )
{
    const int   lprm = 80;
    int         prvlen;
    int         code = MPI_SUCCESS;
    char        prvbuf[KBYTE], *p;
    static char errbuf[KBYTE] = {0};

    p = prvbuf;
    *( int* )p = lprm; p += sizeof( int );
    *( int* )p = length; p += sizeof( int );

    //bcast length first
    code = MPI_Bcast( prvbuf, int( p - prvbuf ), MPI_CHAR, mMPIMasterRank(), MPI_COMM_WORLD );

    if( code == MPI_SUCCESS ) {
        p = prvbuf;
        if( *( int* )p != lprm )
            sprintf( errbuf, "MessageDispatcher(%d): BcastlenMPIMessage: Invalid header.", mMPIMyRank());
        else {
            p += sizeof( int );
            prvlen = *( int* )p;

            if( max_len < prvlen )
                sprintf( errbuf, "MessageDispatcher(%d): BcastlenMPIMessage: "
                                 "Length exceeds maximum value.", mMPIMyRank());
            else {
                code = MPI_Bcast( msg, prvlen, MPI_CHAR, mMPIMasterRank(), MPI_COMM_WORLD );

                if( code == MPI_SUCCESS )
                    return MOK;
            }
        }
    }

    if( code != MPI_SUCCESS )
        sprintf( errbuf, "MessageDispatcher(%d): BcastlenMPIMessage: "
                         "MPI broadcast failed - %d.", mMPIMyRank(), code );
    if( throw_on_error )
        throw myruntime_error( mystring( errbuf ));
    else if( mMPIIAmMaster())
            error( errbuf );

    return code;
}

// -------------------------------------------------------------------------
// SendMPIMessage: sends message from a slave node to the master via MPI;
//     returns code obtained from the result of the operation; on success
//     it is MPI_SUCCESS
// -------------------------------------------------------------------------

int MessageDispatcher::SendMPIMessage( char* msg, size_t length, bool throw_on_error )
{
    int         code = MPI_SUCCESS;
    static char errbuf[KBYTE];

    if( mMPIIAmMaster()) {
        //master cannot send a message to the slave nodes this way
        code = MPI_ERR_RANK;
    } else {
        code = MPI_Send( msg, length, MPI_CHAR, mMPIMasterRank(), mMPIDataTag(), MPI_COMM_WORLD );
    }

    if( code == MPI_SUCCESS )
        return MOK;

    sprintf( errbuf, "MessageDispatcher(%d): Erroneous MPI send - %d.", mMPIMyRank(), code );

    if( throw_on_error )
        throw myruntime_error( mystring( errbuf ));
    else if( mMPIIAmMaster())
            error( errbuf );

    return code;
}

// -------------------------------------------------------------------------
// ReceiveMPIMessage: tries to receive a message from the slave node via MPI
//     returns code obtained from the result of the operation; on success
//     it is MPI_SUCCESS
//
// NOTE: this method is blocking itself until a message or error is received
// -------------------------------------------------------------------------

int MessageDispatcher::ReceiveMPIMessage( char* msg, size_t max_len, bool throw_on_error )
{
    MPI_Status  status;
    int         code = MPI_SUCCESS;
    int         recv = 0;   //number of bytes received
    int         cancelled = 0;

    static const char*  cnlstr = "Cancelled";
    static const int    cnllen = strlen( cnlstr );
    static char         errbuf[KBYTE];
    static char         buffer[MSG_MAX];
    int                 szbuffer = MSG_MAX;

    *buffer = 0;

    if( !GetFMinMsgSize() || !GetFIsMsgValid()) {
        code = ERR_PARAMS;
        if( throw_on_error )
            throw myruntime_error( mystring( "MessageDispatcher: Null function pointers." ));
        else
            error( "MessageDispatcher: Null function pointers." );
        return code;
    }

    if( ! mMPIIAmMaster()) {
        //slave node cannot wait to receive a message this way
        code = MPI_ERR_RANK;
    } else {
        code = MPI_Recv( msg, max_len, MPI_CHAR, MPI_ANY_SOURCE, mMPIDataTag(), MPI_COMM_WORLD, &status );
        if( code == MPI_SUCCESS ) {
            code = MPI_Test_cancelled( &status, &cancelled );
            if( cancelled ) {
                code = ERR_CANCELLED;
                memcpy( buffer, cnlstr, cnllen );
            }
            else if( code == MPI_SUCCESS ) {
                //determine the actual number of bytes received
                code = MPI_Get_count( &status, MPI_CHAR, &recv );

                if( code == MPI_SUCCESS ) {
                    if( recv <  ( *GetFMinMsgSize())() ||
                             !  ( *GetFIsMsgValid())( msg )) {
                        code = ERR_MESSAGE;
                        //assume we have an error message from a slave process
                        if( 0 < recv ) {
                            szbuffer = ( recv < MSG_MAX )? recv: ( MSG_MAX - 1 );
                            memcpy( buffer, msg, szbuffer );
                            buffer[szbuffer] = 0;
                            code = ERR_SLAVEMSG; //mark this message as one obtained from a slave process
                        }
                    }
                }
            }
        }
    }

    if( code == MPI_SUCCESS )
        return MOK;

    if( code == ERR_SLAVEMSG )
        sprintf( errbuf, "Transmitted" );
    else
        sprintf( errbuf, "MessageDispatcher(%d): Erroneous MPI receive - %d", mMPIMyRank(), code );

    if( *buffer )   sprintf( errbuf + strlen( errbuf ), ": %s", buffer );
        else        sprintf( errbuf + strlen( errbuf ), "." );

    if( throw_on_error )
        throw myruntime_error( mystring( errbuf ));
    else if( mMPIIAmMaster())
            error( errbuf );

    return code;
}

// -------------------------------------------------------------------------
// BlockUntilAllReady: the method blocks until all processes in the MPI ring
//     are ready for the ongoing job
//
// NOTE: this method is blocking itself until all processes are ready
// -------------------------------------------------------------------------

int MessageDispatcher::BlockUntilAllReady( bool throw_on_error )
{
    int         code = MPI_SUCCESS;
    static char errbuf[KBYTE];

    code = MPI_Barrier( MPI_COMM_WORLD );

    if( code == MPI_SUCCESS )
        return MOK;

    sprintf( errbuf, "MessageDispatcher(%d): MPI barrier with errors - %d.", mMPIMyRank(), code );

    if( throw_on_error )
        throw myruntime_error( mystring( errbuf ));
    else if( mMPIIAmMaster())
            error( errbuf );

    return code;
}
