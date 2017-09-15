/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __MessageDispatcher__
#define __MessageDispatcher__

#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "rc.h"
#include "data.h"

#include "mystring.h"
#include "myexcept.h"

#include "msgcodes.h"

// _________________________________________________________________________
// Class MessageDispatcher
//

class MessageDispatcher
{
public:
    //typedefs...
    typedef int ( *TBcastFunction )( char* msg, size_t, bool throw_on_error );
    typedef int ( *TSendFunction )( char* msg, size_t, bool throw_on_error );
    typedef int ( *TRecvFunction )( char* msg, size_t, bool throw_on_error );
    typedef int ( *TBlockFunction )( bool throw_on_error );

public:
    typedef size_t ( *TMinMsgSizeDet )();
    typedef bool   ( *TMsgValidator )( const char* );

public:
//     MessageDispatcher( TMinMsgSizeDet, TMsgValidator );
    virtual ~MessageDispatcher();

    virtual long    Run() = 0;

protected:
    explicit MessageDispatcher();
                                                                //put a message if running as master
    void                        MasterMessage( const char*, bool = true );

    //{{MPI logic methods
    static void                 PrivateMPIInit();               //initialization of MPI environment
    static void                 PrivateMPITerminate();          //correct leave of MPI environment

    static int                  mMPIDataTag()                   { return 99;            }
    static int                  mMPIMasterRank()                { return 0;             }
    static bool                 mMPIIAmMaster()                 { return mpi_rank_ == mMPIMasterRank(); }
    static int                  mMPIMyRank()                    { return mpi_rank_;      }
    static int                  mMPIRingSize()                  { return mpi_size_;      }

    static TMinMsgSizeDet       GetFMinMsgSize()                { return fGetMinMsgSize_; }
    static TMsgValidator        GetFIsMsgValid()                { return fIsMsgValid_; }

    static void                 SetFMinMsgSize( TMinMsgSizeDet addr )   { fGetMinMsgSize_ = addr; }
    static void                 SetFIsMsgValid( TMsgValidator addr )    { fIsMsgValid_ = addr; }

    static int                  BcastMPIMessage( char*, size_t, bool throw_on_error = true );
    static int                  BcastlenMPIMessage( char*, size_t, size_t, bool throw_on_error = true );
    static int                  SendMPIMessage( char*, size_t, bool throw_on_error = true );
    static int                  ReceiveMPIMessage( char*, size_t, bool throw_on_error = true );
    static int                  BlockUntilAllReady( bool throw_on_error = true );
    //}}

private:
    static int                  mpi_rank_;      //rank of this process in the MPI ring
    static int                  mpi_size_;      //number of processes running in the MPI ring
    static TMinMsgSizeDet       fGetMinMsgSize_;//min message size determinant
    static TMsgValidator        fIsMsgValid_;   //message validator
};


////////////////////////////////////////////////////////////////////////////
// INLINES
//
// CLASS MessageDispatcher
//
// MasterMessage: prints message if running as master
//
inline
void MessageDispatcher::MasterMessage( const char* msg, bool putnl )
{
    if( mMPIIAmMaster())
        message( msg, putnl );
}

#endif//__MessageDispatcher__
