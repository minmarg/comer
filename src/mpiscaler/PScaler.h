/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __PScaler__
#define __PScaler__

#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "rc.h"
#include "data.h"
#include "mystring.h"
#include "myexcept.h"

#include "lmpi/MessageDispatcher.h"
#include "libpro/srcpro/Database.h"
#include "libpro/srcsco/AbstractUniversalScoreMatrix.h"
#include "ParallelUniversalScoreMatrix.h"
#include "PrlLSOUniversalScoreMatrix.h"
#include "PrlHDPUniversalScoreMatrix.h"

// _________________________________________________________________________
// Class PScaler
//

class PScaler: public MessageDispatcher
{
public:
    PScaler(
            const char* paramconfigfile,
            const char* database,
            const char* output,
            bool no_scaling,
            double infrm_threshold,
            double lambda,
            AbstractScoreMatrix::TScaling a_scaling,
            TMask c_masking,
            AbstractScoreMatrix::TType type
    );
    ~PScaler();

    virtual long    Run();                                      //the beginning of scaling

    const char*     GetParamConfigFile() const  { return paramConfigFile; }
    const char*     GetDatabase() const         { return database_name; }
    const char*     GetOutput() const           { return output_name; }

    bool            NotToPerformScaling() const     { return not_to_scale;  }
    double          GetInformationThreshold() const { return infcontent; }
    double          GetLambdaSpecific() const       { return lambda_spec;   }
    bool            UsingLambdaSpecific() const     { return 0.0 < lambda_spec; }

    void            PrintMethodName( FILE* fp ) const;          //printing of the method name used in scoring alignments
    void            PrintParameterTable( FILE* ) const;         //printing of parameter table

protected:
    explicit PScaler();
                                                    //configuration access routines
    const Configuration&        GetConfiguration( TProcomSchemes ) const;
    Configuration&              GetConfiguration( TProcomSchemes );
    Configuration* const        GetConfiguration()              { return configuration; }

    void                        CreateScoreSystem();            //create member score system
    void                        DestroyScoreSystem();           //destroy score system
    void                        ScaleScoreSystem();             //scale score system
    void                        Terminate();                    //terminate the slave processes

    const AbstractUniversalScoreMatrix* GetScoreSystem() const  { return scoreSystem; }

    void                        ComputationLogicWithProfiles( FrequencyMatrix&, LogOddsMatrix&, GapScheme& );
    void                        PostComputationLogic();

    size_t          GetNoSequences() const      { return no_sequences; }    //number of sequences in database
    Uint64          GetDbSize() const           { return db_size; }         //size of database

    void                        SetNoSequences( size_t value )  { no_sequences = value; }
    void                        SetDbSize( Uint64 value )       { db_size = value; }

    AbstractScoreMatrix::TScaling   GetScaling() const          { return scaling; }
    TMask                           GetMasking() const          { return masking; }
    AbstractScoreMatrix::TType      GetType() const             { return type_; }

    void                        PrintHeader( FILE* );           //print of header
    void                        PrintFooter( FILE* );           //print of footer

private:
    const char*                     paramConfigFile;//full pathname to parameter configuration file
    const char*                     database_name;  //profile database name
    const char*                     output_name;    //output file name, null if standard output
    AbstractUniversalScoreMatrix*   scoreSystem;    //score system used to be scaled
    Configuration                   configuration[NoSchemes];   //parameter configuration

    bool                            not_to_scale;   //whether or not to perform scaling of matrix
    const double                    infcontent;     //information content threshold
    double                          lambda_spec;    //specific reference lambda value

    Database                        profile_db;     //profile database
    size_t                          no_sequences;   //number of sequences within the database
    Uint64                          db_size;        //size of the database

    AbstractScoreMatrix::TScaling   scaling;        //strategy of precision
    TMask                           masking;        //masking approach
    AbstractScoreMatrix::TType      type_;          //type of scoring
};


////////////////////////////////////////////////////////////////////////////
// INLINES
//
// CLASS PScaler
//
// GetConfiguration: returns reference to parameter configuration object
//
inline
const Configuration& PScaler::GetConfiguration( TProcomSchemes ps ) const
{
    if( NoSchemes <= ps )
        throw myruntime_error( mystring( "PScaler: Wrong argument while accessing configuration." ));

    return configuration[ps];
}

// GetConfiguration: overloaded
//
inline
Configuration& PScaler::GetConfiguration( TProcomSchemes ps )
{
    if( NoSchemes <= ps )
        throw myruntime_error( mystring( "PScaler: Wrong argument while accessing configuration." ));

    return configuration[ps];
}

// -------------------------------------------------------------------------
// CreateScoreSystem: creates score system used to be scaled
//
inline
void PScaler::CreateScoreSystem()
{
    switch( GetType()) {
      case AbstractScoreMatrix::ParallelUniversal:
          scoreSystem = new ParallelUniversalScoreMatrix(
                  profile_db.GetStore(),
                  GetConfiguration(),
                  NotToPerformScaling(),
                  UsingLambdaSpecific(),
                  GetInformationThreshold(),
                  Database::GetDistributionType(),
                  mMPIIAmMaster(),
                  mMPIMyRank(),
                  mMPIRingSize(),
                  &BcastMPIMessage,
                  &SendMPIMessage,
                  &ReceiveMPIMessage,
                  &BlockUntilAllReady,
                  GetScaling(),
                  GetMasking()
          );
          break;
      case AbstractScoreMatrix::ParallelLSOUniversal:
          scoreSystem = new PrlLSOUniversalScoreMatrix(
                  profile_db.GetStore(),
                  GetConfiguration(),
                  NotToPerformScaling(),
                  UsingLambdaSpecific(),
                  GetInformationThreshold(),
                  Database::GetDistributionType(),
                  mMPIIAmMaster(),
                  mMPIMyRank(),
                  mMPIRingSize(),
                  &BcastMPIMessage,
                  &SendMPIMessage,
                  &ReceiveMPIMessage,
                  &BlockUntilAllReady,
                  GetScaling(),
                  GetMasking()
          );
          break;
      case AbstractScoreMatrix::ParallelHDPUniversal:
          scoreSystem = new PrlHDPUniversalScoreMatrix(
                  profile_db.GetStore(),
                  GetConfiguration(),
                  NotToPerformScaling(),
                  UsingLambdaSpecific(),
                  GetInformationThreshold(),
                  Database::GetDistributionType(),
                  mMPIIAmMaster(),
                  mMPIMyRank(),
                  mMPIRingSize(),
                  &BcastMPIMessage,
                  &SendMPIMessage,
                  &ReceiveMPIMessage,
                  &BlockUntilAllReady,
                  GetScaling(),
                  GetMasking()
          );
          break;
      default:
          throw myruntime_error("PScaler: Invalid type of scoring method.");
    }

    if( scoreSystem == NULL )
        throw myruntime_error("PScaler: Not enough memory.");
}

// -------------------------------------------------------------------------
// ScaleScoreSystem: scale score system
//
inline
void PScaler::ScaleScoreSystem()
{
#ifdef __DEBUG__
    if( !scoreSystem )
        throw myruntime_error( mystring( "PScaler: Unable to scale score matrix." ));
#endif

    if( scoreSystem->GetType() == AbstractScoreMatrix::ParallelUniversal ||
        scoreSystem->GetType() == AbstractScoreMatrix::ParallelLSOUniversal ||
        scoreSystem->GetType() == AbstractScoreMatrix::ParallelHDPUniversal )
        scoreSystem->ScaleScoringMatrix();
}

// -------------------------------------------------------------------------
// DestroyScoreSystem: destroys score system
//
inline
void PScaler::DestroyScoreSystem()
{
    if( scoreSystem )
        delete scoreSystem;
}

// -------------------------------------------------------------------------
// Terminate: terminate the slave processes in the MPI ring
//
inline
void PScaler::Terminate()
{
    if( ! scoreSystem )
        return;

    if( scoreSystem->GetType() != AbstractScoreMatrix::ParallelUniversal &&
        scoreSystem->GetType() != AbstractScoreMatrix::ParallelLSOUniversal &&
        scoreSystem->GetType() != AbstractScoreMatrix::ParallelHDPUniversal )
        return;

    dynamic_cast<ParallelUniversalScoreMatrix*>( scoreSystem )->TellSlavesToTerminate( false/*throw_on_error*/);
}

// -------------------------------------------------------------------------
// PrintMethodName: prints name of method used to score profile positions
//
inline
void PScaler::PrintMethodName( FILE* fp ) const
{
    const size_t    padding = OUTPUTINDENT;
    char            blanks[padding+1];
    size_t          n = 0;

    for( n = 0; n < padding; blanks[n++] = 32 );
    blanks[n] = 0;

    if( scoreSystem ) {
        fprintf( fp, "Scoring method: %s\n", scoreSystem->GetMethodName());
        if( scoreSystem->GetAutoScaling())
            fprintf( fp, "%s(Increased precision)\n", blanks );
        else if( scoreSystem->GetFPScaling())
                fprintf( fp, "%s(Floating-point precision)\n", blanks );

        fprintf( fp, "Information threshold: %.2f\n", GetInformationThreshold());
        fprintf( fp, "\n");
        return;
    }
}

// -------------------------------------------------------------------------
// PrintParameterTable: prints parameter table
//
inline
void PScaler::PrintParameterTable( FILE* fp ) const
{
    if( !scoreSystem )
        return;

    scoreSystem->PrintFinal( fp );
}

#endif//__PScaler__
