/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "lmpi/msgcodes.h"

#include "rc.h"
#include "data.h"
#include "libpro/srcpro/FrequencyStore.h"
#include "ParallelUniversalScoreMatrix.h"




// =========================================================================
// CLASS ParallelUniversalScoreMatrix
//
// constructor
//
ParallelUniversalScoreMatrix::ParallelUniversalScoreMatrix(
        const FrequencyStore*   store,
        Configuration           config[NoSchemes],
        bool                    no_scaling,
        bool                    using_spec_lambda,
        double                  infrm_threshold,
        TFVectorProbabilities   distrib,
        bool                    master_flag,
        int                     rank,
        int                     ring_size,
        TBcastFunction          func_bcast,
        TSendFunction           func_send,
        TRecvFunction           func_recv,
        TBlockFunction          func_block,
        TScaling                a_scaling,
        TMask                   c_masking
    )
:
    AbstractUniversalScoreMatrix(
            store,
            ParallelUniversal,
            config,
            ComputeStatistics,
            a_scaling,
            c_masking,
            distrib,
            true/*cpu*/),

    master( master_flag ),
    mpi_rank( rank ),
    mpi_size( ring_size ),

    bcast_function( func_bcast ),
    send_function( func_send ),
    recv_function( func_recv ),
    block_function( func_block ),

    srm_buffer( NULL ),
    not_to_scale( no_scaling ),
    use_spec_lambda( using_spec_lambda )
{
    srm_buffer = ( char* )malloc( GetMaxSizeOfProbCalculator());
    if(  srm_buffer == NULL )
            USM_THROW( "ParallelUniversalScoreMatrix: Not enough memory." );

    Init( 0, 0 );   //zeros to indicate that no large matrix in memory will be used
    SetInformationThreshold( infrm_threshold );
}

// constructor overloaded
//
ParallelUniversalScoreMatrix::ParallelUniversalScoreMatrix(
        const FrequencyStore*   store,
        Configuration           config[NoSchemes],
        bool                    no_scaling,
        bool                    using_spec_lambda,
        double                  infrm_threshold,
        TFVectorProbabilities   distrib,
        bool                    master_flag,
        int                     rank,
        int                     ring_size,
        TBcastFunction          func_bcast,
        TSendFunction           func_send,
        TRecvFunction           func_recv,
        TBlockFunction          func_block,
        TScaling                a_scaling,
        TMask                   c_masking,
        TType                   type
    )
:
    AbstractUniversalScoreMatrix(
            store,
            type,
            config,
            ComputeStatistics,
            a_scaling,
            c_masking,
            distrib,
            true/*cpu*/),

    master( master_flag ),
    mpi_rank( rank ),
    mpi_size( ring_size ),

    bcast_function( func_bcast ),
    send_function( func_send ),
    recv_function( func_recv ),
    block_function( func_block ),

    srm_buffer( NULL ),
    not_to_scale( no_scaling ),
    use_spec_lambda( using_spec_lambda )
{
    srm_buffer = ( char* )malloc( GetMaxSizeOfProbCalculator());
    if(  srm_buffer == NULL )
            USM_THROW( "ParallelUniversalScoreMatrix: Not enough memory." );

    Init( 0, 0 );   //zeros to indicate that no large matrix in memory will be used
    SetInformationThreshold( infrm_threshold );
}

// -------------------------------------------------------------------------
// default constructor is invalid for initialization
//
ParallelUniversalScoreMatrix::ParallelUniversalScoreMatrix()
:
    AbstractUniversalScoreMatrix(),
    master( false ),
    mpi_rank( -1 ),
    mpi_size( -1 ),

    bcast_function( NULL ),
    send_function( NULL ),
    recv_function( NULL ),
    block_function( NULL ),

    srm_buffer( NULL ),
    not_to_scale( true ),
    use_spec_lambda( false )
{
}

// -------------------------------------------------------------------------
// destructor: deallocation of resources
//
ParallelUniversalScoreMatrix::~ParallelUniversalScoreMatrix()
{
    if( srm_buffer )
        free( srm_buffer );
}

// -------------------------------------------------------------------------
// TellSlavesToTerminate: broadcasts term. message to the slave processes
//
void ParallelUniversalScoreMatrix::TellSlavesToTerminate( bool throw_on_error )
{
    int code;

    if( ! IsMasterProcess())
        return;

    FormatTermMessage( GetSRMBuffer());

    code = ( *bcast_function )( GetSRMBuffer(), GetSizeOfScaleMessage(), throw_on_error );

    if( code != MOK )
        if( throw_on_error )
            USM_THROW( "ParallelUniversalScoreMatrix: Final MPI broadcast failed." );
}

// -------------------------------------------------------------------------
// ScaleScoringMatrix: master process scales scoring matrix using the
//     slave processes
//
void ParallelUniversalScoreMatrix::ScaleScoringMatrix()
{
    int     code;

    if( IsMasterProcess()) {
        //master process
        if( NotToPerformScaling()) {
//             if( GetAutoScaling())
//                 SetMultiplier( GetAutoScalingFactor());
            ComputeStatisticalParameters();
        } else {
            AbstractScoreMatrix::ScaleScoringMatrix();
        }
        TellSlavesToTerminate( true );
        SetMultiplier( 1.0 );

    } else {
        //slave process
        while( ComputeScoreProbabilitiesCPU( NULL ));
    }
}

// -------------------------------------------------------------------------
// ComputeScore: compute score given two vectors of frequencies
//
double ParallelUniversalScoreMatrix::ComputeScore(
    const FrequencyVector& vect_1st,
    const FrequencyVector& vect_2nd, size_t, size_t, bool )
{
    double  score = 0.0;

    size_t  fst_thckn = vect_1st.GetThickness();
    size_t  sec_thckn = vect_2nd.GetThickness();

    double  log_fst_thckn = 0.0;
    double  log_sec_thckn = 0.0;
    double  fst_weight = 1.0;
    double  sec_weight = 1.0;

#if 0
    if( fst_thckn ) log_fst_thckn = 1.0 - ( 1.0 + log( fst_thckn )) / ( double ) fst_thckn; //log( fst_thckn );
    if( sec_thckn ) log_sec_thckn = 1.0 - ( 1.0 + log( sec_thckn )) / ( double ) sec_thckn; //log( sec_thckn );

    if( log_fst_thckn <= 0.0 && log_sec_thckn <= 0.0 ) {
        log_fst_thckn = 1.0;
        log_sec_thckn = 1.0;
    }

    double  sum_thckn = log_fst_thckn + log_sec_thckn;

    if( sec_thckn == 0.0 )
        sum_thckn = 1.0;

    fst_weight = ( log_fst_thckn + log_fst_thckn ) / sum_thckn;
    sec_weight = ( log_sec_thckn + log_sec_thckn ) / sum_thckn;
#endif

    for( int a = 0; a < NUMAA; a++ ) {
        //computed as frequencies of query multiplied by scores of subject plus
        //  frequencies of subject multiplied by scores of query
        //lambda is implicitly incorporated into the log-odds values of the profiles
        score +=  fst_weight * ( double )vect_2nd.GetScoreAt( a ) / SCALE_CONSTANT * vect_1st.GetValueAt( a ) / FREQUENCY_SUM +
                  sec_weight * ( double )vect_2nd.GetValueAt( a ) / FREQUENCY_SUM  * vect_1st.GetScoreAt( a ) / SCALE_CONSTANT;
    }

    score *= GetMultiplier();

    return score;
}

// -------------------------------------------------------------------------
// VectorScoreProbability: compute score probability given two vectors of 
//  frequencies
//
double ParallelUniversalScoreMatrix::VectorScoreProbability(
    const FrequencyVector& rowvect,
    const FrequencyVector& colvect, size_t, size_t, bool ) const
{
    switch( GetDistributionType()) {
        case DISCRETE:
            //return column vector observed frequency
            return colvect.GetProbability();

        case PROVECTOR:
            break;

        case MULTINOMIAL:
            //return column vector multinomial probability
            return colvect.GetProbability();

        default:
            throw myruntime_error( mystring( "ParallelUniversalScoreMatrix: Unknown vector distribution type." ));
    }

    //compute probability given profile vector distribution type;
    //this is a product of two vector probabilities
    return rowvect.GetProbability() * colvect.GetProbability();
}

// -------------------------------------------------------------------------
// MultiplierProgress: method which is invoked just before saving of
//     multiplier; print if master
//
void ParallelUniversalScoreMatrix::MultiplierProgress( AttributableScores* PATTR_SCORES )
{
    static char buffer[MSG_MAX];
    static int  round = 0;

    if( PATTR_SCORES == NULL )
        throw myruntime_error( mystring( "ParallelUniversalScoreMatrix: MultiplierProgress: Wrong argument." ));

    AbstractScoreMatrix::MultiplierProgress( PATTR_SCORES );

    double  loc_lambda = GetLambda();

    if( GetAutoScaling())
        loc_lambda *= PATTR_SCORES->GetAutoScalingFactor();

    if( IsMasterProcess()) {
        sprintf( buffer, "Pass %2d: lambda = %.4f, scale factor = %.4f", round++, loc_lambda, GetMultiplier());
        message( buffer );
    }
}

// -------------------------------------------------------------------------
// ComputeScoreProbabilitiesCPU: computes probabilities to observe scores at
//     each position (i,j). Probability is computed as
//
//      1     _
//    ----   \   Pfj
//   length  /_ 
//           i,j:
//         sij=sk
//
//  where Pfj are probabilities of frequency vector fj to occur; sij are
//  scores at (i,j), sk is a discrete score value; length is for query.
//  Method iterates through the large table of scores only once saving
//  probability distribution on the fly;
//  returns FALSE only if the process should terminate
//
//  NOTE: This method allocates required space for probabilities!
//
bool ParallelUniversalScoreMatrix::ComputeScoreProbabilitiesCPU( AttributableScores* PATTR_SCORES )
{
    int         code;

    if( IsMasterProcess()) {
        if( PATTR_SCORES )
            SetMultiplier( PATTR_SCORES->GetAutoScalingFactor());
        FormatScaleFactor( GetSRMBuffer());
    }

    //commmon code for the master and slave processes
    code = ( *bcast_function )( GetSRMBuffer(), GetSizeOfScaleMessage(), true/*throw on error*/);

    if( code != MOK )
        USM_THROW( "ParallelUniversalScoreMatrix: MPI broadcast failed." );

    //obtain scale factor for the slave processes
    if( ! IsMasterProcess()) {
        if( IsThatTermMessage( GetSRMBuffer()))
            //the process must terminate on the master's order
            return false;
        DeformatScaleFactor( GetSRMBuffer());
    }


    if( IsMasterProcess()) {
        MasterComputations( PATTR_SCORES );
    } else {
        SlaveComputations();
    }


    //WAIT here until other processes pass this point!
    //Actually broadcast for the slave processes will block them, so
    //there's no need in blocking this way
//     ( *block_function )( true/*throw on error*/);

    return true;
}

// -------------------------------------------------------------------------
// MasterComputations: this code is for the master only; it tells slaves
//     when computations are to be invoked, receives all partially computed
//     probabilities, and performs final calculations
//
void ParallelUniversalScoreMatrix::MasterComputations( AttributableScores* PATTR_SCORES )
{
    if( ! IsMasterProcess())
        return;

    bool    all_negat = true;       //all negative scores
    bool    loc_all_negat = true;   //local variable to check for negativity

    //get storage for scores and probabilities...
    BinarySearchStructure*  loc_probabilities = GetProbCalculator();

#ifdef __DEBUG__
    if( GetMyMPIRank() < 0 || GetMPIRingSize() < 1 )
        USM_THROW( "ParallelUniversalScoreMatrix: Wrong MPI rank or ring size." );

//     if( !GetStore() || !GetStore()->GetFrequencies() ||
//         !loc_probabilities )
//         USM_THROW( "ParallelUniversalScoreMatrix: Unable to compute probabilities." );
#endif

    if( loc_probabilities->GetSize())
        USM_THROW( "ParallelUniversalScoreMatrix: Error in computation of probabilities." );


//     const SimpleVector*     frequencies = GetStore()->GetFrequencies();
//     size_t                  no_freqs = 0;
    double                  probsum = 0.0;      //sum of probabilities of all scores received
    double                  partialsum = 0.0;   //partial sum of probabilities of scores received
    mystring                throwstr;

//     if( frequencies )
//         no_freqs = frequencies->GetSize();

    bool        mpierror = false;
    int         code;

    //begin iteration with 1 since one of the rank ids is of the master
    for( size_t r = 1; r < GetMPIRingSize(); r++ ) {
        //we must receive messages from all but the master node in the ring
        code = ( *recv_function )( GetSRMBuffer(), GetMaxSizeOfProbCalculator(), false/*not to throw on error*/);

        if( code != MOK || mpierror ) {
            //we have error but we must go on waiting for the rest of messages
            mpierror = true;
            continue;
        }

        BinarySearchStructure*  probs_obtained = NULL;

        try {
            probs_obtained = DeformatProbCalculator( GetSRMBuffer());
            loc_all_negat = InfuseProbCalculator( probs_obtained, &partialsum );

            probsum += partialsum;

            if( all_negat && ! loc_all_negat )
                all_negat = false;

        } catch( myexception const& ex ) {
            //error! deallocate resources, set error flag, and continue receiving of messages
            mpierror = true;
            throwstr = ex.what();
        }

        if( probs_obtained )
            //members of probs_obtained are redirected to ProbCalculator of this class,
            //destroy only vector used to contain them
            delete probs_obtained;
    }

    if( mpierror ) {
        if( throwstr.empty())
            USM_THROW( "ParallelUniversalScoreMatrix: Error while receiving data from slave processes." );
        else
            USM_THROW( throwstr.c_str());
    }

    //OK, perform further computations on the master's own
    SetAllNegatives( all_negat );
    ProcessScoreProbabilities( PATTR_SCORES, probsum );
}

// -------------------------------------------------------------------------
// SlaveComputations: this is code for the slave processes; it waits until
//     master tells to restart computations, performs partial computations
//
void ParallelUniversalScoreMatrix::SlaveComputations()
{
    if( IsMasterProcess())
        return;

    bool    all_negat = true;   //all negative scores

    int     loc_multiplier = 1;

    if( GetAutoScaling() && !NotToPerformScaling())
        loc_multiplier = GetAutoScalingFactor();


    //get storage for scores and probabilities...
    BinarySearchStructure*  loc_probabilities = GetProbCalculator();

#ifdef __DEBUG__
    if( GetMyMPIRank() < 0 || GetMPIRingSize() < 1 )
        USM_THROW( "ParallelUniversalScoreMatrix: Wrong MPI rank or ring size." );

    if( !GetStore() || !GetStore()->GetFrequencies() ||
        !loc_probabilities )
        USM_THROW( "ParallelUniversalScoreMatrix: Unable to compute probabilities." );
#endif

    if( loc_probabilities->GetSize())
        USM_THROW( "ParallelUniversalScoreMatrix: Error in computation of probabilities." );


    const SimpleVector*     frequencies = GetStore()->GetFrequencies();
    size_t                  no_freqs = 0;

    bool                    error = false;
    double                  score = 0.0;
    double                  scoreprob = 0.0;
    ScoreProbability*       scp = ScoreProbability::NewScoreProbability();

    if( frequencies )
        no_freqs = frequencies->GetSize();

    //one MPI process is aassigned to master and does not take part in matrix decomposition
    //
//     for( size_t N = 0 + GetMyMPIRank() - 1; N < no_freqs && !error; N += GetMPIRingSize() - 1 )
    for( size_t N = 0; N < no_freqs && !error; N ++ )
    {
        const FrequencyVector   rowvector(( const char* )frequencies->GetValueAt( N ));

        //omit rows for which positional vectors have less information than the threshold
        if(( double ) rowvector.GetInfContent() / INFO_SCALE_CONSTANT < GetInformationThreshold())
            continue;

        //The large matrix is SYMMETRIC with respect to its diagonal, so
        //iterate over the upper triangle only...
//         for( size_t M = N; M < no_freqs && !error; M ++ )
        for( size_t M = N + GetMyMPIRank() - 1; M < no_freqs && !error; M += GetMPIRingSize() - 1 )
        {
            const FrequencyVector   colvector(( const char* )frequencies->GetValueAt( M ));

            //omit columns for which positional vectors have less information than the threshold
            if(( double ) colvector.GetInfContent() / INFO_SCALE_CONSTANT < GetInformationThreshold())
                continue;

            score = ComputeScore( rowvector, colvector, N, M );
            scoreprob = VectorScoreProbability( rowvector, colvector, N, M );

            scp->SetScores( score * loc_multiplier );
            scp->SetProbability( scoreprob );

            if( M != N ) {
                //change order of the vectors...
                scoreprob = VectorScoreProbability( colvector, rowvector, M, N, true );
                //if this isn't a score obtained via the diagonal elements,
                //increase the score probability since it occurs twice: in the upper and lower triangles
                scp->IncProbability( scoreprob );
            }

            if( all_negat && 0 < scp->GetScore())
                all_negat = false;

            if( scp->GetScore() <= SCORE_MIN || scp->GetProbability() <= 0.0 )
                continue;

            if( Push( scp )) {
                //score with probability has been inserted: construct new object
                scp = ScoreProbability::NewScoreProbability();
            }
        }
    }
// PrintProbCalculator( stderr );
    ScoreProbability::Destroy( scp );

    //format data and send them to the master
    size_t  written = FormatProbCalculator( GetSRMBuffer());

    //the function below will throw an exception if an error occurs
    ( *send_function  )( GetSRMBuffer(), written, true );

    //destroy members of loc_probabilities and re-initialize it
    ClearProbCalculator( loc_probabilities );
}

// -------------------------------------------------------------------------
// ComputeScoreProbabilitiesMem: for this class, this method cannot be
//     called
//
bool ParallelUniversalScoreMatrix::ComputeScoreProbabilitiesMem( AttributableScores* )
{
    USM_THROW( "ParallelUniversalScoreMatrix: Illegal method call." );
    return false;
}

// =========================================================================
// PRINT ROUTINES
//
// PrintParameterTable: nothing to do since output of the parameter
//     table is in the final
//
void ParallelUniversalScoreMatrix::PrintParameterTable( TPrintFunction, void* ) const
{
}

// PrintFinal: print table of computed statistical parameter
//     values to a stream
//
void ParallelUniversalScoreMatrix::PrintFinal( TPrintFunction print_func, void* vpn ) const
{
    if( vpn == NULL )
        return;

    if( 0.0 <= GetExpectedScore()) {
        print_func( vpn, "Expected score per position is non-negative, %.4f!\n\n", GetExpectedScore());
        return;
    }

    char    _sk[BUF_MAX], _slambda[BUF_MAX];
    if( 0.0 < GetK())   sprintf( _sk, "%6.4f", GetK());
        else            sprintf( _sk, "n/a" );
    if( 0.0 < GetLambda())  sprintf( _slambda, "%6.4f", GetLambda());
        else                sprintf( _slambda, "n/a" );

    print_func( vpn, "%-25s  %6.4f\n\n", "Scale factor,", GetMultiplier());

    print_func( vpn, "%-25s  %-6s   %-6s\n", " ", "K", "Lambda" );
    if( UsingSpecificLambda())
        print_func( vpn, "%-25s  %6c   %6.4f\n", "Reference user-specified,", 32, GetRefLambda());
    else
        print_func( vpn, "%-25s  %6.4f   %6.4f\n", "Reference ungapped,", GetRefK(), GetRefLambda());
    print_func( vpn, "%-25s  %6s   %6s\n", "Computed  ungapped,", _sk, _slambda );

    print_func( vpn, "Entropy, %6.4f; Expected, %6.4f\n\n", GetEntropy(), GetExpectedScore());

    print_func( vpn, "%-25s  %3d/%-3d\n\n", "min/max scores,", GetMinScore(), GetMaxScore());
}

// -------------------------------------------------------------------------
// PrintScoringMatrix: output the computed scoring system
//
void ParallelUniversalScoreMatrix::PrintScoringMatrix( FILE* fp )
{
    if( !GetStore() || !GetStore()->GetFrequencies())
        return;

    if( fp == NULL )
        return;

    fprintf( fp,"%12c Position-specific profile scoring matrix\n", 32 );

    fprintf( fp, "%9c", 32 );

    const SimpleVector*     frequencies = GetStore()->GetFrequencies();
    size_t                  no_freqs = frequencies->GetSize();
    int                     l = 0;

    for( size_t m = 0; m < no_freqs; m++ )
        fprintf( fp, "%4d", m );


    for( size_t N = 0; N < no_freqs; N++ ) {
        const FrequencyVector   rowvector(( const char* )frequencies->GetValueAt( N ));
        fprintf( fp, "\n%5d %c   ", ++l, 32 );

        for( size_t M = 0; M < no_freqs; M++ )
        {
            const FrequencyVector   colvector(( const char* )frequencies->GetValueAt( M ));
            fprintf( fp, "%3d ", ( int )rint( ComputeScore( rowvector, colvector, N, M )));
        }
    }

    fprintf( fp, "\n" );
}

