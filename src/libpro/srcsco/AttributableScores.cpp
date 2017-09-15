/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rc.h"
#include "data.h"
#include "libpro/srcpro/DistributionMatrix.h"
#include "AbstractScoreMatrix.h"
#include "AttributableScores.h"


// -------------------------------------------------------------------------
// constructor: initialize class members with initial values
// -------------------------------------------------------------------------

AttributableScores::AttributableScores(
        AbstractScoreMatrix*    prnt,
        TProbabilityFunction    probfun,
        TScoreReductFunction    redtfun,
        TrogressFunction        progressfun,
        double                  reflambda,
        double                  refH,
        double                  refK,
        const double**          im,
        const TMask**           msk,
        bool                    keep,
        int                     as_factor )
:
    parent( prnt ),
    parent_probfunction( probfun ),
    parent_redtfunction( redtfun ),
    parent_progressfunction( progressfun ),

    private_multiplier( 1.0 ),
    auto_scaling_factor( as_factor ),

    image( im ),
    mask( msk ),

    queryinfcontent( NULL ),
    sbjctinfcontent( NULL ),
    queryinfprobabs( NULL ),
    sbjctinfprobabs( NULL ),
    queryminmaxscores( NULL ),
    sbjctminmaxscores( NULL ),
    probabilities( NULL ),

    min_score( 0 ),
    max_score( 0 ),
    score_gcd( 1 ),

    priv_prob_vector( NULL ),
    prob_vector_size( 0 ),

    expscore( 0.0 ),
    allnegatives( false ),
    lambda( -1.0 ),
    entropy( 0.0 ),
    parameterK( -1.0 ),

    constant_c( -1.0 ),

    querySize( 0 ),
    subjectSize( 0 ),

    keepinmem( keep )
{
    SetRefLambda( reflambda );
    SetRefH( refH );
    SetRefK( refK );
    NewPrivateProbVector( PARAMETER_K_MAXIT, MAX_RANGE_OF_SCORES );
}

// -------------------------------------------------------------------------
// default construction
// -------------------------------------------------------------------------

AttributableScores::AttributableScores()
:
    parent( NULL ),
    parent_probfunction( NULL ),
    parent_redtfunction( NULL ),
    parent_progressfunction( NULL ),

    private_multiplier( 1.0 ),
    auto_scaling_factor( 1 ),

    image( NULL ),
    mask( NULL ),

    queryinfcontent( NULL ),
    sbjctinfcontent( NULL ),
    queryinfprobabs( NULL ),
    sbjctinfprobabs( NULL ),
    queryminmaxscores( NULL ),
    sbjctminmaxscores( NULL ),
    probabilities( NULL ),

    min_score( 0 ),
    max_score( 0 ),
    score_gcd( 1 ),

    priv_prob_vector( NULL ),
    prob_vector_size( 0 ),

    expscore( 0.0 ),
    allnegatives( false ),
    lambda( -1.0 ),
    entropy( 0.0 ),
    parameterK( -1.0 ),

    constant_c( -1.0 ),

    querySize( 0 ),
    subjectSize( 0 ),

    keepinmem( true )
{
    throw myruntime_error(
            mystring( "AttributableScores: Default initialization is not allowed." ));
}

// -------------------------------------------------------------------------
// destructor: deallocate memory used by this class
// -------------------------------------------------------------------------

AttributableScores::~AttributableScores()
{
    DestroyProbabilities();
    DestroyPrivateProbVector();
    DestroyQueryInfContent();
    DestroySbjctInfContent();
    DestroyQueryInfProbs();
    DestroySbjctInfProbs();
    DestroyQueryMinMaxScores();
    DestroySbjctMinMaxScores();
}

// -------------------------------------------------------------------------
// NewQueryInfProbsAt: allocate memory for probabilities for information
//     content at query position n
// -------------------------------------------------------------------------

void AttributableScores::NewQueryInfProbsAt( int size, int npos )
{
    if( size <= 0 ) {
        warning( "AttributableScores: NewQueryInfProbsAt: Wrong size." );
        return;
    }
#ifdef __DEBUG__
    if( !queryinfprobabs || GetQuerySize() <= npos || npos < 0 )
        throw myruntime_error(
            mystring( "AttributableScores: Memory access error." ));
#endif
    if( queryinfprobabs[npos] )
        free( queryinfprobabs[npos] );

    //NOTE: +1: reserve additional space for norm. term of probabilities
    queryinfprobabs[npos] = ( double* )malloc( sizeof( double ) * ( size + 1 ));

    if( !queryinfprobabs[npos] ) 
        throw myruntime_error( mystring( "AttributableScores: Not enough memory." ));

    for( int k = 0; k < size + 1; k++ )
        queryinfprobabs[npos][k] = 0.0;
}

// -------------------------------------------------------------------------
// NewQueryInfProbs: allocate memory for probabilities for information
//     content of query positions
// -------------------------------------------------------------------------

void AttributableScores::NewQueryInfProbs()
{
    int size = GetQuerySize();
    if( size <= 0 ) {
        warning( "AttributableScores: NewQueryInfProbs: Wrong size." );
        return;
    }
    DestroyQueryInfProbs();

    queryinfprobabs = ( double** )malloc( sizeof( double* ) * size );

    if( !queryinfprobabs ) 
        throw myruntime_error( mystring( "AttributableScores: Not enough memory." ));

    for( int n = 0; n < size; n++ )
        queryinfprobabs[n] = NULL;
}

// -------------------------------------------------------------------------
// DestroyQueryInfProbs: deallocate memory of probabilities for information
//     content of query positions
// -------------------------------------------------------------------------

void AttributableScores::DestroyQueryInfProbs()
{
    int size = GetQuerySize();
    if( queryinfprobabs ) {
        for( int n = 0; n < size; n++ )
            if( queryinfprobabs[n] )
                free( queryinfprobabs[n] );
        free( queryinfprobabs );
        queryinfprobabs = NULL;
    }
}

// -------------------------------------------------------------------------
// NewSbjctInfProbsAt: allocate memory for probabilities for information
//     content at subject position m
// -------------------------------------------------------------------------

void AttributableScores::NewSbjctInfProbsAt( int size, int mpos )
{
    if( size <= 0 ) {
        warning( "AttributableScores: NewSbjctInfProbsAt: Wrong size." );
        return;
    }
#ifdef __DEBUG__
    if( !sbjctinfprobabs || GetSubjectSize() <= mpos || mpos < 0 )
        throw myruntime_error(
            mystring( "AttributableScores: Memory access error." ));
#endif
    if( sbjctinfprobabs[mpos] )
        free( sbjctinfprobabs[mpos] );

    //NOTE: +1: reserve additional space for norm. term of probabilities
    sbjctinfprobabs[mpos] = ( double* )malloc( sizeof( double ) * ( size + 1 ));

    if( !sbjctinfprobabs[mpos] )
        throw myruntime_error( mystring( "AttributableScores: Not enough memory." ));

    for( int k = 0; k < size + 1; k++ )
        sbjctinfprobabs[mpos][k] = 0.0;
}

// -------------------------------------------------------------------------
// NewSbjctInfProbs: allocate memory for probabilities for information
//     content of subject positions
// -------------------------------------------------------------------------

void AttributableScores::NewSbjctInfProbs()
{
    int size = GetSubjectSize();
    if( size <= 0 ) {
        warning( "AttributableScores: NewQueryInfProbs: Wrong size." );
        return;
    }
    DestroySbjctInfProbs();

    sbjctinfprobabs = ( double** )malloc( sizeof( double* ) * size );

    if( !sbjctinfprobabs ) 
        throw myruntime_error( mystring( "AttributableScores: Not enough memory." ));

    for( int m = 0; m < size; m++ )
        sbjctinfprobabs[m] = NULL;
}

// -------------------------------------------------------------------------
// DestroySbjctInfProbs: deallocate memory of probabilities for information
//     content of subject positions
// -------------------------------------------------------------------------

void AttributableScores::DestroySbjctInfProbs()
{
    int size = GetSubjectSize();
    if( sbjctinfprobabs ) {
        for( int m = 0; m < size; m++ )
            if( sbjctinfprobabs[m] )
                free( sbjctinfprobabs[m] );
        free( sbjctinfprobabs );
        sbjctinfprobabs = NULL;
    }
}

// -------------------------------------------------------------------------
// Init: initializes class's private data given query and subject sizes
// IMPORTANT: query and subject sizes are supposed to be initialized
// -------------------------------------------------------------------------

void AttributableScores::Init()
{
    NewQueryInfContent( GetQuerySize());
    NewSbjctInfContent( GetSubjectSize());
    NewQueryInfProbs();
    NewSbjctInfProbs();
    NewQueryMinMaxScores();
    NewSbjctMinMaxScores();
}

// -------------------------------------------------------------------------
// SearchForHSPs: perform search of multiple high-scoring pairs (hsp)
//     in the same diagonal of the score system
// minhspscore,- minimum score of hsp
// hsplen,- length of hsp
// nohsps,- minimum number of hsps required
// maxdist,- max distance between hsps
// possbjct, posquery,- subject and query positions of hsps found in the
//     same diagonal
// returns true if such hsps have been found
// -------------------------------------------------------------------------

bool AttributableScores::SearchForHSPs(
    double minhspscore, int hsplen, int nohsps, int maxdist,
    int* possbjct, int* posquery )
{
    bool            bret = false;
    const TScore    hmin = ( TScore )minhspscore;
    TScore          hmax = 0;
    const int       diaglen = nohsps * 3; //nohsps*2 for indices
    TScore**        diagonals = NULL;   //diagonals keeping best-scoring hsps
    TScore**        hsps = NULL;
    int**           lens = NULL;
    TScore          loc_score, tvalue;
    mystring        errstr;
    int             m, n, ind, sum, k, l;
    int             msb, nqu, mind;

    if( hmin <= 0 ) {
        warning( "Parameter: Non-positive HSP score threshold." );
        return false;
    }
    if( hsplen <= 0 ) {
        warning( "Parameter: Non-positive HSP length." );
        return false;
    }
    if( nohsps <= 0 ) {
        warning( "Parameter: Non-positive number of HSPs." );
        return false;
    }
    if( maxdist <= 0 ) {
        warning( "Parameter: Non-positive distance between HSPs." );
        return false;
    }
    if( GetSubjectSize() <= 0 || GetQuerySize() <= 0 )
        return false;


    diagonals = ( TScore** )malloc( sizeof( TScore* ) * ( GetSubjectSize() + GetQuerySize() - 1 ));
    hsps = ( TScore** )malloc( sizeof( TScore* ) * GetSubjectSize());
    lens = ( int** )malloc( sizeof( int* ) * GetSubjectSize());

    if( !hsps || !lens )
        throw myruntime_error( mystring( "AttributableScoresII: Not enough memory." ));

    for( m = GetSubjectSize(); m < GetSubjectSize() + GetQuerySize() - 1; m++ ) {
        diagonals[m] = ( TScore* )malloc( sizeof( TScore ) * diaglen );
        if( !diagonals[m])
            throw myruntime_error( mystring( "AttributableScoresII: Not enough memory." ));

        memset( diagonals[m], 0, sizeof( TScore ) * diaglen );
    }
    for( m = 0; m < GetSubjectSize(); m++ ) {
        diagonals[m] = ( TScore* )malloc( sizeof( TScore ) * diaglen );
        hsps[m] = ( TScore* )malloc( sizeof( TScore ) * GetQuerySize());
        lens[m] = ( int* )malloc( sizeof( int ) * GetQuerySize());

        if( !hsps[m] || !lens[m] || !diagonals[m] )
            throw myruntime_error( mystring( "AttributableScoresII: Not enough memory." ));

        memset( diagonals[m], 0, sizeof( TScore ) * diaglen );
        memset( hsps[m], 0, sizeof( TScore ) * GetQuerySize());
        memset( lens[m], 0, sizeof( int ) * GetQuerySize());
    }

    //heuristics of several high-scoring pairs in the diagonal
    for( m = 0; m < GetSubjectSize() && !bret; m++ ) {
        for( n = 0; n < GetQuerySize() && !bret; n++ )
        {
//             if( GetMaskedToIgnore( m, n ))
//                 continue;

            loc_score = ( TScore )GetScore( m, n );

            if( loc_score <= SCORE_MIN ) {
                hsps[m][n] = 0;
                continue;
            }

            hsps[m][n] = loc_score;
            lens[m][n] = 0;

            if( 0 < loc_score ) {
                //begin or continue hsp
                hsps[m][n] = loc_score;
                lens[m][n] = 1;
            }

            if( 0 < m && 0 < n )
                //boundary hsp positions should be positive
//               (( ( lens[m-1][n-1] <= 1 || lens[m-1][n-1] == hsplen - 1 ) && 0 < hsps[m-1][n-1] ) ||
//               ( 1 < lens[m-1][n-1] && lens[m-1][n-1] < hsplen - 1 ) ))
            {
                tvalue = hsps[m][n] + hsps[m-1][n-1];
                if( hsps[m][n] < tvalue && 0 < hsps[m][n]) {
                    if( lens[m-1][n-1] < hsplen ) {
                        hsps[m][n] = tvalue;
                        lens[m][n] = lens[m-1][n-1] + 1;
                    }
                }
            }

            //record hsp score in the diagonal
            if( hsplen <= lens[m][n] && hmin <= hsps[m][n] )
            {
                ind = GetQuerySize() - 1 - n + m;

                sum = 0;
                for( k = 0; k < nohsps; k++ )
                {
                    if( diagonals[ind][k+nohsps] &&
                        maxdist < m - diagonals[ind][k+nohsps] - hsplen ) {
                        //distance between two hsps is too small- omit the hsp
                        //all hsps are considered in turn
                        sum = -1;//to break
                        break;
                    }
                    if( diagonals[ind][k] < hsps[m][n] )
                    {
                        for( l = nohsps - 1; k < l && diagonals[ind][l-1] <= 0; l-- );
                        for( ; k < l; l-- ) {
                            if( maxdist < m - diagonals[ind][l-1+nohsps] - hsplen ) {
                                //distance between two hsps is too small- omit the hsp
                                sum = -1;//to break
                                break;
                            }
                            diagonals[ind][l] = diagonals[ind][l-1];//save hsp score
                            diagonals[ind][l+nohsps] = diagonals[ind][l-1+nohsps];//save coordinates..
                            diagonals[ind][l+nohsps+nohsps] = diagonals[ind][l-1+nohsps+nohsps];
                            sum += diagonals[ind][l];
                        }
                        if( sum < 0 )
                            break;
                        diagonals[ind][k] = hsps[m][n];
                        diagonals[ind][k+nohsps] = m;
                        diagonals[ind][k+nohsps+nohsps] = n;
                        sum += diagonals[ind][k];
                        break;
                    }
                    else
                        sum += diagonals[ind][k];
                }
                if( diagonals[ind][nohsps-1] && hmax < sum ) {
                    bret = true;    //at least one group of hsps is found
                    hmax = sum;
                    mind = ind;
                    msb = diagonals[ind][nohsps];//positions of max of hsps
                    nqu = diagonals[ind][nohsps+nohsps];
                    if( possbjct )  *possbjct = msb;
                    if( posquery )  *posquery = nqu;
                }
            }
        }
    }

//{{TESTING BLOCK
// GetParent()->PrintScoringMatrix( stderr );
// fprintf( stderr, "\n\n");
// for( n = 0; n < GetQuerySize(); n++ ) {
//     if( !n ) {
//         fprintf( stderr, "****");
//         for( m = 0; m < GetSubjectSize(); m++ )
//             fprintf( stderr, " %5d", m );
//         fprintf( stderr, "\n");
//         for( m = 0; m < GetSubjectSize(); m++ )
//             fprintf( stderr, "------" );
//         fprintf( stderr, "\n");
//     }
//     fprintf( stderr, "%3d:", n );
//     for( m = 0; m < GetSubjectSize(); m++ ) {
//         fprintf( stderr, " %5.1f", ( double )hsps[m][n] / GetAutoScalingFactor());
//     }
//     fprintf( stderr, "\n");
// }
// fprintf( stderr, "\n%d %5.1f: %d %d\n\n", bret, ( double )hmax/* / GetAutoScalingFactor()*/, nqu, msb );
// if( bret ) {
//     for( k = 0; k < nohsps; k++ ) {
//         fprintf( stderr, "  score=%-.1f at (%d,%d)\n",
//             ( double )diagonals[mind][k]/*/GetAutoScalingFactor()*/,
//             diagonals[mind][k+nohsps+nohsps], diagonals[mind][k+nohsps]);
//     }
//     fprintf( stderr, "--------\n");
//     for( ind = 0; ind < GetSubjectSize() + GetQuerySize() - 1; ind++ ) {
//         if( diagonals[ind][nohsps-1] <= 0 )
//             continue;
//         for( k = 0; k < nohsps; k++ ) {
//             fprintf( stderr, "  score=%-.1f at (%d,%d)\n",
//                 ( double )diagonals[ind][k]/*/GetAutoScalingFactor()*/,
//                 diagonals[ind][k+nohsps+nohsps], diagonals[ind][k+nohsps]);
//         }
//         fprintf( stderr, "\n");
//     }
// }
//}}

    for( m = GetSubjectSize(); m < GetSubjectSize() + GetQuerySize() - 1; m++ ) {
        if( diagonals[m]) { free( diagonals[m]); diagonals[m] = NULL; }
    }
    for( m = 0; m < GetSubjectSize(); m++ ) {
        if( diagonals[m]) { free( diagonals[m]); diagonals[m] = NULL; }
        if( hsps[m]) { free( hsps[m]); hsps[m] = NULL; }
        if( lens[m]) { free( lens[m]); lens[m] = NULL; }
    }
    if( diagonals ) { free( diagonals ); diagonals = NULL; }
    if( hsps ) { free( hsps ); hsps = NULL; }
    if( lens ) { free( lens ); lens = NULL; }

    if( !errstr.empty())
        throw myruntime_error( errstr );

    return bret;
}

// -------------------------------------------------------------------------
// MakeQueryMask: makes query and subject masks according to the given Mask
//     matrix
// -------------------------------------------------------------------------

void AttributableScores::MakeQuerySbjctMasks( bool* queryvect, int qlen, bool* sbjctvect, int slen ) const
{
    int     szquery = GetQuerySize();
    int     szsubject = GetSubjectSize();
    int     m, n;

    if( queryvect == NULL || sbjctvect == NULL ||
      ( GetKeepInMemory() && ( qlen != szquery || slen != szsubject )))
        throw myruntime_error( mystring( "AttributableScores: Unable to make mask vectors." ));

    if( qlen == 0 || slen == 0 )
        return;

//     memset( queryvect, 0, sizeof( *queryvect ) * qlen );
//     memset( sbjctvect, 0, sizeof( *sbjctvect ) * slen );

    for( m = 0; m < szsubject; m++ )    sbjctvect[m] = false;
    for( n = 0; n < szquery; n++ )      queryvect[n] = false;

    for( m = 0; m < szsubject; m++ )    if( !GetMaskedUnmasked( m, 0 )) sbjctvect[m] = true;
    for( n = 0; n < szquery; n++ )      if( !GetMaskedUnmasked( 0, n )) queryvect[n] = true;

//     for( m = 0; m < szsubject; m++ )
//         for( n = 0; n < szquery; n++ )
//             if( !GetMaskedUnmasked( m, n )) {
//                 sbjctvect[m] = true;
//                 queryvect[n] = true;
//             }

}

// -------------------------------------------------------------------------
// ComputeHConstant: computes adjustment constant for relative entropy so
//     that transformation of scores by using it ensures relative entropy
//     desired to achieve
// -------------------------------------------------------------------------

void AttributableScores::ComputeHConstant()
{
    try {

        ComputeScoreProbabilities();
        if( GetExpectedScore() < 0.0 ) {
            ComputeScalingLambda();
            if( 0.0 < GetLambda()) {
//                 ComputeEntropyGivenLambda();
                ComputeConstantToAdjustH();
            }
        }

    } catch( myexception const& ex ) {
        DestroyProbabilities();
        throw( myruntime_error( ex.what(), ex.eclass()));
    }

    DestroyProbabilities();
}

// -------------------------------------------------------------------------
// AdjustGaps: calls apropriate methods to adjust gap costs to enable them
//     in run of the alignment algorithm
// -------------------------------------------------------------------------

void AttributableScores::AdjustGaps(
    const LogOddsMatrix& querylogo, const LogOddsMatrix& sbjctlogo,
    GapScheme& querygaps,           GapScheme& sbjctgaps,
    bool autoc,
    int acwindow ) const
{
    //make gap costs compatible with scores used
//     querygaps.Prepare( GetAutoScalingFactor());
//     sbjctgaps.Prepare( GetAutoScalingFactor());

    querygaps.SetScoresMultiplier( GetPrivateMultiplier());
    sbjctgaps.SetScoresMultiplier( GetPrivateMultiplier());

    querygaps.Prepare( GetAutoScalingFactor(), GetQueryInfContent());
    sbjctgaps.Prepare( GetAutoScalingFactor(), GetSbjctInfContent());

    if( autoc )
//         AdjustGapsByAutocorr( querylogo, sbjctlogo, querygaps, sbjctgaps, acwindow );
        AdjustGapsByAutocorr2( querylogo, sbjctlogo, querygaps, sbjctgaps, acwindow );
}

// -------------------------------------------------------------------------
// AdjustGapsByAutocorr: makes maximum score vectors of the scoring
//     matrix and further calls an appropriate method to adjusts gaps by
//     using these vectors
// -------------------------------------------------------------------------

void AttributableScores::AdjustGapsByAutocorr(
    const LogOddsMatrix& querylogo, const LogOddsMatrix& sbjctlogo,
    GapScheme& querygaps,           GapScheme& sbjctgaps,
    int acwindow ) const
{
#ifdef __DEBUG__
    if( GetParent() == NULL )
        throw myruntime_error( mystring( "AttributableScores: Unable to process for gap penalties." ));
#endif

    int     szquery = GetQuerySize();
    int     szsubject = GetSubjectSize();

    if( ! GetKeepInMemory()) {
        szquery = GetParent()->GetQuerySize();
        szsubject = GetParent()->GetSubjectSize();
    }

    if( szquery <= 0 || szsubject <= 0 )
        throw myruntime_error( mystring( "AttributableScores: Unable to compute gap costs." ));

    bool*   sbjctmasks = ( bool* )malloc( sizeof( bool ) * szsubject );
    bool*   querymasks = ( bool* )malloc( sizeof( bool ) * szquery );

    double* maxsbjctscores = ( double* )malloc( sizeof( double ) * szsubject );
    double* maxqueryscores = ( double* )malloc( sizeof( double ) * szquery );

    mystring    error;
    double      score;
    int         m, n;

    if( maxsbjctscores == NULL  || maxqueryscores == NULL ||
        sbjctmasks == NULL      || querymasks == NULL )
        throw myruntime_error( mystring( "AttributableScores: Not enough memory." ));

    //initialize vectors of the max scores
    for( m = 0; m < szsubject; m++ )    maxsbjctscores[m] = 0.0;
    for( n = 0; n < szquery; n++ )      maxqueryscores[n] = 0.0;

    for( m = 0; m < szsubject; m++ ) {
        for( n = 0; n < szquery; n++ )
        {
            ///verification of scoring scheme used
            if( GetKeepInMemory())
                score = GetScore( m, n );
            else
                score = GetParent()->GetScore( m, n );


            //consider maximum postive scores only
            if( maxsbjctscores[m] < score )
                maxsbjctscores[m] = score;

            if( maxqueryscores[n] < score )
                maxqueryscores[n] = score;
        }
    }

    try {
        MakeQuerySbjctMasks( querymasks, szquery, sbjctmasks, szsubject );
        querygaps.AdjustCosts( maxqueryscores, szquery, querylogo, acwindow, querymasks, szquery, GetParent()->GetEntropy());
        sbjctgaps.AdjustCosts( maxsbjctscores, szsubject, sbjctlogo, acwindow, sbjctmasks, szsubject, GetParent()->GetEntropy());
    } catch( myexception const& ex ) {
        error = ex.what();
    }

    free( sbjctmasks );
    free( querymasks );

    free( maxqueryscores );
    free( maxsbjctscores );

    if( !error.empty())
        throw myruntime_error( error );
}

// -------------------------------------------------------------------------
// Autocorrelation: for each vector slot, computes autocorrelation of
//     values found in the slot window.
// -------------------------------------------------------------------------

void AttributableScores::Autocorrelation( double ( *maxs )[ACWINDOW], int length, double* output, int outlen ) const
{
    if( maxs == NULL || output == NULL )
        throw myruntime_error( mystring( "AttributableScores: Memory access error." ));

    if( length != outlen )
        throw myruntime_error( mystring( "AttributableScores: Wrong lengths in autocorrelation computing." ));

    int midwin = ( ACWINDOW >> 1 ) - (( ACWINDOW & 1 ) ^ 1 );

    if( midwin < 0 )
        throw myruntime_error( mystring( "AttributableScores: Invalid autocorrelation window parameters." ));

    double  ro = 0.0;       //sum of autocorrelations
    double  aci, acj;
    double  icorrect = 0.00;
    size_t  no_pairs = 0;   //number of pairs, which is wsize! / ((wsize-2)! 2! )

    aci = acj = 0.0;

    for( int n = 0; n < length; n++ )
    {
        ro = 0.0;
        no_pairs = 0;

        if( ACWINDOW == 1 ) {
//         aci = PP_SCALE_CONSTANT * 1.5 * ( icorrect - info[0] ) * LN2 / LOSCORES.StatisParam( Ungapped, Lambda );
            ro += ( maxs[n][0] + aci )*( maxs[n][0] + aci );
            no_pairs++;
        }
        else
            for( int i = 0; i < ACWINDOW - 1; i++ ) {
                int j = i + midwin;
                if( ACWINDOW <= j )
                    j -= ACWINDOW;
                ro += ( maxs[n][i] + aci )*( maxs[n][j] + acj );
                no_pairs++;

//                 for( int j = i + 1; j < ACWINDOW; j++ ) {
//                 aci = PP_SCALE_CONSTANT * 1.5 * ( icorrect - info[i] ) * LN2 / LOSCORES.StatisParam( Ungapped, Lambda );
//                 acj = PP_SCALE_CONSTANT * 1.5 * ( icorrect - info[j] ) * LN2 / LOSCORES.StatisParam( Ungapped, Lambda );
//                     ro += ( maxs[n][i] + aci )*( maxs[n][j] + acj );
//                     no_pairs++;
//                 }
            }

        if( no_pairs )
            ro /= ( double ) no_pairs;
        if( 0.0 < ro )
            ro = sqrt( ro );
        else
            if( ro < 0.0 )
                ro = -sqrt( -ro );
        output[n] = ro;
    }
}

// -------------------------------------------------------------------------
// AdjustGapsByAutocorr2: makes autocorrelated vectors of maximum
//     scores of the scoring matrix and further calls an appropriate
//     method to adjusts gaps by using these vectors
// -------------------------------------------------------------------------

void AttributableScores::AdjustGapsByAutocorr2(
    const LogOddsMatrix& querylogo, const LogOddsMatrix& sbjctlogo,
    GapScheme& querygaps,           GapScheme& sbjctgaps,
    int acwindow ) const
{
#ifdef __DEBUG__
    if( GetParent() == NULL )
        throw myruntime_error( mystring( "AttributableScores: Unable to process for gap penalties." ));
#endif

    int     szquery = GetQuerySize();
    int     szsubject = GetSubjectSize();

    if( ! GetKeepInMemory()) {
        szquery = GetParent()->GetQuerySize();
        szsubject = GetParent()->GetSubjectSize();
    }

    if( szquery <= 0 || szsubject <= 0 )
        throw myruntime_error( mystring( "AttributableScores: Unable to compute gap costs." ));

    bool*   sbjctmasks = ( bool* )malloc( sizeof( bool ) * szsubject );
    bool*   querymasks = ( bool* )malloc( sizeof( bool ) * szquery );

    double  ( *maxsbjctscores )[ACWINDOW] = ( double (*)[ACWINDOW] )malloc( sizeof( double ) * ACWINDOW * szsubject );
    double  ( *maxqueryscores )[ACWINDOW] = ( double (*)[ACWINDOW] )malloc( sizeof( double ) * ACWINDOW * szquery );

    double* corrsbjctscores = ( double* )malloc( sizeof( double ) * szsubject );
    double* corrqueryscores = ( double* )malloc( sizeof( double ) * szquery );

    mystring    error;
    double  initial = 0.0; //SCORE_MIN; //initial minimum score
    double  score;
    int     m, n;
    size_t  i, j, k;

    if( maxsbjctscores == NULL  || maxqueryscores == NULL  ||
        corrsbjctscores == NULL || corrqueryscores == NULL ||
        sbjctmasks == NULL      || querymasks == NULL )
        throw myruntime_error( mystring( "AttributableScores: Not enough memory." ));


    //initialize vectors of the max scores
    for( m = 0; m < szsubject; m++ ) {
        for( i = 0; i < ACWINDOW; i++ )
            maxsbjctscores[m][i] = initial;
        corrsbjctscores[m] = initial;
    }

    for( n = 0; n < szquery; n++ ) {
        for( j = 0; j < ACWINDOW; j++ )
            maxqueryscores[n][j] = initial;
        corrqueryscores[n] = initial;
    }


    for( m = 0; m < szsubject; m++ ) {
        for( n = 0; n < szquery; n++ )
        {
            ///verification of scoring scheme used
            if( GetKeepInMemory())
                score = GetScore( m, n );
            else
                score = GetParent()->GetScore( m, n );


            //save the <ACWINDOW> largest values for subject and query
            for( i = 0; i < ACWINDOW; i++ )
                if( maxsbjctscores[m][i] < score )
                {
                    for( k = ACWINDOW - 1; k > i; k-- )
                        maxsbjctscores[m][k] = maxsbjctscores[m][k-1];
                    maxsbjctscores[m][i] = score;
                    break;
                }

            for( j = 0; j < ACWINDOW; j++ )
                if( maxqueryscores[n][j] < score )
                {
                    for( k = ACWINDOW - 1; k > j; k-- )
                        maxqueryscores[n][k] = maxqueryscores[n][k-1];
                    maxqueryscores[n][j] = score;
                    break;
                }
        }
    }


    try {
        MakeQuerySbjctMasks( querymasks, szquery, sbjctmasks, szsubject );

        Autocorrelation( maxqueryscores, szquery, corrqueryscores, szquery );
        Autocorrelation( maxsbjctscores, szsubject, corrsbjctscores, szsubject );

// fprintf( stderr, "\n" );
// for( m = 0; m < szsubject; m++ ) {
//     for( i = 0; i < ACWINDOW; i++ ) fprintf( stderr, "%7.2f ", maxsbjctscores[m][i] );
//     fprintf( stderr, " --- %7.2f\n", corrsbjctscores[m] );
// }
// fprintf( stderr, "\n\n" );
// for( n = 0; n < szquery; n++ ) {
//     for( j = 0; j < ACWINDOW; j++ ) fprintf( stderr, "%7.2f ", maxqueryscores[n][j] );
//     fprintf( stderr, " --- %7.2f\n", corrqueryscores[n] );
// }
// fprintf( stderr, "\n\n" );

        querygaps.AdjustCosts( corrqueryscores, szquery, querylogo, acwindow, querymasks, szquery, GetParent()->GetEntropy());
        sbjctgaps.AdjustCosts( corrsbjctscores, szsubject, sbjctlogo, acwindow, sbjctmasks, szsubject, GetParent()->GetEntropy());

    } catch( myexception const& ex ) {
        error = ex.what();
    }

    free( sbjctmasks );
    free( querymasks );

    free( maxsbjctscores );
    free( maxqueryscores );

    free( corrsbjctscores );
    free( corrqueryscores );

    if( !error.empty())
        throw myruntime_error( error );
}

// -------------------------------------------------------------------------
// ScaleScoringMatrix: scales scoring matrix so that its scores conform a 
//     reference scoring system.
// -------------------------------------------------------------------------

void AttributableScores::ScaleScoringMatrix()
{
    if( GetKeepInMemory()) {
        if( GetQuerySize() <= 0 || GetSubjectSize() <= 0 )
            throw myruntime_error( mystring( "Failed to scale matrix: Data are not initialized." ));

        if( GetAllNegatives())
            throw myruntime_error( mystring( "Scaling failed: Negative scores." ), SCALING );
    }

    ScaleToAttainLambda();
}

// -------------------------------------------------------------------------
// ComputeStatisticalParameters: compute all statistical parameters
//     required to calculate statistical significance of alignments
// -------------------------------------------------------------------------

void AttributableScores::ComputeStatisticalParameters( bool computelambda, bool wrn )
{
    if( computelambda )
        ComputeProbabilitiesAndLambda( wrn );
    ComputeKHGivenLambda();

// GetParent()->PrintScoringMatrix( stderr );
// fprintf( stderr, "  >  lambda=%f, entropy=%f, exp.score=%f\n", GetPrivLambda(), GetEntropy(), GetPrivExpectedScore());

}

// -------------------------------------------------------------------------
// ComputeProbabilitiesAndLambda: compute score probabilities and
//     statistical parameter lambda
//
void AttributableScores::ComputeProbabilitiesAndLambda( bool wrn )
{
    ComputeScoreProbabilities();
    if( GetExpectedScore() < 0.0 )
        ComputeScalingLambda( wrn );
}
// ComputeKHGivenLambda: compute entropy H and Karlin's parameter K
//
void AttributableScores::ComputeKHGivenLambda()
{
    try {
        if( GetExpectedScore() < 0.0 && 0.0 < GetLambda()) {
            ComputeEntropyGivenLambda();
            ComputeKarlinsK();
        }
    } catch( myexception const& ex ) {
        DestroyProbabilities();
        throw( myruntime_error( ex.what(), ex.eclass()));
    }
    //this method is to be last-called in the procedure of computation of staistical parameters
    //we don't need probabilities any more
    DestroyProbabilities();
}

// -------------------------------------------------------------------------
// ComputePositionalInfContents: compute entropies at each position of score
//     system
// -------------------------------------------------------------------------

void AttributableScores::ComputePositionalInfContents()
{
    double  lambda = LOSCORES.StatisParam( Ungapped, Lambda ); //reference lambda to make single profile scores
    int     szquery = GetQuerySize();
    int     szsbjct = GetSubjectSize();
    int     m, n, gcd, sc;
    int     minscore, maxscore;
    double  value, expected, probab;
    bool    lambdafailed = false;

    if( GetAutoScalingFactor() <= 0 )
        throw myruntime_error( mystring( "AttributableScores: ComputePositionalInfContents: "
                    "Wrong scale factor." ));

    lambda /= GetAutoScalingFactor();
    ComputePositionalScoreProbs();

    for( n = 0; n < szquery; n++ ) {
        minscore = GetQueryMinScoreAt( n );
        maxscore = GetQueryMaxScoreAt( n );
        value = expected = 0.0;
        gcd = -minscore;    //the greatest common divisor of the scores
        for( sc = minscore + 1; sc <= maxscore && 1 < gcd; sc++ ) {
            probab = GetQueryInfProbOf( sc, n );
            if( 0.0 < probab ) {
                gcd = PC_Gcd( gcd, sc - minscore );
                expected += probab * sc;
            }
        }
        //do not warn if expected score >=0, set entropy to 0 silently instead
        if( 0.0 <= expected ) {
            SetQueryInfContentAt( n, 0.0 );
            continue;
        }

        if( gcd == 0 ) {
            SetQueryInfContentAt( n, 0.0 );
            continue;
        }

        lambda = FindLambdaRootAtPosition(
            &AttributableScores::lambdaconspos, gcd, 0.0 + EPSILON, 10.0, LAMBDA_ACCURACY, MAXIT,
            minscore,
            maxscore,
            n,
            &AttributableScores::GetQueryInfProbOf
        );

        if( lambda < 0.0 + EPSILON ) {
            //be silent since many warnings may be produced
//             warning( "AttributableScores: ComputePositionalInfContents: Finding root of lambda failed." );
            lambdafailed = true;
            SetQueryInfContentAt( n, 0.0 );
            continue;
        }
        value = ComputePositionalEntropy(
                lambda,
                minscore,
                maxscore,
                n,
                &AttributableScores::GetQueryInfProbOf,
                &AttributableScores::GetQueryInfProbNorm
        );
        if( value < 0.0 )
            //it cannot be so but nevertheless avoid systematic errors
            value = 0.0;
        SetQueryInfContentAt( n, value );
    }
    //same with subject positions...
    //
    for( m = 0; m < szsbjct; m++ ) {
        minscore = GetSbjctMinScoreAt( m );
        maxscore = GetSbjctMaxScoreAt( m );
        value = expected = 0.0;
        gcd = -minscore;    //the greatest common divisor of the scores
        for( sc = minscore + 1; sc <= maxscore && 1 < gcd; sc++ ) {
            probab = GetSbjctInfProbOf( sc, m );
            if( 0.0 < probab ) {
                gcd = PC_Gcd( gcd, sc - minscore );
                expected += probab * sc;
            }
        }

        if( 0.0 <= expected ) {
            SetSbjctInfContentAt( m, 0.0 );
            continue;
        }

        if( gcd == 0 ) {
            SetSbjctInfContentAt( m, 0.0 );
            continue;
        }

        lambda = FindLambdaRootAtPosition(
            &AttributableScores::lambdaconspos, gcd, 0.0 + EPSILON, 10.0, LAMBDA_ACCURACY, MAXIT,
            minscore,
            maxscore,
            m,
            &AttributableScores::GetSbjctInfProbOf
        );

        if( lambda < 0.0 + EPSILON ) {
            //be silent: may be many warnings
//             warning( "AttributableScores: ComputePositionalInfContents: Finding root of lambda failed." );
            lambdafailed = true;
            SetSbjctInfContentAt( m, 0.0 );
            continue;
        }
        value = ComputePositionalEntropy(
                lambda,
                minscore,
                maxscore,
                m,
                &AttributableScores::GetSbjctInfProbOf,
                &AttributableScores::GetSbjctInfProbNorm
        );
        if( value < 0.0 )
            //it cannot be so but nevertheless avoid systematic errors
            value = 0.0;
        SetSbjctInfContentAt( m, value );
    }

    if( lambdafailed )
        warning( "AttributableScores: Unable to find lambda root at some positions." );
}

// -------------------------------------------------------------------------
// ComputePositionalScoreProbs: compute probabilities of scores at each
//     query and subject position
// -------------------------------------------------------------------------

void AttributableScores::ComputePositionalScoreProbs()
{
#ifdef __DEBUG__
    if( GetParent() == NULL )
        throw myruntime_error( mystring( "AttributableScores: ComputePositionalScoreProbs: "
                "Unable to compute probabilities." ));
#endif

    if( GetKeepInMemory()) {
        if( GetQuerySize() <= 0 || GetSubjectSize() <= 0 )
            throw myruntime_error( mystring( "AttributableScores: ComputePositionalScoreProbs: "
                    "Unable to compute probabilities." ));

        if( !GetScores() || !GetMask())
            throw myruntime_error( mystring( "AttributableScores: ComputePositionalScoreProbs: "
                    "Unable to compute score probabilities." ));

        TScore  l_minscore = 0,
                l_maxscore = 0;
        TScore  loc_score  = 0;
        TScore* querymins = NULL;
        TScore* querymaxs = NULL;
        int     m, n;

        querymins = ( TScore* )malloc( sizeof( TScore ) * GetQuerySize());
        querymaxs = ( TScore* )malloc( sizeof( TScore ) * GetQuerySize());

        if( querymins == NULL || querymaxs == NULL )
            throw myruntime_error( mystring( "AttributableScores: ComputePositionalScoreProbs: "
                    "Not enough memory." ));

        for( n = 0; n < GetQuerySize(); n++ ) {
            querymins[n] = querymaxs[n] = 0;
        }
        try {
            for( m = 0; m < GetSubjectSize(); m++ ) {
                l_minscore = l_maxscore = 0;
                for( n = 0; n < GetQuerySize(); n++ )
                {
                    if( GetMaskedToIgnore( m, n ))
                        continue;

                    loc_score = ( TScore )rint( GetScore( m, n ));

                    if( loc_score <= SCORE_MIN )
                        continue;

                    if( loc_score < l_minscore ) l_minscore = loc_score;
                    if( l_maxscore < loc_score ) l_maxscore = loc_score;

                    if( loc_score < querymins[n] ) querymins[n] = loc_score;
                    if( querymaxs[n] < loc_score ) querymaxs[n] = loc_score;
                }

                SetSbjctInfMinMaxScoresAt( l_minscore, l_maxscore, m );
            }

            for( n = 0; n < GetQuerySize(); n++ )
                SetQueryInfMinMaxScoresAt( querymins[n], querymaxs[n], n );

        } catch( myexception const& ex ) {
            if( querymins ) free( querymins );
            if( querymaxs ) free( querymaxs );
            throw myruntime_error( ex.what(), ex.eclass());
        }
        if( querymins ) free( querymins );
        if( querymaxs ) free( querymaxs );

        GetParent()->ComputePositionalScoreProbs( this );
    }
}





// -------------------------------------------------------------------------
// CheckforExpectedScore: compute score probabilities and indirectly adjust 
//  scores by expected value
//
void AttributableScores::CheckforExpectedScore()
{
    const int       maxit = 1;
    bool            reded;
    int i;

    if( !GetParentRedtFunction())
        return;
    for( i = 0; i <= maxit && 0.0 <= GetExpectedScore(); i++ ) {
        if( i ) {
            reded = ( GetParent()->*GetParentRedtFunction())( this );
            if( !reded )
                break;
        }
        ComputeScoreProbabilities();
    }//for(;0<GetExpectedScore();)
}

// -------------------------------------------------------------------------
// ComputeScoreProbabilities: compute probabilities to observe scores at
//     each position of the scoring system
//
void AttributableScores::ComputeScoreProbabilities()
{
    const double    accuracy = 1.e-5;
    char            strbuf[KBYTE];
    double          prob, consv;
    int m, n, sc;

#ifdef __DEBUG__
    if( GetParent() == NULL || GetParentProbFunction() == NULL )
        throw myruntime_error( mystring( "AttributableScores: Unable to compute score probabilities." ));
#endif

    if( GetKeepInMemory()) {
        if( GetQuerySize() <= 0 || GetSubjectSize() <= 0 )
            throw myruntime_error( mystring( "AttributableScores: Unable to compute score probabilities." ));

        if( !GetScores() || !GetMask())
            throw myruntime_error( mystring( "AttributableScores: Unable to compute score probabilities." ));

        TScore  l_minscore = 0,
                l_maxscore = 0;
        TScore  loc_score  = 0;

        for( m = 0; m < GetSubjectSize(); m++ ) {
            for( n = 0; n < GetQuerySize(); n++ ) {
                if( GetMaskedToIgnore( m, n ))
                    continue;

                loc_score = ( TScore )rint( GetScore( m, n ));

                if( loc_score <= SCORE_MIN )
                    continue;

                if( loc_score < l_minscore ) l_minscore = loc_score;
                if( l_maxscore < loc_score ) l_maxscore = loc_score;
            }
        }
        SetMinMaxScores( l_minscore, l_maxscore );
    }
    ( GetParent()->*GetParentProbFunction())( this );

    consv = 0.0;
    for( sc = GetMinScore(); sc <= GetMaxScore(); sc ++ ) {
        prob = GetProbabilityOf( sc );
        if( !prob )
            continue;
        consv += prob;
    }
    if( consv < 1.0 - accuracy || consv > 1.0 + accuracy ) {
        sprintf( strbuf, "Score probabilities not conserved: %f", consv );
        throw myruntime_error( strbuf );
    }
}

// -------------------------------------------------------------------------
// AdjustScoresInMatrix: adjusts scores in the matrix given previously
//     computed scaling parameter lambda
// -------------------------------------------------------------------------

void AttributableScores::AdjustScoresInMatrix()
{
    if( GetLambda() <= 0.0 ) {
        warning( "AttributableScores: Unable to scale matrix." );
        return;
    }

    double multiplier = GetLambda() / GetRefLambda();
    AdjustScoresInMatrix( multiplier, ImageTable );
}

// -------------------------------------------------------------------------
// AdjustScoresInMatrix: adjusts scores in the matrix by multiplying them 
//     by a given factor
// -------------------------------------------------------------------------

void AttributableScores::AdjustScoresInMatrix( double multiplier, TTable table_type, bool direct )
{
#ifdef __DEBUG__
    if( GetParent() == NULL || GetParentProgressFunction() == NULL )
        throw myruntime_error( mystring( "AttributableScores: No parent method address given." ));

    if( GetKeepInMemory() && ( !GetImage() || !GetScores() || !GetMask()))
        throw myruntime_error( mystring( "AttributableScores: Unable to adjust scores." ));
#endif

    double  local_scale_factor = GetAutoScalingFactor();


    if( GetKeepInMemory()) {
        for( int m = 0; m < GetSubjectSize(); m++ ) {
            for( int n = 0; n < GetQuerySize(); n++ )
            {
                if( SCORE_MIN < image[m][n] /*&& ! GetMaskedToIgnore( m, n )*/)
                    if( table_type == ImageTable )
                        SetScore( image[m][n] * multiplier, m, n );
                    else
                        MultiplyScoreBy( multiplier, m, n );
            }
        }
    } else {
        if( direct )
            ( GetParent()->*GetParentProgressFunction())( this );
        SetPrivateMultiplier( multiplier );
    }
}

// -------------------------------------------------------------------------
// ComputeGCD: compute the greatest common divisor of all scores
// -------------------------------------------------------------------------

void AttributableScores::ComputeGCD()
{
    int     gcd = 1;    //the greatest common divisor of the scores

    gcd = -GetMinScore();
    for( int sc = GetMinScore() + 1; sc <= GetMaxScore() && gcd > 1; sc++ )
        if( 0.0 < GetProbabilityOf( sc ))
            gcd = PC_Gcd( gcd, sc - GetMinScore());

    if( gcd == 0 )
        return;

    SetGCD( gcd );
}

// -------------------------------------------------------------------------
// ComputeScalingLambda: compute statistical (scaling) parameter lambda
// -------------------------------------------------------------------------

void AttributableScores::ComputeScalingLambda( bool wrn )
{
    mystring    msgbuf;
    double  newlambda = -1.0;

    if( 0.0 <= GetExpectedScore()) {
        warning( "Expected score per position is positive." );
    } else {

        //rare will be cases when interval between adjacent scores is >1;
        //nevertheless, compute the greatest common divisor of all scores;
        //change of variables then will be made via substitution y=exp(-lambda gcd)
        ComputeGCD();

        newlambda = findLambdaRoot(
            &AttributableScores::conservation, GetGCD(), 0.0 + EPSILON, 10.0,
            LAMBDA_ACCURACY, MAXIT
        );

        if( newlambda <= 0.0 + EPSILON ) {
            //newlambda = LOSCORES.StatisParam( Ungapped, Lambda );   //to leave score values unchanged
            if( wrn ) {
                msgbuf = "Computation of statistical parameter lambda failed";
                if( GetParent() && GetParent()->GetName())
                        msgbuf += mystring(": ")+ GetParent()->GetName();
                else    msgbuf += ".";
                warning( msgbuf.c_str());
//                 warning("*", false );
            }
            //return;
        }
    }

    SetLambda( newlambda );
}

// -------------------------------------------------------------------------
// ComputeConstantToAdjustH: compute constant c to adjust relative entropy
// -------------------------------------------------------------------------

void AttributableScores::ComputeConstantToAdjustH()
{
    double  c = -1.0;       //constant c to be found
    double  lb = EPSILON_C; //EPSILON_C; lower bound for c
    double  ub = 1000.0;    //upper bound for c

    if( 0.0 <= GetExpectedScore()) {
        warning( "Unable to adjust H: Expected score is positive." );
    } else {

        //compute the greatest common divisor of all scores
        ComputeGCD();

        c = findConstantRoot(
            &AttributableScores::relentropy_conservation, GetGCD(), lb, ub,
            ACCURACY_OF_CONSTANT_FOR_H, MAXIT
        );

        if( c < lb ) {
            warning( "Failed to find a root for constant to adjust H." );
            //return;
        }
    }

    SetConstantForHAdjustment( c );
}

// -------------------------------------------------------------------------
// ComputeKarlinsK: computes parameter K analitically as given by
//     Karlin et al. PNAS 87 (1990)
// K is an extreme value distribution parameter related to its location
// parameter mu by relation
//   mu = ln (Kmn) / lambda
// K is computed by the formula
//
//       gcd lambda exp( -2 sigma )
//   K = --------------------------
//       H (1 - exp( -gcd lambda ))
//
// where gcd is the greatest common divisor of the scores, H is relative
// entropy, and sigma is roughly a probability to obtain an arbitrary
// alignment of any length >=1 given scoring scheme
//
//            _   1 (  _           i lambda     _          )
//   sigma = \   -- ( \   P[j](i) e          + \   P[j](i) )
//           /_   j ( /_                       /_          )
//           j>0      i<0                     i>=0
//
// P[j](i) is a probability to obtain an alignment of score i and length j
//              _
//   P[j](i) = \  P[1](k) P[j-1](i-k)
//             /_
//              k
//
// There are three simplifications: if high score = gcd and low score = -gcd,
//
//   K = squared( P[1](gcd) - P[1](-gcd) )  / p[1](-gcd)
//
// if high score = gcd and low score <> -gcd,
//
//   K = H (1 - exp(-gcd lambda )) / (gcd lambda)
//
// if high score <> gcd and low score = -gcd,
//                                                          _
//   K = lambda (1 - exp(-gcd lambda )) / (gcd H) squared( \  i P[1](i) )
//                                                         /_
//                                                        i=l:u
// -------------------------------------------------------------------------

void AttributableScores::ComputeKarlinsK()
{
    double  loc_lambda = GetLambda();
    double  loc_expected_score = GetExpectedScore();

    if( 0.0 <= loc_expected_score ) {
        warning( "Unable to compute statistical parameter K." );
        return;
    }
    if( loc_lambda < 0.0 || GetH() <= 0.0 ) {
        warning( "Unable to compute statistical parameter K." );
        return;
    }

    double  K = -1.0;

    const int       maxrange    = MAX_RANGE_OF_SCORES;
    const int       maxit       = PARAMETER_K_MAXIT;
    const double    precision   = PARAMETER_K_ACCURACY;

    TScore  gcd = GetGCD();
    TScore  range = GetMaxScore() - GetMinScore();
    double  y = exp( -gcd * loc_lambda );

    if( maxrange < range )
        throw myruntime_error( mystring( "AttributableScores: Range of scores exceeds the maximum allowed." ), SCALING );

    if( gcd <= 0 || range <= 0 )
        throw myruntime_error( mystring( "AttributableScores: Wrong range of scores obtained." ), SCALING );


    if( GetMinScore() == -gcd && GetMaxScore() == gcd ) {
        K = SQUARE( GetProbabilityOf( GetMinScore()) - GetProbabilityOf( GetMaxScore()) ) /
                    GetProbabilityOf( GetMinScore());
        SetK( K );
        return;
    }

    if( GetMinScore() == -gcd || GetMaxScore() == gcd ) {
        if( GetMaxScore() != gcd )
            K = loc_lambda / ( gcd * GetH()) * ( 1.0 - y ) * SQUARE( loc_expected_score );
        else
            K = GetH() / ( gcd * loc_lambda ) * ( 1.0 - y );
        SetK( K );
        return;
    }

    double  sigma = 0.0;    //probability of alignment of any length
    double  Pj  = 1.0;      //probability to obtain arbitrary alignment of length j
    double  Pji = 0.0;      //probability to obtain alignment of length j with score i
    double* Pjivector =                 //vector of alignment score probabilities P[j](i),
        GetPrivateProbVector();         // where j is an alignment length, and i is a score
    TScore  min_j_score = 0;    //minimum score of alignment of length j
    TScore  max_j_score = 0;    //minimum score of alignment of length j
    TScore  lower = 0;          //lower bound of scores
    TScore  upper = 0;          //upper bound of scores

    size_t  maxPjisize = maxit * range + 1; //plus to include zero

    if( !Pjivector ) {
        error( "Probability vector is not allocated." );
        return;
    }

    ResetPrivateProbVector( maxPjisize );
    Pjivector[0] = 1.0;     //in order to compute the lowest score probability of alignment with length=1


    //compute sigma, probability of alignment of any length
    // Pjivector is reused in each iteration to restore probabilities of
    // alignment of length one smaller
    for( int j = 1; j <= maxit && precision < Pj; j++ )
    {
        TScore  i;
        min_j_score += GetMinScore();
        max_j_score += GetMaxScore();
        lower = GetMaxScore();
        upper = GetMaxScore();

        //compute Pj for all possible scores i
        for( i = max_j_score; min_j_score <= i; i -= gcd )
        {
            Pji = 0.0;
            //P[j](i) = SUM P[j-1](i-k) P[1](k)
            for( int k = lower; k <= upper; k += gcd ) {
                Pji += Pjivector[ ( i - min_j_score ) - ( k - GetMinScore()) ] * GetProbabilityOf( k );
            }

            if( GetMinScore() < lower )
                lower -= gcd;
            if( i - min_j_score <= range )
                upper -= gcd;

            Pjivector[ i - min_j_score ] = Pji;
        }

        Pj = 0.0;
        //compute SUM (i<0) ( P[j](i) exp(i lambda) );
        //actually we need to compute polynomial of exp(lambda) with
        // coefficients P[j](i)
        for( i = min_j_score; i < 0; i += gcd )
            Pj = Pj * y + Pjivector[ i - min_j_score ];
        Pj *= y;

        //compute SUM (i>=0) ( P[j](i) );
        //just sum previoulsy computed probabilities P[j](i)
        for( ; i <= max_j_score; i += gcd )
            Pj += Pjivector[ i - min_j_score ];

        sigma += Pj / j;
    }

    //we have all variables found to compute K
    K = ( gcd * loc_lambda * exp( -2.0 * sigma )) / ( GetH() * ( 1.0 - y ));
    SetK( K );
}

// -------------------------------------------------------------------------
// rootByNRandBisection: finds the root of a function bracketed between x1
//     and x2 using a combination of Newton-Raphson and bisection methods.
//     A Bisection step is taken whenever Newton-Raphson would take the solution
//     out of bounds, or whenever it is not reducing the size of the brackets
//     rapidly enough.
//     The root will be refined until its accuracy is reached within ±xacc or
//     maximum number of iteration, maxit, has been passed. The method must be
//     supplied with a supplied routine that returns both the function value
//     and the first derivative of the function.
// Computation logic is as found in
//     Numerical Recipes in C by W.H.Press et al.
// -------------------------------------------------------------------------

double AttributableScores::rootByNRandBisection(
        TConservationFunction fdfunction,
        int gcd,
        double x1,
        double x2,
        double xacc,
        int maxit,
        void* params,
        bool warn )
{
    double  df, dx, dxold, f, fh, fl;
    double  temp, xh, xl, rts;
    int     n;

    ( this->*fdfunction )( x1, &fl, &df, gcd, params );
    ( this->*fdfunction )( x2, &fh, &df, gcd, params );

    if(( fl > 0.0 && fh > 0.0 ) || ( fl < 0.0 && fh < 0.0 )) {
        if( warn )
            warning( "Root evaluation: Root must be bracketed in interval supplied." );
        return -1.0;                //indicates error
    }


    if( fl < 0.0 ) {                //Orient the search so that f (xl) < 0.
        xl = x1;
        xh = x2;
    } else {
        xh = x1;
        xl = x2;
    }
    rts = 0.5 * ( x1 + x2 );        //Initialize the guess for root,
    dxold = fabs( x2 - x1 );        //the "stepsize before last,"
    dx = dxold;                     //and the last step.

    ( this->*fdfunction )( rts, &f, &df, gcd, params );

    for( n = 0; n < maxit; n++ )                            //Loop over allowed iterations.
    {                                                       //Bisect if Newton out of range,
        if( ((( rts - xh ) * df - f ) * (( rts - xl ) * df - f ) > 0.0 ) ||
               ( fabs( 2.0 * f ) > fabs( dxold * df ))) {   //or not decreasing fast enough.
            dxold = dx;
            dx = 0.5 * ( xh - xl );
            rts = xl + dx;
            if( xl == rts ) break; //return rts;            //Change in root is negligible.
        } else {                                            //Newton step acceptable. Take it.
            dxold = dx;
            dx = f / df;
            temp = rts;
            rts -= dx;
            if( temp == rts ) break; //return rts;
        }
        if( fabs( dx ) < xacc ) break; //return rts;        //Convergence criterion.
        ( this->*fdfunction )( rts, &f, &df, gcd, params );
        //The one new function evaluation per iteration Maintain the bracket on the root.
        if( f < 0.0 )
            xl = rts;
        else
            xh = rts;
    }

    if( n == maxit ) {
        error( "Root evaluation: maximum number of iterations reached." );
        return -1.0;
    }

    return rts;
}

// -------------------------------------------------------------------------
// testRoot: tests root found by the NR searching algorithm
// -------------------------------------------------------------------------

void AttributableScores::testRoot()
{
    double  f, df;
    conservation( lambda, &f, &df );
    if( 1.e-5 < fabs( f ))
            fprintf( stdout, "Error: Function value is not zero, " );
    else    fprintf( stdout,  "Function value is " );
    fprintf( stdout, "%f\n\n", f );
}

// -------------------------------------------------------------------------
// PrintProbabilities: print score probabilities if available
// -------------------------------------------------------------------------

void AttributableScores::PrintProbabilities( FILE* fp )
{
    if( fp == NULL )
        return;

    if( !probabilities )
        return;

    fprintf( fp,"\n%5c Score probabilities\n", 32 );

    fprintf( fp, "%9c", 32 );

    for( int sc = GetMinScore(); sc <= GetMaxScore(); sc++ ) {
        fprintf( fp, "\n%5d %6.4g", sc, GetProbabilityOf( sc ));
    }
    fprintf( fp, "\n" );
}


