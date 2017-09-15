/***************************************************************************
 *   Copyright (C) 2012 by Mindaugas Margelevicius                         *
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
#include "PrlHDPUniversalScoreMatrix.h"



//if true, background probabilities are multivariate t-dist.
bool        PrlHDPUniversalScoreMatrix::hdpbprobs_ = false;

// =========================================================================
// CLASS PrlHDPUniversalScoreMatrix
//
// constructor
//
PrlHDPUniversalScoreMatrix::PrlHDPUniversalScoreMatrix(
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
    ParallelUniversalScoreMatrix(
        store,
        config,
        no_scaling,
        using_spec_lambda,
        infrm_threshold,
        distrib,
        master_flag,
        rank,
        ring_size,
        func_bcast,
        func_send,
        func_recv,
        func_block,
        a_scaling,
        c_masking,
        ParallelHDPUniversal ),
    basin_( NULL ),
    hdpbase_( NULL ),
    prnormterm_( 0.0 ),
    nosupports_( 1 )
{
}

// -------------------------------------------------------------------------
// default constructor is invalid
//
PrlHDPUniversalScoreMatrix::PrlHDPUniversalScoreMatrix()
:   ParallelUniversalScoreMatrix(),
    basin_( NULL ),
    hdpbase_( NULL ),
    prnormterm_( 0.0 ),
    nosupports_( 1 )
{
}

// -------------------------------------------------------------------------
// destructor
//
PrlHDPUniversalScoreMatrix::~PrlHDPUniversalScoreMatrix()
{
    DestroyHDPbase();
    DestroyBasin();
}

// -------------------------------------------------------------------------
// ScaleScoringMatrix: call InitHDP and parent method
//
void PrlHDPUniversalScoreMatrix::ScaleScoringMatrix()
{
    if( IsMasterProcess())
        message("  Initializing HDP...", false );

    InitHDP();

    if( IsMasterProcess())
        message("  done.");

    ParallelUniversalScoreMatrix::ScaleScoringMatrix();
}

// -------------------------------------------------------------------------
// InitHDP: initialize hdp-related structures;
//  split all vectors to two set; one of which is support set, another one
//  comprise vectors to compare with
//
void PrlHDPUniversalScoreMatrix::InitHDP()
{
    if( !GetStore() || !GetStore()->GetFrequencies())
        USM_THROW( "PrlHDPUniversalScoreMatrix: InitHDP: Null frequencies store." );

    const SimpleVector* frequencies = GetStore()->GetFrequencies();
    size_t              no_freqs = frequencies->GetSize();

    const double    lambda = LOSCORES.StatisParam( Ungapped, Lambda );
    const double    tol = 0.1;//error tolerance for target probabilities
    const int       dim = NUMAA - 1;
    double  tprob;//target probability from profile
    double  pscor;//score from profile
    double  lp;//log-probability
    Dish*   dsh;
    Pslvector*  ps;
    int propsize;//subject proposal size
    int bnpos, npos, kpos, dpos;
    int N, a;

    if( no_freqs < 1 )
        return;

    SetNoSupPoss( SLC_MAX( 1, int( 0.5 *( double )no_freqs )));
    if( GetNoSupPoss() < 1 )
        USM_THROW("PrlHDPUniversalScoreMatrix: InitHDP: Too small number of support positions.");

    InitBasin( no_freqs + 1 ),
    InitHDPbase();

    //one position is for background prob. vector
    GetHDPbase()->ReserveMenu( GetNoSupPoss() + 1 );
    GetHDPbase()->ReserveBasin( GetNoSupPoss() + 1 );
    GetHDPbase()->GetBasin()->SetDestroy( false );

    for( N = 0; N <= no_freqs; N++ ) {
        const FrequencyVector   frqvec(( const char* )frequencies->GetValueAt((N<no_freqs)? N: N-1 ));

        dsh = new Dish( 1 );
        ps = new Pslvector( NUMAA+0 );
        if( ps == NULL )
            throw myruntime_error("PrlHDPUniversalScoreMatrix: InitHDP: Not enough memory.");
        dsh->SetBasin( GetBasin());
        for( a = 0; a < NUMAA; a++ ) {
            tprob = LOSCORES.PROBABility( a );
            //the last position is for background probability vector
            if( N < no_freqs ) {
                pscor = ( double )frqvec.GetScoreAt( a ) /( double )SCALE_CONSTANT;
                tprob *= exp( pscor * lambda );
            }
            ps->SetValueAt( a, tprob );
        }
//         ps->SetValueAt( NUMAA, 1.e-3 );
        HDPbase::LogitNormal2Normal( ps, tol );
        //reduce dimensionality
        ps->DecDim();
        bnpos = GetBasin()->NewValue( ps );
        if( N < GetNoSupPoss() || no_freqs <= N ) {
            if( GetNoSupPoss() <= 3000 || !( N % 1000 ))
                //when large number of various vecs is placed,
                //prior distribution encompasses nearly all space, implying low scores
                npos = GetHDPbase()->AddToBasin( ps );
            kpos = GetHDPbase()->GetMenu()->NewDish( dsh );
//             propsize = 1;
            dpos = dsh->NewVectorNInd( bnpos );
//             dsh->SetReadSize( propsize );
        }
    }

    GetHDPbase()->SetDegFAdjustment( 0.0 );
    GetHDPbase()->CalcPriorParams();
    //reduce weight of prior parameters after calculating them
    GetHDPbase()->GetMenu()->SetKappa0( 10.0 );
    GetHDPbase()->GetMenu()->SetNu0( 10.0 );
    //calculate parameters for each dish
    GetHDPbase()->RecalcMenuParams();

    //precalculate log-probabilities for conditional probabilities
    GetQueValues().Reserve( no_freqs );
    for( N = 0; N < no_freqs; N++ ) {
        if( N < GetNoSupPoss())
            continue;
        ps = GetBasin()->GetValueAt( N );
        if( ps == NULL )
            throw myruntime_error("PrlHDPUniversalScoreMatrix: InitHDP: Memory access error.");
        GetHDPbase()->SetDegFAdjustment( 21.0 );
//         GetHDPbase()->PriorProbVec( ps, &lp );//P(Q)-- prior probability
        GetHDPbase()->ProbVecOfDish( ps, GetNoSupPoss(), &lp );//P(Q|B), where B is background vector
        GetQueValues().SetValueAt( N, lp );
    }

    if( GetUseHDPBProbs()) {
        //precalculate log-probabilities for prior probabilities
        GetSubValues().Reserve( GetNoSupPoss());
        for( N = 0; N < GetNoSupPoss(); N++ ) {
            ps = GetBasin()->GetValueAt( N );
            if( ps == NULL )
                throw myruntime_error("PrlHDPUniversalScoreMatrix: InitHDP: Memory access error.");
            GetHDPbase()->SetDegFAdjustment( 0.0 );
            GetHDPbase()->PriorProbVec( ps, &lp );
            GetSubValues().SetValueAt( N, lp );
        }

        SetProbNormTerm( SLC_MAX( GetQueValues().Max(), GetSubValues().Max()));
    }
}

// -------------------------------------------------------------------------
// ComputeScore: compute score given two vectors of frequencies
//
double PrlHDPUniversalScoreMatrix::ComputeScore(
    const FrequencyVector& vect_1st,
    const FrequencyVector& vect_2nd, 
    size_t N, size_t M, bool revsd )
{
    double  score = SCORE_MIN;
    size_t  fst_enos = vect_1st.GetThickness();
    size_t  sec_enos = vect_2nd.GetThickness();
    size_t  r;

    if(( N < GetNoSupPoss() && M < GetNoSupPoss()) ||
       ( GetNoSupPoss() <= N && GetNoSupPoss() <= M ))
        //this score will be omitted
        return score;

    double  qlp, slp;
    double  qslp, qblp;
    int dim = NUMAA - 1;
    int d;
    Pslvector*  qv, *sv;
//     Pslmatrix   qsv( dim, 2 );

    if( GetHDPbase() == NULL || GetBasin() == NULL )
        throw myruntime_error("PrlHDPUniversalScoreMatrix: ComputeScore: Null HDP structures.");

    //always N < M if !revsd
    if( !revsd ) {
        r = M; M = N; N = r;
        r = sec_enos; sec_enos = fst_enos; fst_enos = r;
    }
    if( fst_enos < sec_enos ) {
        //vectors should be sorted according to eff. nos
        throw myruntime_error("PrlHDPUniversalScoreMatrix: ComputeScore: Vectors are not sorted.");
    }

    qv = GetBasin()->GetValueAt( N );
    if( qv == NULL )
        throw myruntime_error("PrlHDPUniversalScoreMatrix: ComputeScore: Memory access error.");

    GetHDPbase()->SetDegFAdjustment( 21.0 );

    GetHDPbase()->ProbVecOfDish( qv, M, &qslp );//P(Q|S)
    qblp = qlp = GetQueValues().GetValueAt( N );//P(Q|B)

    score = qslp - qblp;
    score = SLC_MIN( 10.0, score );

    score *= GetMultiplier();
    return score;
}

// -------------------------------------------------------------------------
// VectorScoreProbability: compute score probability given two vectors of
//  frequencies
//
double PrlHDPUniversalScoreMatrix::VectorScoreProbability(
    const FrequencyVector& rowvect,
    const FrequencyVector& colvect, 
    size_t N, size_t M, bool revsd ) const
{
    double  pro1, pro2;
    size_t  r;

    if(( N < GetNoSupPoss() && M < GetNoSupPoss()) ||
       ( GetNoSupPoss() <= N && GetNoSupPoss() <= M ))
        //score of this probability will be omitted
        return 0.0;

    if( !GetUseHDPBProbs())
        return rowvect.GetProbability() * colvect.GetProbability();


    if( GetHDPbase() == NULL )
        throw myruntime_error("PrlHDPUniversalScoreMatrix: VectorScoreProbability: Null HDP structure.");

    if( !revsd ) {
        r = M; M = N; N = r;
    }

    pro1 = GetQueValues().GetValueAt( N );
    pro2 = GetSubValues().GetValueAt( M );
    
    pro1 = exp( pro1 - GetProbNormTerm());
    pro2 = exp( pro2 - GetProbNormTerm());

    return pro1 * pro2;
}

