/***************************************************************************
 *   Copyright (C) 2014 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include "liblib/logitnormal.h"
#include "libpro/srcpro/Serializer.h"
#include "libpro/srcpro/GapScheme.h"
#include "libpro/srcpro/DistributionMatrix.h"
#include "RndMSAGenerator.h"


const char* slg_MSAEXT = ".fa";

// -------------------------------------------------------------------------
// CLASS RndMSAGenerator
//
// Constructors
//
RndMSAGenerator::RndMSAGenerator( const char* outdir, const char* namepat )
:
    outdir_( outdir ),
    nampat_( namepat ),
    gen_lenmsa_( 0 ),
    gen_nomsas_( 0 ),
    gen_msawidth_( 0 ),
    gen_msaeffnos_( 0 ),
    gen_trprobs_( NULL ),
    gen_reprobs_( NULL ),
    gen_model_( false ),
    gen_noise_( 0.0 ),
    msa_( NULL ),
    seqs_( NULL ),
    stts_( NULL ),
    maxi_( NULL ),
    options_( NULL ),
    drv_em_( NULL ),
    drv_in_( NULL ),
    drv_tr_m_( NULL ),
    drv_tr_i_( NULL ),
    drv_tr_d_( NULL ),
    mnrv_tr_( NULL ),
    mnrv_em_( NULL ),
    trmv_( 2 ),
    trco_( 2 ),
    remv_( NUMAA-1 ),
    reco_( NUMAA-1 )
{
    Init();
}

// Default constructor
//
RndMSAGenerator::RndMSAGenerator()
:
    outdir_( NULL ),
    nampat_( NULL ),
    gen_lenmsa_( 0 ),
    gen_nomsas_( 0 ),
    gen_msawidth_( 0 ),
    gen_msaeffnos_( 0 ),
    gen_trprobs_( NULL ),
    gen_reprobs_( NULL ),
    gen_model_( false ),
    gen_noise_( 0.0 ),
    msa_( NULL ),
    seqs_( NULL ),
    stts_( NULL ),
    maxi_( NULL ),
    options_( NULL ),
    drv_em_( NULL ),
    drv_in_( NULL ),
    drv_tr_m_( NULL ),
    drv_tr_i_( NULL ),
    drv_tr_d_( NULL ),
    mnrv_tr_( NULL ),
    mnrv_em_( NULL ),
    trmv_( 2 ),
    trco_( 2 ),
    remv_( NUMAA-1 ),
    reco_( NUMAA-1 )
{
    throw( myruntime_error("RndMSAGenerator: Default initialization is not allowed."));
}

// Destructor
//
RndMSAGenerator::~RndMSAGenerator()
{
    DestroyMSA();
    DestroySStructures();
    if( drv_em_ ) { delete drv_em_; drv_em_ = NULL; }
    if( drv_in_ ) { delete drv_in_; drv_in_ = NULL; }
    if( drv_tr_m_ ) { delete drv_tr_m_; drv_tr_m_ = NULL; }
    if( drv_tr_i_ ) { delete drv_tr_i_; drv_tr_i_ = NULL; }
    if( drv_tr_d_ ) { delete drv_tr_d_; drv_tr_d_ = NULL; }
    if( mnrv_tr_ ) { delete mnrv_tr_; mnrv_tr_ = NULL; }
    if( mnrv_em_ ) { delete mnrv_em_; mnrv_em_ = NULL; }
}

// -------------------------------------------------------------------------
// Init: initialization
//
void RndMSAGenerator::Init()
{
    int n;
    time_t tm;
    time( &tm );
    rng_em_.Set((unsigned long)(size_t)( &rng_em_ ) +(unsigned long)tm );
    rng_in_.Set((unsigned long)(size_t)( &rng_in_ ) +(unsigned long)tm );
    rng_tr_.Set((unsigned long)(size_t)( &rng_tr_ ) +(unsigned long)tm );
    drv_em_ = new RVDisc( rng_em_, RVDisc::TRVDsc_Inverse );
    drv_in_ = new RVDisc( rng_in_, RVDisc::TRVDsc_Inverse );
    drv_tr_m_ = new RVDisc( rng_tr_, RVDisc::TRVDsc_Inverse );
    drv_tr_i_ = new RVDisc( rng_tr_, RVDisc::TRVDsc_Inverse );
    drv_tr_d_ = new RVDisc( rng_tr_, RVDisc::TRVDsc_Inverse );
    if( drv_em_ == NULL || drv_in_ == NULL || 
        drv_tr_m_ == NULL || drv_tr_i_ == NULL || drv_tr_d_ == NULL )
        throw( myruntime_error("RndMSAGenerator: Init: Not enough memory."));
    for( n = 0; n < NUMAA; n++ )
        gen_inprobs_[n] = LOSCORES.PROBABility(n);
    //reset original probabilities
    memset( org_trprobs_, 0, P_NSTATES * sizeof(double));
    memset( org_reprobs_, 0, NUMAA * sizeof(double));
    //set mean vectors and covariance matrices for 
    // multivariate normal random variable
    mnrv_tr_ = new RMVNorm( rng_tr_ );
    mnrv_em_ = new RMVNorm( rng_em_ );
    if( mnrv_tr_ == NULL || mnrv_em_ == NULL )
        throw( myruntime_error("RndMSAGenerator: Init: Not enough memory."));
    trmv_.Zero();
    remv_.Zero();
    trco_.Zero(); trco_.Add( 0.5 ); trco_.AddToDiag( 0.5 );
    reco_.Zero(); reco_.Add( 0.5 ); reco_.AddToDiag( 0.5 );
    mnrv_tr_->SetMu( &trmv_ );
    mnrv_tr_->SetSigma( &trco_ );
    mnrv_em_->SetMu( &remv_ );
    mnrv_em_->SetSigma( &reco_ );
}

// -------------------------------------------------------------------------
// Generate: generate random profiles
//
void RndMSAGenerator::Generate()
{
    mystring    outfname, name;
    char        namebuf[BUF_MAX];
    int n;

    try {
        if(! gen_profile_.empty())
            ReadProfile();
        VerifyProbabilities();
        SetRandomVariables();

        if( 0 < GetGenNumber()) {
            sprintf( namebuf, "Generating %d MSA(s)...", GetGenNumber());
            message( namebuf );
        }
        if( GetGenModel() && 0 < gen_pssm_.GetColumns() && 
            gen_pssm_.GetColumns() < GetGenLength()) {
            sprintf( namebuf, "Sampling sequences of a match state lenth of %d", gen_pssm_.GetColumns());
            warning( namebuf );
        }

        NewMSA();
        NewSStructures();
        for( n = 0; n < GetGenNumber(); n++ ) {
            AddNoiseToProbs();
            sprintf( namebuf, "_%05d", n );
            outfname = mystring( GetOutputDir()) + DIRSEP + GetNamePattern() + namebuf + slg_MSAEXT;
            name = mystring( GetNamePattern()) + namebuf;
            GenerateMSA( outfname.c_str(), name.c_str());
        }

    } catch( myexception const& ex )
    {
        throw myruntime_error( ex.what(), ex.eclass());
    }

    message( "Finished." );
}

// -------------------------------------------------------------------------
// ReadProfile: read profile
//
void RndMSAGenerator::ReadProfile()
{
    Serializer  serializer;
    serializer.ReadProfile( gen_profile_.c_str(), gen_freq_, gen_pssm_, gen_gaps_ );
    gen_freq_.CheckForAllZeros();//to ensure valid scores
    //keep original profile as well
    serializer.ReadProfile( gen_profile_.c_str(), org_freq_, org_pssm_, org_gaps_ );
    org_freq_.CheckForAllZeros();//to ensure valid scores
    SetGenModel( true );
}

// -------------------------------------------------------------------------
// VerifyProbabilities: verify probabilities
//
void RndMSAGenerator::VerifyProbabilities()
{
    double*     trprobs = GetGenTrProbs();
    double*     reprobs = GetGenReProbs();
    mystring    preamb = "RndMSAGenerator: VerifyProbabilities: ";
    int m;

    if( GetGenModel()) {
        if( gen_pssm_.GetColumns() < 1 || gen_pssm_.GetColumns() != gen_gaps_.GetColumns() ||
            org_pssm_.GetColumns() != gen_pssm_.GetColumns())
            throw myruntime_error( preamb + "Invalid profile.");

        for( m = -1; m < gen_pssm_.GetColumns(); m++ ) {
            if( gen_gaps_.GetOrgTrProbsAt(m) == NULL || org_gaps_.GetOrgTrProbsAt(m) == NULL )
                throw myruntime_error( preamb + "Null transition probabilities at the model's position.");
            if( 0 <= m &&( gen_pssm_.GetTrgFreqsAt(m) == NULL || org_pssm_.GetTrgFreqsAt(m) == NULL ))
                throw myruntime_error( preamb + "Null target probabilities at the model's position.");

            VerifyTransitionProbabilities((double*)*gen_gaps_.GetOrgTrProbsAt(m), 1.e-3, true );
            VerifyTransitionProbabilities((double*)*org_gaps_.GetOrgTrProbsAt(m), 1.e-3, true );
            if( 0 <= m ) {
                VerifyEmissionProbabilities((double*)*gen_pssm_.GetTrgFreqsAt(m), 1.e-1, true );
                VerifyEmissionProbabilities((double*)*org_pssm_.GetTrgFreqsAt(m), 1.e-1, true );
            }
        }
    }
    else {
        VerifyTransitionProbabilities( trprobs );
        VerifyTransitionProbabilities( org_trprobs_ );
        VerifyEmissionProbabilities( reprobs );
        VerifyEmissionProbabilities( org_reprobs_ );
    }
}

// VerifyTransitionProbabilities: verify transition probabilities
//
void RndMSAGenerator::VerifyTransitionProbabilities( double* trprobs, const double acc, bool adjust )
{
    const int   maxtr = gTPTRANS_NTPS;
    char        infobuf[KBYTE];
    int         n, st, states;
    double      sum;
    if( trprobs == NULL )
        throw myruntime_error("RndMSAGenerator: VerifyTransitionProbabilities: Null probabilities.");
    if( 0.1 < acc )
        throw myruntime_error("RndMSAGenerator: VerifyTransitionProbabilities: Invalid accuracy threshold.");
    for( n = 0, states = 0; states < P_NSTATES; states += maxtr, n++ ) {
        sum = 0.0;
        for( st = states; st < states + maxtr && st < P_NSTATES; st++ )
            sum += trprobs[st];
        if( sum < 1.0 - acc || sum > 1.0 + acc ) {
            sprintf( infobuf, "RndMSAGenerator: VerifyTransitionProbabilities: "
                    "Invalid transition probabilities: state %d: %.6g.", states, sum );
            throw myruntime_error( infobuf );
        }
        if( adjust && sum < 1.0 )
            trprobs[st-1] += 1.0 - sum;
    }
}

// VerifyTransitionProbabilities: verify emission probabilities
//
void RndMSAGenerator::VerifyEmissionProbabilities( double* reprobs, const double acc, bool adjust )
{
    char        infobuf[KBYTE];
    int         n;
    double      sum;
    if( reprobs == NULL )
        throw myruntime_error("RndMSAGenerator: VerifyEmissionProbabilities: Null probabilities.");
    if( 0.1 < acc )
        throw myruntime_error("RndMSAGenerator: VerifyEmissionProbabilities: Invalid accuracy threshold.");
    sum = 0.0;
    for( n = 0; n < NUMAA; n++ )
        sum += reprobs[n];
    if( sum < 1.0 - acc || sum > 1.0 + acc ) {
        sprintf( infobuf, "RndMSAGenerator: VerifyEmissionProbabilities: "
                "Invalid residue probabilities: %.6g.", sum );
        throw myruntime_error( infobuf );
    }
    if( adjust && sum < 1.0 )
        reprobs[n-1] += 1.0 - sum;
}

// -------------------------------------------------------------------------
// AddNoiseToProbs: add noise to probabilities
//
void RndMSAGenerator::AddNoiseToProbs()
{
    double*     trprobs = GetGenTrProbs();
    double*     reprobs = GetGenReProbs();
    mystring    preamb = "RndMSAGenerator: AddNoiseToProbs: ";
    int m;

    if( GetGenNoise() <= 0.0 )
        return;

    if( GetGenModel()) {
        if( gen_pssm_.GetColumns() < 1 || gen_pssm_.GetColumns() != gen_gaps_.GetColumns() ||
            org_pssm_.GetColumns() != gen_pssm_.GetColumns())
            throw myruntime_error( preamb + "Invalid profile.");

        for( m = -1; m < gen_pssm_.GetColumns(); m++ ) {
            if( gen_gaps_.GetOrgTrProbsAt(m) == NULL || org_gaps_.GetOrgTrProbsAt(m) == NULL )
                throw myruntime_error( preamb + "Null transition probabilities at the model's position.");
            if( 0 <= m &&( gen_pssm_.GetTrgFreqsAt(m) == NULL || org_pssm_.GetTrgFreqsAt(m) == NULL))
                throw myruntime_error( preamb + "Null target probabilities at the model's position.");

            AddNoiseToTransitionProbs((double*)*gen_gaps_.GetOrgTrProbsAt(m), *org_gaps_.GetOrgTrProbsAt(m));
            if( 0 <= m )
                AddNoiseToEmissionProbs((double*)*gen_pssm_.GetTrgFreqsAt(m), *org_pssm_.GetTrgFreqsAt(m));
        }
    }
    else {
        AddNoiseToTransitionProbs( trprobs, org_trprobs_ );
        AddNoiseToEmissionProbs( reprobs, org_reprobs_ );
    }
}

// AddNoiseToTransitionProbs: add noise to transition probabilities
//
void RndMSAGenerator::AddNoiseToTransitionProbs( double* trprobs, const double* source )
{
    mystring    preamb = "RndMSAGenerator: AddNoiseToTransitionProbs: ";
    const double nw = GetGenNoise();//noise weight
    const double sw = 1.0 - nw;//signal weight
    const double acc = 1.e-6;
    const int   maxtr = gTPTRANS_NTPS;
    char        infobuf[KBYTE];
    Pslvector   rmv(maxtr-1);
    Pslvector   mlgrv(maxtr);//multivariate logistic normal random variable
    double      sum;
    int         n, st, states;
    int         err;

    if( mnrv_tr_ == NULL )
        throw myruntime_error( preamb + "Null normal random variable.");
    if( trprobs == NULL || source == NULL )
        throw myruntime_error( preamb + "Null probabilities.");
    if( nw < 0.0 || 1.0 < nw )
        throw myruntime_error( preamb + "Invalid noise level.");

    for( states = 0; states < P_NSTATES; states += maxtr )
    {
        if(( err = mnrv_tr_->Gen( &rmv )) !=0 )
            throw myruntime_error( preamb + TranslatePSLError( err ));

        Normal2LogitNormal( rmv.GetVector(), rmv.GetSize(), mlgrv.GetVector(), mlgrv.GetSize());
        sum = mlgrv.Sum();
        if( sum < 1.-acc || 1.+acc < sum )
            throw myruntime_error( preamb + "Invalid logistic normal random variable.");

        for( st = states, n = 0; st < states + maxtr && st < P_NSTATES; n++, st++ )
            trprobs[st] = sw * source[st] + nw * mlgrv.GetValueAt(n);
    }
}

// AddNoiseToEmissionProbs: add noise to emission probabilities
//
void RndMSAGenerator::AddNoiseToEmissionProbs( double* reprobs, const double* source )
{
    mystring    preamb = "RndMSAGenerator: AddNoiseToEmissionProbs: ";
    const double nw = GetGenNoise();//noise weight
    const double sw = 1.0 - nw;//signal weight
    const double acc = 1.e-6;
    char        infobuf[KBYTE];
    Pslvector   rmv(NUMAA-1);
    Pslvector   mlgrv(NUMAA);//multivariate logistic normal random variable
    double      sum;
    int         n, err;

    if( mnrv_em_ == NULL )
        throw myruntime_error( preamb + "Null normal random variable.");
    if( reprobs == NULL || source == NULL )
        throw myruntime_error( preamb + "Null probabilities.");
    if( nw < 0.0 || 1.0 < nw )
        throw myruntime_error( preamb + "Invalid noise level.");

    if(( err = mnrv_em_->Gen( &rmv )) !=0 )
        throw myruntime_error( preamb + TranslatePSLError( err ));

    Normal2LogitNormal( rmv.GetVector(), rmv.GetSize(), mlgrv.GetVector(), mlgrv.GetSize());
    sum = mlgrv.Sum();
    if( sum < 1.-acc || 1.+acc < sum )
        throw myruntime_error( preamb + "Invalid logistic normal random variable.");

    for( n = 0; n < NUMAA; n++ )
        reprobs[n] = sw * source[n] + nw * mlgrv.GetValueAt(n);
}

// -------------------------------------------------------------------------
// SetRandomVariables: set random variables
//
void RndMSAGenerator::SetRandomVariables()
{
    if( gen_reprobs_ == NULL || gen_trprobs_ == NULL  )
        throw( myruntime_error("RndMSAGenerator: SetRandomVariables: Null probabilities."));
    if( drv_em_ == NULL || drv_in_ == NULL || 
        drv_tr_m_ == NULL || drv_tr_i_ == NULL || drv_tr_d_ == NULL )
        throw( myruntime_error("RndMSAGenerator: SetRandomVariables: Null random variables."));

    const int   maxtr = gTPTRANS_NTPS;

    drv_em_->SetProbs( gen_reprobs_, NUMAA );
    drv_in_->SetProbs( gen_inprobs_, NUMAA );

    drv_tr_m_->SetProbs( gen_trprobs_, maxtr );
    drv_tr_i_->SetProbs( gen_trprobs_ + maxtr, maxtr );
    drv_tr_d_->SetProbs( gen_trprobs_ + maxtr + maxtr, maxtr );
}

// -------------------------------------------------------------------------
// GenerateMSA: generate MSA and write it to file
//
void RndMSAGenerator::GenerateMSA( const char* outfilename, const char* name )
{
    if( 0 < GetGenEffnos())
        GenerateMSAofEffnos( name );
    else
        GenerateMSAHelper( name, GetGenWidth());
    ;
    msa_->PutAlignment( outfilename );
}

// -------------------------------------------------------------------------
// GenerateMSAHelper: generate MSA of the given eff. number of observations
//
void RndMSAGenerator::GenerateMSAofEffnos( const char* name )
{
    if( msa_ == NULL )
        throw myruntime_error("RndMSAGenerator: GenerateMSAofEffnos: Null MSA object.");

    const int effns = GetGenEffnos();
    if( effns < 1 )
        return;

    int     lw = 1;
    int     rw = 2000;
    int     frw = 100;//floating right bound
    int     mw = effns;
    char    msgbuf[BUF_MAX];
    mystring msg = name;
    int     peff;

    for( mw = effns; lw <= rw; mw=(lw+frw)>>1 )
    {
        GenerateMSAHelper( name, mw );
        //msa_->ConstructProfile();
        CalculateMSAEffnos();
        peff = (int)msa_->GetEffNoSequences();

        if( peff < effns ) {
            lw = mw + 1;
            if( frw <= lw )
                frw = rw;
        }
        else if( effns < peff ) {
            rw = mw - 1;
            frw = mw - 1;
        } else
            break;
    }
    if( peff != effns ) {
        sprintf( msgbuf, ": MSA of %d eff. number of sequences.", peff );
        msg += msgbuf;
        warning( msg.c_str());
    }
}

// -------------------------------------------------------------------------
// GenerateMSAHelper: generate MSA of the given width
//
void RndMSAGenerator::GenerateMSAHelper( const char* name, int width )
{
    if( width < 1 )
        return;

    PosDescriptionVector* sq = NULL;
    const mystring* seq = NULL;
    const mystring* stt = NULL;
    char    namebuf[BUF_MAX];
    int     n, m, p, i, len;
    int     mposlen = GetGenLength();

    if( GetGenModel() && gen_pssm_.GetColumns() < mposlen )
        mposlen = gen_pssm_.GetColumns();

    NewMSA();
    NewSStructures();

    msa_->clear();
    msa_->SetTitle( name );
    msa_->AppendTitle(" Random multiple sequence alignment");
    if( GetGenModel()) {
        msa_->AppendTitle(" (model: ");
        msa_->AppendTitle( my_basename( GetGenProfile().c_str()));
        msa_->AppendTitle(")");
    }

    for( n = 0; n < width; n++ )
        GenerateSequence();

    //translate from the interior format and fill the alignment structure
    if( seqs_ == NULL || stts_ == NULL || maxi_ == NULL )
        throw myruntime_error("RndMSAGenerator: GenerateMSAHelper: Null helper structures.");

    for( n = 0, len = 0; n < seqs_->GetSize(); n++ ) {
        seq = reinterpret_cast<const mystring*>(seqs_->GetValueAt(n));
        stt = reinterpret_cast<const mystring*>(stts_->GetValueAt(n));
        if( seq == NULL || stt == NULL || seq->length() != stt->length())
            throw myruntime_error("RndMSAGenerator: GenerateMSAHelper: Memory access error.");
        sq = new PosDescriptionVector( len? len: TIMES2(GetGenLength()));
        if( sq == NULL )
            throw myruntime_error("RndMSAGenerator: GenerateMSAHelper: Not enough memory.");
        sprintf( namebuf, "SRND_%06d", n+1 );
        sq->AppendDescription( namebuf, 0, strlen( namebuf ));
        for( p = 0, m = -1; m < mposlen; ) {
            for( i = 0; i < maxi_->GetValueAt(m+1); i++ ) {
                if( p < stt->length() && (*stt)[p] == PS_I )
                    sq->push((*seq)[p++]);
                else
                    sq->push(GAP);
            }
            m++;
            if( mposlen <= m ) {
                if( p < seq->length())
                    throw myruntime_error("RndMSAGenerator: GenerateMSAHelper: Translating failed.");
                break;
            }
            if( seq->length() <= p )
                throw myruntime_error("RndMSAGenerator: GenerateMSAHelper: Invalid sequence.");
            if((*stt)[p] != PS_M && (*stt)[p] != PS_D )
                throw myruntime_error("RndMSAGenerator: GenerateMSAHelper: Invalid state in translation.");
            sq->push((*seq)[p++]);
        }
        msa_->push( sq );
        if( !len )
            len = sq->size();
    }
}

// -------------------------------------------------------------------------
// GenerateSequence: generate a random sequence 
//
void RndMSAGenerator::GenerateSequence()
{
    if( GetGenModel()) {
        SampleSequenceFromModel();
        return;
    }

    if( drv_em_ == NULL || drv_in_ == NULL || 
        drv_tr_m_ == NULL || drv_tr_i_ == NULL || drv_tr_d_ == NULL )
        throw( myruntime_error("RndMSAGenerator: GenerateSequence: Null random variables."));
    if( seqs_ == NULL || stts_ == NULL || maxi_ == NULL )
        throw myruntime_error("RndMSAGenerator: GenerateSequence: Memory access error.");

    mystring    preamb = "RndMSAGenerator: GenerateSequence: ";
    const int   maxlen = 10000;//to prevent livelock
    mystring*   seq = NULL;//sequence sampled
    mystring*   stt = NULL;//state vector for the sequence
    int     n, m, i, st, err;
    int     pst = PS_M;
    char    sym;
    char    infobuf[KBYTE];
    int     mposlen = GetGenLength();

    seq = new mystring;
    stt = new mystring;
    if( seq == NULL || stt == NULL )
        throw myruntime_error( preamb + "Not enough memory.");
    seq->reserve( TIMES2(mposlen));
    stt->reserve( TIMES2(mposlen));
    seqs_->Push( seq );
    stts_->Push( stt );

    for( n = 0, i = 0, m = -1; m < mposlen && n < maxlen; n++ ) {
        if( pst == PS_M )
            SampleTransition( drv_tr_m_, &st );
        else if( pst == PS_I )
            SampleTransition( drv_tr_i_, &st );
        else if( pst == PS_D )
            SampleTransition( drv_tr_d_, &st );
        else {
            sprintf( infobuf, "RndMSAGenerator: GenerateSequence: Invalid state at pos. %d", n );
            throw myruntime_error( infobuf );
        }

        if( seqs_->GetSize() == 1 ) {
            //consider delete states for the first sequence as inserts
            sym = GAP;
            if( st== PS_D )
                st = PS_I;
        }

        pst = st;
        if( st == PS_I ) {
            i++;
            if( maxi_->GetValueAt(m+1) < i )
                maxi_->SetValueAt( m+1, i );
        }
        else {//if( st == PS_M || st == PS_D )
            m++;
            i = 0;
        }

        if( mposlen <= m )
            break;

        if( st == PS_M )
            SampleResidue( drv_em_, &sym );
        else if( st == PS_I && seqs_->GetSize() != 1 )
            SampleResidue( drv_in_, &sym );
        else if( st == PS_D )
            sym = GAP;

        stt->append((char)st);
        seq->append(sym);
    }
}

// -------------------------------------------------------------------------
// SampleSequenceFromModel: sample sequence from the model
//
void RndMSAGenerator::SampleSequenceFromModel()
{
    if( drv_em_ == NULL || drv_in_ == NULL || 
        drv_tr_m_ == NULL || drv_tr_i_ == NULL || drv_tr_d_ == NULL )
        throw( myruntime_error("RndMSAGenerator: SampleSequenceFromModel: Null random variables."));
    if( seqs_ == NULL || stts_ == NULL || maxi_ == NULL )
        throw myruntime_error("RndMSAGenerator: SampleSequenceFromModel: Memory access error.");
    if( !GetGenModel() || gen_pssm_.GetColumns() < 1 )
        throw myruntime_error("RndMSAGenerator: SampleSequenceFromModel: No model.");
    if( gen_pssm_.GetColumns() != gen_gaps_.GetColumns())
        throw myruntime_error("RndMSAGenerator: SampleSequenceFromModel: Invalid profile.");

    mystring    preamb = "RndMSAGenerator: SampleSequenceFromModel: ";
    const int   maxtr = gTPTRANS_NTPS;
    const int   maxlen = 10000;//to prevent livelock
    mystring*   seq = NULL;//sequence sampled
    mystring*   stt = NULL;//state vector for the sequence
    int     trs = 0;
    int     n, m, i, st, err;
    int     pst = PS_M;
    char    sym;
    char    infobuf[KBYTE];
    int     mposlen = GetGenLength();
    RVDisc* drv_tr = NULL;

    if( gen_pssm_.GetColumns() < mposlen )
        mposlen = gen_pssm_.GetColumns();
    seq = new mystring;
    stt = new mystring;
    if( seq == NULL || stt == NULL )
        throw myruntime_error( preamb + "Not enough memory.");
    seq->reserve( TIMES2(mposlen));
    stt->reserve( TIMES2(mposlen));
    seqs_->Push( seq );
    stts_->Push( stt );

    for( n = 0, i = 0, m = -1; m < mposlen && n < maxlen; n++ ) {
        if( gen_gaps_.GetOrgTrProbsAt(m) == NULL )
            throw myruntime_error( preamb + "Null transition probabilities at the model's position.");

        if( pst == PS_M ) {
            drv_tr = drv_tr_m_; trs = 0;
        } else if( pst == PS_I ) {
            drv_tr = drv_tr_i_; trs = maxtr;
        } else if( pst == PS_D ) {
            drv_tr = drv_tr_d_; trs = TIMES2(maxtr);
        } else {
            sprintf( infobuf, "RndMSAGenerator: SampleSequenceFromModel: Invalid state at pos. %d", m );
            throw myruntime_error( infobuf );
        }
        drv_tr->SetProbs( *gen_gaps_.GetOrgTrProbsAt(m)+trs, maxtr );
        SampleTransition( drv_tr, &st );

        if( seqs_->GetSize() == 1 ) {
            //consider delete states for the first sequence as inserts
            sym = GAP;
            if( st== PS_D )
                st = PS_I;
        }

        pst = st;
        if( st == PS_I ) {
            i++;
            if( maxi_->GetValueAt(m+1) < i )
                maxi_->SetValueAt( m+1, i );
        }
        else {//if( st == PS_M || st == PS_D )
            m++;
            i = 0;
        }

        if( mposlen <= m )
            break;

        if( st == PS_M ) {
            if( gen_pssm_.GetTrgFreqsAt(m) == NULL )
                throw myruntime_error( preamb + "Null target probabilities at the model's position.");
            drv_em_->SetProbs( *gen_pssm_.GetTrgFreqsAt(m), NUMAA );
            SampleResidue( drv_em_, &sym );
        }
        else if( st == PS_I && seqs_->GetSize() != 1 )
            SampleResidue( drv_in_, &sym );
        else if( st == PS_D )
            sym = GAP;

        stt->append((char)st);
        seq->append(sym);
    }
}

// -------------------------------------------------------------------------
// SampleResidue: sample a transition
//
void RndMSAGenerator::SampleTransition( RVDisc* rvd, int* st )
{
    mystring    preamb = "RndMSAGenerator: SampleTransition: ";
    int err, ind;
    if( rvd == NULL || st == NULL )
        throw myruntime_error( preamb + "Memory access error.");
    if(( err = rvd->GenI(&ind)) != PSL_OK )
        throw myruntime_error( preamb + TranslatePSLError( err ));
    if( ind < 0 || PS_NSTATES <= ind || rvd->GetProbs()[ind] <= 0.0 )
        throw myruntime_error( preamb + "Invalid transition.");
    *st = ind;
}

// -------------------------------------------------------------------------
// SampleResidue: sample a single residue
//
void RndMSAGenerator::SampleResidue( RVDisc* rvd, char* res )
{
    mystring    preamb = "RndMSAGenerator: SampleResidue: ";
    int err, ind;
    if( rvd == NULL || res == NULL )
        throw myruntime_error( preamb + "Memory access error.");
    if(( err = rvd->GenI(&ind)) != PSL_OK )
        throw myruntime_error( preamb + TranslatePSLError( err ));
    if( ind < 0 || NUMAA <= ind || rvd->GetProbs()[ind] <= 0.0 )
        throw myruntime_error( preamb + "Invalid emission.");
    *res = (char)ind;
}
