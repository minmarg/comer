/***************************************************************************
 *   Copyright (C) 2014 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __RndMSAGenerator__
#define __RndMSAGenerator__

#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "rc.h"
#include "mystring.h"
#include "myexcept.h"
#include "ext/rng.h"
#include "ext/ivector.h"
#include "ext/pslvector.h"
#include "ext/spdmatrix.h"
#include "ext/rv/rvdisc.h"
#include "ext/rv/rmvnorm.h"
#include "liblib/BinarySearchStructure.h"
#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/GapScheme.h"
#include "libpro/srcpro/DistributionMatrix.h"
#include "libpro/srcaln/InputMultipleAlignment.h"


// _________________________________________________________________________
// Class RndMSAGenerator
//

class RndMSAGenerator {
public:
    //
    RndMSAGenerator( const char* outdir, const char* namepat );
    ~RndMSAGenerator();

    void        Generate(); //generate random profiles

    const MOptions* GetOptions() const { return options_; }
    void        SetOptions( MOptions* value ) { options_ = value; }

    int         GetGenLength() const { return gen_lenmsa_; }
    void        SetGenLength( int value ) { gen_lenmsa_ = value; }

    int         GetGenNumber() const { return gen_nomsas_; }
    void        SetGenNumber( int value ) { gen_nomsas_ = value; }

    int         GetGenWidth() const { return gen_msawidth_; }
    void        SetGenWidth( int value ) { gen_msawidth_ = value; }

    int         GetGenEffnos() const { return gen_msaeffnos_; }
    void        SetGenEffnos( int value ) { gen_msaeffnos_ = value; }

    const double* GetGenTrProbs() const { return gen_trprobs_; }
    void        SetGenTrProbs( double* value );

    const double* GetGenReProbs() const { return gen_reprobs_; }
    void        SetGenReProbs( double* value );

    mystring    GetGenProfile() const { return gen_profile_; }
    void        SetGenProfile( mystring value ) { gen_profile_ = value; }

    double      GetGenNoise() const { return gen_noise_; }
    void        SetGenNoise( double value ) { gen_noise_ = value; }

protected:
    explicit RndMSAGenerator();

    void    Init();//initialization

    void    ReadProfile();

    double* GetGenTrProbs() { return gen_trprobs_; }
    double* GetGenReProbs() { return gen_reprobs_; }

    bool    GetGenModel() const { return gen_model_; }
    void    SetGenModel( bool value ) { gen_model_ = value; }

    void    VerifyProbabilities();
    void    VerifyTransitionProbabilities( double*, const double acc = 1.e-6, bool adjust = false );
    void    VerifyEmissionProbabilities( double*, const double acc = 1.e-6, bool adjust = false );
    void    AddNoiseToProbs();
    void    AddNoiseToTransitionProbs( double*, const double* );
    void    AddNoiseToEmissionProbs( double*, const double* );
    void    SetRandomVariables();
    void    GenerateMSA( const char* fname, const char* name );
    void    GenerateMSAofEffnos( const char* name );
    void    GenerateMSAHelper( const char* name, int width );
    void    GenerateSequence();
    void    SampleSequenceFromModel();

    void    SampleTransition( RVDisc*, int* st );
    void    SampleResidue( RVDisc*, char* res );

    void    CalculateMSAEffnos();

    const char* GetOutputDir() const { return outdir_; }
    const char* GetNamePattern() const { return nampat_; }

private:
    void    DestroyMSA();
    void    NewMSA();

    void    DestroySStructures();
    void    NewSStructures();

private:
    const char*         outdir_;            //output directory name
    const char*         nampat_;            //filename pattern for new generated MSAs

    int                 gen_lenmsa_;        //first sequence's length
    int                 gen_nomsas_;        //number of MSAs to generate
    int                 gen_msawidth_;      //number of sequences within MSA
    int                 gen_msaeffnos_;     //effective number of sequences of MSA
    double*             gen_trprobs_;       //transition probabilities
    double*             gen_reprobs_;       //residue probabilities
    double              gen_inprobs_[NUMAA];//insertion probabilities
    double              org_trprobs_[P_NSTATES];//original transition probabilities
    double              org_reprobs_[NUMAA];    //original residue probabilities
    mystring            gen_profile_;       //model to sample sequences from
    FrequencyMatrix     gen_freq_, org_freq_; //frequency matrix
    LogOddsMatrix       gen_pssm_, org_pssm_; //PSSM matrix
    GapScheme           gen_gaps_, org_gaps_; //gap costs
    bool                gen_model_;         //sample from model
    double              gen_noise_;         //noise level to add to probabilities

    MTRng               rng_em_;//RNG for emission probabilities
    MTRng               rng_in_;//RNG for insertion probabilities
    MTRng               rng_tr_;//RNG for transition probabilities
    RVDisc*             drv_em_;//random variable for emission probabilities
    RVDisc*             drv_in_;//random variable for insertion probabilities
    RVDisc*             drv_tr_m_;//random variable for transition probabilities given M state
    RVDisc*             drv_tr_i_;//random variable for transition probabilities given I state
    RVDisc*             drv_tr_d_;//random variable for transition probabilities given D state

    RMVNorm*            mnrv_tr_;//multivariate normal rv for transition noise
    RMVNorm*            mnrv_em_;//multivariate normal rv for emission noise
    Pslvector           trmv_, remv_;
    SPDmatrix           trco_, reco_;

    InputMultipleAlignment* msa_;//alignmemt
    SimpleVector*           seqs_;//sequences in the mystring format
    SimpleVector*           stts_;//states in the mystring format
    Ivector*                maxi_;//maximums of insertion lengths at each model position
    MOptions*           options_;//options
};


////////////////////////////////////////////////////////////////////////////
// INLINES
//
inline
void RndMSAGenerator::DestroyMSA()
{
    if( msa_ )
        delete msa_;
    msa_ = NULL;
}
inline
void RndMSAGenerator::NewMSA()
{
    if( options_ == NULL )
        throw myruntime_error("RndMSAGenerator: NewMSA: Null options.");

    DestroyMSA();
    msa_ = new InputMultipleAlignment();
    if( msa_ == NULL )
        throw myruntime_error("RndMSAGenerator: NewMSA: Not enough memory.");
    msa_->SetKeepTitles( true );

    if( options_->GetHCFILTER()) {
        msa_->SetSEGWindow( options_->GetHCWINDOW());
        msa_->SetSEGLowEntropy( options_->GetHCLOWENT());
        msa_->SetSEGHighEntropy( options_->GetHCHIGHENT());
    }
    else
        msa_->SetUsingSEGFilter( false );

    if( options_->GetLCFILTEREACH()) {
        msa_->SetSeqSEGWindow( options_->GetLCWINDOW());
        msa_->SetSeqSEGLowEntropy( options_->GetLCLOWENT());
        msa_->SetSeqSEGHighEntropy( options_->GetLCHIGHENT());
    }
    else
        msa_->SetUsingSeqSEGFilter( false );

    msa_->SetIdentityLevel( options_->GetIDENTITY());
    msa_->SetComputeDELETEstates( options_->GetDELSTATE());
    msa_->SetExtentMinWindow( options_->GetMINALNPOS());
    msa_->SetExtentMinSeqPercentage( options_->GetMINALNFRN());
    msa_->SetPseudoCountWeight( options_->GetPCFWEIGHT());
}

// -------------------------------------------------------------------------
// CalculateMSAEffnos: construct profile and calculate eff. no. sequences
//
inline
void RndMSAGenerator::CalculateMSAEffnos()
{
    if( msa_ == NULL )
        throw myruntime_error("RndMSAGenerator: GetMSAEffnos: Memory access error.");

    double  avgpeseq;

    msa_->SelectSequences();
    msa_->ComputeGlobSequenceWeights( &avgpeseq );

    msa_->ComputeExtents();
    msa_->ComputeMIDstateSequenceWeights();
    msa_->ComputeTransitionFrequencies( true/*gwghts*/, false );
    msa_->CalculateEffNoSequences();
}

// -------------------------------------------------------------------------
// Destroy and create strcutures related to sampling sequences from the 
//  model
//
inline
void RndMSAGenerator::DestroySStructures()
{
    size_t  s;
    if( seqs_ ) {
        for( s = 0; s < seqs_->GetSize(); s++ )
            if( seqs_->GetValueAt(s)) {
                delete (mystring*)seqs_->GetValueAt(s);
                seqs_->SetValueAt( s, NULL );
            }
        delete seqs_;
        seqs_ = NULL;
    }
    if( stts_ ) {
        for( s = 0; s < stts_->GetSize(); s++ )
            if( stts_->GetValueAt(s)) {
                delete (mystring*)stts_->GetValueAt(s);
                stts_->SetValueAt( s, NULL );
            }
        delete stts_;
        stts_ = NULL;
    }
    if( maxi_ ) {
        delete maxi_;
        maxi_ = NULL;
    }
}
inline
void RndMSAGenerator::NewSStructures()
{
    int len = GetGenLength();
    if( GetGenModel())
        len = gen_pssm_.GetColumns();
    DestroySStructures();
    seqs_ = new SimpleVector( 1000 );
    stts_ = new SimpleVector( 1000 );
    maxi_ = new Ivector( len+1 );
    if( seqs_ == NULL || stts_ == NULL || maxi_ == NULL )
        throw myruntime_error("RndMSAGenerator: NewSStructures: Not enough memory.");
}

// -------------------------------------------------------------------------
// SetGenTrProbs: set transition probabilities
//
inline
void RndMSAGenerator::SetGenTrProbs( double* value )
{
    int n;
    if( value == NULL )
        throw myruntime_error("RndMSAGenerator: SetGenTrProbs: Null argument.");
    gen_trprobs_ = value;
    for( n = 0; n < P_NSTATES; n++ )
        org_trprobs_[n] = gen_trprobs_[n];
}

// SetGenTrProbs: set emission probabilities
//
inline
void RndMSAGenerator::SetGenReProbs( double* value )
{
    int n;
    if( value == NULL )
        throw myruntime_error("RndMSAGenerator: SetGenReProbs: Null argument.");
    gen_reprobs_ = value;
    for( n = 0; n < NUMAA; n++ )
        org_reprobs_[n] = gen_reprobs_[n];
}

// -------------------------------------------------------------------------

#endif//__RndMSAGenerator__
