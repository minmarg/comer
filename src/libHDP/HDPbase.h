/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __HDPbase__
#define __HDPbase__

#include <stdio.h>
#include <stdlib.h>

#include "ext/psl.h"
#include "ext/ivector.h"
#include "rc.h"
#include "myexcept.h"
#include "hdplib.h"
#include "HDPscores.h"
#include "HDP_SSSScores.h"
#include "iHDP_SSSScores.h"
#include "basin.h"
#include "dish.h"
#include "menu.h"
#include "rntchain.h"


//global object for application only
class HDPbase;
extern HDPbase  HDPBASE;
extern HDPbase  HDPctBASE;
void SetHDPBASE( const char* filename );
void SetHDPctBASE( const char* filename );

//define deg.of freedom nu of student's t distribution
// #define NU_t_o2 ( nu_t_o2 - do2 ) //original deg. of freedom
#define NU_t_o2 ( nu_t_o2 + do2 * GetDegFAdjustment()) //adjusted deg.of freedom

// -------------------------------------------------------------------------
// class HDPbase: base class for Hierarchical Dirichlet Process
//
class HDPbase
{
public:
    enum TResType {
        TRT_MHUpdate,//results obtained from MH update
        TRT_GibbsUpdate,//results obtained from Gibbs sampling update
        TRT_novalues
    };

public:
    HDPbase();
    virtual ~HDPbase();

    bool        GetUninfPrior() const { return uninfprior_; }
    void        SetUninfPrior( bool value ) { uninfprior_ = value; }

    bool        GetCluster4Each() const { return clst4each; }
    void        SetCluster4Each( bool value ) { clst4each = value; }

    double      GetS0ScaleFac() const { return S0scalefac_; }
    void        SetS0ScaleFac( double value ) { S0scalefac_ = value; }

    double      GetDegFAdjustment() const { return adjdegf_; }
    void        SetDegFAdjustment( double value ) { adjdegf_ = value; }


    double      GetMixWeight() const { return mixweight_; }
    void        SetMixWeight( double value ) { mixweight_ = value; }

    const HDPscores* GetScores() const {return scores_; }
    void        SetScores( const HDPscores* value );


    void        SetPriorProbs();

    int         GetReadNoGroups() const { return rdnorsts_; }
    void        SetReadNoGroups( int value ) { rdnorsts_ = value; }

    int         GetNoSupClusters() const { return nosupclsts_; }
    void        SetNoSupClusters( int value ) { nosupclsts_ = value + 1; }

    double      GetAdjWeight() const { return adjweight_; }
    void        SetAdjWeight( double value ) { adjweight_ = value; }


    void        ReserveMenu( int );
    void        ReserveBasin( int );
    int         AddToBasin( Pslvector* nv );

    void        CalcPPProbs( Pslvector& lnvar, double* bppr, Pslvector* ppprobs, Ivector* cindcs, 
                             double lpfact = 0.02, bool usepriors = true, bool tonormal = true ) const;
    void        MixCLNVar( Pslvector& lnvar, Pslvector* lnmixed ) const;

    bool        CalcPriorParams( bool addnoise = false, double std = 1.0 );
    int         CalcPriorParams( const Pslvector& stds );
    int         SetUninfPriorParamsS0I( int dim, int ctx, const Pslvector& mvec, double = -1.0, double = -1.0 );
    void        CalcPriorProbFact();
    void        RecalcMenuParams();
    void        RecalcDishParams( int k );
    void        CalcProbOfDish( int k, double* lprob );
    void        CalcProbOfData( double* lprob );

    void            ReadGroups( const char*, Ivector* dids );
    void            PrintGroups( const char*, double* = NULL, int* = NULL );

    void            ReadParameters( const char*, const Ivector* dids = NULL );
    virtual void    PrintParameters( const char* );

    void            PrintDishes( const char*, double* = NULL, int* = NULL );

    static void LogitNormalErrCorr( Pslvector* f, double tol );
    static void LogitNormal2Normal( Pslvector* f, double tol );

    const Basin*    GetBasin() const { return basin_; }
    Basin*          GetBasin() { return basin_; }

    const Menu*     GetMenu() const { return menu_; }
    Menu*           GetMenu() { return menu_; }
    void            SetMenu( Menu* value ) { DestroyMenu(); menu_ = value; }
    void            ReplaceMenu( Menu* value ) { menu_ = value; }

    const RntChain* GetChain() const { return chain_; }
    RntChain*       GetChain() { return chain_; }

    int         GetCtxtSize() const { return ctxtsize_; }
    double      GetDPMTau() const { return tau_; }
    double      GetDPMGamma() const { return gamma_; }

public:
    ///{{M-VARIATE PROBABILITIES
    //prior
    void        PriorProbVec( const Pslvector*, double* lprob ) const;
    void        PriorProbVecHlp( const Pslvector*, double* lprob ) const;
    void        PriorProbVecHlpObs( const Pslvector*, double* lprob ) const;
    //posterior predictive
    void        ProbVecOfDish( const Pslvector*, int k, double* lprob ) const;
    void        ProbVecOfDishHlp( const Pslvector*, int k, double* lprob ) const;
#ifdef PROBMTXINVSM
    void        ProbVecOfDishHlpObs( const Pslvector*, int k, double* lprob ) const;
#endif
    //exclusive posterior predictive
    void        ProbVecOfDishExc( const Pslvector*, int k, double* lprob ) const;
    void        ProbVecOfDishExcHlp( const Pslvector*, int k, double* lprob ) const;
#ifdef PROBMTXINVSM
    void        ProbVecOfDishExcHlpObs( const Pslvector*, int k, double* lprob ) const;
#endif
    ///}}

    ///{{MATRIX PROBABILITIES
    //prior matrix probabilities
    void        PriorProbMtx( const Table*, double* lprob ) const;
    void        PriorProbMtx( const Pslmatrix*, double* lprob ) const;
    void        PriorProbMtxHlp( const Table*, double* lprob ) const;
    void        PriorProbMtxHlp( const Pslmatrix*, double* lprob ) const;
    void        PriorProbMtxHlpObs( const Pslmatrix*, double* lprob ) const;
    //posterior predictive matrix probabilities
    void        ProbMtxOfDish( const Table*, int k, double* lprob ) const;
    void        ProbMtxOfDish( const Pslmatrix*, int k, double* lprob ) const;
    void        ProbMtxOfDishHlp( const Table*, int k, double* lprob ) const;
    void        ProbMtxOfDishHlp( const Pslmatrix*, int k, double* lprob ) const;
#ifdef PROBMTXINVSM
    void        ProbMtxOfDishHlpObs( const Pslmatrix*, int k, double* lprob ) const;
#endif
    //exclusive posterior predictive matrix probabilities
    void        ProbMtxOfDishExc( const Table*, int k, double* lprob ) const;
    void        ProbMtxOfDishExc( const Pslmatrix*, int k, double* lprob ) const;
    void        ProbMtxOfDishExcHlp( const Table*, int k, double* lprob ) const;
    void        ProbMtxOfDishExcHlp( const Pslmatrix*, int k, double* lprob ) const;
#ifdef PROBMTXINVSM
    void        ProbMtxOfDishExcHlpObs( const Pslmatrix*, int k, double* lprob ) const;
#endif
    ///}}

protected:
    int         GetIterationRead() const { return iterread_; }
    void        SetIterationRead( int value ) { iterread_ = value; }

    double      GetMaxLProbData() const { return mlpd_; }
    void        SetMaxLProbData( double value ) { mlpd_ = value; }

    double      GetLastLProbData() const { return lastlpd_; }
    void        SetLastLProbData( double value ) { lastlpd_ = value; }

    TResType    GetResType() const { return restype_; }
    void        SetResType( TResType value ) { restype_ = value; }

    int         GetMHUpdateNum() const { return mhupdate_; }
    void        SetMHUpdateNum( int value ) { mhupdate_ = value; }

    int         GetGibbsIt() const { return gibbsit_; }
    void        SetGibbsIt( int value ) { gibbsit_ = value; }

    bool        GetRestarted() const { return restarted_; }
    void        SetRestarted( bool value ) { restarted_ = value; }

    void        SetDPMTau( double value ) { tau_ = value; }
    void        SetDPMGamma( double value ) { gamma_ = value; }

    void            ReadParameters( FILE*, const Ivector* );
    bool            ReadSample( FILE*, Pslvector**, int dim = 0 );
    virtual void    PrintParameters( FILE* );

    void            ReadPriorParams( FILE* );
    void            PrintPriorParams( FILE* );

    void            ReadDishParams( FILE*, const Ivector* );
    void            PrintDishParams( FILE* );

    void            PrintGroups( FILE*, double* = NULL, int* = NULL );
    void            PrintDishes( FILE*, double* = NULL, int* = NULL );
    void            PrintSummary( const char*, double* = NULL, int* = NULL );
    void            PrintSummary( FILE*, double* = NULL, int* = NULL );
    void            PrintBasin( FILE* );//TEST

    void        Table2Mtx( const Table*, Pslmatrix* ) const;

    void        DestroyBasin();
    void        DestroyMenu();
    void        DestroyChain();

    void        InitBasin( int size );
    void        InitMenu( int size );
    void        InitChain( int size );

    bool        TestMenu() const;

    void        ResetTotalNoSamples()       { totalnos_ = 0; }
    void        IncTotalNoSamples()         { totalnos_++; }
    int         GetTotalNoSamples() const   { return totalnos_; }

    void        SetCtxtSize( int value ) { ctxtsize_ = value; }

    static double   GetDefDPMTau() { return s_defdpmtau_; }
    static double   GetDefDPMGamma() { return s_defdpmgamma_; }

    static int  GetDefSampleDim() { return s_defsampledim_; }
    static int  GetDefDishSize() { return s_defdishsize_; }
    static int  GetDefTableSize() { return s_deftablesize_; }
    static int  GetDefBasinSize() { return s_defbasinsize_; }

    static double   GetDefKappa_pp_a() {return s_defkappa_pp_a_; }
    static double   GetDefKappa_pp_b() {return s_defkappa_pp_b_; }
    static double   GetDefKappaParam() {return s_defkappa_; }
    static double   GetDefNuParam() { return s_defnu_; }

private:
    bool        uninfprior_; //uninformative prior
    bool        clst4each;  //assign each sample to new cluster
    Basin*      basin_; //basin of entire set of vectors
    Menu*       menu_;  //menu of dishes
    RntChain*   chain_; //chain of restaurants
    int         ctxtsize_;  //size of context
    int         rdnorsts_; //proposal (read) number of restaurants (groups)
    int         totalnos_;  //total number of samples
    TResType    restype_;  //results type
    int         mhupdate_;  //MH update number
    int         gibbsit_;  //Gibbs iteration
    bool        restarted_;  //just restarted

    //{{application parameters
    double      mixweight_; //PME weight of mixing
    int         nosupclsts_;//number of support clusters (dishes)
    double      adjweight_; //weight of score adjustment under HDP framework
    const HDPscores* scores_;//dish id-id scores
    //}}

    double      tau_;       //DPM concentration parameter
    double      gamma_;     //DPM concentration parameter

    static double   s_defdpmtau_;//default value of DPM concentration parameter tau_
    static double   s_defdpmgamma_;//default value of DPM concentration parameter gamma_

    static int  s_defsampledim_; //default dimensionality of samples
    static int  s_defdishsize_; //default dish size
    static int  s_deftablesize_;//default table size
    static int  s_defbasinsize_;//default basin size

    static double   s_defkappa_pp_a_;//default value of the prior param. a for kappa
    static double   s_defkappa_pp_b_;//default value of the prior param. b for kappa
    static double   s_defkappa_;//default value of the NIW kappa parameter (scale to matrix)
    static double   s_defnu_;//default value of the NIW nu parameter (deg. of freedom)

    int         iterread_;//iteration number read
    double      mlpd_;//maximum value of log probability of data
    double      lastlpd_;//last obtained value of log probability of data

    double      S0scalefac_;//S0 scale factor
    double      adjdegf_;//adjustment to deg. of freedom in terms of dim. over 2
};


// -------------------------------------------------------------------------
// DestroyBasin: Destroy basin of vectors
//
inline
void HDPbase::DestroyBasin()
{
    if( basin_ )
        delete basin_;
    basin_ = NULL;
}

// -------------------------------------------------------------------------
// DestroyMenu: Destroy menu of dishes
//
inline
void HDPbase::DestroyMenu()
{
    if( menu_ )
        delete menu_;
    menu_ = NULL;
}

// -------------------------------------------------------------------------
// DestroyChain: Destroy chain of restaurants
//
inline
void HDPbase::DestroyChain()
{
    if( chain_ )
        delete chain_;
    chain_ = NULL;
}

// -------------------------------------------------------------------------
// InitBasin: Initialize basin of vectors
//
inline
void HDPbase::InitBasin( int size )
{
    DestroyBasin();

    basin_ = new Basin( size );
    if( basin_ == NULL )
        throw myruntime_error( "HDPbase: Not enough memory." );
    basin_->SetDestroy( true );
}

// -------------------------------------------------------------------------
// InitMenu: Initialize menu of dishes
//
inline
void HDPbase::InitMenu( int size )
{
    DestroyMenu();

    menu_ = new Menu( size );
    if( menu_ == NULL )
        throw myruntime_error( "HDPbase: Not enough memory." );
    menu_->SetDestroy( true );
}

// -------------------------------------------------------------------------
// InitChain: Initialize chain of restaurants
//
inline
void HDPbase::InitChain( int size )
{
    DestroyChain();

    chain_ = new RntChain( size );
    if( chain_ == NULL )
        throw myruntime_error( "HDPbase: Not enough memory." );
    chain_->SetDestroy( true );
}

// -------------------------------------------------------------------------
// ReserveBasin: reserve basin size
//
inline
void HDPbase::ReserveMenu( int size )
{
    InitMenu( size );
}

// -------------------------------------------------------------------------
// ReserveBasin: reserve basin size
//
inline
void HDPbase::ReserveBasin( int size )
{
    InitBasin( size );
}

// =========================================================================
// AddToBasin: add normal multivariate variable to basin
//
inline
int HDPbase::AddToBasin( Pslvector* nv )
{
    if( GetBasin() == NULL )
        InitBasin( GetDefBasinSize());
    return GetBasin()->NewValue( nv );
}

// =========================================================================
// AddToBasin: add normal multivariate variable to basin
//
inline
void HDPbase::SetScores( const HDPscores* value )
{
    scores_ = NULL;
    if( value == NULL )
        return;
    if( !GetMenu() || GetMenu()->GetSize() != value->GetCardinality())
        throw myruntime_error("HDPbase: SetScores: Inconsistent score table.");
    scores_ = value;
}


#endif//__HDPbase__
