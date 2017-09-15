/***************************************************************************
 *   Copyright (C) 2010 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __rvdisc__
#define __rvdisc__

#include "ext/pslcodes.h"
#include "ext/pslvector.h"
#include "ext/ivector.h"
#include "ext/rng.h"
#include "rv.h"

// -------------------------------------------------------------------------
// General discrete random variable
//
class RVDisc: public RVar
{
public:
    enum TRVDsc {
        TRVDsc_Alias,
        TRVDsc_Inverse
    };
    RVDisc( Rng& rng, TRVDsc met = TRVDsc_Alias );
    virtual ~RVDisc();

    //this method should not be called for discrete r.v.
    virtual int     Gen( double* res ) { return PSL_ERR_NOPROG; }
    int             GenI( int* ind );//MAIN method to get index associated with probability

    TRVDsc  GetMet() const { return met_; }
    void    SetMet( TRVDsc value ) { met_ = value; }

    const double*   GetProbs() const { return probs_; }
    void            SetProbs( const double* ps, int sz ) { probs_ = ps; psize_ = sz; SetChanged(); }

    int             GetPSize() const { return psize_; }
    bool            GetChanged() const { return changed_; }

protected:
    int             GenIInvLin( int* ind );

    int             GenIAlias( int* ind );
    int             MakeAliasTables();

    void            SetChanged() { changed_ = true; }
    void            ResetChanged() { changed_ = false; }

private:
    TRVDsc          met_;
    const double*   probs_;//probabilities
    int             psize_;//length of prob. vector
    bool            changed_;//data has changed
    //alias method's data
    Pslvector       q_;//alias probability vector
    Ivector         j_;//alias table index
};

// -------------------------------------------------------------------------
//
inline
int RVDisc::GenI( int* ind )
{
    if( ind == NULL )
        return PSL_ERR_ADDRESS;
    int err;

    switch( GetMet()) {
      case TRVDsc_Alias: err = GenIAlias( ind ); break;
      case TRVDsc_Inverse: err = GenIInvLin( ind ); break;
      default:  return PSL_ERR_INVALID;
    };

    if( err != PSL_OK )
        return err;
    return err;
}

#endif//__rvdisc__
