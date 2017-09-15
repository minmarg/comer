/***************************************************************************
 *   Copyright (C) 2012 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __HybridAlignment__
#define __HybridAlignment__

#include "debug.h"
#include "compdef.h"
#include "rc.h"
#include "mystring.h"
#include "myexcept.h"
#include "libpro/srcaln/ProfileAlignment.h"


// _________________________________________________________________________
// Class HybridAlignment
//  implementation of semi-probabilistic alignment
//
class HybridAlignment: public ProfileAlignment
{
public:
    HybridAlignment(
            const FrequencyMatrix& freq_fst, const LogOddsMatrix& logo_fst, const GapScheme& gaps_fst,
            const FrequencyMatrix& freq_sec, const LogOddsMatrix& logo_sec, const GapScheme& gaps_sec,
            const AbstractScoreMatrix*  usc_system,
            bool ungapped = false
    );

    virtual ~HybridAlignment();
    virtual void    Run();

protected:
    explicit HybridAlignment();
    virtual void    ClearF();
    virtual void    AlignProfiles();
    virtual void    PostProcess();
    virtual void    SetFinalScore( double value );
    virtual void    MakeAlignmentPath();

    static void     GetMax( TPAScore* zss, int noes, TPAScore* maxv, int* indx );

private:

};


////////////////////////////////////////////////////////////////////////////
// INLINES
//
// SetFinalScore: save final alignment score
// 
inline
void HybridAlignment::SetFinalScore( double value )
{
    finscore_ = value;
}



#endif//__HybridAlignment__
