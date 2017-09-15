/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


// #include <time.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "rc.h"
#include "mystring.h"
#include "myexcept.h"
#include "liblib/logitnormal.h"
#include "ext/psl.h"
#include "ext/rng.h"
#include "sort/sort.h"
#include "libpro/srcsco/ProfileMatrix.h"
#include "libseg/SEGSequence.h"
#include "InputMultipleAlignment.h"


//global file variables
//
const bool      gbTOLXs = false;//include X into calculation of sequence weights
const double    gdWCPC = 0.5;//constant power coefficient in weight calculation
const size_t    gszMINEXTINTRV = 10;//minimum extent interval in calculation of seq. weights
const double    gszMINPRCSEQNS = 0.5;//minimum percentage of sequences in extent at a position
const bool      gbUSEXTENTS = true;//use extents in calculation of seq. weights


////////////////////////////////////////////////////////////////////////////
// CLASS InputMultipleAlignment
//
// Initialization constructor
//

InputMultipleAlignment::InputMultipleAlignment( unsigned reservation )
:   alignmentMatrix( NULL ),
    queryDescription( NULL ),
    length_( 0 ),
    capacity_( 0 ),
    identity_level( IDENTITY_LEVEL ),
    effnoseqs_( 0 ),
    name( NULL ),
    titletext( NULL ),
    keeptitles_( false ),
    ignoregapsinquery_( false/*true*/ ),
    deletestateson_( true/*false*/ ),

    tfrmix_( tfrmixNo ),
    scoadj_( scoadjNo ),
    HDPbase_( NULL ),
    HDPctbase_( NULL ),

    usingsegfilt_( true ),
    segfiltwinlenval_( MAC_SEGSEQ_WIN_LENGTH ),
    segfiltlowentval_( MAC_SEGSEQ_LOW_ENTROPY ),
    segfilthighentval_( MAC_SEGSEQ_HIGH_ENTROPY ),

    usingseqseg_( false ),
    seqsegwinlenval_( DEFAULT_SEGSEQ_WIN_LENGTH ),
    seqseglowentval_( DEFAULT_SEGSEQ_LOW_ENTROPY ),
    seqseghighentval_( DEFAULT_SEGSEQ_HIGH_ENTROPY ),

    extminwindow_( MAC_EXTENT_MINWIN ),
    extminseqperc_( MAC_EXTENT_SEQ_COVER ),
    pseudocntweight_( MAC_WEIGHT_PSEUDO_COUNTS )
{
    Realloc( reservation );
}

// -------------------------------------------------------------------------
// Destructor
// -------------------------------------------------------------------------

InputMultipleAlignment::~InputMultipleAlignment()
{
    DeleteExtendedDescriptionVector( queryDescription );

    clear();

    if( alignmentMatrix )
        free( alignmentMatrix );
}

// -------------------------------------------------------------------------
// PlainPreprocess: plainly preprocess input multiple sequence alignment
//
void InputMultipleAlignment::PlainPreprocess()
{
    InitQueryDescription();
    SetStates();
}

// -------------------------------------------------------------------------
// SelectSequences: select non-redundant set of sequences from multiple 
//  alignment
//
void InputMultipleAlignment::SelectSequences()
{
    PreprocessAlignment();
    if( !size())
        throw myruntime_error( mystring( "No sequences: Unable to construct profile." ));

    PurgeAtSequenceIdentity();

    InitQueryDescription();
    SetStates();

    if( GetUsingSeqSEGFilter())
        FilterSequencesWithSEG();

    if( GetUsingSEGFilter())
        RefineWithSEG();

// //{{TEST
// if( GetUsingSEGFilter() || GetUsingSeqSEGFilter())
//     PurgeAtSequenceIdentity();
// //}}

    SetBackgroundProbabilities();
    SetPosNoSequences();
}

// -------------------------------------------------------------------------
// ConstructProfile: main procedure describing steps of profile construction
//
void InputMultipleAlignment::ConstructProfile()
{
    double  avgpeseq;

    SelectSequences();

    ComputeGlobSequenceWeights( &avgpeseq );

#if 1
    if( gbUSEXTENTS )
        ComputeExtents();

// // //     ComputeSequenceWeights();
// //     ComputeMstateSequenceWeights();
// //     ComputeIstateSequenceWeights();
// //     ComputeDstateSequenceWeights();
//     ComputeMIDstateSequenceWeights();
//     ComputeTransitionFrequencies( true );
//     DeriveExpectedNoObservations();

    if( gbUSEXTENTS )
          ComputeMIDstateSequenceWeights();
    else  ComputeMIDstateSequenceWeightsNoExtents();
    ComputeTransitionFrequencies( true/*gwghts*/, false );
#else
    ComputeGWMIDstateSequenceWeights( avgpeseq );
    ComputeTransitionFrequencies( true/*gwghts*/, false );
#endif

    CalculateEffNoSequences();

    AdjustWeights();
    ComputeTargetTransFrequencies();

//     ComputeTargetFrequencies();
    ComputeTargetFrequenciesMDLVar();
    //do not mix target frequencies here; 
    //they will be mixed, if applied, for query just before searching
    if( GetScoAdjmentHDPCtx())
        //NOTE:calculate posterior predictives after profile scaling
        ;//CalcTFPosteriorPredictives();
    else if( GetTarFrMixHDPCtx())
        MixTargetFrequenciesHDPCtx();
    RecalcBackgroundProbabilities();
    CalcPosteriorProbabilities();
    ComputePSSM();
}

// -------------------------------------------------------------------------
// InitQueryDescription: Initialize query description structure
//
void InputMultipleAlignment::InitQueryDescription()
{
    PosDescriptionVector*   sequence = SequenceAt( 0 );
    if( !sequence )
        throw myruntime_error("InitQueryDescription: Null sequence." );

    DeleteExtendedDescriptionVector( queryDescription );
    queryDescription = NewExtendedDescriptionVector( *sequence );
}

// -------------------------------------------------------------------------
// Realloc: reallocate memory
//
void InputMultipleAlignment::Realloc( int newcap )
{
    if( capacity_ == 0 ) {
        alignmentMatrix = ( PosDescriptionVector** )malloc( sizeof( void* ) * newcap );
    } else {
        alignmentMatrix = ( PosDescriptionVector** )realloc( alignmentMatrix, sizeof( void* ) * newcap );
    }

    if( !alignmentMatrix )
        throw myruntime_error( mystring( "Not enough memory." ));

    PosDescriptionVector**  taligns = alignmentMatrix;

    if( capacity_ != 0 ) {
        taligns = alignmentMatrix + capacity_;
    }

    memset( taligns, 0, sizeof( void* ) * ( newcap - capacity_ ));

    capacity_ = newcap;
}

// -------------------------------------------------------------------------
// push: push sequence into the multiple alignment matrix
// -------------------------------------------------------------------------

void InputMultipleAlignment::push( PosDescriptionVector* seq )
{
    if( capacity_ < length_ + 1 ) {
        Realloc(  capacity_ * 2 );
    }

    alignmentMatrix[length_] = seq;

    length_++;
}

// -------------------------------------------------------------------------
// push: clears all the sequences in the multiple alignment matrix
// -------------------------------------------------------------------------

void InputMultipleAlignment::clear()
{
    if( name )      { free( name );      name = NULL; }
    if( titletext ) { free( titletext ); titletext = NULL; }

    if( alignmentMatrix == NULL )
        return;

    for( size_t n = 0; n < length_; n++ )
        DeletePositionVector( alignmentMatrix[n] );

    memset( alignmentMatrix, 0, sizeof( void* ) * capacity_ );

    length_ = 0;
    effnoseqs_ = 0;
}

// -------------------------------------------------------------------------
// PreprocessAlignment: makes terminal gaps of the sequences unused
// -------------------------------------------------------------------------

void InputMultipleAlignment::PreprocessAlignment()
{
    PosDescriptionVector*   one = NULL;
    size_t  i, p;

//     for( Reset(); !Eof(); Inc() ){
//         one = Sequence();
// 
//         if( !one || !one->GetUsed())
//             continue;
// 
//         for( one->Reset(); !one->Eos() && one->Residue() == GAP; one->Inc() )
//             one->SetUnused();
// 
//         for( one->BackReset(); one->Geb() && one->Residue() == GAP; one->Dec() )
//             one->SetUnused();
//     }

    for( i = 0; i < size(); i++ )
    {
        one = SequenceAt( i );
        if( !one || !one->GetUsed())
            continue;

        for( p = 0; p < one->size() && !IsValidResSym( one->ResidueAt( p )); p++ )
            if( i )
                one->SetUnusedAt( p );

        if( p < one->size() && IsValidResSym( one->ResidueAt( p )))
            one->SetFirstUsed( p );

        for( p = one->size(); 1 <= p && !IsValidResSym( one->ResidueAt( p - 1 )); p-- )
            if( i )
                one->SetUnusedAt( p - 1 );

        if( 1 <= p && IsValidResSym( one->ResidueAt( p - 1 )))
            one->SetLastUsed( p - 1 );
    }
}

// -------------------------------------------------------------------------
// LenCompare: compare lenths of sequences
static inline
int LenCompare( const void* lenvec, size_t n1, size_t n2 )
{
    const size_t* vec = ( const size_t* )lenvec;
    const size_t  nn1 = TIMES2(n1);
    const size_t  nn2 = TIMES2(n2);
    if( vec == NULL )
        return 0;
    //[n2]-[n1] -- to sort in descending order
//     return vec[nn2] - vec[nn1];
    return ( int )(( vec[nn2] == vec[nn1])? vec[nn1+1] - vec[nn2+1]: vec[nn2] - vec[nn1] );
}

// -------------------------------------------------------------------------
// IdlCompare: compare identity levels of sequences
static inline
int IdlCompare( const void* key1, const void* key2 )
{   //key2-key1 -- to sort in descending order
    return ( int )(( ssize_t )key1 - ( ssize_t )key2 );
}

// -------------------------------------------------------------------------
// PurgeAtSequenceIdentity: removes sequences that share a certain 
//  percentage of sequence identity
//
void InputMultipleAlignment::PurgeAtSequenceIdentityObs()
{
    const bool   bPASS1 = true; 
    const bool   bPASS2 = false; 
    const double dIDTHR2 = 0.9; //identity threshold in 2nd pass
    const size_t sMINSEG = 5;  //minimum segment length in 2nd pass
    bool    simthr;     //similar under threshold
    size_t  segment;    //length of aligned segment
    size_t  matches;    //number of residue matches in the sequences
    size_t  n, q;
    size_t  i, j, e, p, ith, jth;
    bool            u1, u2;
    unsigned char   r1, r2;
    PosDescriptionVector*   fst, *sec;

    if( !size())
        throw myruntime_error( "InputMultipleAlignment: No sequences." );

    PosDescriptionVector*   query = SequenceAt( 0 );

    if( query == NULL )
        throw myruntime_error( mystring( "InputMultipleAlignment: No query sequence obtained." ));

    const size_t effsize = query->GetEffectiveSize();

    if( effsize == 0 )
        throw myruntime_error( mystring( "InputMultipleAlignment: Query sequence contains no residues." ));

    size_t  locindices[effsize];

    size_t  lens[size()];
    size_t  lindsrt[size()];

    //{{sort lengths of sequences
    for( i = 0; i < size(); i++ ) {
        PosDescriptionVector* sqn = SequenceAt( i );
        if( sqn )
              lens[i] = sqn->GetEffectiveSize();
        else  lens[i] = 0;
    }
    lindsrt[0] = 0;
    //sort with reservation of the first place for query
    HeapSortInd( lindsrt+1, lens+1, size()-1, LenCompare );
    for( i = 1; i < size(); i++ )
        lindsrt[i]++;
    //}}

    // collect indices of match positions in query
    for( q = 0, n = 0; q < query->size(); q++ ) {
        if( query->ResidueAt( q ) == GAP )
            continue;
        if( effsize <= n ) {
            n = 0;
            break;
        }
        locindices[n++] = q;
    }

    if( effsize != n )
        throw myruntime_error( "InputMultipleAlignment: Wrong query effective size." );

    if( bPASS1 ) {
        // 1st pass...
        for( i = 0; i + 1 < size(); i++ )
        {
            fst = SequenceAt( ith=lindsrt[i] );
            if( fst == NULL || !fst->GetUsed())
                continue;

            for( j = i + 1; j < size(); j++ )
            {
                sec = SequenceAt( jth=lindsrt[j] );
                if( sec == NULL || !sec->GetUsed())
                    continue;

                segment = 0;
                matches = 0;

                for( e = 0, p = 0; e < effsize; e++ ) {
                    // iterate over match positions
                    p = locindices[e];

                    u1 = fst->IsUsedAt( p );
                    u2 = sec->IsUsedAt( p );
                    r1 = fst->ResidueAt( p );
                    r2 = sec->ResidueAt( p );

                    if( !u1 && !u2 )
                        continue;
                    if( r1 == GAP || r2 == GAP )
                        continue;
                    if( r1 == X || r2 == X )
                        continue;

                    segment++;
                    if( u1 && u2 && r1 == r2 )
                        matches++;
                }
                if( segment && GetIdentityLevel() <= ( double )matches /( double )segment )
                    sec->SetUsed( false );//set this sequence not to be used
            }
        }
    }

    if( bPASS2 ) {
        // 2nd pass...
        for( i = 0; i + 1 < size(); i++ )
        {
            fst = SequenceAt( ith=lindsrt[i] );
            if( fst == NULL || !fst->GetUsed())
                continue;

            for( j = i + 1; j < size(); j++ )
            {
                sec = SequenceAt( jth=lindsrt[j] );
                if( sec == NULL || !sec->GetUsed())
                    continue;

                simthr = true;
                segment = 0;
                matches = 0;

                for( e = 0, p = 0; e < effsize; e++ ) {
                    // iterate over match positions
                    p = locindices[e];

                    u1 = fst->IsUsedAt( p );
                    u2 = sec->IsUsedAt( p );
                    r1 = fst->ResidueAt( p );
                    r2 = sec->ResidueAt( p );

                    if( !u1 && !u2 )
                        continue;
                    if( r1 == GAP || r2 == GAP ) {
                        if( segment && sMINSEG <= segment && 
                          ( double )matches /( double )segment < dIDTHR2 ) {
                            simthr = false;
                            break;
                        }
                        segment = 0;
                        matches = 0;
                    }
                    if( r1 == X || r2 == X )
                        continue;

                    segment++;
                    if( u1 && u2 && r1 == r2 )
                        matches++;
                }
                if( simthr )
                    sec->SetUsed( false );
            }
        }
    }

// for( i = 0; i < size(); i++ )
//     if( !SequenceAt( i )->GetUsed())
//         fprintf( stderr, "s %d\n", i );
}

// -------------------------------------------------------------------------
// PurgeAtSequenceIdentity: retain diverse subset of sequences
//
void InputMultipleAlignment::PurgeAtSequenceIdentity()
{
    const bool      bRESORT = false;
    const ssize_t   nDETNOSEQS = 100;
    const size_t    szMINNOSEQS = 50;
    double          dLOSID = 0.2;
    const double    dHISID = GetIdentityLevel();
    double          dININC = 0.01;//initial increment
    const double    dUPINC = 0.05;//upper limit of initial increment
    const double    dCPC = 20.0;//power coefficient in map function
    const size_t    szFRAGLEN = 55;//length for fragment-based checking
    const size_t    szFRAGLENh = szFRAGLEN >> 1;
    bool    simthr;     //similar under threshold
    size_t  segment;    //length of aligned segment
    size_t  matches;    //number of residue matches in the sequences
    size_t  n, q;
    size_t  efst, eset, elst;
    size_t  i, ii, j, c, e, p, ith, jth;
    bool            u1, u2;
    unsigned char   r1, r2;
    PosDescriptionVector*   fst, *sec;

    if( !size())
        throw myruntime_error( "InputMultipleAlignment: No sequences." );

    PosDescriptionVector*   query = SequenceAt( 0 );

    if( query == NULL )
        throw myruntime_error( mystring( "InputMultipleAlignment: No query sequence obtained." ));

    const size_t noseqs = size();
    const size_t effsize = query->GetEffectiveSize();

    if( effsize == 0 )
        throw myruntime_error( mystring( "InputMultipleAlignment: Query sequence contains no residues." ));

    size_t  sum, nomps = 0;
    size_t  locindices[effsize];
    size_t  lens[TIMES2(noseqs)];
    size_t  lindsrt[noseqs];
    size_t  rdist[NUMALPH];
    size_t  prefnoseqs = noseqs;
    BinarySearchStructure   idlevs( IdlCompare, noseqs, true/*keep*/);
    const size_t    idlsc = 1000;
    size_t  noids;
    ssize_t idl = 0;
    double  qdist[NUMALPH];
    double  re, val, sfrq;//fraction of set of sequences
    double  medidl = 0.0;
    double  dsid, dthr, dinc;
    double  idthrs[effsize];
    size_t  nosqns[effsize];
    bool    passed, seld;

    //{{collect indices of match positions in query
    for( q = 0, nomps = 0; q < query->size(); q++ ) {
        r1 = query->ResidueAt(q);
        if( !IsValidResSym(r1) || r1 == X )
            continue;
        if( effsize <= nomps )
            throw myruntime_error( "InputMultipleAlignment: PurgeAtSequenceIdentity: Invalid query size." );
        locindices[nomps++] = q;
    }
//     lens[0] = nomps;
//     lens[1] = idlsc;
    lens[0] = idlsc;
    lens[1] = 0;
    //}}

    //{{calculate sequence identity levels w/r to query
    for( i = 1, ii = 2; i < noseqs; i++, ii+=2 )
    {
        lens[ii] = lens[ii+1] = 0;
        if(( sec = SequenceAt(i)) == NULL )
            continue;
        segment = 0;
        matches = 0;
        for( e = 0, p = 0; e < nomps; e++ ) {
            // iterate over match positions
            p = locindices[e];

            u1 = query->IsUsedAt( p );
            u2 = sec->IsUsedAt( p );
            r1 = query->ResidueAt( p );
            r2 = sec->ResidueAt( p );

            if( !u1 || !u2 )
                continue;
            if( r1 == GAP || r2 == GAP )
                continue;
            if( r1 == X || r2 == X )
                continue;

            segment++;
            if( u1 && u2 && r1 == r2 )
                matches++;
        }
        if( segment ) {
            idl = ( ssize_t )rint(( double )( matches * idlsc )/( double )segment );
            idlevs.Push(( const void* )idl );
//             lens[ii] = segment;
//             lens[ii+1] = ( size_t )idl;
            lens[ii] = ( size_t )idl;
            lens[ii+1] = i;
        }
    }
    //}}

    if( bRESORT ) {
        //{{calculate rel. entropy for each sequence
        for( i = 0, ii = 0; i < noseqs; i++, ii+=2 )
        {
            lens[ii+1] = 0;
            if(( sec = SequenceAt(i)) == NULL )
                continue;
            sum = 0;
            memset( rdist, 0, NUMALPH * sizeof( size_t ));
            for( e = 0, p = 0; e < nomps; e++ ) {
                // iterate over match positions
                p = locindices[e];

                u2 = sec->IsUsedAt( p );
                r2 = sec->ResidueAt( p );

                if(!( u2 && IsValidResSym(r2) && r2 != X ))
                    continue;
                rdist[r2]++;
                sum++;
            }
            re = 0.0;
            for( e = 0; e < NUMAA; e++ ) {
                val = (( double )rdist[e] + LOSCORES.PROBABility(e))/( double )( sum + 1 );
                if( !i )
                    qdist[e] = val;
                else if( val && qdist[e])
                    re += val * log( val / qdist[e]);
            }
            lens[ii+1] = ( size_t )rint(( double )idlsc * re );
        }
        //}}
    }

    //{{sort lengths of sequences
    lindsrt[0] = 0;
    //sort with reservation of the first place for query
    HeapSortInd( lindsrt+1, lens+2, noseqs-1, LenCompare );
    for( i = 1; i < noseqs; i++ )
        lindsrt[i]++;
    //}}


    if( 0 < nDETNOSEQS ) {
        prefnoseqs = ( size_t )nDETNOSEQS;
    }
    else if( !nDETNOSEQS ) {
        //get median value of identity levels
        noids = idlevs.GetSize();
        if( noids ) {
            if( noids & 1 )
                medidl = ( double )( ssize_t )idlevs.GetValueAt( noids>>1 )/( double )idlsc;
            else {
                medidl = ( double )((ssize_t)idlevs.GetValueAt((noids>>1)-1) + (ssize_t)idlevs.GetValueAt(noids>>1));
                medidl = medidl /( double )( TIMES2( idlsc ));
            }
            medidl = SLC_MIN( 1.0, medidl );
            //map sid level
    //         sfrq = exp( dCPC * medidl * log( 1.001 - medidl ));
            sfrq = exp(( dCPC*(1.0-SQUARE(medidl)) + medidl )* log(1.001-medidl));
    //         sfrq = 1.001 - medidl;
            prefnoseqs = SLC_MAX( szMINNOSEQS, ( size_t )rint( sfrq *( double )( noids + 1 )));
            if( noids + 1 < prefnoseqs )
                prefnoseqs = noids + 1;
        }
    }

    //{{incrementaly process mutual sequence identity levels
    if( noseqs <= prefnoseqs )
        dLOSID = dHISID;
    dININC = SLC_MAX( dININC, SLC_MIN( dUPINC, ( double )prefnoseqs * dININC * 0.01 ));
    for( i = 1; i < noseqs; i++ ) {
        sec = SequenceAt( ith=lindsrt[i] );
        if( sec )
            sec->SetUsed( false );//NOTE:initially unused
    }
    for( e = 0; e < nomps; e++ ) {
        idthrs[e] = 0.0;
        nosqns[e] = 1;
    }
    //
    for( dsid = dLOSID, dinc = dININC, seld = false; dsid <= dHISID && !seld; dsid += dinc )
    {
        seld = true;
        for( e = 0, p = 0; e < nomps; e++ ) {
            passed = false;
            //borders
            i = ( e <= szFRAGLENh )? 0: e - szFRAGLENh;
            j = i + szFRAGLEN;
            if( nomps <= j ) {
                j = nomps - 1;
                i = 0;
                if( szFRAGLEN < j )
                    i = j - szFRAGLEN;
            }
            for( c = i; c <= j; c++ ) {
                if( prefnoseqs <= nosqns[c]) {
                    passed = true;
                    break;
                }
            }
            if( !passed ) {
                seld = false;
                idthrs[e] = dsid;
            }
        }
//         for( e = 0, p = 0; e < nomps; e++ ) {
//             p = locindices[e];
//             if( nosqns[e] < prefnoseqs ) {
//                 seld = false;
//                 idthrs[e] = dsid;
//             }
//         }
        if( seld )
            break;
        for( i = 0; i < noseqs; i++ )
        {
            fst = SequenceAt( ith=lindsrt[i] );
            if( fst == NULL || fst->GetUsed())
                continue;//cont. if processed

            passed = true;
            dthr = 0.0;
            for( e = 0, p = 0; e < nomps; e++ ) {
                p = locindices[e];
                if( p < fst->GetFirstUsed()) continue;
                if( fst->GetLastUsed() < p ) break;
                if( dthr < idthrs[e])
                    dthr = idthrs[e];
                if( nosqns[e] < prefnoseqs )
                    passed = false;
            }
            if( dthr <= 0.0 )
                continue;
//             if( passed && dLOSID < dthr )
//                 continue;

            seld = false;
            passed = true;
            for( j = 0; j < i; j++ )
            {
                sec = SequenceAt( jth=lindsrt[j] );
                if( sec == NULL || !sec->GetUsed())
                    continue;//cont. if unused

                segment = 0;
                matches = 0;

                for( e = 0, p = 0; e < nomps; e++ ) {
                    //iterate over match positions
                    p = locindices[e];

                    u1 = fst->GetFirstUsed() <= p && p <= fst->GetLastUsed();//fst->IsUsedAt( p );
                    u2 = sec->GetFirstUsed() <= p && p <= sec->GetLastUsed();//sec->IsUsedAt( p );
                    r1 = fst->ResidueAt( p );
                    r2 = sec->ResidueAt( p );

                    if( !u1 || !u2 )
                       continue;
                    if( r1 == GAP || r2 == GAP )
                        continue;
                    if( r1 == X || r2 == X )
                        continue;

                    segment++;
                    if( u1 && u2 && r1 == r2 )
                        matches++;
                }
                if( segment && dthr <= ( double )matches /( double )segment ) {
                    passed = false;
                    break;
                }
            }
            if( !passed )
                continue;
            fst->SetUsed( true );//set used
            //use fragment-based instead of pos.-specific mode
            for( e = efst = eset = elst = 0, p = 0; e < nomps; e++ ) {
                p = locindices[e];
                if( p < fst->GetFirstUsed()) continue;
                if( fst->GetLastUsed() < p ) break;
                if( !eset ) { efst = e; eset = 1; }
                elst = e;
                nosqns[e]++;
            }
//             for( e = efst-1, c = 0; e + 1 && c < szFRAGLENh; e--, c++ )
//                 nosqns[e]++;
//             if( elst )
//                 for( e = elst+1, c = 0; e < nomps && c < szFRAGLENh; e++, c++ )
//                     nosqns[e]++;
        }
    }
    //}}
}

// -------------------------------------------------------------------------
// RefineWithSEG: purifies multiple alignment by applying SEG algorithm on
//     each match position
// -------------------------------------------------------------------------

void InputMultipleAlignment::RefineWithSEG()
{
    if( !queryDescription )
        throw myruntime_error( mystring( "InputMultipleAlignment: Unable to refine." ));

    const size_t    szFRAGLEN = 40;//application length from both ends
    const size_t    nomstates = queryDescription->GetEffectiveSize();
    const size_t    noseqs = size();
    size_t  segwinlenval    = GetSEGWindow();
    double  seglowentval    = GetSEGLowEntropy();
    double  seghighentval   = GetSEGHighEntropy();

    if( noseqs <= 1 || noseqs < segwinlenval )
        return;

    mystring        errstr;
    unsigned char*  validseq = ( unsigned char* )malloc( noseqs * sizeof( unsigned char ));
    unsigned char*  wholeseq = ( unsigned char* )malloc( noseqs * sizeof( unsigned char ));
    char*           omitmask = ( char* )malloc( noseqs * sizeof( char ));

    size_t          validlen = 0;
    size_t          wholelen = 0;
    size_t          p, pm, i;

    PosDescriptionVector*   seqn = NULL;
    unsigned char           res;
    int                     pstate;

    if( validseq == NULL || wholeseq == NULL || omitmask == NULL )
        throw myruntime_error( mystring( "InputMultipleAlignment: Not enough memory." ));

//     memset( validseq, 0, sizeof( unsigned char ) * size());
//     memset( wholeseq, 0, sizeof( unsigned char ) * size());
//     memset( omitmask, 0, sizeof( char ) * size());

    for( p = 0, pm =( size_t )-1; p < queryDescription->size() && errstr.empty(); p++ )
    {
        pstate = queryDescription->GetStateAt( p );
        if( !pstate )
            continue;

        pm++;
        if( szFRAGLEN <= pm && pm + szFRAGLEN < nomstates )
            continue;

        validlen = 0;
        wholelen = 0;

        //iterate over all sequences but the first (query)
        for( i = 1; i < noseqs; i++ )
        {
            seqn = SequenceAt( i );
            if( seqn == NULL ) {
                errstr =  "InputMultipleAlignment: RefineWithSEG: Memory access error.";
                break;
            }
            res = seqn->ResidueAt( p );

            wholelen++;
            omitmask[wholelen-1] = 0;
            wholeseq[wholelen-1] = res;

            if( !seqn->GetUsed() || !seqn->IsUsedAt( p ) || !IsValidResSym( res ) || res == X ) {
                omitmask[wholelen-1] = 1;
                continue;
            }

            validlen++;
            validseq[validlen-1] = res;
        }

        if( validlen < segwinlenval )
            continue;

        try {
            //SEG logic
            SEGSequence segseq(
                validseq,
                validlen,
                true/*hashed*/,
                segwinlenval,
                seglowentval,
                seghighentval
            );
            segseq.SetHighCSearch();//Important!
            segseq.Run();
            segseq.MaskSequence(( char* )wholeseq, wholelen, X, omitmask, 255/*symmask1*/, 255/*symmask2*/ );

        } catch( myexception const& ex ) {
            errstr = ex.what();
            break;
        }

        wholelen = 0;
        //set segged symbols as unused
        for( i = 1; i < noseqs; i++ ) {
            seqn = SequenceAt( i );

            wholelen++;
            if( omitmask[wholelen-1] )
                continue;

            if( wholeseq[wholelen-1] == X )
//                 seqn->SetUnusedAt( p );
                seqn->SetResidueAt( X, p );//sequence modification!
        }
    }

    free( validseq );
    free( wholeseq );
    free( omitmask );

    if( !errstr.empty())
        throw myruntime_error( errstr );
}

// -------------------------------------------------------------------------
// FilterSequencesWithSEG: filters sequences in multiple alignment with SEG
// -------------------------------------------------------------------------

void InputMultipleAlignment::FilterSequencesWithSEG()
{
    if( !queryDescription )
        throw myruntime_error( mystring( "InputMultipleAlignment: Unable to refine." ));

    size_t  seqsize = queryDescription->size();
    size_t  segwinlen   = GetSeqSEGWindow();
    double  seglowent   = GetSeqSEGLowEntropy();
    double  seghighent  = GetSeqSEGHighEntropy();

    if( seqsize <= 1 || seqsize < segwinlen )
        return;

    mystring        errstr;
    unsigned char*  validseq = ( unsigned char* )malloc( sizeof( unsigned char ) * seqsize );
    unsigned char*  wholeseq = ( unsigned char* )malloc( sizeof( unsigned char ) * seqsize );
    char*           omitmask = ( char* )malloc( sizeof( char ) * seqsize );

    size_t          validlen = 0;
    size_t          wholelen = 0;

    PosDescriptionVector*   sequence = NULL;
    unsigned char           residue;

    if( validseq == NULL || wholeseq == NULL || omitmask == NULL )
        throw myruntime_error( mystring( "InputMultipleAlignment: FilterSequencesWithSEG: Not enough memory." ));

//     memset( validseq, 0, sizeof( unsigned char ) * seqsize );
//     memset( wholeseq, 0, sizeof( unsigned char ) * seqsize );
//     memset( omitmask, 0, sizeof( char ) * seqsize );

    //iterate over all sequences
    for( size_t i = 0; i < size() && errstr.empty(); i++ ) {
        sequence = SequenceAt( i );
        if( sequence == NULL ) {
            errstr =  "InputMultipleAlignment: Memory access error in refinement.";
            break;
        }
        if( !sequence->GetUsed())
            continue;

        validlen = 0;
        wholelen = 0;

        //iterate over all positions
        for( size_t p = 0; p < seqsize; p++ )
        {
            residue = sequence->ResidueAt( p );

            wholelen++;
            omitmask[wholelen-1] = 0;
            wholeseq[wholelen-1] = residue;

            if( !sequence->IsUsedAt( p ) || residue == GAP ) {
                omitmask[wholelen-1] = 1;
                continue;
            }

            validlen++;
            validseq[validlen-1] = residue;
        }

        if( validlen < segwinlen )
            continue;

        try {
            //SEG logic
            SEGSequence segseq(
                validseq,
                validlen,
                true/*hashed*/,
                segwinlen,
                seglowent,
                seghighent
            );
            segseq.Run();
            segseq.MaskSequence(( char* )wholeseq, wholelen, X, omitmask, 255/*symmask1*/, 255/*symmask2*/ );

        } catch( myexception const& ex ) {
            errstr = ex.what();
            break;
        }

        wholelen = 0;
        //set segged symbols as unused
        for( size_t p = 0; p < seqsize; p++ )
        {
            wholelen++;
            if( omitmask[wholelen-1] )
                continue;

            if( wholeseq[wholelen-1] == X )
//                 sequence->SetUnusedAt( p );
                sequence->SetResidueAt( X, p );//sequence modification!
        }
    }

    free( validseq );
    free( wholeseq );
    free( omitmask );

    if( !errstr.empty())
        throw myruntime_error( errstr );
}

// -------------------------------------------------------------------------
// SetPosNoSequences: set number of sequences in column contributing to a 
//     position; overall number of sequences is also set
// -------------------------------------------------------------------------

void InputMultipleAlignment::SetPosNoSequences()
{
    if( !queryDescription )
        throw myruntime_error( mystring( "InputMultipleAlignment: Unable to set pos. number of sequences." ));

    size_t  posnos = 0;     //number of sequences per column
    size_t  nos = 0;        //number of sequences
    bool    allgaps;


    for( size_t p = 0; p < queryDescription->size(); p++ ) {
        posnos = 0;
        allgaps = true;

        //iterate over all sequences
        for( size_t i = 0; i < size(); i++ )
        {
            PosDescriptionVector*   sequence = SequenceAt( i );

            if( !sequence->GetUsed())
                continue;

            unsigned char   residue = sequence->ResidueAt( p );

            if( allgaps && IsValidResSym( residue ))
                allgaps = false;
            if(!( sequence->IsUsedAt( p ) && IsValidResSym( residue )))
                continue;

            posnos++;
        }
        queryDescription->SetNoSequencesAt( posnos, p );
        if( allgaps )
            queryDescription->SetUnusedAt(p);
    }

    //set number of sequences...
    for( size_t i = 0; i < size(); i++ )
        if( SequenceAt( i )->GetUsed())
            nos++;

    SetNoSequences( nos );
}

// -------------------------------------------------------------------------
// SetStates: set multiple alignment states
//
void InputMultipleAlignment::SetStates()
{
    if( !queryDescription )
        throw myruntime_error( "SetStates: No extended description vector." );

    unsigned char   residue;
    size_t          p, nn;

    for( p = 0, nn = 0; p < queryDescription->size(); p++ ) {
        residue = queryDescription->ResidueAt( p );
        if( !queryDescription->IsUsedAt( p ) || !IsValidResSym( residue ))
            continue;
        //match states are determined by query  
        queryDescription->SetStateAt( p );
        nn++;
    }
    queryDescription->SetEffectiveSize( nn );
    return;
}

// -------------------------------------------------------------------------
// SetBackgroundProbabilities: compute background probabilities for query
//     sequence with application of pseudo counts
//
void InputMultipleAlignment::SetBackgroundProbabilities()
{
    if( !queryDescription )
        throw myruntime_error( "SetBackgroundProbabilities: No extended description vector." );

    const double    backpseudocounts = 120.0;

    const double    accuracy = 1.0e-4;
    double          prob, consv;
    double          weight;
    unsigned char   residue;
    const size_t    effobslen = NUMAA;
    const size_t    obslen = NUMALPH;
    size_t          observs[obslen];
    size_t          noobs;
    size_t          p;

    static int  symB = HashAlphSymbol('B');
    static int  symZ = HashAlphSymbol('Z');
    static int  resN = HashAlphSymbol('N');
    static int  resD = HashAlphSymbol('D');
    static int  resQ = HashAlphSymbol('Q');
    static int  resE = HashAlphSymbol('E');

    noobs = 0;
    memset( observs, 0, obslen * sizeof( size_t ));

    for( p = 0; p < queryDescription->size(); p++ ) {
        residue = queryDescription->ResidueAt( p );
        if( !queryDescription->IsUsedAt( p ) ||
            residue == GAP || residue == X || residue == ASTERISK )
            continue;

        noobs += 2;

        if( residue == symB ) {
            observs[resN]++;
            observs[resD]++;
            continue;
        }
        if( residue == symZ ) {
            observs[resQ]++;
            observs[resE]++;
            continue;
        }
        if( NUMAA <= residue )
            throw myruntime_error( "SetBackgroundProbabilities: Unrecognized residue." );

        observs[residue] += 2;
    }

    if( backpseudocounts <= 0.0 && noobs < 1 )
        throw myruntime_error( "SetBackgroundProbabilities: Invalid counts." );

    consv = 0.0;
    weight = backpseudocounts /( noobs + backpseudocounts );

    if( noobs < 1 )
        noobs = 1;

    for( p = 0; p < obslen; p++ ) {
        prob = ( 1.0 - weight ) * ( double )observs[p] / ( double )noobs +
            weight * LOSCORES.PROBABility( p );
        consv += prob;

        queryDescription->SetBackProbsAt( p, prob );
    }

    if( consv < 1.0 - accuracy || consv > 1.0 + accuracy )
        throw myruntime_error( mystring( "SetBackgroundProbabilities: Probabilities are not conserved." ));

#ifdef USEPROFBACKPROBS
    LOSCORES.StoreProbabilities( queryDescription->GetBackProbs());
#endif

    return;
}

// -------------------------------------------------------------------------
// ComputeExtents: computes extents (left-right boundaries) for each match
//     position
//   Extent is a reduced multiple alignment constructed for each position
//     separately so that sequences in the constructed alignment 
//     contribute a residue or internal gap symbol
//
void InputMultipleAlignment::ComputeExtents()
{
    if( !queryDescription )
        throw myruntime_error( mystring( "Unable to compute extents: no extended description vector." ));

    unsigned char   residue;
    ssize_t matchbeg = -1;//first match state position
    ssize_t matchend = -1;//last match state position

    size_t  extwinlen = GetExtentMinWindow();
    size_t  extcentre =   extwinlen >> 1;
    size_t  extadjust = ( extwinlen &  1 ) ^ 1; //for even window lengths, adjust
    size_t  i, p;

//     //do not use percenatages as this does not anyway relates to domain boundaries
//     extwinlen = ( size_t )(( double )queryDescription->GetEffectiveSize() * GetExtentMinSeqPercentage());
    if( extwinlen < GetExtentMinWindow())
        extwinlen = GetExtentMinWindow();
    extcentre =   extwinlen >> 1;
    extadjust = ( extwinlen &  1 ) ^ 1; //for even window lengths, adjust

// fprintf(stderr,"W=%d P=%f w=%d c=%d a=%d\n",GetExtentMinWindow(),GetExtentMinSeqPercentage(),extwinlen,extcentre,extadjust);
    size_t*     matchpos = NULL;    //accumulated number of match positions
    size_t*     msvector = NULL;    //accumulated match state vector

    matchpos = ( size_t* )malloc( sizeof( size_t ) * queryDescription->size());
    msvector = ( size_t* )malloc( sizeof( size_t ) * queryDescription->size());
    if( !matchpos || !msvector )
        throw myruntime_error( "ComputeExtents: Not enough memory." );
    memset( matchpos, 0, sizeof( size_t ) * queryDescription->size());
    memset( msvector, 0, sizeof( size_t ) * queryDescription->size());

    for( p = 0; p < queryDescription->size(); p++ ) {
        if( IsValidResSym( queryDescription->ResidueAt( p ))) {
            if( matchbeg < 0 ) matchbeg = p;
            matchend = p;
        }
    }
    if( matchbeg < 0 || matchend < 0 )
        throw myruntime_error( mystring( "ComputeExtents: Invalid match state sequence boundaries." ));

    //iterate over all sequences in multiple alignment
    for( i = 0; i < size(); i++ )
    {
        PosDescriptionVector*   sequence = SequenceAt(i);
        if( sequence == NULL || !sequence->GetUsed())
            continue;

        size_t  seqsize = queryDescription->size();
        size_t  left = SIZE_MAX,
                right = SIZE_MAX;

        //set match pos.: cumulative values
        memset( matchpos, 0, seqsize * sizeof( size_t ));
        for( p = 0; p < seqsize; p++ ) {
            if( p )
                matchpos[p] = matchpos[p-1];
            if( sequence->IsUsedAt( p ) &&
                IsValidResSym( sequence->ResidueAt( p )) && 
                IsValidResSym( queryDescription->ResidueAt( p ))) {
                    if( p ) matchpos[p] = matchpos[p-1] + 1;
                    else    matchpos[p] = 1; 
            }
        }
        //set support match pos.
        if( i == 0 )
            memcpy( msvector, matchpos, seqsize * sizeof( size_t ));


        for( p = 0; p < seqsize; p++ ) {
            if( !sequence->IsUsedAt( p )) {
                continue;
            }
            residue = sequence->ResidueAt( p );

            queryDescription->IncCountAt( p );
            queryDescription->IncDistributionAt( residue, p );

            if( residue != GAP ) {
                if( left == SIZE_MAX ) left = p;
                right = p;
            }
        }

        if( left == SIZE_MAX || right == SIZE_MAX )
            continue;

        size_t  sres, mres, pos, lres, rres;
        ssize_t lmsbder = 0, rmsbder = 0;
        ssize_t lborder = 0, rborder = 0;
        ssize_t leftadj = 0, rightadj = 0;

        if( left < matchbeg ) left = matchbeg;
        if( matchend < right ) right = matchend;

        lres = 0; if(( left? matchpos[left-1]: 0 ) < matchpos[left]) lres = 1;
        rres = 0; if(( right? matchpos[right-1]: 0 ) < matchpos[right]) rres = 1;

        for( p = 0; p < seqsize; p++ ) {
            if( !sequence->IsUsedAt( p )) {
                continue;
            }

            leftadj = extcentre - extadjust;
            rightadj = extcentre;

            lborder = matchpos[p] - 1 - leftadj;
            rborder = matchpos[p] - 1 + rightadj;

            lmsbder = msvector[p] - 1 - leftadj;
            rmsbder = msvector[p] - 1 + rightadj;

            //process cases when window is out of boundaries of match-state sequences
            if( lmsbder < 0 ) { 
                pos = p;
                if( p < matchbeg ) pos = matchbeg;
                sres = 0;
                mres = 0;
                if( 0 < matchbeg ) {
                    if( matchpos[matchbeg-1] < matchpos[matchbeg]) sres = 1;
                    if( msvector[matchbeg-1] < msvector[matchbeg] ) mres = 1;
                } else {
                    if( 0 < matchpos[matchbeg]) sres = 1;
                    if( 0 < msvector[matchbeg] ) mres = 1;
                }
//                 lborder = matchpos[pos] - msvector[pos]; //exactly the same as below
                lborder = matchpos[matchbeg] - sres + 
                        //difference of beginnings of this and match sequences
                        ( matchpos[pos] - matchpos[matchbeg] + sres ) - 
                        ( msvector[pos] - msvector[matchbeg] + mres );
                rborder = lborder + extwinlen - 1; 
            }
            if(( ssize_t )( msvector[matchend] - 1 ) < rmsbder ) {
                pos = p;
                if( matchend < p ) pos = matchend;
                sres = 0;
                mres = 0;
                if( 0 < pos ) {
                    if( matchpos[pos-1] < matchpos[pos]) sres = 1;
                    if( msvector[pos-1] < msvector[pos] ) mres = 1;
                } else {
                    if( 0 < matchpos[pos]) sres = 1;
                    if( 0 < msvector[pos] ) mres = 1;
                }

                rborder = matchpos[matchend] - 1 + //for -1 assume there is always at least 1 res. in sequence
                        //difference of tails of this and match sequences
                        ( msvector[matchend] - msvector[pos] + mres ) - 
                        ( matchpos[matchend] - matchpos[pos] + sres ); 
                lborder = rborder - extwinlen + 1;
            }

            //include sequence in the extent
            if( left <= p && p <= right ) {
                //queryDescription->PushIndexAt( i, p ); //avoid because of required additional memory
                sequence->SetUsedInExtentAt( p ); //this makes computations faster
                queryDescription->IncNoSequencesInExtentAt( p ); //one more sequence in the extent at the pos.
                queryDescription->IncNoSequencesInMSExtentAt( p );
            }

// fprintf(stderr,"%d,%d:: %d - %d <= %d  &&  %d <= %d - %d::   %d <= %d && %d <= %d:: %d\n",
// i,p,matchpos[left],lres,lborder,rborder,matchpos[right],1,
// (ssize_t)(matchpos[left]-lres),lborder,rborder,(ssize_t)(matchpos[right]-1),
// (ssize_t)(matchpos[left]-lres)<=lborder && rborder<=(ssize_t)(matchpos[right]-1));

            //set extent boundaries
            if((( ssize_t ) extwinlen <= ( ssize_t )( msvector[matchend] - msvector[matchbeg] + 1 )) ?
               (( ssize_t )( matchpos[left] - lres )<= lborder && 
                rborder <= ( ssize_t )( matchpos[right] - 1 ))
                :
               ( left <= p && p <= right ))
            {
                // omit positions which are unsued or gaps in query
                if( !queryDescription->IsUsedAt( p ))
                    continue;

                if( queryDescription->GetCountAt( p ) == 1 ||
///                    left < queryDescription->GetLeftExtentAt( p ))       //was worth to verify but worked a bit worse
                    queryDescription->GetLeftExtentAt( p ) < left ) {
                    queryDescription->SetLeftExtentAt( left, p );
                    queryDescription->SetLeftMSExtentAt( left, p );
                }

                if( queryDescription->GetCountAt( p ) == 1 ||
///                    queryDescription->GetRightExtentAt( p ) < right )    //was worth to verify but worked a bit worse
                    right < queryDescription->GetRightExtentAt( p )) {
                    queryDescription->SetRightExtentAt( right, p );
                    queryDescription->SetRightMSExtentAt( right, p );
                }
            }
        }
    }
    //save intervals of extents for each position
    for( p = 0; p < queryDescription->size(); p++ ) {
        // omit positions which are unsued or gaps in query
        if( !queryDescription->IsUsedAt( p ))
            continue;

        size_t  left = queryDescription->GetLeftExtentAt( p );
        size_t  right = queryDescription->GetRightExtentAt( p );
        size_t  interval;

        if( right < left || right == SIZE_MAX || left == SIZE_MAX )
            continue;

//         // adjust left boundary values to properly compute intervals
//         size_t  leftaccpos = left? notused[0][left-1]: 0;
//         size_t  rightaccpos = notused[0][right];
//         interval =  right - left + 1 - ( rightaccpos - leftaccpos );

        left = queryDescription->GetLeftMSExtentAt( p );
        right = queryDescription->GetRightMSExtentAt( p );
        interval =  msvector[right] - msvector[left] + 1;
        queryDescription->SetExtentIntervalAt( interval, p );
        queryDescription->SetMSExtentIntervalAt( interval, p );
// fprintf( stderr, "%d %c %d,%d  %d\n", p, DehashCode( queryDescription->ResidueAt( p )), left, right, interval );
    }

    free( matchpos );
    free( msvector );
}

// -------------------------------------------------------------------------
// ComputeSequenceWeights: compute general sequence weights
//
void InputMultipleAlignment::ComputeSequenceWeights()
{
    if( !queryDescription )
        throw myruntime_error( "Unable to compute sequence weights." );

    size_t  noress = NUMALPH;   // number of residues
    size_t  noeffress = NUMAA;  // effective number of residues
    bool    newset = true;      // indicates a new set of sequences in extent
    size_t  numseq = 1;         // number of sequences in extent
    size_t  extentN = 0, ms_extentN = 0;    // number of different symbols occuring in extent w/o identical columns
    size_t  extentIdN = 0, ms_extentIdN = 0;// number of different symbols occuring in extent (identical columns inc.)
    size_t  column[NUMALPH];    // residue distribution in column
    size_t  diffsyms = 0;       // number of different symbols in column
    size_t  nadiffsyms = 0;     // non-abstract different symbols in column
    double  wghtsum = 0.0;      // sum of weights used to normalize sequence weights
    double*         Mweights = NULL;// sequence weights computed over all extents
    double*         Iweights = NULL;
    double*         Dweights = NULL;
    double          w;
    unsigned char   residue;
    mystring        merror;
    size_t  p, prep, i;         //position indices
    ssize_t pp;


    // iterate over all query positions
    for( p = 0, prep = 0; p < queryDescription->size(); p++ ) {
        // omit positions which are unsued or gaps in query
        if( !queryDescription->IsUsedAt( p ))
                continue;

        // do not compute weights for one sequence
        if( //queryDescription->GetCountAt( p ) <= 1 ||
            //compute it anyway in order to obtain reasonable profile-pair scores
            //queryDescription->GetNoSequencesInExtentAt( p ) <= 1 ||
            !queryDescription->GetExtentIntervalAt( p ))
            continue;

        // determine if it is the same set of sequences;
        // the extent will always be the same if sequences comprise the same set!
        newset = true;
        if( p ) {
            newset = true;
            pp = prep;
//             for( pp = p-1; 0 <= pp && newset; pp-- ) {
            if( pp < p ) {
                for( i = 0; i < size(); i++ )
                    if( SequenceAt( i )->IsUsedInExtentAt( p ) !=
                        SequenceAt( i )->IsUsedInExtentAt( pp ))
                            break;
                if( i == size()) {
                    newset = false;
                    queryDescription->SetSqnWeightsAt( p, PS_M, pp );
                    queryDescription->SetSqnWeightsAt( p, PS_I, pp );
                    queryDescription->SetSqnWeightsAt( p, PS_D, pp );
                    Mweights = queryDescription->GetSqnWeightsAt( p, PS_M );
                    Iweights = queryDescription->GetSqnWeightsAt( p, PS_I );
                    Dweights = queryDescription->GetSqnWeightsAt( p, PS_D );
                    if( queryDescription->GetLeftMSExtentAt( p ) != queryDescription->GetLeftMSExtentAt( pp ) ||
                        queryDescription->GetRightMSExtentAt( p ) != queryDescription->GetRightMSExtentAt( pp ))
                    {
                        newset = true;
                        Mweights = Iweights = Dweights = NULL;
                    }
                }
            }
            if( !newset && Mweights == NULL ) {
                // it can be so that weights are NULL in case when e.g. the
                // beginning of alignment consists of one sequence and it's not processed
                newset = true;
//                 merror = "Memory access error.";
//                 break;
            }
            if( !newset ) {
                if( p < prep + 1 ) {
                    merror = "Unable to compute sequence weights: Invalid position.";
                    break;
                }
                pp = prep;
                for( nadiffsyms = 0; nadiffsyms < noress; nadiffsyms++ ) {
                    queryDescription->SetDistinctHistAt( 
                        queryDescription->GetDistinctHistAt( nadiffsyms, pp ), nadiffsyms, p );
                    queryDescription->SetMSDistinctHistAt( 
                        queryDescription->GetMSDistinctHistAt( nadiffsyms, pp ), nadiffsyms, p );
                }
            }
        }
        if( newset ) {
            prep = p;
            numseq = 0;
            ms_extentN = extentN = 0;
            ms_extentIdN = extentIdN = 0;
            wghtsum = 0.0;
            queryDescription->NewSqnWeightsAt( p, PS_M, size());
            queryDescription->NewSqnWeightsAt( p, PS_I, size());
            queryDescription->NewSqnWeightsAt( p, PS_D, size());
            Mweights = queryDescription->GetSqnWeightsAt( p, PS_M );
            Iweights = queryDescription->GetSqnWeightsAt( p, PS_I );
            Dweights = queryDescription->GetSqnWeightsAt( p, PS_D );

            if( Mweights == NULL || Iweights == NULL || Dweights == NULL ) {
                merror = "Null weights.";
                break;
            }

            size_t  left = queryDescription->GetLeftMSExtentAt( p );
            size_t  right = queryDescription->GetRightMSExtentAt( p );

            for( size_t k = left; 0 <=( ssize_t )right && k <= right; k++ ) {
                // omit positions which are unused or gaps in query
                if( !queryDescription->IsUsedAt( k ))
                    continue;

                diffsyms = nadiffsyms = 0;
                memset( column, 0, sizeof( size_t ) * NUMALPH );

                for( i = 0; i < size(); i++ ) {
                    // omit sequences not in the extent
                    if( !SequenceAt( i )->IsUsedInExtentAt( p ))
                        continue;
                    // compute statistics
                    residue = SequenceAt( i )->ResidueAt( k );
                    if( column[residue]++ == 0 ) {
                        diffsyms++;
                        if( IsValidResSym( residue ) && residue != X )
                            nadiffsyms++;
                    }
                }

                extentIdN += diffsyms;
                if( 1 < diffsyms )
                    extentN += diffsyms;
                if( noeffress < nadiffsyms )
                    nadiffsyms = noeffress;

                queryDescription->IncDistinctHistAt( nadiffsyms, p );

                if( queryDescription->GetStateAt( k )) {
                    ms_extentIdN += diffsyms;
                    if( 1 < diffsyms )
                        ms_extentN += diffsyms;
                    queryDescription->IncMSDistinctHistAt( nadiffsyms, p );
                }

                if( queryDescription->GetStateAt( k )) ///
                    for( i = 0; i < size(); i++ ) {
                        // omit sequences not in the extent
                        if( !SequenceAt( i )->IsUsedInExtentAt( p ))
                            continue;
                        residue = SequenceAt( i )->ResidueAt( k );
                        w = 1.0 / ( column[residue] * diffsyms );
                        Mweights[i] += w;
                        Iweights[i] += w;
                        Dweights[i] += w;
                        wghtsum += w;
                    }
            }
            // normalize sequence weights computed for the extent
            if( wghtsum )
                for( i = 0; i < size(); i++ ) {
                    if( !SequenceAt( i )->IsUsedInExtentAt( p ))
                        continue;
                    Mweights[i] /= wghtsum;
                    Iweights[i] /= wghtsum;
                    Dweights[i] /= wghtsum;
                }
            else
                // actually it cannot pass this way but...
                for( i = 0; i < size(); i++ ) {
                    if( !SequenceAt( i )->IsUsedInExtentAt( p ))
                        continue;
                    if( queryDescription->GetNoSequencesInMSExtentAt( p ))
                        Mweights[i] = Iweights[i] = Dweights[i] = 
                            1.0 / queryDescription->GetNoSequencesInMSExtentAt( p );
                    else
                        Mweights[i] = Iweights[i] = Dweights[i] = 0.0;
                }
        }

        // if the extent contains at least one non-identical column
//      if( extentN ) {
//      }

        queryDescription->SetNoSymbolsInExtentAt( extentIdN, p );
        queryDescription->SetNoSymbolsInMSExtentAt( ms_extentIdN, p );

        for( i = 0; i < size(); i++ ) {
            if( !SequenceAt( i )->IsUsedInExtentAt( p ))
                continue;
            residue = SequenceAt( i )->ResidueAt( p );
            queryDescription->IncMatchWeightsAt( Mweights[i], residue, p );
        }
    }

    //check for error
    if( !merror.empty())
        throw myruntime_error( merror );

    return;
}





// =========================================================================
// ComputeMstateSequenceWeights: compute M state sequence weights
//
void InputMultipleAlignment::ComputeMstateSequenceWeights()
{
    if( !queryDescription )
        throw myruntime_error("ComputeMstateSequenceWeights: Null description data.");

    size_t  noress = NUMALPH;   // number of residues
    size_t  noeffress = NUMAA;  // effective number of residues
    bool    newset = true;      // indicates a new set of sequences in extent
    size_t  numseq = 1;         // number of sequences in extent
    size_t  extentIdN = 0;      // number of different symbols occuring in extent (identical columns inc.)
    size_t  column[NUMALPH];    // residue distribution in column
    size_t  diffsyms = 0;       // number of different symbols in column
    size_t  nadiffsyms = 0;     // non-abstract different symbols in column
    size_t  histval;            // histogram value for residue
    double  wghtsum = 0.0;      // sum of weights used to normalize sequence weights
    double*         weights = NULL;// sequence weights computed over all extents
    double          w;
    unsigned char   residue, pres;
    char            errbuf[KBYTE];
    const mystring  preamble = "ComputeMstateSequenceWeights: ";
    mystring        merror;
    size_t  p, prep, k, i; //position indices
    ssize_t pp;


    // iterate over all query positions
    for( p = 0, prep = 0; p < queryDescription->size() && merror.empty(); p++ ) {
        // omit unsued and not match-state positions
        if( !queryDescription->IsUsedAt( p ) || 
            !queryDescription->GetStateAt( p ))
            continue;

        // do compute weights even for one sequence
        if( !queryDescription->GetExtentIntervalAt( p ))
            continue;

        // determine if it is the same set of sequences;
        // the extent will always be the same if sequences comprise the same set!
        newset = true;
        if( p ) {
            newset = true;
            pp = prep;
            if( pp < p ) {
                for( i = 0; i < size(); i++ ) {
                    //TODO: use state extents
                    if( SequenceAt( i )->IsUsedInExtentAt( p ) !=
                        SequenceAt( i )->IsUsedInExtentAt( pp ))
                            break;
                    pres = SequenceAt( i )->ResidueAt( p );
                    residue = SequenceAt( i )->ResidueAt( pp );
                    if( IsValidResSym( pres ) != IsValidResSym( residue ))
                        break;
                }
                if( i == size()) {
                    newset = false;
                    queryDescription->SetSqnWeightsAt( p, PS_M, pp );
                    weights = queryDescription->GetSqnWeightsAt( p, PS_M );
                    if( queryDescription->GetLeftMSExtentAt( p ) != queryDescription->GetLeftMSExtentAt( pp ) ||
                        queryDescription->GetRightMSExtentAt( p ) != queryDescription->GetRightMSExtentAt( pp ))
                    {
                        newset = true;
                        weights = NULL;
                    }
                }
            }
            if( !newset && weights == NULL ) {
                // it can be so that weights are NULL in case when e.g. the
                // beginning of alignment consists of one sequence and it's not processed;
                // should not be the case when even a single sequence is processed
//                 newset = true;
                sprintf( errbuf, "Unexpected null weights for pos. %d.", pp );
                merror = errbuf;
                break;
            }
            if( !newset ) {
                if( p < prep + 1 ) {
                    merror = "Invalid index of previous position.";
                    break;
                }
                pp = prep;
                for( nadiffsyms = 0; nadiffsyms < noress; nadiffsyms++ ) {
                    histval = queryDescription->GetDistinctHistAt( pp, PS_M, nadiffsyms );
                    queryDescription->SetDistinctHistAt( p, PS_M, nadiffsyms, histval );
                }
            }
        }
        if( newset ) {
            prep = p;
            extentIdN = 0;
            wghtsum = 0.0;
            queryDescription->NewSqnWeightsAt( p, PS_M, size());
            weights = queryDescription->GetSqnWeightsAt( p, PS_M );

            if( weights == NULL ) {
                merror = "Null weights.";
                break;
            }

            size_t  left = queryDescription->GetLeftMSExtentAt( p );
            size_t  right = queryDescription->GetRightMSExtentAt( p );

            for( k = left; 0 <=( ssize_t )right && k <= right; k++ ) {
                // omit unsued and unmatched positions in the extent
                if( !queryDescription->IsUsedAt( k ))
                    continue;
//                 if( !queryDescription->GetStateAt( k ))
//                     continue;

                diffsyms = nadiffsyms = 0;
                memset( column, 0, sizeof( size_t ) * NUMALPH );

                for( i = 0; i < size(); i++ ) {
                    // omit sequences not in the extent
                    if( !SequenceAt( i )->IsUsedInExtentAt( p ))
                        continue;
                    pres = SequenceAt( i )->ResidueAt( p );
                    residue = SequenceAt( i )->ResidueAt( k );
                    // leave only sequences which are in match state at position `p'
                    if( pres == GAP || !IsValidResSym( pres ))
                        continue;
                    if( column[residue]++ == 0 ) {
                        diffsyms++;
                        if( IsValidResSym( residue ))
                            nadiffsyms++;
                    }
                }

                if( noeffress < nadiffsyms )
                    nadiffsyms = noeffress;

                //update histogram only in match state
                if( queryDescription->GetStateAt( k )) {
                    extentIdN += nadiffsyms;
                    queryDescription->IncDistinctHistAt( p, PS_M, nadiffsyms );
                }

                if( nadiffsyms ) {
                    for( i = 0; i < size(); i++ ) {
                        // omit sequences not in the extent
                        if( !SequenceAt( i )->IsUsedInExtentAt( p ))
                            continue;
                        pres = SequenceAt( i )->ResidueAt( p );
                        residue = SequenceAt( i )->ResidueAt( k );
                        // use only sequences which are in match state at position `p'
                        if( pres == GAP || !IsValidResSym( pres ))
                            continue;
                        if( !IsValidResSym( residue ))
                            continue;
                        w = 1.0 /( double )( column[residue] * nadiffsyms );
                        weights[i] += w;
                        wghtsum += w;
                    }
                }
            }

            // normalize sequence weights computed for the extent
            if( !wghtsum ) {
                numseq = 0;
                for( i = 0; i < size(); i++ ) {
                    if( !SequenceAt( i )->IsUsedInExtentAt( p ))
                        continue;
                    pres = SequenceAt( i )->ResidueAt( p );
                    if( pres == GAP || !IsValidResSym( pres ))
                        continue;
                    weights[i] = 1.0;
                    numseq++;
                }
                wghtsum = numseq;
            }
            if( !wghtsum ) {
                sprintf( errbuf, "No contribution to sequence weights at pos. %d.", p );
                merror = errbuf;
                break;
            }
            for( i = 0; i < size(); i++ ) {
                if( !SequenceAt( i )->IsUsedInExtentAt( p ))
                    continue;
                pres = SequenceAt( i )->ResidueAt( p );
                if( pres == GAP || !IsValidResSym( pres ))
                    continue;
                weights[i] /= wghtsum;
            }
        }

        queryDescription->SetNoSymbolsInExtentAt( p, PS_M, extentIdN );

        for( i = 0; i < size(); i++ ) {
            if( !SequenceAt( i )->IsUsedInExtentAt( p ))
                continue;
            pres = SequenceAt( i )->ResidueAt( p );
            if( pres == GAP || !IsValidResSym( pres ))
                continue;
            queryDescription->IncMatchWeightsAt( weights[i], pres, p );
        }
    }

    //check for error
    if( !merror.empty())
        throw myruntime_error( preamble + merror );
}

// -------------------------------------------------------------------------
// ComputeIstateSequenceWeights: compute I state sequence weights
//
void InputMultipleAlignment::ComputeIstateSequenceWeights()
{
    if( !queryDescription )
        throw myruntime_error("ComputeIstateSequenceWeights: Null description data.");

    size_t  noress = NUMALPH;   // number of residues
    size_t  noeffress = NUMAA;  // effective number of residues
    bool    newset = true;      // indicates a new set of sequences in extent
    size_t  numseq = 1;         // number of sequences in extent
    size_t  extentIdN = 0;      // number of different symbols occuring in extent (identical columns inc.)
    size_t  Ihist[NUMALPH];     // accumulated histogram data
    size_t  column[NUMALPH];    // residue distribution in column
    size_t  diffsyms = 0;       // number of different symbols in column
    size_t  nadiffsyms = 0;     // non-abstract different symbols in column
    size_t  histval;            // histogram value for residue
    size_t  rdif;
    int     pstate;
    double  wghtsum = 0.0;      // sum of weights used to normalize sequence weights
    double*         weights = NULL;// sequence weights computed over all extents
    double          w;
    unsigned char   residue;
    char            errbuf[KBYTE];
    const mystring  preamble = "ComputeIstateSequenceWeights: ";
    mystring        merror;
    size_t  p, pM, prep, k, i; //position indices
    ssize_t pp;

    memset( Ihist, 0, sizeof( size_t ) * NUMALPH );


    // iterate over all query positions
    for( p = 0, pM = prep = 0; p < queryDescription->size() && merror.empty(); p++ ) {
        // omit unsued and not insert-state positions
        if( !queryDescription->IsUsedAt( p ))
            continue;

        // do compute weights even for one sequence
        if( !queryDescription->GetExtentIntervalAt( p ))
            continue;

        pstate = queryDescription->GetStateAt( p );

        if( pstate ) {
            if( pM )
                //save accumulated histogram data of I state
                for( rdif = 0; rdif < noress; rdif++ ) {
                    queryDescription->SetDistinctHistAt( pM, PS_I, rdif, Ihist[rdif]);
                    Ihist[rdif] = 0;
                }
            pM = p;
            continue;
        }

        // determine if it is the same set of sequences;
        // the extent will always be the same if sequences comprise the same set!
        newset = true;
        if( p ) {
            newset = true;
            pp = prep;
            if( pp < p ) {
                for( i = 0; i < size(); i++ )
                    //TODO: use state extents
                    if( SequenceAt( i )->IsUsedInExtentAt( p ) !=
                        SequenceAt( i )->IsUsedInExtentAt( pp ))
                            break;
                if( i == size()) {
                    newset = false;
                    queryDescription->SetSqnWeightsAt( p, PS_I, pp );
                    weights = queryDescription->GetSqnWeightsAt( p, PS_I );
                    if( queryDescription->GetLeftMSExtentAt( p ) != queryDescription->GetLeftMSExtentAt( pp ) ||
                        queryDescription->GetRightMSExtentAt( p ) != queryDescription->GetRightMSExtentAt( pp ))
                    {
                        newset = true;
                        weights = NULL;
                    }
                }
            }
//             if( !newset && weights == NULL ) {
// //                 newset = true;
//                 sprintf( errbuf, "Unexpected null weights for pos. %d.", pp );
//                 merror = errbuf;
//                 break;
//             }
            if( !newset && weights ) {
                if( p < prep + 1 ) {
                    merror = "Invalid index of previous position.";
                    break;
                }
                pp = prep;
                for( nadiffsyms = 0; nadiffsyms < noress; nadiffsyms++ ) {
                    histval = queryDescription->GetDistinctHistAt( pp, PS_I, nadiffsyms );
                    queryDescription->SetDistinctHistAt( p, PS_I, nadiffsyms, histval );
                }
            }
        }
        if( newset ) {
            prep = p;
            extentIdN = 0;
            wghtsum = 0.0;
            queryDescription->NewSqnWeightsAt( p, PS_I, size());
            weights = queryDescription->GetSqnWeightsAt( p, PS_I );

            if( weights == NULL ) {
                merror = "Null weights.";
                break;
            }

            size_t  left = queryDescription->GetLeftMSExtentAt( p );
            size_t  right = queryDescription->GetRightMSExtentAt( p );

            for( k = left; 0 <=( ssize_t )right && k <= right; k++ ) {
                // omit unsued and unmatched positions in the extent
                if( !queryDescription->IsUsedAt( k ))
                    continue;
//                 if( queryDescription->GetStateAt( k ))
//                     continue;

                diffsyms = nadiffsyms = 0;
                memset( column, 0, sizeof( size_t ) * NUMALPH );

                for( i = 0; i < size(); i++ ) {
                    // omit sequences not in the extent
                    if( !SequenceAt( i )->IsUsedInExtentAt( p ))
                        continue;
                    residue = SequenceAt( i )->ResidueAt( k );
                    if( column[residue]++ == 0 ) {
                        diffsyms++;
                        if( IsValidResSym( residue ))
                            nadiffsyms++;
                    }
                }

                if( noeffress < nadiffsyms )
                    nadiffsyms = noeffress;

                //update histogram only in insert state
                if( !queryDescription->GetStateAt( k )) {
                    extentIdN += nadiffsyms;
                    queryDescription->IncDistinctHistAt( p, PS_I, nadiffsyms );
                }

                if( nadiffsyms ) {
                    for( i = 0; i < size(); i++ ) {
                        // omit sequences not in the extent
                        if( !SequenceAt( i )->IsUsedInExtentAt( p ))
                            continue;
                        residue = SequenceAt( i )->ResidueAt( k );
                        if( !IsValidResSym( residue ))
                            continue;
                        w = 1.0 /( double )( column[residue] * nadiffsyms );
                        weights[i] += w;
                        wghtsum += w;
                    }
                }
            }

            //normalize if applicable
            if( wghtsum )
                //if inserts are present in an extent, normalize;
                //otherwise, all weights will be 0
                for( i = 0; i < size(); i++ ) {
                    if( !SequenceAt( i )->IsUsedInExtentAt( p ))
                        continue;
                    weights[i] /= wghtsum;
                }
        }

        queryDescription->SetNoSymbolsInExtentAt( p, PS_I, extentIdN );
        //accumulate histogram data of I state
        for( rdif = 0; rdif < noress; rdif++ ) {
            histval = queryDescription->GetDistinctHistAt( p, PS_I, rdif );
            Ihist[rdif] += histval;
        }
    }

    //check for error
    if( !merror.empty())
        throw myruntime_error( preamble + merror );
}

// -------------------------------------------------------------------------
// ComputeMstateSequenceWeights: compute D state sequence weights
//
void InputMultipleAlignment::ComputeDstateSequenceWeights()
{
    if( !queryDescription )
        throw myruntime_error("ComputeDstateSequenceWeights: Null description data.");

    size_t  noress = NUMALPH;   // number of residues
    size_t  noeffress = NUMAA;  // effective number of residues
    bool    newset = true;      // indicates a new set of sequences in extent
    size_t  numseq = 1;         // number of sequences in extent
    size_t  extentIdN = 0;      // number of different symbols occuring in extent (identical columns inc.)
    size_t  column[NUMALPH];    // residue distribution in column
    size_t  diffsyms = 0;       // number of different symbols in column
    size_t  nadiffsyms = 0;     // non-abstract different symbols in column
    size_t  histval;            // histogram value for residue
    double  wghtsum = 0.0;      // sum of weights used to normalize sequence weights
    double*         weights = NULL;// sequence weights computed over all extents
    double          w;
    unsigned char   residue, pres;
    char            errbuf[KBYTE];
    const mystring  preamble = "ComputeDstateSequenceWeights: ";
    mystring        merror;
    size_t  p, prep, k, i; //position indices
    ssize_t pp;


    // iterate over all query positions
    for( p = 0, prep = 0; p < queryDescription->size() && merror.empty(); p++ ) {
        // omit unsued and not match-state positions
        if( !queryDescription->IsUsedAt( p ) || 
            !queryDescription->GetStateAt( p ))
            continue;

        // do compute weights even for one sequence
        if( !queryDescription->GetExtentIntervalAt( p ))
            continue;

        // determine if it is the same set of sequences;
        // the extent will always be the same if sequences comprise the same set!
        newset = true;
        if( p ) {
            newset = true;
            pp = prep;
            if( pp < p ) {
                for( i = 0; i < size(); i++ ) {
                    //TODO: use state extents
                    if( SequenceAt( i )->IsUsedInExtentAt( p ) !=
                        SequenceAt( i )->IsUsedInExtentAt( pp ))
                            break;
                    pres = SequenceAt( i )->ResidueAt( p );
                    residue = SequenceAt( i )->ResidueAt( pp );
                    if( IsValidResSym( pres ) != IsValidResSym( residue ))
                        break;
                }
                if( i == size()) {
                    newset = false;
                    queryDescription->SetSqnWeightsAt( p, PS_D, pp );
                    weights = queryDescription->GetSqnWeightsAt( p, PS_D );
                    if( queryDescription->GetLeftMSExtentAt( p ) != queryDescription->GetLeftMSExtentAt( pp ) ||
                        queryDescription->GetRightMSExtentAt( p ) != queryDescription->GetRightMSExtentAt( pp ))
                    {
                        newset = true;
                        weights = NULL;
                    }
                }
            }
            if( !newset && weights == NULL ) {
//                 newset = true;
                sprintf( errbuf, "Unexpected null weights for pos. %d.", pp );
                merror = errbuf;
                break;
            }
            if( !newset ) {
                if( p < prep + 1 ) {
                    merror = "Invalid index of previous position.";
                    break;
                }
                pp = prep;
                for( nadiffsyms = 0; nadiffsyms < noress; nadiffsyms++ ) {
                    histval = queryDescription->GetDistinctHistAt( pp, PS_D, nadiffsyms );
                    queryDescription->SetDistinctHistAt( p, PS_D, nadiffsyms, histval );
                }
            }
        }
        if( newset ) {
            prep = p;
            extentIdN = 0;
            wghtsum = 0.0;
            queryDescription->NewSqnWeightsAt( p, PS_D, size());
            weights = queryDescription->GetSqnWeightsAt( p, PS_D );

            if( weights == NULL ) {
                merror = "Null weights.";
                break;
            }

            size_t  left = queryDescription->GetLeftMSExtentAt( p );
            size_t  right = queryDescription->GetRightMSExtentAt( p );

            for( k = left; 0 <=( ssize_t )right && k <= right; k++ ) {
                // omit unsued and unmatched positions in the extent
                if( !queryDescription->IsUsedAt( k ))
                    continue;
//                 if( !queryDescription->GetStateAt( k ))
//                     continue;

                diffsyms = nadiffsyms = 0;
                memset( column, 0, sizeof( size_t ) * NUMALPH );

                for( i = 0; i < size(); i++ ) {
                    // omit sequences not in the extent
                    if( !SequenceAt( i )->IsUsedInExtentAt( p ))
                        continue;
                    pres = SequenceAt( i )->ResidueAt( p );
                    residue = SequenceAt( i )->ResidueAt( k );
                    // leave only sequences which are in delete state at position `p'
                    if( pres != GAP )
                        continue;
                    if( column[residue]++ == 0 ) {
                        diffsyms++;
                        if( IsValidResSym( residue ))
                            nadiffsyms++;
                    }
                }

                if( noeffress < nadiffsyms )
                    nadiffsyms = noeffress;

                //update histogram only in delete state
                if( queryDescription->GetStateAt( k )) {
                    extentIdN += nadiffsyms;
                    queryDescription->IncDistinctHistAt( p, PS_D, nadiffsyms );
                }

                if( nadiffsyms ) {
                    //if gaps are present in a column
                    for( i = 0; i < size(); i++ ) {
                        // omit sequences not in the extent
                        if( !SequenceAt( i )->IsUsedInExtentAt( p ))
                            continue;
                        pres = SequenceAt( i )->ResidueAt( p );
                        residue = SequenceAt( i )->ResidueAt( k );
                        // use only sequences which are in delete state at position `p'
                        if( pres != GAP )
                            continue;
                        if( !IsValidResSym( residue ))
                            continue;
                        w = 1.0 /( double )( column[residue] * nadiffsyms );
                        weights[i] += w;
                        wghtsum += w;
                    }
                }
            }

            //normalize if applicable
            if( wghtsum )
                //if gaps are present in a column, normalize;
                //otherwise, all weights will be 0
                for( i = 0; i < size(); i++ ) {
                    if( !SequenceAt( i )->IsUsedInExtentAt( p ))
                        continue;
                    pres = SequenceAt( i )->ResidueAt( p );
                    if( pres != GAP )
                        continue;
                    weights[i] /= wghtsum;
                }
        }

        queryDescription->SetNoSymbolsInExtentAt( p, PS_D, extentIdN );
    }

    //check for error
    if( !merror.empty())
        throw myruntime_error( preamble + merror );
}





// =========================================================================
// ComputeGlobSequenceWeights: calculate globally sequence weights
//  avgpeseq, average fraction of distinct residue per one matched position
//
void InputMultipleAlignment::ComputeGlobSequenceWeights( double* pavgsq )
{
    if( !queryDescription )
        throw myruntime_error("ComputeGlobSequenceWeights: Null description data.");
    if( !size())
        throw myruntime_error("ComputeGlobSequenceWeights: No data.");

    const size_t    MINUSDLEN = 30;
    const size_t    noress = NUMALPH;
    const size_t    noeffress = NUMAA;

    const int       noseqs = size();
    const size_t    nomstates = queryDescription->GetEffectiveSize();
    const size_t    left = 0;
    const size_t    right = nomstates;
    const size_t    interval = 1;//no avgs: histograms of one column

    size_t  mss;//number of sequences in M state
    size_t  column[NUMALPH];
    size_t  diffsyms;
    size_t  nadiffsyms;
    size_t  usdlen[noseqs];//length used in weight calculation
    double  avgpeseq;//average fraction of distinct residue per one matched position
    double  w;
    double  gwtsum;//sum of global weights
    double  gwghts[noseqs];

    unsigned char   residue;
    int             pstate;
    size_t  p, i; //position indices


    gwtsum = 0.0;
    avgpeseq = 0.0;
    memset( usdlen, 0, sizeof( size_t ) * noseqs );
    memset( gwghts, 0, sizeof( double ) * noseqs );

    //calculate global weights
    for( p = 0; p < queryDescription->size(); p++ ) {
        // omit unsued positions
        if( !queryDescription->IsUsedAt( p ))
            continue;

        //reset boundaries
        queryDescription->SetLeftMSExtentAt( left, p );
        queryDescription->SetRightMSExtentAt( right, p );
        queryDescription->SetMSExtentIntervalAt( interval, p );

        pstate = queryDescription->GetStateAt( p );
        if( !pstate )
            continue;

        mss = 0;
        diffsyms = 0;
        nadiffsyms = 0;
        memset( column, 0, sizeof( size_t ) * NUMALPH );

        for( i = 0; i < noseqs; i++ ) {
            // check for valid sequence position
            if( !SequenceAt(i)->GetUsed())
                continue;
            if( !SequenceAt( i )->IsUsedAndInfAt( p ))
                continue;
            residue = SequenceAt( i )->ResidueAt( p );
            if( !gbTOLXs && residue == X )
                continue;
            if( IsValidResSym( residue )) {
                usdlen[i]++;
                mss++;
            }
            if( column[residue]++ == 0 ) {
                diffsyms++;
                if( IsValidResSym( residue ))
                    nadiffsyms++;
            }
        }

        if( noeffress < nadiffsyms )
            nadiffsyms = noeffress;

        if( pstate )
            if( mss && nadiffsyms )
                //fraction of distinct residue per one matched position
                avgpeseq += ( double )nadiffsyms/( double )mss;

        if( nadiffsyms )
            for( i = 0; i < noseqs; i++ ) {
                // check for valid sequence position
                if( !SequenceAt(i)->GetUsed())
                    continue;
                if( !SequenceAt( i )->IsUsedAndInfAt( p ))
                    continue;
                residue = SequenceAt( i )->ResidueAt( p );
                if( !gbTOLXs && residue == X )
                    continue;
                if( IsValidResSym( residue )) {
                    w = 1.0 /( double )( column[residue] * nadiffsyms );
                    gwghts[i] += w;
//                     gwtsum += w;
                }
            }
    }//1:for(p...)

    //average no. distinct residues per one matched position
    if( avgpeseq && nomstates )
        avgpeseq /= ( double )nomstates;
    if( pavgsq )
        *pavgsq = avgpeseq;

    //adjust weights by sequence lengths;
    // each column has total weight =1
    for( i = 0; i < noseqs; i++ )
        if( gwghts[i]) {
            gwghts[i] /= ( double )SLC_MAX( usdlen[i] + 1, MINUSDLEN );
            gwtsum += gwghts[i];
        }
    if( !gwtsum )
        throw myruntime_error("ComputeGlobSequenceWeights: Null weights.");
    //normalize weights
    for( i = 0; i < noseqs; i++ ) {
        if( gwghts[i])
            gwghts[i] /= gwtsum;
        if( SequenceAt(i))
            SequenceAt(i)->SetGlbWeight( gwghts[i]);
    }
}

// -------------------------------------------------------------------------
// CalculateExpNoResiduesAt: calculate expected number of distinct 
//  residues at a given position over the extent
//
void InputMultipleAlignment::CalculateExpNoResiduesAt( 
    size_t p, bool usdgwght, size_t left, size_t right, 
    const double* wghts, double* express, size_t scale )
{
    if( !queryDescription || !wghts || !express )
        throw myruntime_error("CalculateExpNoResiduesAt: Null addresses.");

    const bool      bMED = true;
    const int       noseqs = size();
    const size_t    noress = NUMALPH;
    const size_t    noeffress = NUMAA;
    const size_t    rscale = SLC_MIN( 1000, SLC_MAX( 10, scale ));
    const size_t    maxres = 20;
    const size_t    scsz = maxres * rscale;
    size_t          expnoress[scsz+1];
    size_t          nosinext;//no. sequences in extent
    double          wghtsum[PS_NSTATES] = {0.0,0.0,0.0};
    double          mtwghts[NUMALPH];//match weights
    double          w, ww = 0.0, wwsum = 0.0;
    size_t          mss;//number of sequences in M state
    unsigned char   residue;
    int             kstate;
    size_t          k, i, r, n;

    if( usdgwght ) {
        left = 0;
        right = queryDescription->size()-1;
    }

    memset( expnoress, 0, ( scsz+1 ) * sizeof( size_t ));
    nosinext = usdgwght?queryDescription->GetNoSequencesAt( p ):
                        queryDescription->GetNoSequencesInMSExtentAt( p );

    //iterate over the extent
    for( n = 0, k = left; 0 <=( ssize_t )right && k <= right; k++ ) {
        if( !queryDescription->IsUsedAt(k))
            continue;
        kstate = queryDescription->GetStateAt(k);
        if( !kstate )
            continue;
        //use sequence weights calculated for this position
        mss = 0;
        wghtsum[PS_M] = 0.0;
        memset( mtwghts, 0, sizeof( double ) * NUMALPH );
        for( i = 0; i < noseqs; i++ ) {
            // omit sequences not in the extent
            if( usdgwght? !SequenceAt(i)->IsUsedAndInfAt( p ):
                          !SequenceAt(i)->IsUsedInExtentAt( p ))
                continue;
            residue = SequenceAt(i)->ResidueAt( k );
            if( !gbTOLXs && residue == X )
                continue;
            if( IsValidResSym( residue )) {
                mtwghts[residue] += wghts[i];
                wghtsum[PS_M] += wghts[i];
                mss++;
            }
        }
        if( mss < ( size_t )(( double )nosinext * gszMINPRCSEQNS ))
            continue;
        if( wghtsum[PS_M]) {
            for( r = 0; r < noress; r++ )
                if( mtwghts[r])
                    mtwghts[r] /= wghtsum[PS_M];
            AdjustWeightsAt(( size_t )-1, &mtwghts );//adjust in `mtwghts'
            ww = 0.0;
            for( r = 0; r < noeffress; r++ ) {
                w = mtwghts[r];
                if( 0.0 < w )
                    ww -= w * log(w);
            }
            ww = exp( ww );//value in [1,20]
            if( bMED )
                expnoress[(size_t)rint( SLC_MIN( maxres, ww ))*rscale]++;
            else {
                wwsum += ww;
                n++; 
            }
        }
    }//for(k...)
    if( bMED )
        MedianValue( expnoress, scsz, &ww, rscale );
    else if( n )
        ww = wwsum /( double )n;
    *express = ww;
}

// -------------------------------------------------------------------------
// ComputeMIDstateSequenceWeights: calculate MID state sequence weights in 
//  one run
//
void InputMultipleAlignment::ComputeMIDstateSequenceWeights()
{
    if( !queryDescription )
        throw myruntime_error("ComputeMstateSequenceWeights: Null description data.");

    const bool      bLGWEIGHTS = true;//use locally global weights
    const int       noseqs = size();
    const size_t    nomstates = queryDescription->GetEffectiveSize();
    const size_t    noress = NUMALPH;
    const size_t    noeffress = NUMAA;
    const double    err = 1.e-6;

    size_t  histval;//histogram data value
    size_t  nosinext;//no. sequences in extent
    bool    newset = true;
    size_t  Ihist[NUMALPH];
    size_t  extentIdN[PS_NSTATES] = {0,0,0};
    size_t  column[PS_NSTATES][NUMALPH];
    size_t  diffsyms[PS_NSTATES] = {0,0,0};
    size_t  nadiffsyms[PS_NSTATES] = {0,0,0};
    bool    gwgflgs[nomstates];//flags of usage of global weights
    double  gexprss = 0.0;//expected number of distinct residues computed globally
    double  expnors[nomstates];//recorded expected number of distinct residues
    double  gwghts[noseqs];
    double  wghtsum[PS_NSTATES] = {0.0,0.0,0.0};
    double* weights[PS_NSTATES] = {NULL,NULL,NULL};
    bool    usdgwght = true;

    unsigned char   press[noseqs];

    double          w, ww;
    double  emss, eiss, edss;//expected numbers of observations
    unsigned char   residue;
    int             pstate, kstate;
    char            errbuf[KBYTE];
    const mystring  preamble = "ComputeMIDstateSequenceWeights: ";
    mystring        merror;
    size_t  mss;//number of sequences in M state
    size_t  iss, ilen;//number of sequences and length in I state
    size_t  dss;//number of sequences in D state
    size_t  p, pM, k, i, r, rdif; //position indices
    ssize_t ppM, pp, ppp;

    memset( Ihist, 0, sizeof( size_t ) * NUMALPH );
    memset( expnors, 0, sizeof( double ) * nomstates );
    memset( gwgflgs, 0, sizeof( bool ) * nomstates );


    for( i = 0; i < noseqs; i++ ) {
        if( SequenceAt(i))
            gwghts[i] = SequenceAt(i)->GetGlbWeight();
        else
            gwghts[i] = 0.0;
    }

//     //{{calculate exp. no. residues along entire length
//     for( p = 0, pp = 0; p < queryDescription->size(); p++ ) {
//         if( !queryDescription->GetStateAt(p))
//             continue;
//         if( nomstates>>1 <= ++pp )
//             break;
//     }
//     CalculateExpNoResiduesAt( p, true, 0, queryDescription->size()-1, gwghts, &gexprss, 100 );
//     //}}

    //1:calc. sequence weights, match weights, expected no. residues
    for( p = 0, pp = ppp = -1, pM = ppM = 0; p < queryDescription->size() && merror.empty(); p++ ) {
        // omit unsued positions
        if( !queryDescription->IsUsedAt( p ))
            continue;

        pstate = queryDescription->GetStateAt( p );
        if( !pstate )
            continue;

        //{{record residues
        for( i = 0; i < noseqs; i++ ) {
            // omit sequences not in the extent
            if( !SequenceAt(i)->IsUsedInExtentAt( p ))
                    press[i] = GAP;
            else    press[i] = SequenceAt(i)->ResidueAt( p );
        }
//             if( pM )
//                 //save accumulated histogram data of I state
//                 for( rdif = 0; rdif < noress; rdif++ ) {
//                     queryDescription->SetDistinctHistAt( pM, PS_I, rdif, Ihist[rdif]);
//                     Ihist[rdif] = 0;
//                 }
        pM = p;
        pp++;
        if( nomstates <= pp ) {
            sprintf( errbuf, "Inconsistent number of match positions.");
            merror = errbuf;
            break;
        }
        //}}
        // determine if it is the same set of sequences;
        // the extent will always be the same if sequences comprise the same set!
        newset = true;
        if( p ) {
            if( p <= ppM ) {
                sprintf( errbuf, "Invalid position index at pos. %d.", p );
                merror = errbuf;
                break;
            }

            for( i = 0; i < noseqs; i++ ) {
                if( SequenceAt(i)->IsUsedInExtentAt( p ) !=
                    SequenceAt(i)->IsUsedInExtentAt( ppM ))
                        break;
                residue = SequenceAt(i)->ResidueAt( ppM );
                if( IsValidResSym( press[i]) != IsValidResSym( residue ))
                    break;
                if( press[i] != residue &&( press[i] == X || residue == X ))
                    break;
            }
            if( i == noseqs ) {
                newset = false;
                queryDescription->SetSqnWeightsAt( p, PS_M, ppM );
//                 queryDescription->SetSqnWeightsAt( p, PS_D, ppM );
//                 queryDescription->SetSqnWeightsAt( p, PS_I, ppM );
                weights[PS_M] = queryDescription->GetSqnWeightsAt( p, PS_M );
//                 weights[PS_D] = queryDescription->GetSqnWeightsAt( p, PS_D );
//                 weights[PS_I] = queryDescription->GetSqnWeightsAt( p, PS_I );
                if( queryDescription->GetLeftMSExtentAt( p ) != queryDescription->GetLeftMSExtentAt( ppM ) ||
                    queryDescription->GetRightMSExtentAt( p ) != queryDescription->GetRightMSExtentAt( ppM ))
                {
                    newset = true;
                    weights[PS_M] = weights[PS_I] = weights[PS_D] = NULL;
                }
            }

            if( !newset && !( weights[PS_M] /*&& weights[PS_D] && weights[PS_I]*/)) {
//                 newset = true;
                sprintf( errbuf, "Unexpected null weights for pos. %d.", ppM );
                merror = errbuf;
                break;
            }
        }

        if( newset ) {
            ppM = p;
            ppp = pp;

            extentIdN[PS_M] = extentIdN[PS_I] = extentIdN[PS_D] = 0;
            wghtsum[PS_M] = wghtsum[PS_I] = wghtsum[PS_D] = 0.0;

            queryDescription->NewSqnWeightsAt( p, PS_M, size());
//             queryDescription->NewSqnWeightsAt( p, PS_D, size());
//             queryDescription->NewSqnWeightsAt( p, PS_I, size());
            weights[PS_M] = queryDescription->GetSqnWeightsAt( p, PS_M );
//             weights[PS_D] = queryDescription->GetSqnWeightsAt( p, PS_D );
//             weights[PS_I] = queryDescription->GetSqnWeightsAt( p, PS_I );
            if( !weights[PS_M] /*|| !weights[PS_D] || !weights[PS_I]*/) {
                merror = "Null weights.";
                break;
            }

            size_t  left = queryDescription->GetLeftMSExtentAt( p );
            size_t  right = queryDescription->GetRightMSExtentAt( p );
            size_t  interval = queryDescription->GetMSExtentIntervalAt( p );
            size_t  nusdposs = 0;

            usdgwght = interval < gszMINEXTINTRV;

            //{{calculate sequence weights
            if( usdgwght ) {
                gwgflgs[pp] = true;
                for( i = 0; i < noseqs; i++ ) {
                      weights[PS_M][i] = gwghts[i];
//                       weights[PS_I][i] = gwghts[i];
//                       weights[PS_D][i] = gwghts[i];
                  }
            }
            else {
                nusdposs = 0;
                nosinext = queryDescription->GetNoSequencesInMSExtentAt( p );
                for( k = left; 0 <=( ssize_t )right && k <= right; k++ ) {
                    //omit unsued and unmatched positions in the extent
                    if( !queryDescription->IsUsedAt( k ))
                        continue;

                    kstate = queryDescription->GetStateAt( k );
                    if( !kstate )
                        continue;

                    mss = 0;
                    diffsyms[PS_M] = diffsyms[PS_I] = diffsyms[PS_D] = 0;
                    nadiffsyms[PS_M] = nadiffsyms[PS_I] = nadiffsyms[PS_D] = 0;

                    memset( column[PS_M], 0, sizeof( size_t ) * NUMALPH );
                    memset( column[PS_I], 0, sizeof( size_t ) * NUMALPH );
                    memset( column[PS_D], 0, sizeof( size_t ) * NUMALPH );

                    for( i = 0; i < noseqs; i++ ) {
                        //omit sequences not in the extent
                        if( !SequenceAt(i)->IsUsedInExtentAt( p ))
                            continue;
                        residue = SequenceAt(i)->ResidueAt( k );
                        if( !gbTOLXs && residue == X )
                            continue;
                        if( IsValidResSym( residue ))
                            mss++;
                        if( column[PS_M][residue]++ == 0 ) {
                            diffsyms[PS_M]++;
                            if( IsValidResSym( residue ))
                                nadiffsyms[PS_M]++;
                        }
                    }

                    if( mss < ( size_t )(( double )nosinext * gszMINPRCSEQNS ))
                        continue;

                    nusdposs++;
                    if( noeffress < nadiffsyms[PS_M] ) nadiffsyms[PS_M] = noeffress;
//                     if( noeffress < nadiffsyms[PS_D] ) nadiffsyms[PS_D] = noeffress;
//                     if( noeffress < nadiffsyms[PS_I] ) nadiffsyms[PS_I] = noeffress;

                    //update histogram
                    extentIdN[PS_M] += nadiffsyms[PS_M];
//                     extentIdN[PS_D] += nadiffsyms[PS_D];
//                     extentIdN[PS_I] += nadiffsyms[PS_I];
    //                 queryDescription->IncDistinctHistAt( p, PS_M, nadiffsyms[PS_M]);
    //                 queryDescription->IncDistinctHistAt( p, PS_D, nadiffsyms[PS_D]);
    //                 queryDescription->IncDistinctHistAt( p, PS_I, nadiffsyms[PS_I]);

                    if( nadiffsyms[PS_M]) {
                        for( i = 0; i < noseqs; i++ ) {
                            // omit sequences not in the extent
                            if( !SequenceAt(i)->IsUsedInExtentAt( p ))
                                continue;
                            residue = SequenceAt(i)->ResidueAt( k );
                            if( !gbTOLXs && residue == X )
                                continue;
                            if( IsValidResSym( residue )) {
                                w = 1.0 /( double )( column[PS_M][residue] * nadiffsyms[PS_M]);
                                weights[PS_M][i] += w;
                                wghtsum[PS_M] += w;
                            }
                        }
                    }
                }//extent: for(k...)

                if( !wghtsum[PS_M] && gbTOLXs ) {
                    sprintf( errbuf, "No contribution to sequence weights at pos. %d.", p );
                    merror = errbuf;
                    break;
                }
                //check again
                usdgwght = nusdposs < gszMINEXTINTRV;
                //copy sequence weights
                if( usdgwght ) {
                    gwgflgs[pp] = true;
                    for( i = 0; i < noseqs; i++ ) {
                          weights[PS_M][i] = gwghts[i];
//                           weights[PS_I][i] = gwghts[i];
//                           weights[PS_D][i] = gwghts[i];
                    }
                } else if( wghtsum[PS_M])
                    for( i = 0; i < noseqs; i++ ) {
                        if( weights[PS_M][i])
                            weights[PS_M][i] /= wghtsum[PS_M];
                    }
            }
            //}}
            //{{calculate expected number of distinct residues
            if( bLGWEIGHTS )
//                 ww = gexprss;
                CalculateExpNoResiduesAt( p, true, left, right, gwghts, &ww, 10 );
            else
                CalculateExpNoResiduesAt( p, usdgwght, left, right, weights[PS_M], &ww, 10 );
            expnors[pp] = SLC_MAX( 1.0, ww );
            //}}
        }//newset

//         if( newset && !pstate )
//             //accumulate histogram data for I state
//             for( rdif = 0; rdif < noress; rdif++ ) {
//                 histval = queryDescription->GetDistinctHistAt( p, PS_I, rdif );
//                 Ihist[rdif] += histval;
//             }

        if( !newset ) {
            //copy expected number of residues
            expnors[pp] = expnors[ppp];
            gwgflgs[pp] = gwgflgs[ppp];
//             //copy histogram data
//             for( rdif = 0; rdif < noress; rdif++ ) {
//                 histval = queryDescription->GetDistinctHistAt( ppM, PS_M, rdif );
//                 queryDescription->SetDistinctHistAt( p, PS_M, rdif, histval );
//                 histval = queryDescription->GetDistinctHistAt( ppM, PS_D, rdif );
//                 queryDescription->SetDistinctHistAt( p, PS_D, rdif, histval );
//                 histval = queryDescription->GetDistinctHistAt( ppM, PS_I, rdif );
//                 queryDescription->SetDistinctHistAt( p, PS_I, rdif, histval );
//             }
        }

        //{{set match weights
        wghtsum[PS_M] = 0.0;
        for( i = 0; i < noseqs; i++ ) {
            if( usdgwght? !SequenceAt(i)->IsUsedAndInfAt( p ):
                          !SequenceAt(i)->IsUsedInExtentAt( p ))
                continue;
            if( !gbTOLXs && press[i] == X )
                continue;
            if( IsValidResSym( press[i] )) {
                wghtsum[PS_M] += weights[PS_M][i];
                queryDescription->IncMatchWeightsAt( weights[PS_M][i], press[i], p );
            }
        }
        if( !wghtsum[PS_M] && gbTOLXs ) {
            sprintf( errbuf, "No contribution to sequence weights at pos. %d.", p );
            merror = errbuf;
            break;
        }
        if( wghtsum[PS_M])
            for( r = 0; r < noress; r++ ) {
                w = queryDescription->GetMatchWeightsAt( r, p );
                if( w )
                    queryDescription->SetMatchWeightsAt( w / wghtsum[PS_M], r, p );
            }
        AdjustWeightsAt( p );
        //}}
    }//1:for(p...)

    iss = ilen = 0;
    eiss = 0.0;
    pM = -1;

    //set expectations at the beginning position
    queryDescription->SetMIDExpNoObservationsAt( 1.0, -1, PS_M );
    queryDescription->SetMIDExpNoObservationsAt( 0.0, -1, PS_D );
    weights[PS_M] = weights[PS_D] = weights[PS_I] = NULL;

    //2:calc. accumulated weights for each state, expected no. observations
    for( p = 0, pp = 0; p < queryDescription->size() && merror.empty(); p++ ) {
        // omit unsued positions
        if( !queryDescription->IsUsedAt( p ))
            continue;

        pstate = queryDescription->GetStateAt( p );

        if( pstate )
            weights[PS_M] = queryDescription->GetSqnWeightsAt( p, PS_M );

        mss = dss = 0;
        wghtsum[PS_M] = wghtsum[PS_I] = wghtsum[PS_D] = 0.0;

        //calc. accumulated weights for each state
        for( i = 0; i < noseqs; i++ ) {
            if( !SequenceAt(i)->GetUsed())
                continue;
            if( !SequenceAt(i)->IsUsedInExtentAt( p ))
                continue;
            //
            residue = SequenceAt(i)->ResidueAt( p );
            if( !gbTOLXs && residue == X )
                continue;
            if( pstate ) {
                if( IsValidResSym( residue )) {
                    wghtsum[PS_M] += gwghts[i];//weights[PS_M][i];
                    mss++;
                }
                else if( residue == GAP ) {
                    wghtsum[PS_D] += gwghts[i];
                    dss++;
                }
            } else
                if( IsValidResSym( residue )) {
                    wghtsum[PS_I] += gwghts[i];
                    eiss += gwghts[i];
                    iss++;
                }
        }

        //calc. expected no. observations (set histogram data)
        if( pstate ) {
            if( !wghtsum[PS_M]) {
                if( !gbTOLXs )
                    wghtsum[PS_M] = 1.0;
                else {
                    sprintf( errbuf, "No contribution to match weights at pos. %d.", p );
                    merror = errbuf;
                    break;
                }
            }

            if( iss && ilen ) {
                //average observations
                eiss /= ( double )iss;
                iss = ( size_t )rint(( double )iss /( double )ilen );
            }

            if( wghtsum[PS_M] < 0.0 || 1.0+err < wghtsum[PS_M] || eiss < 0.0 || 1.0+err < eiss ||
                wghtsum[PS_D] < 0.0 || 1.0+err < wghtsum[PS_D] ) {
                sprintf( errbuf, "Invalid MID state weights at pos. %d.", p );
                merror = errbuf;
                break;
            }

//             emss = SLC_MIN(( double )mss * avgpeseq, ( double )noeffress );
//             eiss = SLC_MIN(( double )iss * avgpeseq, ( double )noeffress );
//             edss = SLC_MIN(( double )dss * avgpeseq, ( double )noeffress );

//             emss = SLC_MIN( expnors[pp], ( double )noeffress );
//             eiss = SLC_MIN( expnors[pp] * eiss / wghtsum[PS_M], ( double )noeffress );
//             edss = SLC_MIN( expnors[pp] * wghtsum[PS_D] / wghtsum[PS_M], ( double )noeffress );

            //convex function to avoid of rapid drop of expectation as weight decreases
            if( bLGWEIGHTS || gwgflgs[pp])
                emss = SLC_MIN(expnors[pp]* exp( gdWCPC*(1.0-wghtsum[PS_M])*log(wghtsum[PS_M])),( double )noeffress);
            else
                emss = SLC_MIN(expnors[pp], ( double )noeffress );
            if( eiss )
                eiss = SLC_MIN(expnors[pp]* exp(gdWCPC*(1.0-eiss)*log(eiss)),( double )noeffress);
            if( wghtsum[PS_D])
                edss = SLC_MIN(expnors[pp]* exp(gdWCPC*(1.0-wghtsum[PS_D])*log(wghtsum[PS_D])),(double)noeffress);
            else
                edss = 0.0;

            //M state no. diff residues is required
            queryDescription->SetNoSymbolsInExtentAt( p, PS_M, SLC_MAX( 1.0, emss ));

///             emss = ExtendedDescriptionVector::GetExpNoObservations( emss );
///             eiss = ExtendedDescriptionVector::GetExpNoObservations( eiss );
///             edss = ExtendedDescriptionVector::GetExpNoObservations( edss );

            if(( double )( mss + dss ) < emss ) emss = ( double )( mss + dss );
            if(( double )( mss + dss ) < eiss ) eiss = ( double )( mss + dss );
            if(( double )( mss + dss ) < edss ) edss = ( double )( mss + dss );

            if( emss < 1.0 ) emss = 1.0;
            if( eiss < 1.0 && eiss ) eiss = 1.0;
            if( edss < 1.0 && edss ) edss = 1.0;

            //set expectations directly, avoiding using of histogram data
            queryDescription->SetMIDExpNoObservationsAt( emss, p, PS_M );
            queryDescription->SetMIDExpNoObservationsAt( eiss, pM, PS_I );
            queryDescription->SetMIDExpNoObservationsAt( edss, p, PS_D );
//             queryDescription->IncDistinctHistAt( p, PS_M, mss );
//             queryDescription->IncDistinctHistAt( pM, PS_I, iss );
//             queryDescription->IncDistinctHistAt( p, PS_D, dss );
            iss = ilen = 0;
            eiss = 0.0;
            pM = p;
            pp++;
        }
        else
            //increase insertion length
            ilen++;

    }//2:for(p...)

    if( iss && ilen ) {
        //average observations
        eiss /= ( double )iss;
        iss = ( size_t )rint(( double )iss /( double )ilen );
        if( /*wghtsum[PS_M] &&*/ eiss && nomstates ) {
//             eiss = SLC_MIN( expnors[nomstates-1] * eiss / wghtsum[PS_M], ( double )noeffress );
            eiss = SLC_MIN(expnors[nomstates-1]* exp( gdWCPC*(1.0-eiss)*log(eiss)),( double )noeffress);
///             eiss = ExtendedDescriptionVector::GetExpNoObservations( eiss );
            if(( double )( mss + dss ) < eiss ) eiss = ( double )( mss + dss );
            queryDescription->SetMIDExpNoObservationsAt( eiss, pM, PS_I );
        }
//         queryDescription->IncDistinctHistAt( pM, PS_I, ( int )rint(( double )iss * avgpeseq ));
    }

    //check for error
    if( !merror.empty())
        throw myruntime_error( preamble + merror );
}





// -------------------------------------------------------------------------
// ComputeMIDstateSequenceWeightsNoExtents: calculate MID state sequence 
//  weights with no using of extents
//
void InputMultipleAlignment::ComputeMIDstateSequenceWeightsNoExtents()
{
    if( !queryDescription )
        throw myruntime_error("ComputeMstateSequenceWeights: Null description data.");

    const int       noseqs = size();
    const size_t    nomstates = queryDescription->GetEffectiveSize();
    const size_t    noress = NUMALPH;
    const size_t    noeffress = NUMAA;
    const double    err = 1.e-6;

    size_t  histval;//histogram data value
    size_t  nosinext;//no. sequences in extent
    bool    newset = true;
    size_t  extentIdN[PS_NSTATES] = {0,0,0};
    size_t  column[PS_NSTATES][NUMALPH];
    size_t  diffsyms[PS_NSTATES] = {0,0,0};
    size_t  nadiffsyms[PS_NSTATES] = {0,0,0};
    double  expnors[nomstates];//expected number of residues
    double  gwghts[noseqs];
    double  wghtsum[PS_NSTATES] = {0.0,0.0,0.0};
    double* weights[PS_NSTATES] = {NULL,NULL,NULL};
    double  mtwghts[NUMALPH];//match weights

    unsigned char   press[noseqs];

    double          w, ww;
    const size_t    rscale = 10;
    const size_t    maxres = 20;
    const size_t    scsz = maxres * rscale;
    size_t          expnoress[scsz+1];
    double  emss, eiss, edss;//expected numbers of observations
    unsigned char   residue;
    int             pstate, kstate;
    char            errbuf[KBYTE];
    const mystring  preamble = "ComputeMIDstateSequenceWeightsNoExtents: ";
    mystring        merror;
    size_t  mss;//number of sequences in M state
    size_t  iss, ilen;//number of sequences and length in I state
    size_t  dss;//number of sequences in D state
    size_t  p, pM, k, i, r, rdif; //position indices
    ssize_t ppM, pp, ppp;

    memset( expnors, 0, sizeof( double ) * nomstates );


    for( i = 0; i < noseqs; i++ ) {
        if( SequenceAt(i))
            gwghts[i] = SequenceAt(i)->GetGlbWeight();
        else
            gwghts[i] = 0.0;
    }

    //1:calc. sequence weights, match weights, expected no. residues
    for( p = 0, pp = ppp = -1, pM = ppM = 0; p < queryDescription->size() && merror.empty(); p++ ) {
        // omit unsued positions
        if( !queryDescription->IsUsedAt( p ))
            continue;

        pstate = queryDescription->GetStateAt( p );
        if( !pstate )
            continue;

        //{{record residues
        for( i = 0; i < noseqs; i++ ) {
            // omit sequences not in the extent
            if( !SequenceAt(i)->IsUsedAndInfAt( p ))
                    press[i] = GAP;
            else    press[i] = SequenceAt(i)->ResidueAt( p );
        }
        pM = p;
        pp++;
        if( nomstates <= pp ) {
            sprintf( errbuf, "Inconsistent number of match positions.");
            merror = errbuf;
            break;
        }
        //}}
        // determine if it is the same set of sequences;
        // the extent will always be the same if sequences comprise the same set!
        newset = true;
        if( p ) {
            if( p <= ppM ) {
                sprintf( errbuf, "Invalid position index at pos. %d.", p );
                merror = errbuf;
                break;
            }

            for( i = 0; i < noseqs; i++ ) {
                if( SequenceAt(i)->IsUsedAndInfAt( p ) != SequenceAt(i)->IsUsedAndInfAt( ppM ))
                        break;
                residue = SequenceAt(i)->ResidueAt( ppM );
                if( IsValidResSym( press[i]) != IsValidResSym( residue ))
                    break;
            }
            if( i == noseqs ) {
                newset = false;
                queryDescription->SetSqnWeightsAt( p, PS_M, ppM );
                weights[PS_M] = queryDescription->GetSqnWeightsAt( p, PS_M );
            }

            if( !newset && !weights[PS_M]) {
                sprintf( errbuf, "Unexpected null weights for pos. %d.", ppM );
                merror = errbuf;
                break;
            }
        }

        //NEWSET...
        if( newset ) {
            ppM = p;
            ppp = pp;

            extentIdN[PS_M] = extentIdN[PS_I] = extentIdN[PS_D] = 0;
            wghtsum[PS_M] = wghtsum[PS_I] = wghtsum[PS_D] = 0.0;

            queryDescription->NewSqnWeightsAt( p, PS_M, size());
            weights[PS_M] = queryDescription->GetSqnWeightsAt( p, PS_M );
            if( !weights[PS_M]) {
                merror = "Null weights.";
                break;
            }

            size_t  nusdposs = 0;
            bool    usdgwght = false;

            //{{calculate sequence weights
            nosinext = queryDescription->GetNoSequencesAt( p );
            for( k = 0; k < queryDescription->size(); k++ ) {
                //omit unsued and unmatched positions in the extent
                if( !queryDescription->IsUsedAt( k ))
                    continue;

                kstate = queryDescription->GetStateAt( k );
                if( !kstate )
                    continue;

                mss = 0;
                diffsyms[PS_M] = diffsyms[PS_I] = diffsyms[PS_D] = 0;
                nadiffsyms[PS_M] = nadiffsyms[PS_I] = nadiffsyms[PS_D] = 0;
                memset( column[PS_M], 0, sizeof( size_t ) * NUMALPH );

                for( i = 0; i < noseqs; i++ ) {
                    //omit sequences not in the extent
                    if( !SequenceAt(i)->IsUsedAndInfAt( p ))
                        continue;
                    residue = SequenceAt(i)->ResidueAt( k );
                    if( !gbTOLXs && residue == X )
                        continue;
                    if( IsValidResSym( residue ))
                        mss++;
                    if( column[PS_M][residue]++ == 0 ) {
                        diffsyms[PS_M]++;
                        if( IsValidResSym( residue ))
                            nadiffsyms[PS_M]++;
                    }
                }

                if( mss < ( size_t )(( double )nosinext * gszMINPRCSEQNS ))
                    continue;

                nusdposs++;
                if( noeffress < nadiffsyms[PS_M] ) nadiffsyms[PS_M] = noeffress;

                extentIdN[PS_M] += nadiffsyms[PS_M];

                if( nadiffsyms[PS_M]) {
                    for( i = 0; i < noseqs; i++ ) {
                        // omit sequences not in the extent
                        if( !SequenceAt(i)->IsUsedAndInfAt( p ))
                            continue;
                        residue = SequenceAt(i)->ResidueAt( k );
                        if( !gbTOLXs && residue == X )
                            continue;
                        if( IsValidResSym( residue )) {
                            w = 1.0 /( double )( column[PS_M][residue] * nadiffsyms[PS_M]);
                            weights[PS_M][i] += w;
                            wghtsum[PS_M] += w;
                        }
                    }
                }
            }//extent: for(k...)

            if( !wghtsum[PS_M] && gbTOLXs ) {
                sprintf( errbuf, "No contribution to sequence weights at pos. %d.", p );
                merror = errbuf;
                break;
            }
            //check again
            usdgwght = nusdposs < gszMINEXTINTRV;
            //copy sequence weights
            if( usdgwght )
                for( i = 0; i < noseqs; i++ ) {
                      weights[PS_M][i] = gwghts[i];
                }
            else if( wghtsum[PS_M])
                for( i = 0; i < noseqs; i++ ) {
                    if( weights[PS_M][i])
                        weights[PS_M][i] /= wghtsum[PS_M];
                }
            //}}
            //{{calculate expected number of distinct residues
            CalculateExpNoResiduesAt( p, true, 0, queryDescription->size()-1, weights[PS_M], &ww, 10 );
            expnors[pp] = SLC_MAX( 1.0, ww );
            //}}
        }//newset

        if( !newset ) {
            //copy expected number of residues
            expnors[pp] = expnors[ppp];
        }

        //{{set match weights
        wghtsum[PS_M] = 0.0;
        for( i = 0; i < noseqs; i++ ) {
            if( !SequenceAt(i)->IsUsedAndInfAt( p ))
                continue;
            if( !gbTOLXs && press[i] == X )
                continue;
            if( IsValidResSym( press[i] )) {
                wghtsum[PS_M] += weights[PS_M][i];
                queryDescription->IncMatchWeightsAt( weights[PS_M][i], press[i], p );
            }
        }
        if( !wghtsum[PS_M] && gbTOLXs ) {
            sprintf( errbuf, "No contribution to sequence weights at pos. %d.", p );
            merror = errbuf;
            break;
        }
        if( wghtsum[PS_M])
            for( r = 0; r < noress; r++ ) {
                w = queryDescription->GetMatchWeightsAt( r, p );
                if( w )
                    queryDescription->SetMatchWeightsAt( w / wghtsum[PS_M], r, p );
            }
        AdjustWeightsAt( p );
        //}}
    }//1:for(p...)

    iss = ilen = 0;
    eiss = 0.0;
    pM = -1;

    //set expectations at the beginning position
    queryDescription->SetMIDExpNoObservationsAt( 1.0, -1, PS_M );
    queryDescription->SetMIDExpNoObservationsAt( 0.0, -1, PS_D );
    weights[PS_M] = weights[PS_D] = weights[PS_I] = NULL;

    //2:calc. accumulated weights for each state, expected no. observations
    for( p = 0, pp = 0; p < queryDescription->size() && merror.empty(); p++ ) {
        // omit unsued positions
        if( !queryDescription->IsUsedAt( p ))
            continue;

        pstate = queryDescription->GetStateAt( p );

        if( pstate )
            weights[PS_M] = queryDescription->GetSqnWeightsAt( p, PS_M );

        mss = dss = 0;
        wghtsum[PS_M] = wghtsum[PS_I] = wghtsum[PS_D] = 0.0;

        //calc. accumulated weights for each state
        for( i = 0; i < noseqs; i++ ) {
            if( !SequenceAt(i)->GetUsed())
                continue;
            if( !SequenceAt(i)->IsUsedAndInfAt( p ))
                continue;
            //
            residue = SequenceAt(i)->ResidueAt( p );
            if( !gbTOLXs && residue == X )
                continue;
            if( pstate ) {
                if( IsValidResSym( residue )) {
                    wghtsum[PS_M] += gwghts[i];//weights[PS_M][i];
                    mss++;
                }
                else if( residue == GAP ) {
                    wghtsum[PS_D] += gwghts[i];
                    dss++;
                }
            } else
                if( IsValidResSym( residue )) {
                    wghtsum[PS_I] += gwghts[i];
                    eiss += gwghts[i];
                    iss++;
                }
        }

        //calc. expected no. observations (set histogram data)
        if( pstate ) {
            if( !wghtsum[PS_M]) {
                if( !gbTOLXs )
                    wghtsum[PS_M] = 1.0;
                else {
                    sprintf( errbuf, "No contribution to match weights at pos. %d.", p );
                    merror = errbuf;
                    break;
                }
            }

            if( iss && ilen ) {
                //average observations
                eiss /= ( double )iss;
                iss = ( size_t )rint(( double )iss /( double )ilen );
            }

            if( wghtsum[PS_M] < 0.0 || 1.0+err < wghtsum[PS_M] || eiss < 0.0 || 1.0+err < eiss ||
                wghtsum[PS_D] < 0.0 || 1.0+err < wghtsum[PS_D] ) {
                sprintf( errbuf, "Invalid MID state weights at pos. %d.", p );
                merror = errbuf;
                break;
            }

//             emss = SLC_MIN(( double )mss * avgpeseq, ( double )noeffress );
//             eiss = SLC_MIN(( double )iss * avgpeseq, ( double )noeffress );
//             edss = SLC_MIN(( double )dss * avgpeseq, ( double )noeffress );

//             emss = SLC_MIN( expnors[pp], ( double )noeffress );
//             eiss = SLC_MIN( expnors[pp] * eiss / wghtsum[PS_M], ( double )noeffress );
//             edss = SLC_MIN( expnors[pp] * wghtsum[PS_D] / wghtsum[PS_M], ( double )noeffress );

            //convex function to avoid of rapid drop of expectation as weight decreases
            emss = SLC_MIN(expnors[pp]* exp( gdWCPC*(1.0-wghtsum[PS_M])*log(wghtsum[PS_M])),( double )noeffress);
            if( eiss )
                eiss = SLC_MIN(expnors[pp]* exp(gdWCPC*(1.0-eiss)*log(eiss)),( double )noeffress);
            if( wghtsum[PS_D])
                edss = SLC_MIN(expnors[pp]* exp(gdWCPC*(1.0-wghtsum[PS_D])*log(wghtsum[PS_D])),(double)noeffress);
            else
                edss = 0.0;

            //M state no. diff residues is required
            queryDescription->SetNoSymbolsInExtentAt( p, PS_M, SLC_MAX( 1.0, emss ));

///             emss = ExtendedDescriptionVector::GetExpNoObservations( emss );
///             eiss = ExtendedDescriptionVector::GetExpNoObservations( eiss );
///             edss = ExtendedDescriptionVector::GetExpNoObservations( edss );

            if(( double )( mss + dss ) < emss ) emss = ( double )( mss + dss );
            if(( double )( mss + dss ) < eiss ) eiss = ( double )( mss + dss );
            if(( double )( mss + dss ) < edss ) edss = ( double )( mss + dss );

            if( emss < 1.0 ) emss = 1.0;
            if( eiss < 1.0 && eiss ) eiss = 1.0;
            if( edss < 1.0 && edss ) edss = 1.0;

            //set expectations directly, avoiding using of histogram data
            queryDescription->SetMIDExpNoObservationsAt( emss, p, PS_M );
            queryDescription->SetMIDExpNoObservationsAt( eiss, pM, PS_I );
            queryDescription->SetMIDExpNoObservationsAt( edss, p, PS_D );
//             queryDescription->IncDistinctHistAt( p, PS_M, mss );
//             queryDescription->IncDistinctHistAt( pM, PS_I, iss );
//             queryDescription->IncDistinctHistAt( p, PS_D, dss );
            iss = ilen = 0;
            eiss = 0.0;
            pM = p;
            pp++;
        }
        else
            //increase insertion length
            ilen++;

    }//2:for(p...)

    if( iss && ilen ) {
        //average observations
        eiss /= ( double )iss;
        iss = ( size_t )rint(( double )iss /( double )ilen );
        if( /*wghtsum[PS_M] &&*/ eiss && nomstates ) {
//             eiss = SLC_MIN( expnors[nomstates-1] * eiss / wghtsum[PS_M], ( double )noeffress );
            eiss = SLC_MIN(expnors[nomstates-1]* exp( gdWCPC*(1.0-eiss)*log(eiss)),( double )noeffress);
///             eiss = ExtendedDescriptionVector::GetExpNoObservations( eiss );
            if(( double )( mss + dss ) < eiss ) eiss = ( double )( mss + dss );
            queryDescription->SetMIDExpNoObservationsAt( eiss, pM, PS_I );
        }
//         queryDescription->IncDistinctHistAt( pM, PS_I, ( int )rint(( double )iss * avgpeseq ));
    }

    //check for error
    if( !merror.empty())
        throw myruntime_error( preamble + merror );
}





// =========================================================================
// ComputeGWMIDstateSequenceWeights: calculate globally sequence weights at 
//  all states
//  avgpeseq, average fraction of distinct residue per one matched position
//
void InputMultipleAlignment::ComputeGWMIDstateSequenceWeights( double avgpeseq )
{
    if( !queryDescription )
        throw myruntime_error("ComputeGWMIDstateSequenceWeights: Null description data.");
    if( !size())
        throw myruntime_error("ComputeGWMIDstateSequenceWeights: No data.");

    //pos. normalization curves weight distribution
    const bool      bPOSNORMAL = false;//positional normalization of seq. weights
    const size_t    NODRLEN = 40;//length for observed no. distinct residues
    const size_t    NODRLENh = NODRLEN >> 1;
    const size_t    noress = NUMALPH;
    const size_t    noeffress = NUMAA;

    const int       noseqs = size();
    const size_t    nomstates = queryDescription->GetEffectiveSize();
    const double    err = 1.e-6;

    size_t  mss;//number of sequences in M state
    size_t  iss, ilen;//number of sequences and length in I state
    size_t  dss;//number of sequences in D state
    size_t  diffsyms;
    size_t  nadiffsyms;
    double  emss, eiss, edss;//expected numbers of observations
    double  w, ww, ndr[nomstates], avglim, avgndr[nomstates];
    double  gwghts[noseqs];
    double  poswsum[PS_NSTATES];//positional weight sums
    double* weights[PS_NSTATES] = {NULL,NULL,NULL};


    unsigned char   residue;
    bool            sused, extused;
    int             pstate;
    char            errbuf[KBYTE];
    const mystring  preamble = "ComputeGWMIDstateSequenceWeights: ";
    mystring        merror;
    size_t  p, pp, i, r, n, nn;//position indices
    ssize_t pM;


    for( i = 0; i < noseqs; i++ ) {
        if( SequenceAt(i))
            gwghts[i] = SequenceAt(i)->GetGlbWeight();
        else
            gwghts[i] = 0.0;
    }

    //distribute global weights at all positions
    if( !bPOSNORMAL ) {
        queryDescription->NewSqnWeightsAt( 0, PS_M, size());
        queryDescription->NewSqnWeightsAt( 0, PS_I, size());
        queryDescription->NewSqnWeightsAt( 0, PS_D, size());
        weights[PS_M] = queryDescription->GetSqnWeightsAt( 0, PS_M );
        weights[PS_I] = queryDescription->GetSqnWeightsAt( 0, PS_I );
        weights[PS_D] = queryDescription->GetSqnWeightsAt( 0, PS_D );
        if( !weights[PS_M] || !weights[PS_I] || !weights[PS_D]) {
            merror = "Null weights.";
        }

        for( i = 0; i < noseqs; i++ ) {
            weights[PS_M][i] = gwghts[i];
            weights[PS_I][i] = gwghts[i];
            weights[PS_D][i] = gwghts[i];
        }
    }


    //check for error
    if( !merror.empty())
        throw myruntime_error( preamble + merror );

    memset( avgndr, 0, sizeof( double ) * nomstates );
    avglim = 0.0;
    nn = NODRLENh;

    //adjust recalculate weights, set MID weights
    for( p = 0, pp = 0; p < queryDescription->size() && merror.empty(); p++ ) {
        // omit unsued positions
        if( !queryDescription->IsUsedAt( p ))
            continue;

        pstate = queryDescription->GetStateAt( p );


        if( bPOSNORMAL ) {
            if( pstate ) {
                queryDescription->NewSqnWeightsAt( p, PS_M, size());
                queryDescription->NewSqnWeightsAt( p, PS_D, size());

                weights[PS_M] = queryDescription->GetSqnWeightsAt( p, PS_M );
                weights[PS_D] = queryDescription->GetSqnWeightsAt( p, PS_D );
                if( !weights[PS_M] || !weights[PS_D]) {
                    merror = "Null weights.";
                    break;
                }
            } else {
                queryDescription->NewSqnWeightsAt( p, PS_I, size());
                weights[PS_I] = queryDescription->GetSqnWeightsAt( p, PS_I );
                if( !weights[PS_I]) {
                    merror = "Null weights.";
                    break;
                }
            }
        }//not bPOSNORMAL
        else {
            if( p ) {
                queryDescription->SetSqnWeightsAt( p, PS_M, 0 );
                queryDescription->SetSqnWeightsAt( p, PS_I, 0 );
                queryDescription->SetSqnWeightsAt( p, PS_D, 0 );
            }
        }


        poswsum[PS_M] = poswsum[PS_I] = poswsum[PS_D] = 0;

        //sequence weights at each position are normalized by definition
        for( i = 0; i < noseqs; i++ ) {
            if( !SequenceAt(i)->GetUsed())
                continue;
            sused = SequenceAt( i )->IsUsedAndInfAt( p );
            extused = SequenceAt( i )->IsUsedInExtentAt( p );
            // check for valid sequence position and 
            //  modify belonging-to-extent flag
            if( !sused ) {
                if( extused )
                    SequenceAt( i )->UnsetUsedInExtentAt( p );
                continue;
            }
            if( !extused )
                SequenceAt( i )->SetUsedInExtentAt( p );
            //
            queryDescription->IncNoSequencesInMSExtentAt( p );
            //
            residue = SequenceAt( i )->ResidueAt( p );
            if( !gbTOLXs && residue == X )
                continue;
            if( pstate ) {
                if( IsValidResSym( residue )) {
                    poswsum[PS_M] += gwghts[i];
                }
                else if( residue == GAP )
                    if( bPOSNORMAL )
                        poswsum[PS_D] += gwghts[i];
            } else
                if( IsValidResSym( residue ))
                    if( bPOSNORMAL )
                        poswsum[PS_I] += gwghts[i];
        }

        if( pstate ) {
            if( !poswsum[PS_M]) {
                if( !gbTOLXs )
                    poswsum[PS_M] = 1.0;
                else {
                    sprintf( errbuf, "No contribution to match weights at pos. %d.", p );
                    merror = errbuf;
                    break;
                }
            }
            //set match weights
            for( i = 0; i < noseqs; i++ ) {
                if( !SequenceAt(i)->GetUsed())
                    continue;
                if( !SequenceAt( i )->IsUsedAndInfAt( p ))
                    continue;
                residue = SequenceAt( i )->ResidueAt( p );
                if( !gbTOLXs && residue == X )
                    continue;
                if( bPOSNORMAL ) {
                    if( IsValidResSym( residue ))
                        weights[PS_M][i] = gwghts[i] / poswsum[PS_M];
                    else if( residue == GAP && poswsum[PS_D])
                        weights[PS_D][i] = gwghts[i] / poswsum[PS_D];
                    if( IsValidResSym( residue ))
                        queryDescription->IncMatchWeightsAt( weights[PS_M][i], residue, p );
                } else // !bPOSNORMAL
                    if( IsValidResSym( residue ))
                        queryDescription->IncMatchWeightsAt( weights[PS_M][i] / poswsum[PS_M], residue, p );
            }
        } else if( bPOSNORMAL ) {
            if( poswsum[PS_I])
                for( i = 0; i < noseqs; i++ ) {
                    if( !SequenceAt(i)->GetUsed())
                        continue;
                    if( !SequenceAt( i )->IsUsedAndInfAt( p ))
                        continue;
                    residue = SequenceAt( i )->ResidueAt( p );
                    if( !gbTOLXs && residue == X )
                        continue;
                    if( IsValidResSym( residue ))
                        weights[PS_I][i] = gwghts[i] / poswsum[PS_I];
                }
        }

        //calculate observed number of distinct residues
        if( pstate ) {
            ww = 0.0;
            for( r = 0; r < noeffress; r++ ) {
                w = queryDescription->GetMatchWeightsAt( r, p );
                if( 0.0 < w )
                    ww -= w * log(w);
            }
            //exp(ww) varies between 1 and 20
            ww = exp( ww );
            ndr[pp] = ww;
            avglim += ww;
            pp++;
            if( NODRLEN <= pp ) {
                avgndr[nn] = avglim /( double )NODRLEN;
                if( NODRLEN == pp )
                    for( n = 0; n < nn; n++ )
                        avgndr[n] = avgndr[nn];
                if( nomstates == pp )
                    for( n = nn+1; n < nomstates; n++ )
                        avgndr[n] = avgndr[nn];
                avglim -= ndr[pp-NODRLEN];
                nn++;
            }
        }
    }//2:for(p...)

    if( pp < NODRLEN && nomstates ) {
        avgndr[0] = avglim /( double )nomstates;
        for( n = 1; n < nomstates; n++ )
            avgndr[n] = avgndr[0];
    }
    iss = ilen = 0;
    eiss = 0.0;
    pM = -1;

    //set expectations at the beginning state
    queryDescription->SetMIDExpNoObservationsAt( 1.0, -1, PS_M );
    queryDescription->SetMIDExpNoObservationsAt( 0.0, -1, PS_D );

    //set histogram data for MID states
    for( p = 0, pp = 0; p < queryDescription->size() && merror.empty(); p++ ) {
        // omit unsued positions
        if( !queryDescription->IsUsedAt( p ))
            continue;

        pstate = queryDescription->GetStateAt( p );

        mss = dss = 0;
        poswsum[PS_M] = poswsum[PS_I] = poswsum[PS_D] = 0;

        //sequence weights at each position are normalized by definition
        for( i = 0; i < noseqs; i++ ) {
            if( !SequenceAt(i)->GetUsed())
                continue;
            if( !SequenceAt( i )->IsUsedAndInfAt( p ))
                continue;
            residue = SequenceAt( i )->ResidueAt( p );
            if( !gbTOLXs && residue == X )
                continue;
            if( pstate ) {
                if( IsValidResSym( residue )) {
                    poswsum[PS_M] += gwghts[i];
                    mss++;
                }
                else if( residue == GAP ) {
                    poswsum[PS_D] += gwghts[i];
                    dss++;
                }
            } else
                if( IsValidResSym( residue )) {
                    poswsum[PS_I] += gwghts[i];
                    eiss += gwghts[i];
                    iss++;
                }
        }

        //set histogram data
        if( pstate ) {
            if( !poswsum[PS_M]) {
                if( !gbTOLXs )
                    poswsum[PS_M] = 1.0;
                else {
                    sprintf( errbuf, "No contribution to match weights at pos. %d.", p );
                    merror = errbuf;
                    break;
                }
            }

            if( iss && ilen ) {
                //average observations
                eiss /= ( double )iss;
                iss = ( size_t )rint(( double )iss /( double )ilen );
            }

            if( poswsum[PS_M] < 0.0 || 1.0+err < poswsum[PS_M] || eiss < 0.0 || 1.0+err < eiss ||
                poswsum[PS_D] < 0.0 || 1.0+err < poswsum[PS_D] ) {
                sprintf( errbuf, "Invalid MID state weights weights at pos. %d.", p );
                merror = errbuf;
                break;
            }

//             emss = SLC_MIN(( double )mss * avgpeseq, ( double )noeffress );
//             eiss = SLC_MIN(( double )iss * avgpeseq, ( double )noeffress );
//             edss = SLC_MIN(( double )dss * avgpeseq, ( double )noeffress );

//             emss = SLC_MIN( avgndr[pp], ( double )noeffress );
//             eiss = SLC_MIN( avgndr[pp] * eiss / poswsum[PS_M], ( double )noeffress );
//             edss = SLC_MIN( avgndr[pp] * poswsum[PS_D] / poswsum[PS_M], ( double )noeffress );

            //convex function to avoid of rapid drop of expectation as weight decreases
            emss = SLC_MIN(avgndr[pp]* exp( gdWCPC*(1.0-poswsum[PS_M])*log(poswsum[PS_M])),( double )noeffress);
            if( eiss )
                eiss = SLC_MIN(avgndr[pp]* exp(gdWCPC*(1.0-eiss)*log(eiss)),( double )noeffress);
            if( poswsum[PS_D])
                edss = SLC_MIN(avgndr[pp]* exp(gdWCPC*(1.0-poswsum[PS_D])*log(poswsum[PS_D])),(double)noeffress);
            else
                edss = 0.0;

            //M state no. diff residues is required
            queryDescription->SetNoSymbolsInExtentAt( p, PS_M, SLC_MAX( 1.0, emss ));

            emss = ExtendedDescriptionVector::GetExpNoObservations( emss );
            eiss = ExtendedDescriptionVector::GetExpNoObservations( eiss );
            edss = ExtendedDescriptionVector::GetExpNoObservations( edss );

            if(( double )( mss + dss ) < emss ) emss = ( double )( mss + dss );
            if(( double )( mss + dss ) < eiss ) eiss = ( double )( mss + dss );
            if(( double )( mss + dss ) < edss ) edss = ( double )( mss + dss );

            if( emss < 1.0 ) emss = 1.0;
            if( eiss < 1.0 && eiss ) eiss = 1.0;
            if( edss < 1.0 && edss ) edss = 1.0;

            //set expectations directly, avoiding using of histogram data
            queryDescription->SetMIDExpNoObservationsAt( emss, p, PS_M );
            queryDescription->SetMIDExpNoObservationsAt( eiss, pM, PS_I );
            queryDescription->SetMIDExpNoObservationsAt( edss, p, PS_D );
//             queryDescription->IncDistinctHistAt( p, PS_M, mss );
//             queryDescription->IncDistinctHistAt( pM, PS_I, iss );
//             queryDescription->IncDistinctHistAt( p, PS_D, dss );
            iss = ilen = 0;
            eiss = 0.0;
            pM = p;
            pp++;
        }
        else
            //increase insertion length
            ilen++;
    }//3:for(p...)


    if( iss && ilen ) {
        //average observations
        eiss /= ( double )iss;
        iss = ( size_t )rint(( double )iss /( double )ilen );
        if( /*poswsum[PS_M] &&*/ eiss && nomstates ) {
//             eiss = SLC_MIN( avgndr[nomstates-1] * eiss / poswsum[PS_M], ( double )noeffress );
            eiss = SLC_MIN(avgndr[nomstates-1]* exp( gdWCPC*(1.0-eiss)*log(eiss)),( double )noeffress);
            eiss = ExtendedDescriptionVector::GetExpNoObservations( eiss );
            if(( double )( mss + dss ) < eiss ) eiss = ( double )( mss + dss );
            queryDescription->SetMIDExpNoObservationsAt( eiss, pM, PS_I );
        }
//         queryDescription->IncDistinctHistAt( pM, PS_I, ( int )rint(( double )iss * avgpeseq ));
    }


    //check for error
    if( !merror.empty())
        throw myruntime_error( preamble + merror );
}





// =========================================================================
// CalcTW: calculate combined transition weight given two 
//  weights of adjacent positions
//
double CalcTW( double wp, double w )
{
    const bool  clbPROD = false;
    double  result = wp;
    if( clbPROD && wp && w )
        result = wp * w;
    return result;
}

// -------------------------------------------------------------------------
// ComputeTransitionFrequencies: compute observed weighted transition 
//     frequencies at each position
//
void InputMultipleAlignment::ComputeTransitionFrequencies( bool usegwgts, bool expmids )
{
    const int       noseqs = size();
    const int       maxtr = gTPTRANS_NTPS;  // number of transitions per state
    int             states, st;         // state values
    double          gwghts[noseqs];     // sequence weights computed globally
    double*         Mweights = NULL;    // M state sequence weights at a position
    double*         Iweights = NULL;    // I state sequence weights
    double*         Dweights = NULL;    // D state sequence weights
    double*         ppmMwghs = NULL;    // M state sequence weights at previous position
    double*         ppmIwghs = NULL;    // I state sequence weights at previous position
    double*         ppmDwghs = NULL;    // D state sequence weights at previous position
    double          wM, wI, wD;         // weights of single sequence at a position
    double          wppmM, wppmI, wppmD;
    ssize_t*        inserts = NULL;     // numbers of insertions between adjacent match states
    unsigned char   residue, ppmres;
    mystring        merror;
    size_t          p, i, nosx;
    ssize_t         ppp, ppm, ppm_t;    //indices of previous positions
    ssize_t         pppstate;

    if( !queryDescription )
        throw myruntime_error("ComputeTransitionFrequencies: Unable to compute frequencies.");

    inserts = ( ssize_t* )malloc( sizeof( ssize_t ) * ( noseqs + 1 ));

     if( inserts == NULL )
        throw myruntime_error("ComputeTransitionFrequencies: Not enough memory.");

    memset( inserts, 0, sizeof( ssize_t ) * ( noseqs + 1 ));

    ppp = ppm = -1;
    pppstate = 1;//match


    if( usegwgts ) {
        for( i = 0; i < noseqs; i++ ) {
            if( SequenceAt(i))
                gwghts[i] = SequenceAt(i)->GetGlbWeight();
            else
                gwghts[i] = 0.0;
        }
        Mweights = Iweights = Dweights = gwghts;
        ppmMwghs = ppmIwghs = ppmDwghs = gwghts;
    }


    // iterate over all query positions
    for( p = 0; p < queryDescription->size(); p++ ) {
        // omit unsued positions
        if( !queryDescription->IsUsedAt( p ))
            continue;

        if( !usegwgts ) {
            Mweights = Iweights = Dweights = queryDescription->GetSqnWeightsAt( p, PS_M );
//             Iweights = queryDescription->GetSqnWeightsAt( p, PS_I );
//             Dweights = queryDescription->GetSqnWeightsAt( p, PS_D );
        }

        // weights may be not computed if there were no extents for the position
//         if( weights == NULL )
//             continue;

        ppm_t = p;
        if( 0 <= ppm )
            ppm_t = ppm;

        if( !usegwgts ) {
            ppmMwghs = ppmIwghs = ppmDwghs = queryDescription->GetSqnWeightsAt( ppm_t, PS_M );
//             ppmIwghs = queryDescription->GetSqnWeightsAt( ppm_t, PS_I );
//             ppmDwghs = queryDescription->GetSqnWeightsAt( ppm_t, PS_D );
        }

//         if( ppmwghs == NULL )
//             continue;
//         nosx = queryDescription->GetNoSequencesInMSExtentAt( ppm_t );
//         if( !nosx )
//             continue;

        for( i = 0; i < noseqs; i++ ) {
            if( !SequenceAt(i)->GetUsed())
                continue;
//             if( !SequenceAt(i)->IsUsedInExtentAt( p ))
            if( !SequenceAt(i)->IsUsedAndInfAt( p ))
                continue;
            residue = SequenceAt(i)->ResidueAt( p );

            wM = wI = wD = 0.0;
            if( Mweights ) wM = Mweights[i];
            if( Iweights ) wI = Iweights[i];
            if( Dweights ) wD = Dweights[i];

            //{{PROCESS TRANSITION WEIGHTS
            ppmres = X;//match

            wppmM = wppmI = wppmD = 0.0;
            if( ppmMwghs ) wppmM = ppmMwghs[i];
            if( ppmIwghs ) wppmI = ppmIwghs[i];
            if( ppmDwghs ) wppmD = ppmDwghs[i];

            if( 0 <= ppm ) {
//                 if( SequenceAt(i)->IsUsedInExtentAt( ppm ))
                if( SequenceAt(i)->IsUsedAndInfAt( ppm ))
                    ppmres = SequenceAt(i)->ResidueAt( ppm );
                else
                    continue;//require both positions to be in use
            }

            if( queryDescription->GetStateAt( p )) {
                if( pppstate ) {
                    if( IsValidResSym( residue )) {
                        if( IsValidResSym( ppmres )) {
                            queryDescription->IncTransWeightsAt( CalcTW( wppmM, wM ), P_MM, ppm );
                            if( expmids ) queryDescription->IncMIDExpNoObservationsAt( wppmM, ppm, PS_M );
                        } else {
                            queryDescription->IncTransWeightsAt( CalcTW( wppmD, wM ), P_DM, ppm );
                            if( expmids ) queryDescription->IncMIDExpNoObservationsAt( wppmD, ppm, PS_D );
                        }
                    }
                    else {
                        if( IsValidResSym( ppmres )) {
                            queryDescription->IncTransWeightsAt( CalcTW( wppmM, wD ), P_MD, ppm );
                            if( expmids ) queryDescription->IncMIDExpNoObservationsAt( wppmM, ppm, PS_M );
                        } else {
                            queryDescription->IncTransWeightsAt( CalcTW( wppmD, wD ), P_DD, ppm );
                            if( expmids ) queryDescription->IncMIDExpNoObservationsAt( wppmD, ppm, PS_D );
                        }
                    }
                }
                else {
                    if( IsValidResSym( residue )) {
                        if( 0 < inserts[i]) {
                            queryDescription->IncTransWeightsAt( CalcTW( wppmI, wM ), P_IM, ppm );
                            if( expmids ) queryDescription->IncMIDExpNoObservationsAt( wppmI, ppm, PS_I );
                        } else {
                            if( IsValidResSym( ppmres )) {
                                queryDescription->IncTransWeightsAt( CalcTW( wppmM, wM ), P_MM, ppm );
                                if( expmids ) queryDescription->IncMIDExpNoObservationsAt( wppmM, ppm, PS_M );
                            } else {
                                queryDescription->IncTransWeightsAt( CalcTW( wppmD, wM ), P_DM, ppm );
                                if( expmids ) queryDescription->IncMIDExpNoObservationsAt( wppmD, ppm, PS_D );
                            }
                        }
                    }
                    else {
                        if( 0 < inserts[i]) {
                            queryDescription->IncTransWeightsAt( CalcTW( wppmI, wD ), P_ID, ppm );
                            if( expmids ) queryDescription->IncMIDExpNoObservationsAt( wppmI, ppm, PS_I );
                        } else {
                            if( IsValidResSym( ppmres )) {
                                queryDescription->IncTransWeightsAt( CalcTW( wppmM, wD ), P_MD, ppm );
                                if( expmids ) queryDescription->IncMIDExpNoObservationsAt( wppmM, ppm, PS_M );
                            } else {
                                queryDescription->IncTransWeightsAt( CalcTW( wppmD, wD ), P_DD, ppm );
                                if( expmids ) queryDescription->IncMIDExpNoObservationsAt( wppmD, ppm, PS_D );
                            }
                        }
                    }
                }
            }
            else {
                if( IsValidResSym( residue ))
                    inserts[i]++;

                if( IsValidResSym( residue )) {
                    if( 1 < inserts[i]) {
                        queryDescription->IncTransWeightsAt( CalcTW( wppmI, wI ), P_II, ppm );
                        if( expmids ) queryDescription->IncMIDExpNoObservationsAt( wppmI, ppm, PS_I );
                    } else {//first time when insertion residue is met
                        if( IsValidResSym( ppmres )) {
                            queryDescription->IncTransWeightsAt( CalcTW( wppmM, wI ), P_MI, ppm );
                            if( expmids ) queryDescription->IncMIDExpNoObservationsAt( wppmM, ppm, PS_M );
                        } else {
                            queryDescription->IncTransWeightsAt( CalcTW( wppmD, wI ), P_DI, ppm );
                            if( expmids ) queryDescription->IncMIDExpNoObservationsAt( wppmD, ppm, PS_D );
                        }
                    }
                }
            }
            //}}PROCESS TW
        }

        ppp = p;//last position
        pppstate = queryDescription->GetStateAt( p );
        if( pppstate ) {
            //reset inserts: this is match state
            memset( inserts, 0, sizeof( ssize_t ) * ( noseqs + 1 ));
            ppm = p;//last match state position
        }
    }

    //{{PROCESS END-STATE
    //to assign transition probabilities to the end state...
    //back-iterate over query positions down to match state
    for( ppp = queryDescription->size() - 1; 0 <= ppp ; ppp-- ) {
        // omit unsued positions
        if( !queryDescription->IsUsedAt( ppp ))
            continue;

        bool state = queryDescription->GetStateAt( ppp );
        if( !state ) {
            if( 0 <= ppm ) {
                if( ppm <= ppp ) {
                    if( !queryDescription->IsUsedAt( ppm ))
                        warning( "End-state transitions: Unused last match state position." );
                    else
                        ppp = ppm;
                } else {
                    warning( "End-state transitions: Invalid match state position." );
                    break;
                }
            }
        }

        if( !usegwgts ) {
            Mweights = Iweights = Dweights = queryDescription->GetSqnWeightsAt( ppp, PS_M );
//             Iweights = queryDescription->GetSqnWeightsAt( ppp, PS_I );
//             Dweights = queryDescription->GetSqnWeightsAt( ppp, PS_D );
        }

        if( Mweights == NULL && Iweights == NULL && Dweights == NULL ) {
            if( state ) {
                warning( "End-state transitions: Null state sequence weights." );
                break;
            }
            continue;
        }

        for( i = 0; i < noseqs; i++ ) {
            if( !SequenceAt(i)->GetUsed())
                continue;
//             if( !SequenceAt( i )->IsUsedInExtentAt( ppp ))
            if( !SequenceAt( i )->IsUsedAndInfAt( ppp ))
                continue;
            residue = SequenceAt( i )->ResidueAt( ppp );

            wM = wI = wD = 0.0;
            if( Mweights ) wM = Mweights[i];
            if( Iweights ) wI = Iweights[i];
            if( Dweights ) wD = Dweights[i];

            //{{PROCESS TRANSITION WEIGHTS
            if( queryDescription->GetStateAt( ppp )) {
                if( IsValidResSym( residue )) {
                    queryDescription->IncTransWeightsAt( CalcTW( wM, wM ), P_MM, ppp );
                    if( expmids ) queryDescription->IncMIDExpNoObservationsAt( wM, ppp, PS_M );
                } else {
                    queryDescription->IncTransWeightsAt( CalcTW( wD, wM ), P_DM, ppp );
                    if( expmids ) queryDescription->IncMIDExpNoObservationsAt( wD, ppp, PS_D );
                }
            }
            else {
                if( 0 < inserts[i]) {
                    queryDescription->IncTransWeightsAt( CalcTW( wI, wM ), P_IM, ppp );
                    if( expmids ) queryDescription->IncMIDExpNoObservationsAt( wI, ppp, PS_I );
                } else {
                    if( IsValidResSym( residue ) || ppm < 0 ) {
                        queryDescription->IncTransWeightsAt( CalcTW( wM, wM ), P_MM, ppm < 0? ppm: ppp );
                        if( expmids ) queryDescription->IncMIDExpNoObservationsAt( wM, ppm < 0? ppm: ppp, PS_M );
                    } else {
                        queryDescription->IncTransWeightsAt( CalcTW( wD, wM ), P_DM, ppp );
                        if( expmids ) queryDescription->IncMIDExpNoObservationsAt( wD, ppp, PS_D );
                    }
                }
            }
            //}}PROCESS TW
        }
        break;//leave loop
    }
    //}}PROCESS END

    //NORMALIZE Transition weights and observed state frequencies
    for( ppp = -1; ppp < ( ssize_t )queryDescription->size(); ppp++ ) {
        // omit unsued positions
        if( 0 <= ppp && !queryDescription->IsUsedAt( ppp ))
            continue;
        if( 0 <= ppp && !queryDescription->GetStateAt( ppp ))
            continue;

        //process MID states
        for( states = 0; states < P_NSTATES; states += maxtr ) {
            wM = 0.0;//for norm.
            for( st = states; st < states + maxtr && st < P_NSTATES; st++ ) {
//!                 if( st == P_ID || st == P_DI ) continue;//IGNORE states P_ID and P_DI
                wM += queryDescription->GetTransWeightsAt( st, ppp );
            }
            if( wM )
                for( st = states; st < states + maxtr && st < P_NSTATES; st++ ) {
                    //if( st == P_ID || st == P_DI ) continue;//IGNORE states P_ID and P_DI
                    wppmM = queryDescription->GetTransWeightsAt( st, ppp ) / wM;
                    queryDescription->SetTransWeightsAt( wppmM, st, ppp );
                }
        }
        if( expmids ) {
            //normalize observed state frequencies
            wM = 0.0;//for norm.
            for( st = 0; st < PS_NSTATES; st++ )
                wM += queryDescription->GetMIDExpNoObservationsAt( ppp, st );
            if( wM )
                for( st = 0; st < PS_NSTATES; st++ ) {
                    wppmM = queryDescription->GetMIDExpNoObservationsAt( ppp, st ) / wM;
                    queryDescription->SetMIDExpNoObservationsAt( wppmM, ppp, st );
                }
        }
    }

// // queryDescription->PrintTransWeights( stderr );
// // queryDescription->PrintMIDExpNoObservations( stderr );
    free( inserts );
    return;
}





// =========================================================================
// ComputeTargetTransFrequencies: compute target transition  frequencies at 
//     each position
//
void InputMultipleAlignment::ComputeTargetTransFrequencies()
{
    const size_t    nomstates = queryDescription->GetEffectiveSize();
    const double ( *ptargs )[P_NSTATES];
    double          trns[P_NSTATES];
    double          expn;
    double          nm, ni, nd;
    ssize_t         ppp;
    int             st;

    if( !queryDescription )
        throw myruntime_error( "ComputeTargetTransFrequencies: Memory access error." );

    if(( ssize_t )queryDescription->size() <= 0 )
        return;

    trns[P_MM] = 1.0; trns[P_MI] = 0.0; trns[P_MD] = 0.0;
    trns[P_IM] = 1.0; trns[P_II] = 0.0; trns[P_ID] = 0.0;
    trns[P_DM] = 1.0; trns[P_DI] = 0.0; trns[P_DD] = 0.0;

    queryDescription->SetTargetTranstAt( &trns, ppp = -1 );

    // iterate over all query positions; beginning state, -1
    for( ppp = 0; ppp < ( ssize_t )queryDescription->size(); ppp++ ) {
        if( 0 <= ppp && !queryDescription->IsUsedAt( ppp ))
            continue;
        if( 0 <= ppp && !queryDescription->GetStateAt( ppp ))
            continue;

//         if( ppp < 0 )
//               expn = queryDescription->GetExpNoObservationsAt( 0 );
//         else  expn = queryDescription->GetExpNoObservationsAt( ppp );

        nm = queryDescription->GetMIDExpNoObservationsAt( ppp, PS_M );
        ni = queryDescription->GetMIDExpNoObservationsAt( ppp, PS_I );
        nd = queryDescription->GetMIDExpNoObservationsAt( ppp, PS_D );

        //one exp. sequence bears no inform.
        nm = SLC_MAX( 0.0, nm - 1.0 );

        TRANSPROBS.SetEffNoSequences( nm, ni, nd );
        TRANSPROBS.PME( queryDescription->GetTransWeightsAt( ppp ));//conserv. there
        ptargs = TRANSPROBS.GetPMEstimators();
        if( ptargs == NULL )
            throw myruntime_error("ComputeTargetTransFrequencies: Failed to compute target trans. frequencies.");
        queryDescription->SetTargetTranstAt( ptargs, ppp );
    }

    if( nomstates )
        queryDescription->SetTargetTranstAt( &trns, ppp = ( ssize_t )nomstates - 1 );
// // queryDescription->PrintTargetTranst( stderr );
    return;
}





// =========================================================================
// ExpectedNoObservationsAt: calculate expected number of independent
//  observations; 
//
double InputMultipleAlignment::ExpectedNoObservationsAt( size_t pos, int st ) const
{
    if( queryDescription == NULL || queryDescription->size() <= pos || 
        st < 0 || PS_NSTATES <= st )
        throw myruntime_error("ExpectedNoObservationsAt: Memory access error.");

    if( !queryDescription->IsUsedAt( pos ))
        return 0.0;

    ssize_t interval = ( ssize_t )queryDescription->GetMSExtentIntervalAt( pos );
    ssize_t noeffres = NUMAA;
    ssize_t halfint = 0;
    ssize_t estnocols;  //number of columns having certain number of distinct residues
    ssize_t sumnocols;  //accumulated number of columns
    ssize_t distinct;   //accumulated number of distinct residues in half of the extent
    double  avgdistinct;//average number of distinct residues
    double  expnobs;    //expected number of observations
    ssize_t d;

    if( interval < 0 || queryDescription->GetEffectiveSize() < interval ) {
        warning( "ExpectedNoObservationsAt: Extent interval out of range." );
        return 0.0;
    }
    if( st == PS_I && interval < 1 )
        return 0.0;
    if( interval < 1 ) {
        warning( "ExpectedNoObservationsAt: Match state extent interval is 0." );
        return 0.0;
    }

    halfint = ( interval + 1 ) >> 1;

    sumnocols = 0;
    distinct = 0;

    //calculate average number of distinct residues in half of the extent
    //starting with the most distinct amino acids
    for( d = noeffres; d && sumnocols < halfint; d-- ) {
        estnocols = queryDescription->GetDistinctHistAt( pos, st, d );
        sumnocols += estnocols;
        distinct += estnocols * d;

        if( halfint < sumnocols ) {
            distinct -= ( sumnocols - halfint ) * d;
            sumnocols = halfint;
        }
    }

    if( sumnocols < 1 )
//         throw myruntime_error( "ExpectedNoObservationsAt: No observed counts." );
        return 0.0;

    avgdistinct = ( double )distinct /( double )sumnocols;

    //to calculate expected number of independent observations
    //corresponding to average number of distinct residues,
    //use the Altschul et al. method as in NAR, 2009

    expnobs = ExtendedDescriptionVector::GetExpNoObservations( avgdistinct );

    if( queryDescription->GetNoSequencesInMSExtentAt( pos ) < expnobs )
        expnobs = queryDescription->GetNoSequencesInMSExtentAt( pos );

    //gaps are not included in calc. of distinct res. hist.
    if( expnobs < 1.0 )
        expnobs = 1.0;

    return expnobs;
}

// -------------------------------------------------------------------------
// DeriveExpectedNoObservations: calculate expected numbers of 
//     observations at each position
//
void InputMultipleAlignment::DeriveExpectedNoObservations()
{
    double  wm, wi, wd; //MID state weights 
    double  expn[PS_NSTATES];
    size_t  p;

    if( !queryDescription )
        throw myruntime_error( "DeriveExpectedNoObservations: Memory access error." );

    //{{NOTE;below:to be removed.
    //set MID expected observations at the beginning position
//     if( queryDescription->size()) {
//         pp = -1;
//         wm = queryDescription->GetMIDExpNoObservationsAt( pp, PS_M );
//         wi = queryDescription->GetMIDExpNoObservationsAt( pp, PS_I );
//         wd = queryDescription->GetMIDExpNoObservationsAt( pp, PS_D );
//         if( wm <= 0.0 ) {
//             warning("DeriveExpectedNoObservations: Null observed match state frequencies.");
//             wm = 1.0;
//         }
//         for( r = 0; r < queryDescription->size(); r++ )
//             if(!queryDescription->IsUsedAt( r ))
//                 continue;
//         if( r < queryDescription->size()) {
//             expn = ExpectedNoObservationsAt( r );
//             //overwrite MID obs. frequencies with MID expected observations
//             queryDescription->SetMIDExpNoObservationsAt( expn, pp, PS_M );
//             queryDescription->SetMIDExpNoObservationsAt( expn * wi / wm, pp, PS_I );
//             queryDescription->SetMIDExpNoObservationsAt( expn * wd / wm, pp, PS_D );
//         }
//     }
    //}}

    //iterate over all query positions
    for( p = 0; p < queryDescription->size(); p++ ) {
        if(!queryDescription->IsUsedAt( p ))
            continue;
        if(!queryDescription->GetStateAt( p ))
            continue;

        //MID observed frequencies should be precalculated
        wm = queryDescription->GetMIDExpNoObservationsAt( p, PS_M );
        wi = queryDescription->GetMIDExpNoObservationsAt( p, PS_I );
        wd = queryDescription->GetMIDExpNoObservationsAt( p, PS_D );
        if( wm <= 0.0 ) {
            if( queryDescription->GetStateAt( p ))
                warning("DeriveExpectedNoObservations: Null observed match state frequencies.");
            wm = 1.0;
        }
        expn[PS_M] = expn[PS_I] = expn[PS_D] = 0.0;
        expn[PS_M] = ExpectedNoObservationsAt( p, PS_M );
        expn[PS_I] = ExpectedNoObservationsAt( p, PS_I );
        expn[PS_D] = ExpectedNoObservationsAt( p, PS_D );
        queryDescription->SetExpNoObservationsAt( expn[PS_M], p );
        //overwrite MID obs. frequencies with MID expected observations
        queryDescription->SetMIDExpNoObservationsAt( expn[PS_M], p, PS_M );
        queryDescription->SetMIDExpNoObservationsAt( expn[PS_I], p, PS_I );
        queryDescription->SetMIDExpNoObservationsAt( expn[PS_D], p, PS_D );
//         queryDescription->SetMIDExpNoObservationsAt( expn, p, PS_M );
//         queryDescription->SetMIDExpNoObservationsAt( expn * wi / wm, p, PS_I );
//         queryDescription->SetMIDExpNoObservationsAt( expn * wd / wm, p, PS_D );
    }
}

// -------------------------------------------------------------------------
// MedianValue: calculate median value of given distribution
//  scale -- histogram values scaled by
//
void MedianValue( size_t* histdata, size_t szdat, double* median, size_t scale )
{
    if( histdata == NULL )
        return;
    size_t  sumoccs = 0; //sum of occurences
    size_t  novals = 0;//number of values
    size_t  hnovals = 0;//half number
    size_t  mass = 0;
    double  medn = 0.0;
    size_t  e;

    for( e = 0; e <= szdat; e++ ) {
        novals += histdata[e];
        mass += histdata[e] * e;
    }
    hnovals = novals >> 1;
    if( hnovals && ( novals & 1 ))
        hnovals++;//1-based further sum

    for( e = 0; e <= szdat; e++ ) {
        sumoccs += histdata[e];
        if( hnovals <= sumoccs ) {
            medn = ( double )e;//median value
            //if length is even and it is boundary of bins or bin contains <2 elms.
            if( e < szdat  && !( novals & 1 ) && ( hnovals == sumoccs || histdata[e] < 2 )) { 
                for( e++; e <= szdat && !histdata[e]; e++ );
                if( e <= szdat )
                    //take average value of adjacent expected values
                    medn = ( medn + ( double )e ) * 0.5;
            }
            break;
        }
    }
    if( medn <= 0.0 && novals )
        medn = ( double )mass /( double )novals;//average then
    medn /= ( double )scale;
    if( median )
        *median = medn;
}

// -------------------------------------------------------------------------
// CalculateEffNoSequences: assign effective number of sequences to
//  expected number of observations calculated over all positions
//
void InputMultipleAlignment::CalculateEffNoSequences()
{
    double  expn[PS_NSTATES];
    double  effnos = 0.0;
    bool    pstate;
    size_t  p;

    if( !queryDescription )
        throw myruntime_error( "CalculateEffNoSequences: Memory access error." );

    const int   maxexpval = ExtendedDescriptionVector::GetSizeOfExpNoDistinctRes();
    size_t      expnhist[maxexpval+1];

    memset( expnhist, 0, sizeof( size_t )*( maxexpval + 1 ));

    // iterate over all query positions
    for( p = 0; p < queryDescription->size(); p++ ) {
        if( !queryDescription->IsUsedAt( p ))
            continue;
        pstate = queryDescription->GetStateAt( p );
        if( !pstate )
            continue;

        expn[PS_M] = expn[PS_I] = expn[PS_D] = 0.0;
        expn[PS_M] = queryDescription->GetMIDExpNoObservationsAt( p, PS_M );
        expn[PS_I] = queryDescription->GetMIDExpNoObservationsAt( p, PS_I );
        expn[PS_D] = queryDescription->GetMIDExpNoObservationsAt( p, PS_D );
        if( expn[PS_M] <= 0.0 ) {
            if( pstate )
                warning("CalculateEffNoSequences: Null observed match state frequencies.");
            expn[PS_M] = 1.0;
        }
        expnhist[( int )rint( SLC_MIN( maxexpval, expn[PS_M]))]++;
    }

    //calculate median value for expected no. observations
    MedianValue( expnhist, maxexpval, &effnos );
    SetEffNoSequences( effnos );
}

    


// =========================================================================
// AdjustWeights: adjusts weights
//
void InputMultipleAlignment::AdjustWeights()
{
    size_t  p;

    if( !queryDescription )
        throw myruntime_error("AdjustWeights: Null description vector." );

    //iterate over all query positions
    for( p = 0; p < queryDescription->size(); p++ ) {
        // omit positions consisting of one symbol or which are unsued in query
        if( !queryDescription->IsUsedAt( p ))
            continue;
        if( !queryDescription->GetStateAt( p ))
            continue;

        //{{NOTE:below:to be removed.
        //if position p was not used to compute extents (NotUsed flag)
//         if( !queryDescription->GetExtentIntervalAt( p ))
//             continue;
        //}}

        AdjustWeightsAt( p );
    }
}

// -------------------------------------------------------------------------
// AdjustWeights: adjusts weights
//
void InputMultipleAlignment::AdjustWeightsAt( size_t p, double (*mtwghts)[NUMALPH] )
{
    const double    accuracy = 1.0e-4;
    const int       efective_number = NUMAA;// effective number of residues
    double  posum = 0.0,                // sum of weights in a column
            gapw = 0.0,                 // gap weight for a position
            Xw = 0.0,                   // weight for symbol X
            Bw = 0.0,                   // weight for symbol B
            Zw = 0.0;                   // weight for symbol Z
    double  mw;
    char    strbuf[KBYTE];
    unsigned char r, e;

    static int  symB = HashAlphSymbol('B');
    static int  symZ = HashAlphSymbol('Z');
    static int  resN = HashAlphSymbol('N');
    static int  resD = HashAlphSymbol('D');
    static int  resQ = HashAlphSymbol('Q');
    static int  resE = HashAlphSymbol('E');

    if( !queryDescription )
        return;

    double  ( *pmw )[NUMALPH] = mtwghts;
    if( pmw == NULL )
        pmw = const_cast<double(*)[NUMALPH]>( queryDescription->GetMatchWeightsAt( p ));

    if( pmw == NULL )
        throw myruntime_error("AdjustWeightsAt: Null match weight vector.");

//     //{{NOTE:below:to be removed.
//     if( queryDescription->ResidueAt( p ) == X ) {
//         for( e = 0; e < NUMALPH; e++ ) {
//             if( 0.0 < LOSCORES.PROBABility( e ))
//                 queryDescription->SetMatchWeightsAt( LOSCORES.PROBABility( e ), e, p );
//             else
//                 queryDescription->SetMatchWeightsAt( 0.0, e, p );
//         }
//         //NOTE:below:to be removed.
// //         queryDescription->SetNoSymbolsInExtentAt( 0, p );
//         continue;
//     }
//     //}}

    posum = 0.0;
//     for( r = 0; r < NUMALPH; r++ )
    for( r = 0; IsValidResSym( r ); r++ )
        posum += (*pmw)[r];

    if( !posum )
        for( r = 0; r < NUMAA; r++ ) {
            mw = LOSCORES.PROBABility( r );
            (*pmw)[r] = mw;
            posum += mw;
        }

    if( posum < 1.0 - accuracy || posum > 1.0 + accuracy ) {
        sprintf( strbuf, "Weights are not properly normalized, %-12g (pos. %d)", posum, p );
        throw myruntime_error( strbuf );
    }

    gapw = (*pmw)[GAP];
    Xw = (*pmw)[X];
    Bw = (*pmw)[symB];
    Zw = (*pmw)[symZ];

    // spread out X weight!
    for( r = 0; r < NUMALPH; r++ ) {
        if( 0.0 < LOSCORES.PROBABility( r ) && r != X )
            (*pmw)[r] += LOSCORES.PROBABility( r ) * Xw;
    }

    (*pmw)[X] = 0.0;
    (*pmw)[GAP] = 0.0;
    if( mtwghts == NULL )
        queryDescription->SetGapWeightsAt( gapw, p );


    static double   Bprob = LOSCORES.PROBABility( resN ) + LOSCORES.PROBABility( resD );
    static double   Zprob = LOSCORES.PROBABility( resQ ) + LOSCORES.PROBABility( resE );

    // spread out B weight: N or D
    if( 0.0 < Bw  ) {
        (*pmw)[resN] += LOSCORES.PROBABility( resN ) * Bw / Bprob;
        (*pmw)[resD] += LOSCORES.PROBABility( resD ) * Bw / Bprob;
        (*pmw)[symB] = 0.0;
    }
    // spread out Z weight: Q or E
    if( 0.0 < Zw  ) {
        (*pmw)[resQ] += LOSCORES.PROBABility( resQ ) * Zw / Zprob;
        (*pmw)[resE] += LOSCORES.PROBABility( resE ) * Zw / Zprob;
        (*pmw)[symZ] = 0.0;
    }


//     //{{NOTE:below:to be removed.
//     posum = 0.0;
//     for( r = 0; r < NUMAA; r++ ) {
//         mweight = queryDescription->GetMatchWeightsAt( r, p );
//         if( 0.0 < mweight ) {
//             posum += mweight;
//             //round all match weights to a floating point number with precision 2;
//             //such a value discretization ensures that for the same distribution of
//             //match weights scores will also be the same
//             queryDescription->SetMatchWeightsAt( rint( FREQUENCY_SUM * mweight ) / FREQUENCY_SUM, r, p );
//         }
//     }
// 
//     if( posum < 1.0 - accuracy || posum > 1.0 + accuracy )
//         throw myruntime_error( mystring( "Weight adjustment failed." ));
//     //}}
}

// -------------------------------------------------------------------------
// ComputeEstimatedFrequencies: compute estimated probabilities for each
//  residue at each query position;
//  estimated probabilities are computed using pseudo frequencies g and 
//  weighted frequecies f;
//
void InputMultipleAlignment::ComputeTargetFrequencies()
{
    if( !queryDescription )
        throw myruntime_error("Unable to estimate target probabilities.");

    double  information = 0.0;  // information content per position
    double  pseudoFreqn = 0.0;  // pseudocount frequency
    double  obsFrequencyWeight = 0.0;   // weight for observed frequencies (weighted), alpha
    double  denominator = 0.0;  // sum of weights of pseudocount and observed frequencies
    double  estimated = 0.0;    // estimated probability
    unsigned char   a, b;
    size_t  p;

    // iterate over all query positions
    for( p = 0; p < queryDescription->size(); p++ ) {
        if( !queryDescription->GetStateAt( p ))
            continue;

        //{{NOTE:below:to be removed.
        // omit positions consisting of Xs and gaps in query
//         if(  queryDescription->ResidueAt( p ) == X ||
//              queryDescription->ResidueAt( p ) == ASTERISK )
//                 // in this case all target frequencies (est. probs) should be equal to
//                 continue;
        //}}

        //{{NOTE:below:to be removed.
        // theoretically interval size can be zero if this position
        //  was not used to compute extents (NotUsed flag)
//         if( !queryDescription->GetExtentIntervalAt( p ))
//             continue;
        //}}

        // weight computed as expected value of different residues per extent column
        obsFrequencyWeight = queryDescription->ComputeObsFrequencyWeightAt( p );
        information = 0.0;
        denominator = obsFrequencyWeight + GetPseudoCountWeight();  // alpha + beta

        for( a = 0; a < NUMAA; a++ ) {
            //unable to estimate probabilities if background probabilities are zero
            if( LOSCORES.PROBABility( a ) <= 0.0 )
                continue;

            pseudoFreqn = 0.0;

            for( b = 0; b < NUMAA; b++ )
                pseudoFreqn += queryDescription->GetMatchWeightsAt( b, p ) * LOSCORES.FreqRatio( a, b );

            pseudoFreqn *= LOSCORES.PROBABility( a );
            estimated = ( obsFrequencyWeight * queryDescription->GetMatchWeightsAt( a, p ) + 
                          GetPseudoCountWeight() * pseudoFreqn ) / denominator;
            queryDescription->SetTargetFreqnsAt( estimated, a, p );

            if( 0.0 < estimated )
                // Q[a] log2( Q[a]/p[a] )
                information += estimated * log( estimated / LOSCORES.PROBABility( a )) / LN2;
        }
        //save information content
        queryDescription->SetInformationAt(( 0.0 < information )? information: 0.0, p );
    }
}

// =========================================================================
// ComputeTargetFrequenciesMDL: pseudocounts are computed by the MDL theory
// (see Altschul et al. theory in NAR, 2009)
//
void InputMultipleAlignment::ComputeTargetFrequenciesMDL()
{
    if( !queryDescription )
        throw myruntime_error("Unable to estimate target probabilities.");

    const double    backcount = 2.0;            //pseudocount for initialy mixing of frequencies
    double          expobscount;                //expected number of observations at certain position
    double          pseudomweights[NUMALPH];    //obs. frequencies mixed with background probabilities
    double          pseudomsum;                 //sum of mixed counts
    double          pseudominfo;                //relative entropy of frequencies to background probabilities
    double          pseudofreqweight;           //weight for pseudo frequencies (alpha)
    const double    psnumerator = 0.0806;       //numerator a in expression a/relent^b
    const double    psexponent = 0.53;          //exponent b of denominator in a/relent^b
    double          psdenominator;              //denominator in expression a/relent^b
    const double    lowvalue = 0.0001;

    double  information = 0.0;          //information content per position
    double  pseudoFreqn = 0.0;          //pseudocount frequency
    double  estimated = 0.0;            //estimated probability
    unsigned char   a, b;
    size_t  p;

    //iterate over all query positions
    for( p = 0; p < queryDescription->size(); p++ ) {
        if( !queryDescription->GetStateAt( p ))
            continue;

    //{{
        memset( pseudomweights, 0, NUMALPH * sizeof( double ));
//         expobscount = queryDescription->GetMIDExpNoObservationsAt( p, PS_M );
        expobscount = queryDescription->GetNoSymbolsInExtentAt( p, PS_M );
        pseudomsum = 0.0;
        pseudominfo = 0.0;
        pseudofreqweight = 1.0;

        for( a = 0; a < NUMAA; a++ ) {
            pseudomweights[a] = ( queryDescription->GetMatchWeightsAt( a, p ) * expobscount ) +
                                ( LOSCORES.PROBABility( a ) * backcount );
            pseudomsum += pseudomweights[a];
        }
        if( pseudomsum )
            for( a = 0; a < NUMAA; a++ ) {
                pseudomweights[a] /= pseudomsum;
                if( 0.0 < pseudomweights[a] && 0.0 < LOSCORES.PROBABility( a ))
                    pseudominfo += pseudomweights[a] * log( pseudomweights[a] / LOSCORES.PROBABility( a ));
            }

        pseudominfo /= LN2;
        if( pseudominfo < 0.0 )
            pseudominfo = 0.0;
        if( lowvalue < pseudominfo ) {
            psdenominator = pow( pseudominfo, psexponent );
            if( psdenominator )
                pseudofreqweight = psnumerator / psdenominator;
        }

        if( pseudofreqweight < 0.0 ) pseudofreqweight = 0.0;
        if( 1.0 < pseudofreqweight ) pseudofreqweight = 1.0;
    //}}

        //pseudocounts do not depend on observations n (m=n*alpha/(1-alpha))
        information = 0.0;
        for( a = 0; a < NUMAA; a++ ) {
            if( LOSCORES.PROBABility( a ) <= 0.0 )
                continue;

            pseudoFreqn = 0.0;
            if( pseudofreqweight )
                for( b = 0; b < NUMAA; b++ )
                    pseudoFreqn += queryDescription->GetMatchWeightsAt( b, p ) * LOSCORES.FreqRatio( a, b );

            pseudoFreqn *= LOSCORES.PROBABility( a );
            estimated = pseudofreqweight * pseudoFreqn + ( 1.0 - pseudofreqweight ) *
                        queryDescription->GetMatchWeightsAt( a, p );
            queryDescription->SetTargetFreqnsAt( estimated, a, p );

            if( 0.0 < estimated )
                information += estimated * log( estimated / LOSCORES.PROBABility( a ));
        }

        information /= LN2;
        if( information < 0.0 )
            information = 0.0;
        // save information content
        queryDescription->SetInformationAt( information, p );
    }
}

// =========================================================================
// ComputeTargetFrequenciesMDLVar: pseudocounts are computed by the MDL theory
//
void InputMultipleAlignment::ComputeTargetFrequenciesMDLVar()
{
    if( !queryDescription )
        throw myruntime_error("Unable to estimate target probabilities.");

    double          expobscount;                //expected number of observations at certain position
    double          pseudofreqweight;           //weight for pseudo frequencies (alpha)
    const double    psnumerator = 0.65;         //numerator a in expression a/n^b
    const double    psexponent = 0.90;          //exponent b of denominator in a/n^b
    double          psdenominator;              //denominator in expression a/n^b
    const double    lowvalue = 0.0001;

    double  information = 0.0;          //information content per position
    double  pseudoFreqn = 0.0;          //pseudocount frequency
    double  estimated = 0.0;            //estimated probability
    unsigned char   a, b;
    size_t  p;

    //iterate over all query positions
    for( p = 0; p < queryDescription->size(); p++ ) {
        if( !queryDescription->GetStateAt( p ))
            continue;

    //{{
//         expobscount = queryDescription->GetMIDExpNoObservationsAt( p, PS_M );
        expobscount = queryDescription->GetNoSymbolsInExtentAt( p, PS_M );
        pseudofreqweight = 1.0;

        psdenominator = pow( expobscount, psexponent );
        if( psdenominator )
            pseudofreqweight = psnumerator / psdenominator;

        if( pseudofreqweight < 0.0 ) pseudofreqweight = 0.0;
        if( 1.0 < pseudofreqweight ) pseudofreqweight = 1.0;
    //}}

        //pseudocounts do not depend on observations n (m=n*alpha/(1-alpha))
        information = 0.0;
        for( a = 0; a < NUMAA; a++ ) {
            if( LOSCORES.PROBABility( a ) <= 0.0 )
                continue;

            pseudoFreqn = 0.0;
            if( pseudofreqweight )
                for( b = 0; b < NUMAA; b++ )
                    pseudoFreqn += queryDescription->GetMatchWeightsAt( b, p ) * LOSCORES.FreqRatio( a, b );

            pseudoFreqn *= LOSCORES.PROBABility( a );
            estimated = pseudofreqweight * pseudoFreqn + ( 1.0 - pseudofreqweight ) *
                        queryDescription->GetMatchWeightsAt( a, p );
            queryDescription->SetTargetFreqnsAt( estimated, a, p );

            if( 0.0 < estimated )
                information += estimated * log( estimated / LOSCORES.PROBABility( a ));
        }

        information /= LN2;
        if( information < 0.0 )
            information = 0.0;
        // save information content
        queryDescription->SetInformationAt( information, p );
    }
}





// =========================================================================
// CalcTFPosteriorPredictives: calculate posterior predictives 
//  probabilities of target frequencies
//
void InputMultipleAlignment::CalcTFPosteriorPredictives()
{
    mystring errstr;
    mystring preamb = "CalcTFPosteriorPredictives: ";
    if( !queryDescription )
        throw myruntime_error( preamb + "Null query description structure.");

    const size_t    noress = NUMALPH;
    const size_t    noeffress = NUMAA;
    int             ndx = 0;//index of cluster
    double          ppr = 0.0;//posterior probability
    double          tfr = 0.0;//target frequency
    const HDPbase*  hdpbase = GetHDPbase();

    if( hdpbase == NULL )
        throw myruntime_error( preamb + "Null HDP structure.");

    int     nosupcls = hdpbase->GetNoSupClusters();
    int     ctxtsz = hdpbase->GetCtxtSize();
    size_t  psize = queryDescription->GetEffectiveSize();
    size_t  p, pp;
    int     left, right, r, n;
    int     parity = ( ctxtsz & 1 ) ^ 1;
    int     hlf = ctxtsz >> 1;
    int     mid = hlf - parity;

    if( nosupcls < 1 && nosupcls != -1 )
        throw myruntime_error( preamb + "Wrong number of HDP support clusters.");
    if( ctxtsz < 1 )
        throw myruntime_error( preamb + "Wrong HDP context size.");
    if(( int )psize < ctxtsz ) {
        warning("Profile length is less than context size. Not using posteriors.");
        return;
    }
    Pslmatrix   promtx(( int)psize, (int)noeffress );//profile matrix
    Pslmatrix   ctxmtx;//context matrix
    Pslvector   ctxsck;//context stack
    Pslvector   ctxnrm(( int)( noeffress*ctxtsz ));//normal transform
    Pslvector   ppprobs;//vector of posterior probabilities
    Ivector     cindcs;//indices of clusters

    //make matrix representation of profile
    for( p = 0, pp = 0; p < queryDescription->size(); p++ ) {
        if( !queryDescription->GetStateAt( p ))
            continue;

        for( r = 0; r < noeffress; r++ ) {
            tfr = queryDescription->GetTargetFreqnsAt( r, p );
            promtx.SetValueAt(( int)pp, (int)r, tfr );
        }
        pp++;
    }

    //iterate over all query positions
    for( p = 0, pp = ( size_t )-1; p < queryDescription->size(); p++ )
    {
        if( !queryDescription->GetStateAt( p ))
            continue;

        pp++;

        right = ( int )SLC_MIN( psize-1, pp+hlf );
        left = SLC_MAX( 0, right-ctxtsz+1 );

        ctxmtx = promtx.SubMatrix( left, 0, ctxtsz, (int)noeffress );
        ctxsck = ctxmtx.Stack();
        ctxnrm.Copy( ctxsck );

        //calculate posterior predictive probabilities;
        //input vector ctxnrm will be transformed;
        hdpbase->CalcPPProbs( ctxnrm, &ppr, &ppprobs, &cindcs );

        if( /*ppprobs.GetSize() < 1 || */ppprobs.GetSize() != cindcs.GetSize())
            throw myruntime_error( preamb + "Invalid number of posterior predictive probabilities.");

        //save posteriors
        queryDescription->SetBckPPProbAt( ppr, p );
        queryDescription->SetPPProbsAt( p, ppprobs.GetVector(), cindcs.GetVector(), ppprobs.GetSize());
    }
}

// =========================================================================
// MixTargetFrequenciesHDPCtx: mix target frequencies under HDP framework
//
void InputMultipleAlignment::MixTargetFrequenciesHDPCtx()
{
    mystring errstr;
    mystring preamb = "MixTargetFrequenciesHDPCtx: ";
    if( !queryDescription )
        throw myruntime_error( preamb + "Null query description structure.");

    const size_t    noress = NUMALPH;
    const size_t    noeffress = NUMAA;
    double          pme = 0.0;//estimated probability
    double          tfr = 0.0;//target frequency
    const HDPbase*  hdpbase = GetHDPbase();

    if( hdpbase == NULL )
        throw myruntime_error( preamb + "Null HDP structure.");

    int     ctxtsz = hdpbase->GetCtxtSize();
    size_t  psize = queryDescription->GetEffectiveSize();
    size_t  p, pp;
    int     left, right, r;
    int     parity = ( ctxtsz & 1 ) ^ 1;
    int     hlf = ctxtsz >> 1;
    int     mid = hlf - parity;

    if( ctxtsz < 1 )
        throw myruntime_error( preamb + "Wrong HDP context size.");
    if(( int )psize < ctxtsz ) {
        warning("Profile length is less than context size. Not mixing.");
        return;
    }
    Pslmatrix   promtx(( int)psize, (int)noeffress );//profile matrix
    Pslmatrix   ctxmtx;//context matrix
    Pslvector   ctxsck;//context stack
    Pslvector   ctxnrm(( int)( noeffress*ctxtsz ));//normal transform
    Pslvector   mixed((int)noeffress );//mixed vector
    double      infrm;

    //make matrix representation of profile
    for( p = 0, pp = 0; p < queryDescription->size(); p++ ) {
        if( !queryDescription->GetStateAt( p ))
            continue;

        for( r = 0; r < noeffress; r++ ) {
            tfr = queryDescription->GetTargetFreqnsAt( r, p );
            promtx.SetValueAt(( int)pp, (int)r, tfr );
        }
        pp++;
    }

    //iterate over all query positions
    for( p = 0, pp = ( size_t )-1; p < queryDescription->size(); p++ )
    {
        if( !queryDescription->GetStateAt( p ))
            continue;

        pp++;

        right = ( int )SLC_MIN( psize-1, pp+hlf );
        left = SLC_MAX( 0, right-ctxtsz+1 );

        ctxmtx = promtx.SubMatrix( left, 0, ctxtsz, (int)noeffress );
        ctxsck = ctxmtx.Stack();
        ctxnrm.Copy( ctxsck );

        //mix central vector of context;
        //input vector ctxnrm will be transformed;
        //mixed will be returned in logit-normal space
        hdpbase->MixCLNVar( ctxnrm, &mixed );

        //write PME back to target frequencies
        infrm = 0.0;
        for( r = 0; r < noeffress; r++ ) {
            pme = mixed.GetValueAt(r);
            queryDescription->SetTargetFreqnsAt( pme, r, p );
            //calculate relative entropy
            if( LOSCORES.PROBABility(r) <= 0.0 )
                continue;
            if( 0.0 < pme )
                infrm += pme * log( pme / LOSCORES.PROBABility(r));
        }
        //save relative entropy
        infrm /= LN2;
        infrm = SLC_MAX( 0.0, infrm );
        queryDescription->SetInformationAt( infrm, p );
    }
}





// -------------------------------------------------------------------------
// RecalcBackgroundProbabilities: recalculate background probabilities by 
//     mixing background and target probabilities
// NOTE: PARAMETERS: tadj
// NOTE: bg probabilities are also recalculated in 
//  ProGenerator::RecalcBgProbs
//
void InputMultipleAlignment::RecalcBackgroundProbabilities()
{
    if( !queryDescription )
        throw myruntime_error("RecalcBackgroundProbabilities: No extended description vector.");

    const size_t    noress = NUMALPH;
    const size_t    noeffress = NUMAA;
    const size_t    nomstates = queryDescription->GetEffectiveSize();
    const double    effnos = GetEffNoSequences();
    const double    maxtgtobss = 20.0;
    //const double    tadj = 1.+4.*exp(-0.5*1.e-6*pow(effnos-18.,6.));//5.0
    const double    tadj = 1.;//5.
    const double    bckpseudocounts = SLC_MAX(1., tadj * effnos );//10.0;
    const double    errtol = 1.0e-3;
    char            errbuf[KBYTE];

    double          bpp[noeffress], mixp[noeffress];
    double          denm, consv;
    double          bckp, tgtp;
    double          expM, tgtpc;
    size_t          p, r, nn;

//     consv = 0.0;
//     for( r = 0; r < noeffress; r++ ) {
//         bckp = LOSCORES.PROBABility( r );
//         bpp[r] = 100.0 / effnos * bckp;
//         mixp[r] = bpp[r];
//         consv += bpp[r];
//     }
//     for( p = 0; p < queryDescription->size(); p++ ) {
//         if( !queryDescription->GetStateAt( p ))
//             continue;
//         for( r = 0; r < noeffress; r++ ) {
//             tgtp = queryDescription->GetTargetFreqnsAt( r, p );
//             mixp[r] += tgtp;
//             consv += tgtp;
//         }
//     }


    for( r = 0; r < noeffress; r++ ) {
        bckp = LOSCORES.PROBABility( r );
        bpp[r] = bckpseudocounts * bckp;
        mixp[r] = 0.0;
    }

    consv = 0.0;
    for( p = 0, nn = 0; p < queryDescription->size(); p++ ) 
    {
        if( !queryDescription->GetStateAt( p ))
            continue;

        nn++;
//         expM = queryDescription->GetMIDExpNoObservationsAt( p, PS_M );
        expM = queryDescription->GetNoSymbolsInExtentAt( p, PS_M );
        tgtpc = SLC_MIN( maxtgtobss, expM );
        denm = 1.0 /( bckpseudocounts + tgtpc );

        for( r = 0; r < noeffress; r++ ) {
            tgtp = queryDescription->GetTargetFreqnsAt( r, p );
            tgtp = ( bpp[r] + tgtpc * tgtp ) * denm;
            mixp[r] += tgtp;
            consv += tgtp;
        }
    }

    denm = 1.0;
    if( consv )
        denm /= consv;

    consv = 0.0;
    for( r = 0; r < noeffress; r++ ) {
        mixp[r] *= denm;
        consv += mixp[r];
        queryDescription->SetBackProbsAt( r, mixp[r]);
    }
    for( ; r < noress; r++ )
        queryDescription->SetBackProbsAt( r, 0.0 );
    if( consv < 1.0 - errtol || consv > 1.0 + errtol ) {
        sprintf( errbuf, "RecalcBackgroundProbabilities: Probabilities are not conserved: %g", consv );
        throw myruntime_error( errbuf );
    }
    return;

//     denm = 1.0;
//     if( nomstates )
//         denm /= ( double )nomstates;
//     for( r = 0; r < noeffress; r++ )
//         mixp[r] *= denm;
// 
//     LogitNormalErrCorr( mixp, noeffress, errtol );//will throw on error
// 
//     for( r = 0; r < noeffress; r++ )
//         queryDescription->SetBackProbsAt( r, mixp[r]);
//     for( ; r < noress; r++ )
//         queryDescription->SetBackProbsAt( r, 0.0 );
//     return;
}

// -------------------------------------------------------------------------
// CalcPosteriorProbabilities: calculate posterior (generalized target) 
//     probabilities
//
void InputMultipleAlignment::CalcPosteriorProbabilities()
{
    if( !queryDescription )
        throw myruntime_error("CalcPosteriorProbabilities: No extended description vector.");

    const size_t    noress = NUMALPH;
    const size_t    noeffress = NUMAA;
    const size_t    nomstates = queryDescription->GetEffectiveSize();
    const double    effnos = GetEffNoSequences();
    const double    errtol = 1.0e-3;
    char            errbuf[KBYTE];

    double          ppp[noeffress];
    double          denm, consv;
    double          tgtp;
    size_t          p, r;

    for( r = 0; r < noeffress; r++ )
        ppp[r] = 0.0;

    consv = 0.0;
    for( p = 0; p < queryDescription->size(); p++ )
    {
        if( !queryDescription->GetStateAt( p ))
            continue;

        for( r = 0; r < noeffress; r++ ) {
            tgtp = queryDescription->GetTargetFreqnsAt( r, p );
            ppp[r] += tgtp;
            consv += tgtp;
        }
    }

    denm = 1.0;
    if( consv )
        denm /= consv;

    consv = 0.0;
    for( r = 0; r < noeffress; r++ ) {
        ppp[r] *= denm;
        consv += ppp[r];
        queryDescription->SetPostProbsAt( r, ppp[r]);
    }
    for( ; r < noress; r++ )
        queryDescription->SetPostProbsAt( r, 0.0 );
    if( consv < 1.0 - errtol || consv > 1.0 + errtol ) {
        sprintf( errbuf, "CalcPosteriorProbabilities: Probabilities are not conserved: %g", consv );
        throw myruntime_error( errbuf );
    }
    return;
}

// -------------------------------------------------------------------------
// ComputePSSM: compute PSSM matrix given estimated target frequencies
//
void InputMultipleAlignment::ComputePSSM()
{
    if( !queryDescription )
        throw myruntime_error("Unable to compute PSSM matrix.");

    double  estimated = 0.0;    // estimated probability
    double  pssmvalue = 0.0;    // PSSM value to be saved
    double  ratiosval = 0.0;    // LOSCORES frequency ratio value
    int     scoresval = 0;      // LOSCORES matrix value
    double  lambda = LOSCORES.StatisParam( Ungapped, Lambda ); //reference lambda
    unsigned char   r, residue;
    size_t  p;

    // iterate over all query positions
    for( p = 0; p < queryDescription->size(); p++ ) 
    {
        if( !queryDescription->GetStateAt( p ))
            continue;

        residue = queryDescription->ResidueAt( p );

        for( r = 0; r < NUMALPH; r++ ) {
            estimated = queryDescription->GetTargetFreqnsAt( r, p );
            if( LOSCORES.PROBABility( r ) <= 0.0 || estimated <= 0.0 ) {
                //unable to estimate probabilities if background probabilities are zero
//                 scoresval = LOSCORES.Entry( residue, r );
                ratiosval = LOSCORES.FreqRatio( residue, r );

                if( /*scoresval == SCORE_MIN || */ratiosval <= 0.0 )
                    //unable to compute PSSM value
                    pssmvalue = SCORE_MIN;
                else
                    // this is how LOSCORES matrix values are obtained;
                    // with frequency ratios LOSCORES matrix values are 1/s * log2(ratio)
                    // (s is scaling constant for LOSCORES values);
                    //pssmvalue = log( ratiosval ) / LN2 * ScaleFactor;
                    // use precomputed values instead...
                    pssmvalue = LOSCORES.PrecomputedEntry( residue, r );

            } else {
                pssmvalue = log( estimated / LOSCORES.PROBABility( r )) / lambda;
            }
            queryDescription->SetPSSMEntryAt( pssmvalue, r, p );
        }
    }
}

// -------------------------------------------------------------------------
// ExportFrequenciesTo: export weighted observed frequencies
//
void InputMultipleAlignment::ExportFrequenciesTo( FrequencyMatrix& frequencies ) const
{
    if( !queryDescription )
        throw myruntime_error("ExportFrequenciesTo: Null query description.");

    if( !queryDescription->size() || !queryDescription->GetEffectiveSize())
        return;

    const double ( *weights )[NUMALPH];
    unsigned char   residue;
    size_t  p;

    frequencies.Reserve( queryDescription->GetEffectiveSize());

    for( p = 0; p < queryDescription->size(); p++ ) {
        if( !queryDescription->GetStateAt( p ))
            continue;

        residue = queryDescription->ResidueAt( p );
        weights = queryDescription->GetMatchWeightsAt( p );
        frequencies.Push( *weights, ( char )residue );
    }
}

// -------------------------------------------------------------------------
// ExportFrequenciesTo: export position-specific scoring matrix
//
void InputMultipleAlignment::ExportPSSMTo( LogOddsMatrix& pssm, bool appfilename ) const
{
    if( !queryDescription )
        throw myruntime_error("ExportPSSMTo: Null query description.");

    if( !queryDescription->size() || !queryDescription->GetEffectiveSize())
        return;

    size_t  noress = NUMALPH;
    size_t  queffsize = queryDescription->GetEffectiveSize();
    double  expMIDs[PS_NSTATES];
    const double ( *scores )[NUMALPH];
    double          freqweight;
    double          information;
    double          expnobs;
    double          bppprob;//background posterior predicitive
    const double*   ppprobs;//posterior predictive probabilities
    const int*      pppndxs;//indices of posteriors
    size_t          noppps;//number of p.p.probability values
    size_t          noseqs;
    unsigned char   residue;
    size_t  p;

    pssm.Reserve( queffsize );

    if( queffsize ) {
        p = -1;
        expMIDs[PS_M] = queryDescription->GetMIDExpNoObservationsAt( p, PS_M );
        expMIDs[PS_I] = queryDescription->GetMIDExpNoObservationsAt( p, PS_I );
        expMIDs[PS_D] = queryDescription->GetMIDExpNoObservationsAt( p, PS_D );
        pssm.SetMIDExpNoObservationsBeg( expMIDs );
    }

    for( p = 0; p <noress; p++ )
        pssm.SetBackProbsAt( p, queryDescription->GetBackProbsAt( p ));

    for( p = 0; p <noress; p++ )
        pssm.SetPostProbsAt( p, queryDescription->GetPostProbsAt( p ));

    for( p = 0; p < queryDescription->size(); p++ ) {
        if( !queryDescription->GetStateAt( p ))
            continue;

        residue = queryDescription->ResidueAt( p );

        //save scaled PSSM scores
        scores = queryDescription->GetPSSMVectorAt( p );
        freqweight = queryDescription->ComputeObsFrequencyWeightAt( p );
        information = queryDescription->GetInformationAt( p );
        expnobs = queryDescription->GetExpNoObservationsAt( p );
        noseqs = queryDescription->GetNoSequencesInExtentAt( p );//queryDescription->GetNoSequencesAt( p );

        if( !scores )
            throw myruntime_error("ExportPSSMTo: Null scores.");

        expMIDs[PS_M] = queryDescription->GetMIDExpNoObservationsAt( p, PS_M );
        expMIDs[PS_I] = queryDescription->GetMIDExpNoObservationsAt( p, PS_I );
        expMIDs[PS_D] = queryDescription->GetMIDExpNoObservationsAt( p, PS_D );

        bppprob = queryDescription->GetBckPPProbAt(p);
        ppprobs = queryDescription->GetPPProbsAt(p);
        pppndxs = queryDescription->GetPPPIndsAt(p);
        noppps = queryDescription->GetNoPPProbsAt(p);

        pssm.Push( *scores, (char)residue, freqweight, information, expMIDs );
    }


    Configuration   config[NoSchemes]; //no setting of filename; read or write attempt raises an error
    SetUngappedParams( config[ProcomUngapped]); //fill values with parameter values of ungapped configuration
    SetUngappedParams( config[ProcomGapped]); //not using gapped configuration, make it identical to ungpd. one

    //scale matrix to make scores distributed according to the known distribution
    ProfileMatrix   matrix( //table to scale scores and compute statistics
            (( const LogOddsMatrix& )pssm).GetVector(),
            pssm.GetResidues(),
            pssm.GetColumns(),
            config
    );

    if( matrix.GetSupportOptimFreq())
        matrix.OptimizeTargetFrequencies();

#ifdef SCALE_PROFILES
    matrix.ScaleScoringMatrix();
#else
    matrix.ComputeStatisticalParameters();
#endif


    for( int c = 0; c < pssm.GetColumns(); c++ ) {
        residue = pssm.GetResidueAt( c );
        scores = matrix.GetVectorAt( c );//scaled scores
        if( !scores )
            throw myruntime_error("ExportPSSMTo: Null scaled scores.");
        pssm.PushAt( *scores, residue, c );
    }

    pssm.SetNoSequences( GetNoSequences());
    pssm.SetEffNoSequences( GetEffNoSequences());

    pssm.SetRefLambda(      matrix.GetRefLambda());
    pssm.SetRefK(           matrix.GetRefK());
    pssm.SetLambda(         matrix.GetLambda());
    pssm.SetEntropy(        matrix.GetEntropy());
    pssm.SetK(              matrix.GetK());
    pssm.SetExpectedScore(  matrix.GetExpectedScore());

    if( appfilename )
        pssm.SetName( GetName());
    pssm.SetDescription( GetTitle());
    pssm.Finalize();
    if( GetScoAdjmentHDPCtx()) {
        //calculate posteriors after pssm.Finalize()
        pssm.CalcTFPostPredictives( GetHDPbase());
        if( GetHDPctBase())
            pssm.Calc_ctTFPostPredictives( GetHDPctBase());
    }
    pssm.CalcCtxVector();
}

// -------------------------------------------------------------------------
// ExportGapWeightsTo: exports gap weights
// -------------------------------------------------------------------------

void InputMultipleAlignment::ExportGapWeightsTo( GapScheme& gaps ) const
{
    if( !queryDescription )
        throw myruntime_error("ExportGapWeights: Null query description.");

    if( !queryDescription->size() || !queryDescription->GetEffectiveSize())
        return;

    const double ( *ttps )[P_NSTATES] = NULL;//target transition probabilities
    unsigned char   residue;
    double          weight = 0.0;
    double          delbeg = 0.0;   //deletion beginning weight
    double          delend = 0.0;   //deletion end weight
    int             delta = 0;      //length of deletion
    size_t          p;

    gaps.Reserve( queryDescription->GetEffectiveSize());

    gaps.SetOrgTrProbsBeg( queryDescription->GetTargetTranstBeg());

    for( p = 0; p < queryDescription->size(); p++ ) {
        if( !queryDescription->GetStateAt( p ))
            continue;

        residue = queryDescription->ResidueAt( p );
        ttps = queryDescription->GetTargetTranstAt( p );
        weight = queryDescription->GetGapWeightsAt( p );

        gaps.Push( ttps, weight, delbeg, delend, delta, ( char )residue );
    }
}

// -------------------------------------------------------------------------
// SetName: set name of the multiple alignment
//
void InputMultipleAlignment::SetName( const char* newname )
{
    size_t  newlength = strlen( newname );
    if( !newlength )
        throw myruntime_error( mystring( "Wrong name argument." ));

    if( name ) free( name );

    name = ( char* )malloc( newlength + 1 );
    if( !name )
        throw myruntime_error( mystring( "Not enough memory." ));

    strncpy( name, newname, newlength + 1 ); // include the terminating null symbol
}

// -------------------------------------------------------------------------
// SetTitle: replaces text of the title to the new one
// -------------------------------------------------------------------------

void InputMultipleAlignment::SetTitle( const char* newtitle )
{
    size_t  newlength = strlen( newtitle );
    if( !newlength )
        throw myruntime_error( mystring( "InputMultipleAlignment: Memory access error." ));

    if( titletext ) free( titletext );

    titletext = ( char* )malloc( newlength + 1 );
    if( !titletext )
        throw myruntime_error( mystring( "InputMultipleAlignment: Not enough memory." ));

    strncpy( titletext, newtitle, newlength + 1 ); // include the terminating null symbol
}

// -------------------------------------------------------------------------
// AppendTitle: appends text to the title if one exists; if title is null,
//     the text is simply copied to the title
// -------------------------------------------------------------------------

void InputMultipleAlignment::AppendTitle( const char* newtitle, size_t beg, size_t end )
{
    if( end < beg )
        return;

    size_t  from = 0;
    size_t  newlength = strlen( newtitle );

    if( !newlength )
        throw myruntime_error( mystring( "InputMultipleAlignment: Memory access error." ));

    if( newlength < beg )
        return;

    if( newlength < end )
        end = newlength;

    newlength = end - beg + 1;

    if( titletext ) {
        size_t  curlength = strlen( titletext );
        titletext = ( char* )realloc( titletext, curlength + newlength + 1 );
        from = curlength;
    } else
        titletext = ( char* )malloc( newlength + 1 );

    if( !titletext )
        throw myruntime_error( mystring( "InputMultipleAlignment: Not enough memory." ));

    strncpy( titletext + from, newtitle + beg, newlength );
    titletext[ from + newlength ] = 0;
}

// -------------------------------------------------------------------------

enum {
    eninfmtSTOCKHOLM,
    eninfmtFASTA,
    eninfmtNO
};
int gnINFORMAT = eninfmtNO;

// -------------------------------------------------------------------------
// ReadAlignment: read multiple sequence alignment
//
void InputMultipleAlignment::ReadAlignment( const char* filename )
{
    mystring    errsto = "STOCKHOLM: ";
    mystring    errfas = "FASTA: ";
    mystring    errstr;

    //try STOCKHOLM
    try { ReadSTOCKHOLM1( filename );
    } catch( myexception const& ex1 ) {
        errsto += ex1.what();
        if( gnINFORMAT == eninfmtSTOCKHOLM )
            throw myruntime_error( errsto );
        //try FASTA
        try { ReadFASTA( filename );
        } catch( myexception const& ex2 ) {
            errfas += ex2.what();
            throw myruntime_error( errfas );
        }
    }
}

// -------------------------------------------------------------------------
// ReadFASTA: reads multiple sequence alignment in fasta from file
//
void InputMultipleAlignment::ReadFASTA( const char* filename )
{
    FILE*   fp = fopen( filename, "r" );

    if( fp == NULL )
        throw myruntime_error("Failed to open file for reading multiple sequence alignment.");

    char                    p = 0;
    char                    buffer[KBYTE];
    int                     toread = KBYTE;
    size_t                  length = 0;     //length of sequences to be read from the file
    size_t                  cread = 1;      //number of characters read with last unformated read operation
    size_t                  begin = 0;      //beginning of the title substring
    size_t                  pos = 0;        //position-of-sequences counter
    PosDescriptionVector*   one = NewPositionVector();
    PosDescriptionVector    fake( ALLOCPOS);//fake description vector to delineate gap positions in the query (first sequence)
    bool                    title = false;

    clear();

#if 1 
    SetName( my_basename( filename ));
#else
    // or use the posix function
    SetName( basename( filename ));
#endif

    while( !feof( fp ) && cread )
    {
        cread = fread( buffer, 1, toread - 1, fp );
        buffer[cread] = 0;
        begin = 0;

        for( size_t n = 0; n < cread; n++, p = buffer[n-1] )
        {
            if(( !p || p == '\n' ) && buffer[n] == '>' ) {
                begin = n + 1;
                title = true;
                pos = 0;

                if( one->size()) {
                    if( !GetTitle())
                        throw myruntime_error("Invalid data format: No profile description.");

                    if( length ) {
                        if( one->size() != length ) {
                            fclose( fp );
                            throw myruntime_error("Invalid data format: Sequences are not equal in alignment length.");
                        }
                    } else
                        length = one->size();

                    push( one );
                    one = NewPositionVector( length );
                }
                continue;
            }

            if( title ) {
                // if about the end of buffer or next line to occur,
                //  save the title text of the first sequence
                if( n + 1 == cread || buffer[n] == '\n' )
                    if( n ) {
                        size_t  bufbeg = ( cread <= begin )? 0: begin;
                        size_t  buflen =(( buffer[n] == '\n' )? n - 1: n ) - bufbeg + 1;
                        if( !length )
                            AppendTitle( buffer, bufbeg, buflen );
                        if( GetKeepTitles() && one )
                            one->AppendDescription( buffer, bufbeg, buflen );
                    }

                if( buffer[n] == '\n' )
                    title = false;
                continue;
            }

            if( !title ) {
                if( buffer[n] == ' ' || buffer[n] == '\t' || buffer[n] == '\r' || buffer[n] == '\n' )
                    continue;

                try {
                    if( GetIgnoreGapsInQuery()) {
                        //if this is the first sequence
                        if( !length ) {
                            fake.push( HashAlphSymbol( buffer[n] ));
                        }
                        if( fake.ResidueAt( pos ) != GAP )
                            //push alignment symbols
                            one->push( HashAlphSymbol( buffer[n] ));
                    }
                    else
                        //push alignment symbols
                        one->push( HashAlphSymbol( buffer[n] ));    
                    pos++;

                } catch( myexception const& ex )
                {
                    fclose( fp );
                    mystring    tothrow( ex.what());
                    tothrow += ". Invalid data format.";
                    throw myruntime_error( tothrow, ex.eclass());
                }
            }
        }
    }

    if( !feof( fp ))
        warning("Not all data from file processed." );
    fclose( fp );

    if( one->size()) {
        if( !GetTitle())
            throw myruntime_error("Invalid data format: No description of profile." );

        if( length && one->size() != length )
            throw myruntime_error("Invalid data format: All aligned sequences must be equal in length.");
        push( one );
    }
}

// -------------------------------------------------------------------------

const char* patstrFMTSTO = "# STOCKHOLM";
const char* patstrSTOTER = "//";
const char* patstrSTOGF = "#=GF";
const char* patstrSTOGFfts[] = {"AC","DE"};
const int   noSTOGFfts = 2;
const char* patstrSTOGS = "#=GS";
const char* patstrSTOGSDE = "DE";
const char* patstrSTOGC = "#=GC";
const char* patstrSTOseqcons = "seq_cons";
const int   lenstrFMTSTO = strlen( patstrFMTSTO );

// -------------------------------------------------------------------------
// ReadSTOCKHOLM1: read multiple sequence alignment in Stockholm 1.0 format
//
void InputMultipleAlignment::ReadSTOCKHOLM1( const char* filename )
{
    FILE* fp = fopen( filename, "r" );

    if( fp == NULL )
        throw myruntime_error("Failed to open file for reading multiple sequence alignment.");

    mystring preamb = "ReadSTOCKHOLM1: ";
    size_t          length, rbts;
    const size_t    defsize = KBYTE;
    const size_t    locsize = KBYTE;
    char            locbuffer[locsize+1] = {0};
    bool            stoter = false;
    bool            descin = false;
    bool            circled = false;
    const char*     p, *pp;
    char            ch;
    int             emsg, n, sind;
    bool            statesempty;
    mystring        seqcons, states;
    mystring        line, errstr, str, *sname;
    const mystring* svn;
    SimpleVector    svhdngs( KBYTE );//headings
    SimpleVector    svnames( KBYTE );//names
    SimpleVector    svaseqs( KBYTE );//alignment sequences
    PosDescriptionVector*   svs, *seq;

    clear();
    states.reserve(TIMES4(KBYTE));

#if 1 
    SetName( my_basename( filename ));
#else
    // or use the posix function
    SetName( basename( filename ));
#endif

    //read header of file
    if(( emsg = skip_comments( fp, locbuffer, SLC_MIN( locsize, lenstrFMTSTO+1 ), &length, 0 )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw myruntime_error("Not the STOCKHOLM format.");

    if(( p = strstr( locbuffer, patstrFMTSTO )) == NULL )
        throw myruntime_error("Not the STOCKHOLM format.");

    //STOCKHOLM format
    gnINFORMAT = eninfmtSTOCKHOLM;
    sind = 0;
    circled = false;

    try {
        while( !feof( fp )) {
            //read full line
            if(( emsg = skip_comments( fp, line, 0 )) != 0 )
                throw myruntime_error( TranslateReadError( emsg ));

            if( feof( fp ) || line.empty())
                continue;

            if( line[0] == '\n' || line[0] == '\r')
                continue;

            for( n = (int)line.length() - 1; 
                 0 <= n && ( line[n] == '\n' || line[n] == '\r' ); n-- ) line[n] = 0;

            //save family description if available 
            if(( p = (char*)strstr( line.c_str(), patstrSTOGF )) != NULL ) {
                p += strlen( patstrSTOGF );
                for( ; *p == ' ' || *p == '\t' ; p++ );
                for( n = 0; n < noSTOGFfts; n++ ) {
                    if( strncmp( p, patstrSTOGFfts[n], strlen( patstrSTOGFfts[n])) == 0 ) {
                        pp = p + strlen( patstrSTOGFfts[n]);
                        for( ; *pp == ' ' || *pp == '\t' ; pp++ );
                        if( strlen( pp ))
                            str = mystring( pp );
                        if( !str.empty()) {
                            if( GetTitle()) 
                                AppendTitle(" ");
                            AppendTitle( str.c_str());
                            descin = true;
                        }
                        break;
                    }
                }
                continue;
            }

            //save description if not saved 
            if(( p = ( char* )strstr( line.c_str(), patstrSTOGS )) != NULL ) {
                pp = p + strlen( patstrSTOGS );
                for( p = pp; *p == ' ' || *p == '\t' ; p++ );
                for( pp = p, n = 0; *pp && *pp != ' ' && *pp != '\t'; pp++, n++ );
                str = mystring( p, n );
                for( ; *pp == ' ' || *pp == '\t' ; pp++ );
                //if(( p = strstr( pp, patstrSTOGSDE )) != NULL ) {
                //    pp = p + strlen( patstrSTOGSDE );
                for( ; *pp && *pp != ' ' && *pp != '\t' ; pp++ );//jump over <feature>
                    for( ; *pp == ' ' || *pp == '\t' ; pp++ );
                    if( strlen( pp ))
                        str += ' ' + mystring( pp );
                    if( !str.empty() && !descin )
                        AppendTitle( str.c_str());
                    svhdngs.Push( new mystring( str ));
                //}
                descin = true;
                continue;
            }

            //if end of MSA exit the loop
            if( strlen( patstrSTOTER ) <= strlen( line.c_str()))
                if( strncmp( line.c_str(), patstrSTOTER, strlen( patstrSTOTER )) == 0 ) {
                    stoter = true;
                    if( !feof( fp ))
                        if(( emsg = skip_comments( fp, line, 0 )) != 0 )
                            throw myruntime_error( TranslateReadError( emsg ));
                    break;
                }

            //read consensus sequence seq_cons
            if(( p = ( char* )strstr( line.c_str(), patstrSTOGC )) != NULL ) {
                p += strlen( patstrSTOGC );
                for( ; *p == ' ' || *p == '\t' ; p++ );
                if( strncmp( p, patstrSTOseqcons, strlen( patstrSTOseqcons )) == 0 ) {
                    p += strlen( patstrSTOseqcons );
                    for( ; *p == ' ' || *p == '\t' ; p++ );
                    for( pp = p, n = 0; *pp && *pp != ' ' && *pp != '\t'; pp++, n++ );
                    if( 0 < n )
                        seqcons += mystring( p, n );
                }
                continue;
            }

            if( line[0] == '#' )
                continue;

            //process alignment sequences
            p = line.c_str();
            for( ; *p == ' ' || *p == '\t' ; p++ );
            for( pp = p, n = 0; *pp && *pp != ' ' && *pp != '\t'; pp++, n++ );
            if( n <= 0 )
                throw myruntime_error("Sequence name expected.");

            //check name
            sname = new mystring( p, n );
            if( sname == NULL )
                throw myruntime_error("Not enough memory.");

            if( svnames.GetSize()) {
                if(( svn = ( mystring* )svnames.GetValueAt( 0 )) == NULL )
                    throw myruntime_error("Memory access error.");
                if( *svn == *sname ) {
                    sind = 0;
                    circled = true;
                }
            }

            if( circled ) {
                if( svnames.GetSize() <= sind || svaseqs.GetSize() <= sind )
                    throw myruntime_error("Unexpected sequences.");
                if(( svn = ( mystring* )svnames.GetValueAt( sind )) == NULL  || 
                   ( svs = ( PosDescriptionVector* )svaseqs.GetValueAt( sind )) == NULL )
                    throw myruntime_error("Memory access error.");
                if( *svn != *sname )
                    throw myruntime_error("Inconsistency.");
                delete sname;
                sname = NULL;
            }
            else {
                svnames.Push( sname );
                sname = NULL;
                if(( svs = NewPositionVector( defsize )) == NULL )
                    throw myruntime_error("Not enough memory.");
                push( svs );
                svaseqs.Push( svs );
                if(( svn = ( mystring* )svhdngs.GetValueAt( sind )) == NULL )
                    throw myruntime_error("Memory access error.");
                if( GetKeepTitles())
                    svs->AppendDescription( svn->c_str(), 0, svn->length());
                delete svn;
                svhdngs.SetValueAt( sind, svn = NULL );
            }

            sind++;
            statesempty = states.empty();

            //save alignment sequence
            for( p = pp; *p == ' ' || *p == '\t' ; p++ );
            for( pp = p; *pp && *pp != ' ' && *pp != '\t'; pp++ ) {
                if( 'A' <= *pp && *pp <= 'Z' )
                    ch = *pp;
                else if( 'a' <= *pp && *pp <= 'z' )
                    ch = *pp;
                else if( *pp == '-' || *pp == '.' || *pp == '_' || *pp == '~' )
                    ch = '-';
                else
                    throw myruntime_error("Unrecognized alignment symbol.");
                svs->push( HashAlphSymbol(ch));
                if( statesempty ) {
                    if(('A' <= *pp && *pp <= 'Z')|| *pp == '-')
                        states.append('1');//match or delete state
                    else
                        states.append('0');
                }
            }
        }
    } catch( myexception const& ex ) {
        errstr = preamb + ex.what();
    }

    if( sname ) {
        delete sname;
        sname = NULL;
    }
    for( n = 0; n < ( int )svhdngs.GetSize(); n++ )
        if(( svn = ( mystring* )svhdngs.GetValueAt( n )) != NULL ) {
            delete svn;
            svhdngs.SetValueAt( n, NULL );
        }
    for( n = 0; n < ( int )svnames.GetSize(); n++ )
        if(( svn = ( mystring* )svnames.GetValueAt( n )) != NULL ) {
            delete svn;
            svnames.SetValueAt( n, NULL );
        }

    if( stoter && !feof( fp ))
        warning(( preamb + "Only a single Multiple sequence alignment processed.").c_str());
    fclose( fp );

    if( !stoter && errstr.empty())
        warning(( preamb + "Unterminated file.").c_str());

    if( !errstr.empty())
        throw myruntime_error( errstr );

    if( !GetTitle())
        warning(( preamb + "No description of Multiple sequence alignment.").c_str());

    if( 0 < size()) {
        seq = SequenceAt(0);
        if( seq == NULL )
            throw myruntime_error( preamb + "Memory access error.");
    }

    //check alignment sequence lengths
    for( n = 1; n < size(); n++ ) {
        svs = SequenceAt(n);
        if( !svs || !seq )
            throw myruntime_error( preamb + "Memory access error.");
        if( svs->size() != seq->size())
            throw myruntime_error( preamb + "Length of aligned sequences differs.");
    }

    //check consensus sequence length, format and save it
//     if( seq && !seqcons.empty()) {
//         if( seq->size() != seqcons.length())
//             throw myruntime_error( preamb + "Inconsistent length of consensus sequence.");
//         //move the first sequence to the end
//         push( seq );
//         //allocate space for consensus sequence in the first place
//         if(( svs = NewPositionVector( seq->size())) == NULL )
//             throw myruntime_error( preamb + "Not enough memory.");
//         SetSequenceAt( 0, svs );
//         if( GetKeepTitles())
//             svs->AppendDescription( GetTitle(), 0, strlen( GetTitle()));
//         //translate and save consensus sequence
//         TranslateSTOConsSequence( seqcons, svs );
//     }

    //check state sequence length, format and save it
    //NOTE: should be consistent with ssp2.pl
    if( seq ) {
        if( seq->size() != states.length())
            throw myruntime_error( preamb + "Inconsistent length of state sequence.");
        //move the first sequence to the end
        push( seq );
        //allocate space for consensus sequence in the first place
        if(( svs = NewPositionVector( seq->size())) == NULL )
            throw myruntime_error( preamb + "Not enough memory.");
        SetSequenceAt( 0, svs );
        if( GetKeepTitles())
            svs->AppendDescription( GetTitle(), 0, strlen( GetTitle()));
        //process state sequence, make it support sequence of profile
        TranslateSTOStates( states, svs );
    }
}

// -------------------------------------------------------------------------
// TranslateSTOConsSequence: translate STOCKHOLM consensus sequence
//
//CntCmp: private function: count comparator
static inline
int CntCmp( const void* vector, size_t n1, size_t n2 )
{
    const size_t* vec = ( const size_t* )vector;
    if( vec == NULL )
        return 0;
    //[n2]-[n1], to sort in descending order
    return (int)(vec[n2] - vec[n1]);
}
//
void InputMultipleAlignment::TranslateSTOConsSequence( 
    const mystring& seqcons, PosDescriptionVector* seq )
{
    int     n, s, slcted;
    char    ch, cd; 
    unsigned char rr;
    char    bufstr[KBYTE];
    size_t  resdst[NUMAA];
    size_t  srtind[NUMAA];
    size_t  chind;
    MTRng   rng;
    PosDescriptionVector* svs;
    mystring preamb = "TranslateSTOConsSequence: ";

    mystring    code("olach-p+sut");
    const int   lencode = code.length();
    mystring    codeRES[lencode]; 
    n = 0;
    codeRES[n++] = "ST";              //o, alcohol
    codeRES[n++] = "ILV";             //l, aliphatic
    codeRES[n++] = "FHWY";            //a, aromatic
    codeRES[n++] = "DEHKR";           //c, charged
    codeRES[n++] = "ACFGHIKLMRTVWY";  //h, hydrophobic
    codeRES[n++] = "DE";              //-, negative
    codeRES[n++] = "CDEHKNQRST";      //p, polar
    codeRES[n++] = "HKR";             //+, positive
    codeRES[n++] = "ACDGNPSTV";       //s, small
    codeRES[n++] = "AGS";             //u, tiny
    codeRES[n++] = "ACDEGHKNQRST";    //t, turnlike

    if( seq == NULL )
        throw myruntime_error( preamb + "Memory access error.");
    seq->clear();

    rng.Set((unsigned long)(size_t)(&rng));// +(unsigned long)tm );

    for( n = 0; n < seqcons.length(); n++ ) {
        ch = seqcons[n];
        if( 'A' <= ch && ch <= 'Z' ) {
            seq->push( HashAlphSymbol(ch));
            continue;
        }
        if( '.' == ch ) {
            //any residue, gap
            seq->push( HashAlphSymbol('-'));
            continue;
        }
        if(( chind = code.find(ch)) == mystring::npos ) {
            sprintf( bufstr, "Unrecognized consensus symbol: %c.", ch );
            throw myruntime_error( preamb + bufstr );
        }
        //find out residue distribution at the position
        memset( resdst, 0, NUMAA*sizeof(size_t));
        //the 0th position reserved for seq!
        for( s = 1; s < size(); s++ ) {
            svs = SequenceAt(s);
            if( svs == NULL )
                throw myruntime_error( preamb + "Memory access error.");
            rr = svs->ResidueAt(n);
            if( NUMAA <= rr )
                continue;
            resdst[rr]++;
        }
        //sort by count
        HeapSortInd( srtind, resdst, NUMAA, CntCmp );
        //select a residue observed max number of times
        slcted = 0;
        for( s = 0; s < NUMAA; s++ ) {
            rr = ( unsigned char )srtind[s];
            if( NUMAA <= rr )
                throw myruntime_error( preamb + "Memory access error: res index.");
            if( resdst[rr] < 1 )
                break;
            if( codeRES[chind].find(DehashCode(rr)) != mystring::npos ) {
                seq->push( rr );
                slcted = 1;
                break;
            }
        }
        //if consensus symbol is inconsistent with maximally observed residues,
        // select randomly from the class
        if( !slcted ) {
            s = int( rng.GetDouble()*(double)(codeRES[chind].length()));
            rr = HashAlphSymbol( codeRES[chind][s]);
            seq->push( rr );
            sprintf( bufstr, "Unmatched consensus symbol '%c' at position %d.", ch, n );
            warning(( preamb + bufstr ).c_str());
        }
    }
}

// -------------------------------------------------------------------------
// TranslateSTOStates: process state sequence for making the support 
//  sequence of profile
//
void InputMultipleAlignment::TranslateSTOStates( 
    const mystring& states, PosDescriptionVector* seq )
{
    //minimum fraction of delete state symbols per position to consider 
    //  position to be in delete state
    const double frcDEL = 0.7;
    //calculate delete states using frcDEL as a criterion for the minimum 
    //  fraction of gaps per column;
    //  if false, all match (or delete) states are considered match states 
    //  and the delete state probability for each column is calculated while 
    //  processing the profile; otherwise (true), delete states will become
    //  insertions during processing the profile
    const bool   cbCALCDELSTATES = false;
    //
    int     n, s, slcted;
    unsigned char rr;
    char    ch; 
    char    bufstr[KBYTE];
    size_t  resdst[NUMALPH];
    size_t  srtind[NUMALPH];
    size_t  restot;
    PosDescriptionVector* svs;
    mystring preamb = "TranslateSTOStates: ";

    if( seq == NULL )
        throw myruntime_error( preamb + "Memory access error.");
    seq->clear();

    for( n = 0; n < states.length(); n++ ) {
        ch = states[n];
        if( '1' != ch ) {
            //any residue, gap
            seq->push( HashAlphSymbol('-'));
            continue;
        }
        //find out residue distribution at the position
        memset( resdst, 0, NUMALPH*sizeof(size_t));
        restot = 0;
        //the 0th position reserved for seq!
        for( s = 1; s < size(); s++ ) {
            svs = SequenceAt(s);
            if( svs == NULL )
                throw myruntime_error( preamb + "Memory access error.");
            rr = svs->ResidueAt(n);
            if( NUMALPH <= rr )
                continue;
            resdst[rr]++;
            restot++;
        }
        //sort by count
        HeapSortInd( srtind, resdst, NUMALPH, CntCmp );
        //if most observed is GAP
        rr = (unsigned char)srtind[0];
        if( GAP == rr && restot && cbCALCDELSTATES ) {
            if( frcDEL < (double)resdst[rr]/(double)restot ) {
                seq->push( HashAlphSymbol('-'));
                continue;
            }
        }
        //select a residue observed max number of times
        slcted = 0;
        for( s = 0; s < NUMALPH; s++ ) {
            rr = ( unsigned char )srtind[s];
            if( GAP == rr )
                continue;
            if( NUMALPH <= rr )
                throw myruntime_error( preamb + "Memory access error: res index.");
            if( resdst[rr] < 1 )
                break;
            if( NUMAA <= rr && s+1 < NUMALPH && 
                srtind[s+1] < NUMAA && resdst[rr] == resdst[srtind[s+1]])
                //count of actual amino acid is the same
                continue;
            seq->push( rr );
            slcted = 1;
            break;
        }
        //if no residue selected
        if( !slcted ) {
            //delete state columns containing only GAP may occur
            seq->push( HashAlphSymbol('-'));
            continue;
            sprintf( bufstr, "Unmatched state at position %d.", n );
            throw myruntime_error( preamb + bufstr );
        }
    }
}

// -------------------------------------------------------------------------
// PutAlignment: writes multiple sequence alignment to the file
// -------------------------------------------------------------------------

void InputMultipleAlignment::PutAlignment( const char* filename )
{
    FILE*   fp = stdout;

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw myruntime_error(
          "Failed to open file for writing multiple alignment." );

    PosDescriptionVector* one = NULL;
    size_t  i, p;

    for( i = 0; i < size(); i++ ){
        one = SequenceAt(i);
        if( i == 0 && queryDescription )
            one = queryDescription;

        if( !one || !one->GetUsed())
            continue;

        if( !i )
            fprintf( fp, ">%s\n", GetTitle());
        else {
            if( GetKeepTitles() && one->GetDescription())
                fprintf( fp, ">%s\n", one->GetDescription());
            else
                fprintf( fp, ">Sequence %d\n", i+1 );
        }

        if( queryDescription ) {
            for( p = 0; p < queryDescription->size(); p++ ) {
                if( !queryDescription->IsUsedAt( p ))
                    continue;
                putc( DehashCode( one->ResidueAt(p)), fp );
            }
        }
        else
            for( p = 0; p < one->size(); p++ ) {
                putc( DehashCode( one->ResidueAt(p)), fp );
            }
        putc( '\n', fp );
    }

    if( fp != stdout )
        fclose( fp );
}

// -------------------------------------------------------------------------
// OutputWeightedFrequencies: outputs observed weighted frequencies
// -------------------------------------------------------------------------

void InputMultipleAlignment::OutputWeightedFrequencies( const char* filename )
{
    FILE*   fp = stdout;

    if( !queryDescription )
        throw myruntime_error( mystring(
        "Unable to print observed weighted frequencies." ));

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw myruntime_error(
            mystring( "Failed to open file for writing." ));

    queryDescription->PrintMatchWeights( fp );

    if( fp != stdout )
        fclose( fp );
}

// -------------------------------------------------------------------------
// OutputPSSM: outputs PSSM matrix
// -------------------------------------------------------------------------

void InputMultipleAlignment::OutputPSSM( const char* filename )
{
    FILE*   fp = stdout;

    if( !queryDescription )
        throw myruntime_error( mystring(
        "Unable to print PSSM matrix." ));

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw myruntime_error(
            mystring( "Failed to open file for writing." ));

    queryDescription->PrintPSSMatrix( fp );

    if( fp != stdout )
        fclose( fp );
}

// -------------------------------------------------------------------------
// OutputSuppressedProfile: outputs PSSM matrix and weighted frequencies
// -------------------------------------------------------------------------

void InputMultipleAlignment::OutputSuppressedProfile( const char* filename )
{
    FILE*   fp = stdout;

    if( !queryDescription )
        throw myruntime_error( mystring(
        "Unable to print Profile information." ));

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw myruntime_error(
            mystring( "Failed to open file for writing." ));

    queryDescription->PrintSuppressedPSSMandWeights( fp );

    if( fp != stdout )
        fclose( fp );
}

// -------------------------------------------------------------------------
// OutputProfile: outputs full profile information in the text format
// -------------------------------------------------------------------------

void InputMultipleAlignment::OutputProfile( const char* filename )
{
    FILE*   fp = stdout;

    if( !queryDescription )
        throw myruntime_error( mystring(
        "Unable to print Profile information." ));

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw myruntime_error(
            mystring( "Failed to open file for writing." ));

    if( GetName()) {
        fprintf( fp, "Multiple alignment, %s\n", GetName());
    }
    if( GetTitle()) {
        fprintf( fp, "First sequence description, %s\n\n", GetTitle());
    }

    queryDescription->PrintProfile( fp );

    if( fp != stdout )
        fclose( fp );
}

