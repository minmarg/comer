/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include <stdlib.h>

#include "segdata.h"
#include "SEGAbstract.h"


// /////////////////////////////////////////////////////////////////////////
// CLASS LinearSearchStructure
//
// constructor
//

LinearSearchStructure::LinearSearchStructure( TComparator compfunc, TEqComparer eqfunc, size_t size )
:   BinarySearchStructure( compfunc, size, false /*do not keep duplicates*/),
    eqcomparer( eqfunc )
{
    if( eqcomparer == NULL )
        throw myruntime_error(
            mystring( "LinearSearchStructure: No comparison function provided." ));
}

// destructor
//

LinearSearchStructure::~LinearSearchStructure()
{
}

// -------------------------------------------------------------------------
// Find: find a key in the structure and set the location of the found key;
//     if no key is found then set the location to point at the last element
// -------------------------------------------------------------------------

bool LinearSearchStructure::Find( const void* key, int* loc ) const
{
    int n = 0;

    for( n = 0; n < ( int )GetSize(); n++ )
        if(( *eqcomparer )( GetValueAt( n ), key )) {
            if( loc )
                *loc = n;
            return true;
        }

    if( loc )
        *loc = n;
    return false;
}

// /////////////////////////////////////////////////////////////////////////
// CLASS BinaryModifier
//
// constructor
//

BinaryModifier::BinaryModifier( TComparator compfunc, size_t size )
:   BinarySearchStructure( compfunc, size, true /*keep_duplicates*/)
{
}

// destructor
//

BinaryModifier::~BinaryModifier()
{
}

// -------------------------------------------------------------------------
// IncValue: increases count value; if several count elements with value
//     exist, modifies the last but less than the next element
// -------------------------------------------------------------------------

bool BinaryModifier::IncValue( tcount value )
{
    size_t  loc;
    tcount  locvalue = 0;
    int     location = -1;

    if( !Find(( void* ) value, &location ) || location < 0 )
        return false;

    //value found, proceed
    for( loc = location; loc < GetSize() && (( tcount )GetValueAt( loc )) == value; loc++ );

    loc--;

    //modification retains the ordinance of structure
    locvalue = ( tcount )GetValueAt( loc ) + 1;
    SetValueAt( loc, ( void* ) locvalue );

    return true;
}

// -------------------------------------------------------------------------
// DecValue: reduces count value; if several count elements with value
//     exist, modifies the first but greater than the preceeding element
// -------------------------------------------------------------------------

bool BinaryModifier::DecValue( tcount value )
{
    size_t  loc;
    tcount  locvalue = 0;
    int     location = -1;

    if( value == 0 )
        return false;

    if( !Find(( void* ) value, &location ) || location < 0 )
        return false;

    //value found, proceed
    for( loc = location; ( locvalue = ( tcount )GetValueAt( loc )) == value && loc; loc-- );

    if( locvalue != value )
        loc++;

    //modification retains the ordinance of structure
    locvalue = ( tcount )GetValueAt( loc ) - 1;
    SetValueAt( loc, ( void* ) locvalue );

    return true;
}

// /////////////////////////////////////////////////////////////////////////
// CLASS Window
//
// constructor:
//
// NOTE: address is supposed to point at a vector of pointers void*
//

Window::Window(
    TVerifier verffunc,
    TComparator comparator,
    TEqComparer eqcomparer,
    const void* address,
    size_t len,
    size_t szalpha )
:
    verifier( verffunc ),
    window( address ),
    length( len ),
    szalphabet( szalpha ),

    entropy( -1.0 ),
    composition( NULL ),
    compaddresses( NULL ),
    state( NULL )
{
    if( verifier == NULL )
        throw myruntime_error( mystring( "Window: No instantiation provided." ));

    if( length == 0 )
        throw myruntime_error( mystring( "Window: Window of length 0." ));

    size_t  reservation = ( szalphabet == 0 )? length: szalphabet;

    composition = new SimpleVector( reservation );

    if( eqcomparer )
        compaddresses = new LinearSearchStructure( comparator, eqcomparer, reservation );
    else
        compaddresses = new BinarySearchStructure( comparator, reservation, false/*do not keep duplicates*/ );

    state = new BinaryModifier( SimpleIntComparer, reservation );

    if( composition == NULL || compaddresses == NULL || state == NULL )
        throw myruntime_error( mystring( "Window: Not enough memory." ));
}

// destructor
//

Window::~Window()
{
    if( composition )
        delete composition;

    if( compaddresses )
        delete compaddresses;

    if( state )
        delete state;
}

// -------------------------------------------------------------------------
// Initialize: computes composition, arranges state vector, and computes
//     entropy
// -------------------------------------------------------------------------

void Window::Initialize()
{
    CompositionOn();
    StateOn();
    ComputeEntropy();
}

// -------------------------------------------------------------------------
// CompositionOn: computes composition of window
// -------------------------------------------------------------------------

void Window::CompositionOn()
{
    for( size_t pos = 0; pos < GetLength(); pos++ )
        RefreshComposition( GetWindowAt( pos ));
}

// -------------------------------------------------------------------------
// StateOn: arranges ordered state vector
// -------------------------------------------------------------------------

void Window::StateOn()
{
    for( size_t n = 0; n < GetCompositionSize(); n++ )
        PushStateValue( GetCompositionAt( n ));
}

// -------------------------------------------------------------------------
// ComputeEntropy: computes entropy of window according to information in
//     state vector;
// there is actually a need to compute sum of
//     P = n[i]/L (log2 n[i]/L),
// where n[i] is a value from state vector; we have precomputed values of
//     p(n[i]) = n[i] log2 n[i];
// so the needed P = 1/L (p(n[i]) - n[i]/L p(L)).
// The sum of P is
//     1/L (SUM p(n[i]) - p(L)), since SUM n[i] == L
// -------------------------------------------------------------------------

void Window::ComputeEntropy()
{
    double  partial = 0.0;
    double  winentr = 0.0;
    tcount  stateval = 0;
    tcount  L = 0;

    for( size_t n = 0; n < GetStateSize(); n++ ) {
        L += stateval = GetStateAt( n );
        partial += PRT_ENTROPIES.GetValueOf(( size_t )stateval );
    }

    if( 0 < L )
        winentr = ( partial - PRT_ENTROPIES.GetValueOf(( size_t )L )) / ( double )L;
    SetEntropy( winentr );
}

// -------------------------------------------------------------------------
// LogProbability: computes logarithm of probability of this window;
//     probability is equal to product of multinomial term and number of
//     compositions divided by total number of available sequences:
//     P = M * F / N^L,
// where M is multinomial term, F is compositional term, N is size of
// alphabet, L is length of window
// -------------------------------------------------------------------------

double Window::LogProbability() const
{
    size_t  N = GetAlphabetSize();
    double  logP = 0.0;

    if( N == 0 )
        //get window length then
        N = GetLength();

    logP = LogMultinomial() + LogNoCompositions( N ) - GetLength() * LOGARITHMS.GetValueOf( N );
    return logP;
}

// -------------------------------------------------------------------------
// LogMultinomial: computes logarithm of multinomial term to final
//     probability of window composition; the result gives the number of
//     distinct sequences given composition vector
// -------------------------------------------------------------------------

double Window::LogMultinomial() const
{
    double  denominator = 0.0;
    tcount  stateval = 0;
    tcount  L = 0;

    for( size_t n = 0; n < GetStateSize(); n++ ) {
        L += stateval = GetStateAt( n );
        denominator += LOG_GAMMA.GetValueOf(( size_t )stateval );
    }

//     return LOG_GAMMA.GetValueOf(( size_t )L ) - denominator;
    return LOG_GAMMA.GetValueOf( GetLength()) - denominator;
}

// -------------------------------------------------------------------------
// LogNoCompositions: computes logarithm of number of distinct compositions
//     the state vector can have; the number is computed as
//     N! / (PRODUCT r[i] (N-K)!),
// N is size of alphabet, PRODUCT is to C, number of distinct counts in
// state vector, r[i] is repetition number of count i in state vector;
//     SUM r[i] == K;
// K is size of state vector (without the 0 elements).
// -------------------------------------------------------------------------

double Window::LogNoCompositions( size_t N ) const
{
    if( N == 0 )
        //get window length then
        N = GetLength();

    double  result = LOG_GAMMA.GetValueOf( N );

    if( GetStateSize() == 0 )
        return result;

    double  denominator = 0.0;
    size_t  repnum = 1;
    size_t  K = 0;
    tcount  stateval = 0;
    tcount  statevalp1 = GetStateAt( 0 );
    tcount  L = 0;

    for( size_t n = 0; n < GetStateSize(); n++ ) {
        L += stateval = statevalp1;

        if( n < GetStateSize() - 1 )
            statevalp1 = GetStateAt( n + 1 );

        if( stateval == 0 )
            continue;
        K++;

        if( n < GetStateSize() - 1 && stateval == statevalp1 ) {
            repnum++;
            continue;
        }

        denominator += LOG_GAMMA.GetValueOf( repnum );
        repnum = 1;
    }

#ifdef __DEBUG__
    if( N < K )
        throw myruntime_error( mystring( "Window: Wrong size of alphabet." ));
#endif

    if( 0 < K ) {
        denominator += LOG_GAMMA.GetValueOf( N - K );
        result -= denominator;
    }
    return result;
}

// -------------------------------------------------------------------------
// SafeShift: safe shift of window by 1 position
// -------------------------------------------------------------------------

void Window::SafeShift()
{
    if( Uninitialized()) {
        Initialize();
        return;
    }

    Shift();
}

// -------------------------------------------------------------------------
// Shift: shifts window by 1 position to the right; refreshes
//     composition and state vectors, recomputes entropy
// -------------------------------------------------------------------------

void Window::Shift()
{
    if( GetLength() == 0 )
        return;

    size_t      address = 0;
    bool        appended = false;
    tcount      begvalue;
    tcount      newvalue;
    const void* value = NULL;
    bool        changed = false;

    value = GetWindowAt( 0 );

    if(( *GetVerifier())( value ))
    {
        //0th element must be found
        if( !FindCompAddress( value, &address ))
            throw myruntime_error( mystring( "Window: Composition value not found." ));

        begvalue = GetCompositionAt( address ); //composition value at the beginning of window

        if( begvalue == 0 )
            throw myruntime_error( mystring( "Window: Illegal composition value found." ));

        DecCompositionAt( address );    //decrease composition value since sliding window will lose its beginning
        DecStateValue( begvalue );      //decrease state of the beginning value
        changed = true;
    }


    //reinitialize window at the next position
    window = GetWindowAt( 1 );


    value = GetWindowAt( GetLength() - 1 );

    if(( *GetVerifier())( value ))
    {
        if( PushCompAddress( value, &address )) {
            //the first occurence of the element
            InsertCompositionAt( address, 1 );
            PushStateValue( 1 );
        } else {
            //element found, get its composition value
            newvalue = GetCompositionAt( address ); //composition value at the boundary of window

            //newvalue can be 0 if it was met before
            IncCompositionAt( address );
            IncStateValue( newvalue );  //increase state of the newly met value
        }
        changed = true;
    }

    if( changed )
        ComputeEntropy();
}

// =========================================================================
// Methods for testing
// Print: prints window composition, state, and entropy
//

void Window::Print( FILE* fp )
{
    if( fp == NULL )
        return;

    fprintf( fp, "\nWindow\n Composition:\n " );
    for( size_t n = 0; n < GetCompositionSize(); n++ )
        fprintf( fp, " %3d", GetCompositionAt( n ));

    fprintf( fp, "\n State:\n " );
    for( size_t n = 0; n < GetStateSize(); n++ )
        fprintf( fp, " %3d", GetStateAt( n ));

    fprintf( fp, "\n Entropy:\n  %.3f", GetEntropy());
    fprintf( fp, "\n Log-multinomial term:\n  %.4f", LogMultinomial());
    fprintf( fp, "\n Log-No-compositions term:\n  %.4f", LogNoCompositions());
    fprintf( fp, "\n Log-probability:\n  %.4f\n\n", LogProbability());
}

// /////////////////////////////////////////////////////////////////////////
// CLASS SEGAbstract
//
// constructor:
//
// NOTE: address is supposed to point at a vector of pointers void*
//

SEGAbstract::SEGAbstract(
    TValidator  valdfunc,
    TVerifier   verffunc,
    TComparator compfunc,
    TEqComparer eqcomfunc,
    void*   address,
    size_t  runlen,
    size_t  winlen,
    double  lowent,
    double  highent,
    size_t  maxdiff,
    size_t  szalpha )
:
    validator( valdfunc ),
    verifier( verffunc ),
    comparator( compfunc ),
    eqcomparer( eqcomfunc ),
    runaddress( NULL ),
    runlength( 0 ),
    winlength( winlen ),
    szalphabet( szalpha ),
    lowentropy( lowent ),
    highentropy( highent ),
    maxdiffextent( maxdiff ),

    inverted( false ),

    entropies( NULL ),
    segments( DEF_NO_SEGMENTS )
{
    if( validator == NULL || verifier == NULL || comparator == NULL )
        throw myruntime_error( mystring( "SEGAbstract: No instantiation provided." ));

    SetRunAddress( address, runlen );
}

// constructor: alternative
//
SEGAbstract::SEGAbstract(
    TValidator  valdfunc,
    TVerifier   verffunc,
    TComparator compfunc,
    TEqComparer eqcomfunc,
    size_t  winlen,
    double  lowent,
    double  highent,
    size_t  maxdiff,
    size_t  szalpha )
:
    validator( valdfunc ),
    verifier( verffunc ),
    comparator( compfunc ),
    eqcomparer( eqcomfunc ),
    runaddress( NULL ),
    runlength( 0 ),
    winlength( winlen ),
    szalphabet( szalpha ),
    lowentropy( lowent ),
    highentropy( highent ),
    maxdiffextent( maxdiff ),

    inverted( false ),

    entropies( NULL ),
    segments( DEF_NO_SEGMENTS )
{
    if( validator == NULL || verifier == NULL || comparator == NULL )
        throw myruntime_error( mystring( "SEGAbstract: No instantiation provided." ));
}

// destructor
//

SEGAbstract::~SEGAbstract()
{
    Destroy();
}

// -------------------------------------------------------------------------
// SetRunAddress: resets run at the address given
// -------------------------------------------------------------------------

void SEGAbstract::SetRunAddress( void* address, size_t runlen )
{
    Allocate( runlen );

    runaddress = address;
    *const_cast<size_t*>( &runlength )= runlen;

    if( runlength < winlength )
        warning( "SEGAbstract: Ineffective: Window length is greater than run length." );
}

// -------------------------------------------------------------------------
// Allocate: allocates memory required for entropy
// -------------------------------------------------------------------------

void SEGAbstract::Allocate( size_t newlen )
{
    if( newlen == 0 )
        throw myruntime_error( mystring( "SEGAbstract: Run of length 0." ));

    entropies = ( double* )malloc( sizeof( double ) * newlen );

    if( entropies == NULL )
        throw myruntime_error( mystring( "SEGAbstract: Not enough memory." ));

    memset( entropies, 0, sizeof( double ) * newlen );
}

// -------------------------------------------------------------------------
// Destroy: destroys allocated memory
// -------------------------------------------------------------------------

void SEGAbstract::Destroy()
{
    if( entropies ) {
        free( entropies );
        entropies = NULL;
    }
}

// -------------------------------------------------------------------------
// Run: implements all steps of algorithm
// -------------------------------------------------------------------------

void SEGAbstract::Run()
{
    Initialize();
    FindSegments();
}

// -------------------------------------------------------------------------
// Initialize: initializes and prepares the run for search of low entropy
//     regions
// -------------------------------------------------------------------------

void SEGAbstract::Initialize()
{
    ComputeEntropy();
}

// -------------------------------------------------------------------------
// Entropy: compute Entropy of segment of the entire run length
// -------------------------------------------------------------------------

double SEGAbstract::Entropy() const
{
    Window  window(
            GetVerifier(),
            GetComparator(),
            GetEqComparer(),
            GetRunAddress(),
            GetRunLength(),// !!- GetWinLength(),
            GetAlphabetSize()
    );
    window.SafeShift();
    return window.GetEntropy();
}

// -------------------------------------------------------------------------
// LogProbability: compute log-probability of segment of the entire run
//     length
// -------------------------------------------------------------------------

double SEGAbstract::LogProbability() const
{
    Window  window(
            GetVerifier(),
            GetComparator(),
            GetEqComparer(),
            GetRunAddress(),
            GetRunLength(),// !!- GetWinLength(),
            GetAlphabetSize()
    );
    window.SafeShift();
    return window.LogProbability();
}

// -------------------------------------------------------------------------
// FindSegments: the main method to be called to find segments of low
//     entropy
// -------------------------------------------------------------------------

void SEGAbstract::FindSegments()
{
    if( GetRunLength() == 0 )
        return;

    RecurrentSegFind( 0, GetRunLength() - 1 );
}

// -------------------------------------------------------------------------
// RecurrentSegFind: finds segments of low entropy within the run bounded by
//     left and right
// -------------------------------------------------------------------------

void SEGAbstract::RecurrentSegFind( size_t left, size_t right )
{
    size_t  n;
    size_t  locleft, locright;
    size_t  optleft, optright;
    size_t  centre =   GetWinLength()>> 1;
    size_t  adjust = ( GetWinLength() & 1 ) ^ 1;    //for even win lengths, adjust
//     size_t  negadj = ( GetWinLength() & 1 );        //negative adjust
    double  entro;

    if( GetRunLength() <= left || GetRunLength() <= right )
        throw myruntime_error( mystring( "SEGAbstract: Wrong extent boundaries." ));

    if( right < left )
        return;

    for( n = left; n <= right; n++ )
    {
        entro = GetEntropyAt( n );
        if( entro < 0.0 )
            continue;

        if( GetInverted()?
            entro < GetHighEntropy() :
            GetLowEntropy() < entro  )
            continue;

        locleft = locright = n;
        FindExtentAt( n, left, right, &locleft, &locright );

        if( locleft + adjust < centre || GetRunLength() <= locright + centre )
            throw myruntime_error( mystring( "SEGAbstract: Wrong extent found." ));

        //draw out to window boundaries
        optleft = locleft = locleft - centre + adjust;
        optright = locright += centre;

        OptimizeExtent( &optleft, &optright );

        //if found optimal extent is beyond the current position,
        //recurse to the left side to find left-hand side segments
        if( n + centre <= optleft && 1 + centre <= optleft ) {
            RecurrentSegFind( locleft + centre - adjust, optleft - 1 - centre );
        }

        //push segment...
        segments.Push( optleft, optright );
        //move to the next's centre from the leftmost window found
        n = PCMIN( locright - centre, optright + centre - adjust );
        left = n + 1;
    }
}

// -------------------------------------------------------------------------
// ComputeEntropy: computes entropy of each window within the run
// -------------------------------------------------------------------------

void SEGAbstract::ComputeEntropy()
{
    size_t  n, p, r, l;
    size_t  centre = GetWinLength() >> 1;
    size_t  adjust = ( GetWinLength() & 1 ) ^ 1;    //for even win lengths, adjust
//     size_t  negadj = ( GetWinLength() & 1 );        //negative adjust
    size_t  invalid = ( size_t ) -1;                //last invalid position

    Window  window(
            GetVerifier(),
            GetComparator(),
            GetEqComparer(),
            GetRunAddress(),
            GetWinLength(),
            GetAlphabetSize()
    );

    //reset entropy to the position of window centre
    for( n = 0; n < GetRunLength() && n < centre - adjust; n++ )
    {
        SetEntropyAt( n, -1.0 );
        //if position is not valid to consider
        if( ! ( *GetValidator())( GetRunAt( n )))
            invalid = n;
    }

    //find the last invalid position within the window
    for( p = n; p < GetRunLength() && p <= n + centre; p++ )
        if( ! ( *GetValidator())( GetRunAt( p )))
            invalid = p;


    for( ;  n < GetRunLength() &&
            n + centre < GetRunLength();
            n++ )
    {
        window.SafeShift();
        l = n - centre + adjust;
        r = n + centre;

        if( ! ( *GetValidator())( GetRunAt( r )))
            invalid = r;

        //mask windows with negative entropy that contain invalid positions
        if( l <= invalid && invalid <= r )
            SetEntropyAt( n, -1.0 );
        else
            SetEntropyAt( n, window.GetEntropy());
    }


    //reset entropy for the half of positions of the last window
    for( ; n < GetRunLength(); n++ )
        SetEntropyAt( n, -1.0 );
}

// -------------------------------------------------------------------------
// FindExtentAt: finds the largest low-entropy extent starting from the
//     position given
// -------------------------------------------------------------------------

void SEGAbstract::FindExtentAt( size_t n, size_t left, size_t right, size_t* locleft, size_t* locright )
{
    if( GetRunLength() <= left || GetRunLength() <= right || right < left )
        throw myruntime_error( mystring( "SEGAbstract: Wrong extent boundaries." ));

    if( locleft == NULL || locright == NULL )
        return;

    double  entro;
    size_t  p;

    for( p = n; p >= left; p-- ) {
        entro = GetEntropyAt( p );
        if( entro < 0.0 )
            break;

        if( GetInverted()?
            entro < GetLowEntropy() :
            GetHighEntropy() < entro )
            break;

        *locleft = p;

        if( p == 0 )
            break;
    }

    for( p = n; p <= right; p++ ) {
        entro = GetEntropyAt( p );
        if( entro < 0.0 )
            break;
        if( GetInverted()?
            entro < GetLowEntropy() :
            GetHighEntropy() < entro )
            break;
        *locright = p;
    }
}

// -------------------------------------------------------------------------
// OptimizeExtents: finds optimal extent within the run segment bounded by
//     left and right
// -------------------------------------------------------------------------

void SEGAbstract::OptimizeExtent( size_t* locleft, size_t* locright )
{
    if( locleft == NULL || locright == NULL )
        return;

    if( GetRunLength() <= *locleft || GetRunLength() <= *locright || *locright < *locleft )
        throw myruntime_error( mystring( "SEGAbstract: Wrong extent boundaries." ));

    double  logprob = 0.0;
    double  minlogprob = 0.0;
    size_t  extlength = *locright - *locleft + 1;
    size_t  minlength = ( GetMaxExtentDifference() < extlength  )? extlength - GetMaxExtentDifference(): 1;
    size_t  length;
    size_t  optleft = *locleft;
    size_t  optright = *locright;

    if( GetInverted())
        minlogprob = 100.0;

    for( length = extlength; length > minlength; length-- )
    {
        Window  window(
                GetVerifier(),
                GetComparator(),
                GetEqComparer(),
                GetRunAt( *locleft ),
                length,
                GetAlphabetSize()
        );

        for( size_t l = *locleft; l + length - 1 <= *locright; l++ )
        {
            window.SafeShift();
            logprob = window.LogProbability();
            if( GetInverted())
                logprob = -logprob;
            if( logprob < minlogprob ) {
                minlogprob = logprob;
                optleft = l;
                optright = l + length - 1;
            }
        }
    }

    *locleft = optleft;
    *locright = optright;
}

// =========================================================================
// PrintSequence: prints formatted sequence 
// -------------------------------------------------------------------------

void SEGAbstract::PrintSequence( FILE* fp, TResAcceder getres, size_t width )
{
    if( fp == NULL || getres == NULL )
        return;

    size_t  n = 0;
    size_t  p = 0;
    char*   buffer = ( char* )malloc( sizeof( char ) * ((( width < GetRunLength())? width: GetRunLength()) + 1 ));

    if( buffer == NULL )
        throw myruntime_error( mystring( "SEGAbstract: Not enough memory." ));

    for( n = 0; n < GetRunLength(); n += p )
    {
        for( p = 0; p < width && n + p < GetRunLength(); p++ )
            buffer[p] = DehashCode(( *getres )( GetRunAt(  n + p )));

        buffer[p] = 0;
        fprintf( fp, "%s\n", buffer );
    }

    fprintf( fp, "\n" );
    free( buffer );
}

// -------------------------------------------------------------------------
// PrintSeggedSequence: prints formatted sequence with segments found by
//     running the algorithm and masked with Xs
// -------------------------------------------------------------------------

void SEGAbstract::PrintSeggedSequence( FILE* fp, TResAcceder getres, size_t width )
{
    if( fp == NULL || getres == NULL )
        return;

    char*   residues = ( char* )malloc( sizeof( char ) * ( GetRunLength() + 1 ));
    char*   buffer = ( char* )malloc( sizeof( char ) * ((( width < GetRunLength())? width: GetRunLength()) + 1 ));

    if( residues == NULL || buffer == NULL )
        throw myruntime_error( mystring( "SEGSequence: Not enough memory." ));

    size_t  n = 0;
    size_t  p = 0;
    size_t  left, right;
    const SegmentStructure&     loc_segments = GetSegments();


    for( n = 0; n < GetRunLength(); n ++ )
        residues[n] = DehashCode(( *getres )( GetRunAt( n )));

    for( n = 0; n < loc_segments.GetSize(); n++ ) {
        left = loc_segments.GetLeftAt( n );
        right = loc_segments.GetRightAt( n );

        for( p = left; p <= right; p++ )
            if(( *getres )( GetRunAt( n )) != GAP )
                residues[p] = 'x';//DehashCode( X );
    }


    for( n = 0; n < GetRunLength(); n += p )
    {
        memcpy( buffer, residues + n, p = ( n + width < GetRunLength())? width: GetRunLength() - n );
        buffer[p] = 0;
        fprintf( fp, "%s\n", buffer );
    }

    fprintf( fp, "\n" );
    free( residues );
    free( buffer );
}

// -------------------------------------------------------------------------
// MaskSequence: masks sequence with Xs with accordance of intervals
//     found by the SEG algorithm;
//     msym, symbol to replace others with in the segged intervals,
//     omitmask, positions not to account for
//     omit, symbols not to replace
// -------------------------------------------------------------------------

void SEGAbstract::MaskSequence(
        char* sequence,
        size_t length,
        char msym,
        const char* omitmask,
        char omit1,
        char omit2 ) const
{
    if( sequence == NULL || omitmask == NULL || length == 0 )
        return;

    size_t  nn = 0;         //sequence iterator
    size_t  n = 0;
    size_t  p = 0;
    size_t  left, right;
    const SegmentStructure& loc_segments = GetSegments();


    for( n = 0; n < loc_segments.GetSize(); n++ ) {
        left = loc_segments.GetLeftAt( n );
        right = loc_segments.GetRightAt( n );

        for( ; p < left && nn < length; nn++ ) {
            if( omitmask[nn] )
                continue;
            p++;
        }

        for( ; p <= right && nn < length; nn++ ) {
            if( omitmask[nn] )
                continue;

            if( sequence[nn] != omit1 && sequence[nn] != omit2 )
                sequence[nn] = msym;
            p++;
        }
    }
}

// =========================================================================
// Methods for testing
// Print: prints vector of computed entropies and segments found
//

void SEGAbstract::Print( FILE* fp )
{
    if( fp == NULL )
        return;

    fprintf( fp, "\nSEGAbstract\n Entropies:\n " );
    for( size_t n = 0; n < GetRunLength(); n++ )
        fprintf( fp, " %6.3f", GetEntropyAt( n ));

    fprintf( fp, "\n" );
    segments.Print( fp );
}

