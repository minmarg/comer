/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __SEGAbstract__
#define __SEGAbstract__

#include <stdio.h>
#include "BinarySearchStructure.h"
#include "SegmentStructure.h"

//default number of segments
#define DEF_NO_SEGMENTS 20
#define DEFAULT_SEGSEQ_PRINT_WIDTH  60

//count type
typedef size_t  tcount;

//Validator function:
// functions of this type will return true if element of run is valid/
//  proper for analysis, false otherwise
typedef bool ( *TValidator )( const void* key );
typedef bool ( *TVerifier )( const void* key );

//residue access function; needed only for subclass of the problems
typedef size_t ( *TResAcceder )( const void* key );

//comparison function for equality only; returns true if two keys are
//equal, false otherwise
typedef bool   ( *TEqComparer )( const void* key1, const void* key2 );
// STATIC FUNCTIONS
//
// SimpleIntComparer: comparison function for state vector values
//

static int SimpleIntComparer( const void* key1, const void* key2 )
{
    return ( int )(( ssize_t )key1 - ( ssize_t )key2 );
}


// =========================================================================
// Class LinearSearchStructure
//

class LinearSearchStructure: public BinarySearchStructure
{
public:
    LinearSearchStructure( TComparator, TEqComparer, size_t size );
    virtual ~LinearSearchStructure();

    virtual bool    Find( const void* key, int* loc = NULL ) const;

private:
    TEqComparer     eqcomparer;     //comparison function for equality only
};

// =========================================================================
// Class BinaryModifier
//

class BinaryModifier: public BinarySearchStructure
{
public:
    BinaryModifier( TComparator, size_t size );
    virtual ~BinaryModifier();

    bool    IncValue( tcount value );
    bool    DecValue( tcount value );

private:
};

// =========================================================================
// Class Window
//

class Window
{
public:
    Window(
        TVerifier   verffunc,
        TComparator comparator,
        TEqComparer eqcomparer,
        const void* address,
        size_t      len,
        size_t      szalpha = 0
    );

    ~Window();

    size_t      GetLength() const           { return length; }
    size_t      GetAlphabetSize() const     { return szalphabet; }
    double      GetEntropy() const          { return entropy; }

    const TVerifier     GetVerifier() const     { return verifier; }
    const void*         GetWindowAt( size_t pos ) const;
    tcount              GetCompositionAt( size_t n ) const;
    tcount              GetStateAt( size_t n ) const;

    bool        Uninitialized() const       { return GetEntropy() < 0.0; }

    void        Initialize();
    void        SafeShift();

    double      LogProbability() const;

    void        Print( FILE* );    //for testing

protected:
    const void*             GetWindow()                 { return window; }
    SimpleVector*           GetComposition()            { return composition; }
    BinarySearchStructure*  GetCompAddresses()          { return compaddresses; }
    BinarySearchStructure*  GetState()                  { return state; }

    size_t                  GetCompositionSize() const;
    size_t                  GetCompAddressesSize() const;
    size_t                  GetStateSize() const;

    bool        RefreshComposition( const void* value, size_t* address = NULL );

    void        IncCompositionAt( size_t n );   //increment composition value by 1
    void        DecCompositionAt( size_t n );   //decrement composition value by 1
    void        SetCompositionAt( size_t n, tcount value );
    void        InsertCompositionAt( size_t n, tcount value );

    bool        FindCompAddress( const void* value, size_t* address );
    bool        PushCompAddress( const void* value, size_t* address );

    void        PushStateValue( tcount value );
    void        IncStateValue( tcount value );
    void        DecStateValue( tcount value );

    void        Shift();

    void        SetEntropy( double value )  { entropy = value; }

    void        StateOn();
    void        CompositionOn();
    void        ComputeEntropy();

    double      LogMultinomial() const;
    double      LogNoCompositions( size_t = 0 ) const;

private:
    const TVerifier         verifier;           //verifier of window positions
    const void*             window;             //beginning of actual window
    const size_t            length;             //length of window
    const size_t            szalphabet;         //size of alphabet which can be given or not
    double                  entropy;            //entropy of window
    SimpleVector*           composition;        //composition vector of window
    BinarySearchStructure*  compaddresses;      //addresses of composition vector
    BinaryModifier*         state;              //state vector (sorted composition) of window
};

// _________________________________________________________________________
// Class SEGAbstract
//

class SEGAbstract
{
public:
    SEGAbstract(
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
        size_t  szalpha = 0
    );

    SEGAbstract(
        TValidator  valdfunc,
        TVerifier   verffunc,
        TComparator compfunc,
        TEqComparer eqcomfunc,
        size_t  winlen,
        double  lowent,
        double  highent,
        size_t  maxdiff,
        size_t  szalpha = 0
    );

    virtual ~SEGAbstract();

    void            SetRunAddress( void* address, size_t runlen );

    void            Run();

    void            Initialize();
    void            FindSegments();

    double          Entropy() const;            //compute entropy of segment of the entire run length
    double          LogProbability() const;     //compute log-probability of segment of the entire run length

    const TValidator    GetValidator() const    { return validator; }
    const TVerifier     GetVerifier() const     { return verifier; }
    const TComparator   GetComparator() const   { return comparator; }
    const TEqComparer   GetEqComparer() const   { return eqcomparer; }
    const void*         GetRunAddress() const   { return runaddress; }
    size_t              GetRunLength() const    { return runlength; }
    size_t              GetWinLength() const    { return winlength; }
    size_t              GetAlphabetSize() const { return szalphabet; }
    double              GetLowEntropy() const   { return lowentropy; }
    double              GetHighEntropy() const  { return highentropy; }
    size_t              GetMaxExtentDifference() const { return maxdiffextent; }

    void            Print( FILE* );    //for testing

    void            SetLowCSearch()             { inverted = false; }
    void            SetHighCSearch()            { inverted = true; }

    void            MaskSequence( char*, size_t length, char msym, const char* omitmask, char, char ) const;

protected:
    const void*     GetRunAt( size_t pos ) const;

    double          GetEntropyAt( size_t pos ) const;
    void            SetEntropyAt( size_t pos, double value );

    const SegmentStructure& GetSegments() const { return segments; }

    bool            GetInverted() const         { return inverted; }

    void            Allocate( size_t newlen );
    void            Destroy();

    void            ComputeEntropy();
    void            RecurrentSegFind( size_t left, size_t right );
    void            FindExtentAt( size_t pos, size_t left, size_t right, size_t* locleft, size_t* locright );
    void            OptimizeExtent( size_t* locleft, size_t* locright );

    //-- semi-abstract methods
    void            PrintSequence( FILE*, TResAcceder, size_t width = DEFAULT_SEGSEQ_PRINT_WIDTH );
    void            PrintSeggedSequence( FILE*, TResAcceder, size_t width = DEFAULT_SEGSEQ_PRINT_WIDTH );

private:
    const TValidator    validator;  //validator function of run elements
    const TVerifier     verifier;   //verifier function for run elements
    const TComparator   comparator; //comparer of run elements
    const TEqComparer   eqcomparer; //comparison function for equality only
    const void*         runaddress; //address of run
    const size_t        runlength;  //length of run
    const size_t        winlength;  //length of window to use
    const size_t        szalphabet; //size of alphabet a run comprises
    const double        lowentropy; //low entropy threshold
    const double        highentropy;//high entropy threshold
    const size_t        maxdiffextent;  //maximum difference between extent (segment) length and minimum extent length

    bool                inverted;   //if true, high complexity segment search is in effect

    double*             entropies;  //entropies computed and attributed to window centres
    SegmentStructure    segments;   //storage for found segments

};

////////////////////////////////////////////////////////////////////////////
// Class Window INLINES
//
// GetWindowAt: returns address of window element at the given position

inline
const void* Window::GetWindowAt( size_t pos ) const
{
#ifdef __DEBUG__
    if( !window || length <= pos )
        throw myruntime_error( mystring( "Window: Memory access error." ));
#endif
    return ( const char* )window + pos * sizeof( void* );
}


// -------------------------------------------------------------------------
// GetCompositionSize: size of vector composition

inline
size_t Window::GetCompositionSize() const
{
#ifdef __DEBUG__
    if( !composition )
        throw myruntime_error( mystring( "Window: Memory access error." ));
#endif
    return composition->GetSize();
}

// -------------------------------------------------------------------------
// GetCompositionSize: size of vector of composition addresses compaddresses

inline
size_t Window::GetCompAddressesSize() const
{
#ifdef __DEBUG__
    if( !compaddresses )
        throw myruntime_error( mystring( "Window: Memory access error." ));
#endif
    return compaddresses->GetSize();
}

// -------------------------------------------------------------------------
// GetStateSize: size of state vector state

inline
size_t Window::GetStateSize() const
{
#ifdef __DEBUG__
    if( !state )
        throw myruntime_error( mystring( "Window: Memory access error." ));
#endif
    return state->GetSize();
}


// -------------------------------------------------------------------------
// GetCompositionAt: gets composition value of element of alphabet

inline
tcount Window::GetCompositionAt( size_t n ) const
{
#ifdef __DEBUG__
    //number of composition values can be greater than length of window
    if( !composition )
        throw myruntime_error( mystring( "Window: Memory access error." ));
#endif
    return ( tcount )composition->GetValueAt( n );
}

// -------------------------------------------------------------------------
// SetCompositionAt: sets composition value for element of alphabet

inline
void Window::InsertCompositionAt( size_t n, tcount value )
{
#ifdef __DEBUG__
    if( !composition )
        throw myruntime_error( mystring( "Window: Memory access error." ));
#endif
    composition->InsertValueAt( n, ( void* )( value ));
}

// -------------------------------------------------------------------------
// SetCompositionAt: sets composition value for element of alphabet

inline
void Window::SetCompositionAt( size_t n, tcount value )
{
#ifdef __DEBUG__
    if( !composition )
        throw myruntime_error( mystring( "Window: Memory access error." ));
#endif
    composition->SetValueAt( n, ( void* )( value ));
}

// -------------------------------------------------------------------------
// IncCompositionAt: increments composition value by 1

inline
void Window::IncCompositionAt( size_t n )
{
    SetCompositionAt( n, GetCompositionAt( n ) + 1 );
}

// -------------------------------------------------------------------------
// DecCompositionAt: decrements composition value by 1

inline
void Window::DecCompositionAt( size_t n )
{
    SetCompositionAt( n, GetCompositionAt( n ) - 1 );
}


// -------------------------------------------------------------------------
// FindCompAddress: finds individual value of alphabet; a position
//     corresponds to the address; returns true if found, false
//     otherwise

inline
bool Window::FindCompAddress( const void* value, size_t* address )
{
#ifdef __DEBUG__
    if( !compaddresses )
        throw myruntime_error( mystring( "Window: Memory access error." ));
#endif
    int     location = -1;
    bool    found;

    found = compaddresses->Find( value, &location );

    if( found && location < 0 )
        throw myruntime_error( mystring( "Window: Failed to determine address of composition value." ));

    if( found && address )
        *address = ( size_t ) location;

    return found;
}

// -------------------------------------------------------------------------
// PushCompAddress: inserts individual value of alphabet; a position will
//     determine the address; returns true if insert was successful, false
//     otherwise but address will be refreshed anyway

inline
bool Window::PushCompAddress( const void* value, size_t* address )
{
#ifdef __DEBUG__
    if( !compaddresses )
        throw myruntime_error( mystring( "Window: Memory access error." ));
#endif
    int     location = -1;
    bool    inserted;

    inserted = compaddresses->Push( value, &location );

    if( location < 0 )
        throw myruntime_error( mystring( "Window: Unable to determine address of composition value." ));

    if( address )
        *address = ( size_t ) location;

    return inserted;
}

// -------------------------------------------------------------------------
// RefreshComposition: refreshes compositional information as this new
//     value appears; returns flag of whether value is new to the
//     composition

inline
bool Window::RefreshComposition( const void* value, size_t* saveaddr )
{
    if( ! ( *GetVerifier())( value ))
        return false;

    size_t  address = 0;
    bool    appended = PushCompAddress( value, &address );
    //compositions will be stored continuously as this is provided by
    //addresses' BinarySearchStructure
    if( appended )
        InsertCompositionAt( address, 0 );

    IncCompositionAt( address );

    if( saveaddr )
        *saveaddr = address;

    return appended;
}


// -------------------------------------------------------------------------
// GetStateAt: gets state value of element of alphabet

inline
tcount Window::GetStateAt( size_t n ) const
{
#ifdef __DEBUG__
    //number of state values can be greater than length of window
    if( !state )
        throw myruntime_error( mystring( "Window: Memory access error." ));
#endif
    return ( tcount )state->GetValueAt( n );
}

// -------------------------------------------------------------------------
// PushStateValue: inserts state value 

inline
void Window::PushStateValue( tcount value )
{
#ifdef __DEBUG__
    if( !state )
        throw myruntime_error( mystring( "Window: Memory access error." ));
#endif
    state->Push(( void* ) value );
}

// -------------------------------------------------------------------------
// IncStateValue: increases state value by 1

inline
void Window::IncStateValue( tcount value )
{
#ifdef __DEBUG__
    if( !state )
        throw myruntime_error( mystring( "Window: Memory access error." ));
#endif

    if( !state->IncValue( value ))
        throw myruntime_error( mystring( "Window: State value not found." ));
}

// -------------------------------------------------------------------------
// DecStateValue: decreases state value by 1

inline
void Window::DecStateValue( tcount value )
{
#ifdef __DEBUG__
    if( !state )
        throw myruntime_error( mystring( "Window: Memory access error." ));
#endif

    if( !state->DecValue( value ))
        throw myruntime_error( mystring( "Window: State value not found." ));
}


////////////////////////////////////////////////////////////////////////////
// Class SEGAbstract INLINES
//
// GetRunAt: returns address of element of the run at the given position

inline
const void* SEGAbstract::GetRunAt( size_t pos ) const
{
#ifdef __DEBUG__
    if( !runaddress || runlength <= pos )
        throw myruntime_error( mystring( "SEGAbstract: Memory access error." ));
#endif
    return ( const char* )runaddress + pos * sizeof( void* );
}

// -------------------------------------------------------------------------
// GetEntropyAt: returns entropy at the specified run position

inline
double SEGAbstract::GetEntropyAt( size_t pos ) const
{
#ifdef __DEBUG__
    if( !entropies || runlength <= pos )
        throw myruntime_error( mystring( "SEGAbstract: Memory access error." ));
#endif
    return entropies[pos];
}

// -------------------------------------------------------------------------
// SetEntropyAt: sets entropy value at the run position

inline
void SEGAbstract::SetEntropyAt( size_t pos, double value )
{
#ifdef __DEBUG__
    if( !entropies || runlength <= pos )
        throw myruntime_error( mystring( "SEGAbstract: Memory access error." ));
#endif
    entropies[pos] = value;
}



#endif//__SEGAbstract__
