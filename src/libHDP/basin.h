/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __basin__
#define __basin__

#include <stdio.h>
#include <stdlib.h>
#include "rc.h"
#include "myexcept.h"
#include "ext/pslvector.h"

// -------------------------------------------------------------------------
// class Basin: reusable vector of values (objects)
//
class Basin
{
public:
    Basin( int size );
    Basin( const Basin& );
    explicit Basin();
    ~Basin();

    Basin&  operator=( const Basin& );

    int         GetSize() const { return length_; }
    int         GetActualSize() const { return actlen_; }

    Pslvector*  GetValueAt( int n ) const;
    void        SetValueAt( int n, Pslvector* value );
    int         NewValue( Pslvector* value );//push value
    void        RemValueAt( int loc, Pslvector* value );//remove value

    bool        GetProcessedAt( int n ) const;
    void        SetProcessedAt( int n, bool value );

    void        Copy( const Basin& vector );     //copy elements
    void        Clear();                        //clear all elements
    void        Print( FILE* ) const;           //print vector

    bool        GetDestroy() const { return destroy_; }
    void        SetDestroy( bool value ) { destroy_ = value; }

    void        Reserve( int size );
    void        ReserveValues( int size );
    void        ReserveVacans( int size );

protected:
    void        Destroy();

    void        Realloc( int newcap );
    void        DestroyValues();
    int         GetCapacity() const { return capacity_; }

    Pslvector** GetVector() const { return values_; }

    void        ReallocVacans( int newcap );
    void        DestroyVacans();
    int         GetCapVacans() const { return capvac_; }

    const int*  GetVacancies() const { return vacancies_; }
    int         GetVacantAt( int n ) const;
    void        SetVacantAt( int n, int index );
    void        PushVacant( int index );

    int         GetNoVacans() const { return novacs_; }
    void        SetNoVacans( int value ) { novacs_ = value; }

    void        SetSize( int value ) { length_ = value; }
    void        SetActualSize( int value ) { actlen_ = value; }

private:
    Pslvector** values_;        //values attributable to the dish
    bool*       processed_;     //flags of processed objects
    int         actlen_;        //actual number of values
    int         length_;        //length of vector
    int         capacity_;      //current capacity of the vector
    bool        destroy_;       //flag of object destroy
    int*        vacancies_;     //vector of indices of vacancies
    int         novacs_;        //number of vacant elements
    int         capvac_;        //capcity of vacancies_
};


// -------------------------------------------------------------------------
// Reserve: Reserve space for vectors
//
inline
void Basin::Reserve( int size )
{
    ReserveValues( size );
    ReserveVacans( size );
}

// -------------------------------------------------------------------------
// ReserveValues: Reserve space for vector of values
//
inline
void Basin::ReserveValues( int size )
{
    if( 0 < size ) {
        Realloc( size );
//         SetSize( size );
    }
}

// -------------------------------------------------------------------------
// ReserveVacans: Reserve space for vector of vacancies
//
inline
void Basin::ReserveVacans( int size )
{
    if( 0 < size ) {
        ReallocVacans( size );
    }
}

// -------------------------------------------------------------------------
// Destroy: Destroy vectors
//
inline
void Basin::Destroy()
{
    DestroyValues();
    DestroyVacans();
}

// -------------------------------------------------------------------------
// DestroyValues: Destroy vector of values
//
inline
void Basin::DestroyValues()
{
    int n;
    if( values_ ) {
        if( GetDestroy())
            for( n = 0; n < length_; n++ )
                if( values_[n] )
                    delete values_[n];
        free( values_ );
    }
    if( processed_ )
        free( processed_ );
    values_ = NULL;
    processed_ = NULL;
    actlen_ = 0;
    length_ = 0;
    capacity_ = 0;
}

// -------------------------------------------------------------------------
// DestroyVacans: Destroy vector of vacancies
//
inline
void Basin::DestroyVacans()
{
    if( vacancies_ )
        free( vacancies_ );
    vacancies_ = NULL;
    novacs_ = 0;
    capvac_ = 0;
}

// -------------------------------------------------------------------------
// GetValueAt: get value at the position
//
inline
Pslvector* Basin::GetValueAt( int loc ) const
{
#ifdef __DEBUG__
    if( !values_ || loc < 0 || length_ <= loc )
        throw myruntime_error( "Memory access error" );
#endif
    return values_[loc];
}

// -------------------------------------------------------------------------
// SetValueAt: set value at the position
//
inline
void Basin::SetValueAt( int loc, Pslvector* value )
{
#ifdef __DEBUG__
    if( !values_ || loc < 0 || length_ <= loc )
        throw myruntime_error( "Basin: Memory access error" );
#endif
    values_[loc] = value;
}

// -------------------------------------------------------------------------
// NewValue: push value
//
inline
int Basin::NewValue( Pslvector* value )
{
    int loc;
    if( 0 < novacs_ ) {
        loc = vacancies_[novacs_-1];
        if( loc < 0 || length_ <= loc )
            throw myruntime_error( "Basin: NewValue: Memory access error" );
        values_[loc] = value;
        novacs_--;
        actlen_++;
        return loc;
    }
    if( capacity_ <= ( loc = length_ ))
        Realloc( TIMES2( capacity_ + 1 ));
    values_[loc] = value;
    length_++;
    actlen_++;
    return loc;
}

// -------------------------------------------------------------------------
// RemValueAt: remove value at the position
//
inline
void Basin::RemValueAt( int loc, Pslvector* value )
{
#ifdef __DEBUG__
    if( !values_ || loc < 0 || length_ <= loc )
        throw myruntime_error( "Basin: Memory access error" );
#endif
    if( values_[loc] != value )
        throw myruntime_error( "Basin: RemValueAt: Memory access error" );
    values_[loc] = NULL;
    processed_[loc] = false;
    PushVacant( loc );
    actlen_--;
}

// -------------------------------------------------------------------------
// GetProcessedAt: get flag of processed object at the position
//
inline
bool Basin::GetProcessedAt( int loc ) const
{
#ifdef __DEBUG__
    if( !processed_ || loc < 0 || length_ <= loc )
        throw myruntime_error( "Table: Memory access error" );
#endif
    return processed_[loc];
}

// -------------------------------------------------------------------------
// SetProcessedAt: set flag of processed object at the position
//
inline
void Basin::SetProcessedAt( int loc, bool value )
{
#ifdef __DEBUG__
    if( !processed_ || loc < 0 || length_ <= loc )
        throw myruntime_error( "Table: Memory access error" );
#endif
    processed_[loc] = value;
}



// -------------------------------------------------------------------------
// GetVacantAt: get index from the vector of vacancies
//
inline
int Basin::GetVacantAt( int loc ) const
{
#ifdef __DEBUG__
    if( !vacancies_ || loc < 0 || novacs_ <= loc )
        throw myruntime_error( "Memory access error" );
#endif
    return vacancies_[loc];
}

// -------------------------------------------------------------------------
// SetVacantAt: write index in the vector of vacancies
//
inline
void Basin::SetVacantAt( int loc, int index )
{
#ifdef __DEBUG__
    if( !vacancies_ || loc < 0 || novacs_ <= loc )
        throw myruntime_error( "Basin: Memory access error" );
#endif
    vacancies_[loc] = index;
}

// -------------------------------------------------------------------------
// PushVacant: push index in the vector of vacancies
//
inline
void Basin::PushVacant( int index )
{
    if( capvac_ <= novacs_ )
        ReallocVacans( TIMES2( capvac_ + 1 ));
    vacancies_[novacs_++] = index;
}

#endif//__basin__
