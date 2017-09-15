/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __table__
#define __table__

#include <stdio.h>
#include <stdlib.h>
#include "rc.h"
#include "myexcept.h"
#include "ext/pslvector.h"
#include "basin.h"
#include "menu.h"

// -------------------------------------------------------------------------
// class Table: reusable vector of values (objects)
//
class Table
{
public:
    Table( int size );
    Table( const Table& );
    explicit Table();
    ~Table();

    Table&  operator=( const Table& );

    int         GetSize() const { return length_; }
    int         GetActualSize() const { return actlen_; }

    int         GetVectorNIndAt( int loc ) const;
    Pslvector*  GetVectorNAt( int loc ) const;
    void        SetVectorNIndAt( int loc, int n );
    int         GetVectorNDishIndAt( int loc ) const;
    void        SetVectorNDishIndAt( int loc, int d );

    int         NewVectorNInd( int n, int d );//push vector with indices in basin and in dish
    void        RemValueAt( int loc, const Pslvector* value );//remove value


    bool        GetProcessedAt( int n ) const;
    void        SetProcessedAt( int n, bool value );

    bool        GetProcessed() const { return tblprocessed_; }
    void        SetProcessed( bool value ) { tblprocessed_ = value; }

    void        Copy( const Table& vector );     //copy elements
    void        Clear();                        //clear all elements
    void        Print( FILE* ) const;           //print vector

    bool        GetDestroy() const { return destroy_; }
    void        SetDestroy( bool value ) { destroy_ = value; }

    const Basin* GetBasin() const { return basin_; }
    void         SetBasin( const Basin* value ) { basin_ = value; }

    Menu*       GetMenu() const { return menu_; }
    void        SetMenu( Menu* value ) { menu_ = value; }

    Dish*       GetDish() const;
    int         GetDishIndex() const { return dishndx_; }
    void        SetDishIndex( int value ) { dishndx_ = value; }

    void        Reserve( int size );
    void        ReserveValues( int size );
    void        ReserveVacans( int size );

protected:
    void        Destroy();

    void        Realloc( int newcap );
    void        DestroyValues();
    int         GetCapacity() const { return capacity_; }

    int*        GetVector() const { return values_; }

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
    const Basin* basin_;        //global basin of vectors
    Menu*       menu_;          //menu of all dishes
    int         dishndx_;       //index of parent structure: dish
    int*        values_;        //indices of vectors attributable to the table and dish
    int*        dshndxvals_;    //dish indices of vectors
    bool*       processed_;     //flags of processed objects
    bool        tblprocessed_;  //flag of this table being processed
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
void Table::Reserve( int size )
{
    ReserveValues( size );
    ReserveVacans( size );
}

// -------------------------------------------------------------------------
// ReserveValues: Reserve space for vector of values
//
inline
void Table::ReserveValues( int size )
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
void Table::ReserveVacans( int size )
{
    if( 0 < size ) {
        ReallocVacans( size );
    }
}

// -------------------------------------------------------------------------
// Destroy: Destroy vectors
//
inline
void Table::Destroy()
{
    DestroyValues();
    DestroyVacans();
}

// -------------------------------------------------------------------------
// DestroyValues: Destroy vector of values
//
inline
void Table::DestroyValues()
{
    int n;
    if( values_ )
        free( values_ );
    if( dshndxvals_ )
        free( dshndxvals_ );
    if( processed_ )
        free( processed_ );
    values_ = NULL;
    dshndxvals_ = NULL;
    processed_ = NULL;
    basin_ = NULL;
    menu_ = NULL;
    actlen_ = 0;
    length_ = 0;
    capacity_ = 0;
}

// -------------------------------------------------------------------------
// DestroyVacans: Destroy vector of vacancies
//
inline
void Table::DestroyVacans()
{
    if( vacancies_ )
        free( vacancies_ );
    vacancies_ = NULL;
    novacs_ = 0;
    capvac_ = 0;
}


// -------------------------------------------------------------------------
// GetDish: get dish object
//
inline
Dish* Table::GetDish() const
{
    if( !menu_ )
        throw myruntime_error( "Table: GetDish: Memory access error" );
    if( dishndx_ < 0 || menu_->GetSize() <= dishndx_ )
        throw myruntime_error( "Table: GetDish: Memory access error" );
    return menu_->GetDishAt( dishndx_ );
}

// -------------------------------------------------------------------------
// GetVectorNIndAt: get index of vector from the basin, stored at `loc'
//
inline
int Table::GetVectorNIndAt( int loc ) const
{
    if( !values_ || loc < 0 || length_ <= loc )
        throw myruntime_error( "Table: GetVectorNIndAt: Memory access error" );
    return values_[loc];
}

// -------------------------------------------------------------------------
// GetVectorNAt: get vector from the basin with index stored at `loc'
//
inline
Pslvector* Table::GetVectorNAt( int loc ) const
{
    if( !basin_ || !values_ || loc < 0 || length_ <= loc )
        throw myruntime_error( "Table: GetVectorNAt: Memory access error" );
    int ndx = values_[loc];
    if( ndx < 0 )
        return NULL;
    if( basin_->GetSize() <= ndx )
        throw myruntime_error( "Table: GetVectorNAt: Memory access error" );
    return basin_->GetValueAt( ndx );
}

// -------------------------------------------------------------------------
// SetVectorNIndAt: set vector's basin index at `loc'
//
inline
void Table::SetVectorNIndAt( int loc, int n )
{
#ifdef __DEBUG__
    if( !values_ || loc < 0 || length_ <= loc )
        throw myruntime_error( "Table: SetVectorNIndAt: Memory access error" );
#endif
    values_[loc] = n;
}

// -------------------------------------------------------------------------
// GetVectorNDishIndAt: get dish index of vector from the basin
//
inline
int Table::GetVectorNDishIndAt( int loc ) const
{
    if( !dshndxvals_ || loc < 0 || length_ <= loc )
        throw myruntime_error( "Table: GetVectorNDishIndAt: Memory access error" );
    return dshndxvals_[loc];
}

// -------------------------------------------------------------------------
// SetVectorNDishIndAt: set vector's dish index at `loc'
//
inline
void Table::SetVectorNDishIndAt( int loc, int d )
{
#ifdef __DEBUG__
    if( !dshndxvals_ || loc < 0 || length_ <= loc )
        throw myruntime_error( "Table: SetVectorNDishIndAt: Memory access error" );
#endif
    dshndxvals_[loc] = d;
}


// -------------------------------------------------------------------------
// NewVectorNInd: push vector's basin and dish indices
//
inline
int Table::NewVectorNInd( int n, int d )
{
    if( !values_ || !dshndxvals_ )
        throw myruntime_error( "Table: NewVectorNInd: Memory access error" );
    if( n < 0 || d < 0 )
        throw myruntime_error( "Table: NewVectorNInd: Invalid indices." );
    int loc;
    if( 0 < novacs_ ) {
        loc = vacancies_[novacs_-1];
        if( loc < 0 || length_ <= loc )
            throw myruntime_error( "Table: NewVectorNInd: Memory access error" );
        values_[loc] = n;
        dshndxvals_[loc] = d;
        novacs_--;
    }
    else {
        if( capacity_ <= ( loc = length_ ))
            Realloc( TIMES2( capacity_ + 1 ));
        values_[loc] = n;
        dshndxvals_[loc] = d;
        length_++;
    }
    actlen_++;
    return loc;
}

// -------------------------------------------------------------------------
// RemValueAt: remove value at the position
//
inline
void Table::RemValueAt( int loc, const Pslvector* value )
{
#ifdef __DEBUG__
    if( !values_ || loc < 0 || length_ <= loc )
        throw myruntime_error( "Table: RemValueAt: Memory access error" );
#endif
    if( GetVectorNAt( loc ) != value )
        throw myruntime_error( "Table: RemValueAt: Memory access error" );
    values_[loc] = -1;
    dshndxvals_[loc] = -1;
    processed_[loc] = false;
    PushVacant( loc );
    actlen_--;
}


// -------------------------------------------------------------------------
// GetProcessedAt: get flag of processed object at the position
//
inline
bool Table::GetProcessedAt( int loc ) const
{
#ifdef __DEBUG__
    if( !processed_ || loc < 0 || length_ <= loc )
        throw myruntime_error( "Table: GetProcessedAt: Memory access error" );
#endif
    return processed_[loc];
}

// -------------------------------------------------------------------------
// SetProcessedAt: set flag of processed object at the position
//
inline
void Table::SetProcessedAt( int loc, bool value )
{
#ifdef __DEBUG__
    if( !processed_ || loc < 0 || length_ <= loc )
        throw myruntime_error( "Table: SetProcessedAt: Memory access error" );
#endif
    processed_[loc] = value;
}



// -------------------------------------------------------------------------
// GetVacantAt: get index from the vector of vacancies
//
inline
int Table::GetVacantAt( int loc ) const
{
#ifdef __DEBUG__
    if( !vacancies_ || loc < 0 || novacs_ <= loc )
        throw myruntime_error( "Table: GetVacantAt: Memory access error" );
#endif
    return vacancies_[loc];
}

// -------------------------------------------------------------------------
// SetVacantAt: write index in the vector of vacancies
//
inline
void Table::SetVacantAt( int loc, int index )
{
#ifdef __DEBUG__
    if( !vacancies_ || loc < 0 || novacs_ <= loc )
        throw myruntime_error( "Table: SetVacantAt: Memory access error" );
#endif
    vacancies_[loc] = index;
}

// -------------------------------------------------------------------------
// PushVacant: push index in the vector of vacancies
//
inline
void Table::PushVacant( int index )
{
    if( capvac_ <= novacs_ )
        ReallocVacans( TIMES2( capvac_ + 1 ));
    vacancies_[novacs_++] = index;
}

#endif//__table__
