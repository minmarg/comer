/***************************************************************************
 *   Copyright (C) 2010 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __pslvector__
#define __pslvector__

#include <stdio.h>
#include <stdlib.h>
#include "pslcodes.h"

class Pslmatrix;

class Pslvector
{
public:
    Pslvector( int size );
    Pslvector( const Pslvector& );
    explicit Pslvector();
    ~Pslvector();

    Pslvector&  operator=( const Pslvector& );

    int         GetSize() const { return length_; }
    int         GetCapacity() const { return capacity_; }

    int         AddGNoise( const Pslvector& stds );
    void        SetAllToValue( double value );
    void        SetUnity() { SetAllToValue( 1.0 ); }

    const Pslvector SubVector( int offset, int n ) const;

    double      GetValueAt( int n ) const;          //get value at the position
    void        SetValueAt( int n, double value );  //set value at the position
    void        AddValueAt( int n, double value );  //add value at the position
    void        MulValueAt( int n, double value );  //multiply by value at the position
    int         DivValueAt( int n, double value );  //divide by value at the position

    void        InsertAt( int loc, double value );  //insert value at the position
    void        Push( double value );               //push value at the end

    void        Copy( const Pslvector& vector );    //copy elements
    void        Zero();                             //assign all elements to zero
    void        Clear();                            //clear all elements
    void        Print( FILE* ) const;               //print vector
    void        Print( FILE*, const char* format ) const;//print vector

    //LINEAR ALGEBRA
    double      Min() const;
    double      Max() const;
    double      Sum() const;                                            //sum of all members
    double      Norm2() const;                                          //norm
    int         DotProduct( const Pslvector& vect2, double* res ) const;//dot product of two vectors
    int         Superposition( double alpha, const Pslvector& vect2 );  //addition of two vectors
    int         MultiplyBy( double alpha );                             //multiplication by scalar
    int         Mul( const Pslmatrix& mt, const Pslvector& v ); //multiplication of matrix by vector
    int         Transpose( Pslmatrix& tr ) const; //transpose of vector
    int         Exp();//exponentiate vector

    void        Allocate( int cap );
    void        Reserve( int size );

    const double*   GetVector() const { return values_; }
    double*         GetVector() { return values_; }
    void            DecDim() { if( !GetMaster()) return; if( length_ ) length_--; }

protected:
    void        Realloc( int newcap );
    void        Destroy();

    void        SetVector( double* vect, int len );

    void        priverror( const char* ) const;


    bool        GetMaster() const { return master_; }
    void        SetMaster() { master_ = true; }
    void        SetSlave() { master_ = false; }

    int         GetStride() const { return stride_; }
    void        SetStride( int value ) { stride_ = value; }

    void        SetSize( int value ) { if( !GetMaster()) return; length_ = value; }

    friend class Pslmatrix;

private:
    double*     values_;        //double values of vector
    int         length_;        //length of vector
    int         capacity_;      //current capacity of the sequence
    int         stride_;        //stride of vector
    bool        master_;        //flag of mastering an object
};


// -------------------------------------------------------------------------
// Allocate: Allocate space for vector
//
inline
void Pslvector::Allocate( int size )
{
    if( !GetMaster())
        priverror( "Only master is allowed to allocate space." );
    if( 0 < size ) {
        SetStride( 1 );
        Realloc( size );
    }
}

// -------------------------------------------------------------------------
// Reserve: Reserve space for vector
//
inline
void Pslvector::Reserve( int size )
{
    if( !GetMaster())
        priverror( "Only master is allowed to reserve space." );
    if( 0 < size ) {
        SetStride( 1 );
        Realloc( size );
        SetSize( size );
    }
}

// -------------------------------------------------------------------------
// Destroy: Destroy vector
//
inline
void Pslvector::Destroy()
{
    if( master_ )
        if( values_ )
            free( values_ );
    values_ = NULL;
    length_ = 0;
    capacity_ = 0;
}

// -------------------------------------------------------------------------
// SetVector: set vector and its size
//
inline
void Pslvector::SetVector( double* vect, int len )
{
    Destroy();
    SetSlave();
    length_ = len;
    values_ = vect;
}

// -------------------------------------------------------------------------
// GetValueAt: get value at the position
//
inline
double Pslvector::GetValueAt( int loc ) const
{
#ifdef __DEBUG__
    if( !values_ || loc < 0 || length_ <= loc || stride_ < 1 )
        priverror( "Memory access error" );
#endif
    if( stride_ == 1 )
        return values_[loc];
    return values_[loc*stride_];
}

// -------------------------------------------------------------------------
// SetValueAt: set value at the position
//
inline
void Pslvector::SetValueAt( int loc, double value )
{
    if( GetMaster()) {
#ifdef __DEBUG__
        if( !values_ || loc < 0 || capacity_ <= loc || stride_ < 1 )
            priverror( "Memory access error" );
//         if( length_ + 1 <= loc )
//             priverror( "Memory access error" );
#endif
        if( length_ <= loc )
            length_ = loc + 1;
    }
    else {//slave
#ifdef __DEBUG__
        if( !values_ || loc < 0 || length_ <= loc || stride_ < 1 )
            priverror( "Memory access error" );
#endif
    }
    if( stride_ == 1 )
        values_[loc] = value;
    else
        values_[loc*stride_] = value;
}

// -------------------------------------------------------------------------
// AddValueAt: add value at the position
//
inline
void Pslvector::AddValueAt( int loc, double value )
{
    if( GetMaster()) {
#ifdef __DEBUG__
        if( !values_ || loc < 0 || capacity_ <= loc || stride_ < 1 )
            priverror( "Memory access error" );
//         if( length_ + 1 <= loc )
//             priverror( "Memory access error" );
#endif
        if( length_ <= loc )
            length_ = loc + 1;
    }
    else {//slave
#ifdef __DEBUG__
        if( !values_ || loc < 0 || length_ <= loc || stride_ < 1 )
            priverror( "Memory access error" );
#endif
    }
    if( stride_ == 1 )
        values_[loc] += value;
    else
        values_[loc*stride_] += value;
}

// -------------------------------------------------------------------------
// MulValueAt: multiply by value at the position
//
inline
void Pslvector::MulValueAt( int loc, double value )
{
    if( GetMaster()) {
#ifdef __DEBUG__
        if( !values_ || loc < 0 || capacity_ <= loc || stride_ < 1 )
            priverror( "Memory access error" );
//         if( length_ + 1 <= loc )
//             priverror( "Memory access error" );
#endif
        if( length_ <= loc )
            length_ = loc + 1;
    }
    else {//slave
#ifdef __DEBUG__
        if( !values_ || loc < 0 || length_ <= loc || stride_ < 1 )
            priverror( "Memory access error" );
#endif
    }
    if( stride_ == 1 )
        values_[loc] *= value;
    else
        values_[loc*stride_] *= value;
}

// -------------------------------------------------------------------------
// DivValueAt: divide by value at the position
//
inline
int Pslvector::DivValueAt( int loc, double value )
{
    if( GetMaster()) {
#ifdef __DEBUG__
        if( !values_ || loc < 0 || capacity_ <= loc || stride_ < 1 )
            priverror( "Memory access error" );
//         if( length_ + 1 <= loc )
//             priverror( "Memory access error" );
#endif
        if( length_ <= loc )
            length_ = loc + 1;
    }
    else {//slave
#ifdef __DEBUG__
        if( !values_ || loc < 0 || length_ <= loc || stride_ < 1 )
            priverror( "Memory access error" );
#endif
    }
    if( value == 0.0 )
        return PSL_ERR_ILLEGAL;

    if( stride_ == 1 )
        values_[loc] /= value;
    else
        values_[loc*stride_] /= value;

    return PSL_SUCCESS;
}

#endif//__pslvector__
