/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __HashTable__
#define __HashTable__

#include "debug.h"
#include "rc.h"
#include "data.h"

#include "mystring.h"
#include "myexcept.h"

//Hash function:
// key is the frequency vector comprising of effective number of residues
typedef size_t ( *THashFunction )( const void* key );


// _________________________________________________________________________
// Class HashTable
//

class HashTable
{
    enum { numPrimes = 28 };

public:
    HashTable( THashFunction func, size_t size );
    ~HashTable();

    size_t  GetAddress( const void* key ) const;                    //get address to this hash table given key

    size_t      GetHashSize() const { return hashsize; }            //hash size
    const void* GetValueAt( size_t address ) const;                 //obtain value at the position
    void        SetValueAt( size_t address, const void* value );    //set value at the position

protected:
    explicit HashTable();
    void    Reserve( size_t );

private:
    THashFunction   function;   //hash function
    const void**    values;     //table of values (void pointers); each value corresponds to a hash key
    size_t          hashsize;   //size of hash table

    static const unsigned long prime_list[numPrimes];

};

// INLINES ...

// -------------------------------------------------------------------------
// GetAddress: Get address to the hash table given key

inline
size_t HashTable::GetAddress( const void* key ) const
{
#ifdef __DEBUG__
    if( !function || !key )
        throw myruntime_error( mystring( "HashTable: Unable to compute address." ));
#endif
    size_t  address = ( *function )( key );
#ifdef __DEBUG__
    if( hashsize <= address )
        throw myruntime_error( mystring( "HashTable: Wrong address computed." ));
#endif
    return address;
}

// -------------------------------------------------------------------------
// GetValueAt: Return value at the position

inline
const void* HashTable::GetValueAt( size_t address ) const
{
#ifdef __DEBUG__
    if( address < hashsize )
#endif
        return values[address];

    throw myruntime_error( mystring( "HashTable: Memory access error." ));
}

// -------------------------------------------------------------------------
// SetValueAt: Set value at the position

inline
void HashTable::SetValueAt( size_t address, const void* value )
{
#ifdef __DEBUG__
    if( hashsize <= address )
        throw myruntime_error( mystring( "HashTable: Memory access error." ));
#endif
    values[address] = value;
}

#endif//__HashTable__
