/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include "myexcept.h"


// -------------------------------------------------------------------------
// CLASS myexception
//
// what: error description
//
const char* myexception::what() const throw()
{
    return "[empty]";
}

// eclass: error class
//
int myexception::eclass() const throw()
{
    return NOCLASS;
}

// -------------------------------------------------------------------------
// CLASS myruntime_error
//
// Constructor
//
myruntime_error::myruntime_error( const mystring& arg, int ecldesc )
:   myexception(),
    _errdsc( arg ),
    _class( ecldesc )
{
}

// Destructor
//
myruntime_error::~myruntime_error() throw()
{
}

// -------------------------------------------------------------------------
// what: returns error description
//
const char* myruntime_error::what() const throw()
{
    return _errdsc.c_str();
}

// -------------------------------------------------------------------------
// eclass: returns class of the error
//
int myruntime_error::eclass() const throw()
{
    return _class;
}
