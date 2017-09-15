/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __myexcept_h__
#define __myexcept_h__

#include "debug.h"
#include "mystring.h"

enum {
    NOCLASS,
    SCALING,
    CRITICAL
};

// _________________________________________________________________________
// CLASS myexception
// for convenient exception handling
//
class myexception
{
public:
    myexception() throw() {}
    virtual ~myexception() throw() {};
    virtual const char* what() const throw();
    virtual int         eclass() const throw();
};

// _________________________________________________________________________
// CLASS myruntime_error
// class implementing runtime error exception
//
class myruntime_error: public myexception
{
public:
    explicit myruntime_error( const mystring& arg, int ecl = NOCLASS );
    virtual ~myruntime_error() throw();

    virtual const char* what() const throw();       //returns cause of the error
    virtual int         eclass() const throw();     //returns class of the error

private:
    mystring    _errdsc;    //error description
    int         _class;     //error class

};


#endif//__myexcept_h__
