/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __HDPCalcPrior__
#define __HDPCalcPrior__

#include <stdio.h>
#include <stdlib.h>

#include "rc.h"
#include "myexcept.h"
#include "libHDP/HDPbase.h"

// -------------------------------------------------------------------------
// class HDPCalcPrior: interface of calculator of prior parameters
//
class HDPCalcPrior: public HDPbase
{
public:
    HDPCalcPrior();
    ~HDPCalcPrior();

    double      GetUninfScaleFactor() { return GetS0ScaleFac(); }
    void        SetUninfScaleFactor( double value ) { SetS0ScaleFac( value ); }

    const char* GetFilename() const { return filename_; }
    void        SetFilename( const char* value ) { filename_ = value; }

    const char* GetOutputFile() const { return outputfile_; }
    void        SetOutputFile( const char* value ) { outputfile_ = value; }

    void        Run();

    void        Read();
//     void        SaveBestAndLast( double*, int, bool last = false );
//     void        Save( mystring& fname, double* = NULL, int* = NULL );
    virtual void    PrintParameters( const char* filename );

protected:
    virtual void    PrintParameters( FILE* );

    int         SetUninfPriorParamsS0Cov( int dim, int ctx, const Pslvector& mvec );

private:
    const char* filename_;//name of file containing the set of vectors
    const char* outputfile_;//name of output file
};


#endif//__HDPCalcPrior__
