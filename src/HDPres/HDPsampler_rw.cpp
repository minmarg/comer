/***************************************************************************
 *   Copyright (C) 2011 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ext/ivector.h"
#include "HDPsampler.h"

// =========================================================================
// Read: read set of vectors from file; vectors are supposed to be divided 
//  into groups
//
void HDPsampler::Read()
{
    if( !GetFilename())
        throw myruntime_error( "No filename pattern." );

    Ivector dids;//dish indices
    mystring grpfile = mystring( GetFilename()) + _sch_ext_grp_;
    mystring parfile = mystring( GetFilename()) + _sch_ext_par_;
    HDPbase::ReadGroups( grpfile.c_str(), &dids );
    if( file_exists( parfile.c_str())) {
        HDPbase::ReadParameters( parfile.c_str(), &dids);
        SetParsRead( true );
    }
}

// =========================================================================
