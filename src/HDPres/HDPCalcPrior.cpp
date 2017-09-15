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
#include "data.h"
#include "libpro/srcpro/datapro.h"
#include "HDPCalcPrior.h"

// -------------------------------------------------------------------------
// constructor: initialization
//
HDPCalcPrior::HDPCalcPrior()
:   filename_( NULL ),
    outputfile_( NULL )
{
}

// -------------------------------------------------------------------------
// destructor:
//
HDPCalcPrior::~HDPCalcPrior()
{
}

// -------------------------------------------------------------------------
// Run: begin sampling procedure
//
void HDPCalcPrior::Run()
{
    if( GetDegFAdjustment() < -1.0 )
        throw myruntime_error("Invalid value of adjustment to deg. of freedom.");

    const double    tol = 0.1;//error tolerance for target probabilities
    int         ctx = 1;
    int         dim = NUMAA;
    Pslvector   mvec( dim );
    int a;

    if( !GetUninfPrior() && GetFilename()) {
        message("Read...");
        Read();
        message("Calculating prior parameters...");
        CalcPriorParams();
    }
    else {
        message("Setting prior parameters...");
        for( a = 0; a < dim; a++ )
            mvec.SetValueAt( a, LOSCORES.PROBABility( a ));
        HDPbase::LogitNormal2Normal( &mvec, tol );
        //reduce dimensionality
        mvec.DecDim();
        SetUninfPriorParamsS0I( mvec.GetSize(), ctx, mvec );
//         SetUninfPriorParamsS0Cov( mvec.GetSize(), ctx, mvec );
    }

    message("Saving...");
    PrintParameters( GetOutputFile());
}


// -------------------------------------------------------------------------
// SetUninfPriorParamsS0Cov: set uninformative prior parameters of NIW 
//     distribution; scale matrix S0 is proportional to Cov matrix obtained 
//     from the default Odds matrix
//
int HDPCalcPrior::SetUninfPriorParamsS0Cov( int dim, int ctx, const Pslvector& mvec )
{
    if( GetS0ScaleFac() <= 0.0 )
        throw myruntime_error("HDPCalcPrior: SetUninfPriorParamsS0Cov: Invalid scale factor.");
    if( dim != NUMAA - 1 || ctx < 1 )
        throw myruntime_error("HDPCalcPrior: SetUninfPriorParamsS0Cov: Invalid dimensionality.");
    if( dim != mvec.GetSize())
        throw myruntime_error("HDPCalcPrior: SetUninfPriorParamsS0Cov: Inconsistent dimensionality.");

    InitMenu( 1 );

    const int           effnores = NUMAA;
    SPDmatrix           O( effnores );//odds
    SPDmatrix           Oinv( effnores );//odds inverse
    Pslvector           bprob( effnores);//background probabilities consistent with the odds
    Pslvector           unit( effnores);//unit vector
    Pslvector*          mu0 = NULL;
    SPDmatrix*          S0 = NULL;
    SPDmatrix*          S0inv = NULL;
    mystring    preamb = "SetUninfPriorParamsS0Cov: ";
    mystring    errstr;
    double      scale = GetS0ScaleFac();
    double      ldet = 0.0;
    double      val;
    int         err;
    int         a, b;

    try {
        unit.SetUnity();
        mu0 = new Pslvector( mvec );
        S0 = new SPDmatrix( dim );
        S0inv = new SPDmatrix( dim );
        if( mu0 == NULL || S0 == NULL || S0inv == NULL )
            throw myruntime_error( "Not enough memory." );

        for( a = 0; a < effnores; a++ )
            for( b = 0; b < effnores; b++ )
                O.SetValueAt( a, b, LOSCORES.FreqRatio( a, b ));

        Oinv = O;

        if(( err = Oinv.CholeskyDecompose()) == 0 ) {
            if(( err = Oinv.CDedLogDet( &ldet )) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
            if(( err = Oinv.CDedInvert()) != 0 )
                throw myruntime_error( TranslatePSLError( err ));
        }
        else
            throw myruntime_error( TranslatePSLError( err ));

        if(( err = bprob.Mul( Oinv, unit )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        for( a = 0; a < dim; a++ )
            for( b = 0; b < dim; b++ ) {
                //<a,b> - <a><b>
                val = O.GetValueAt( a, b ) * bprob.GetValueAt(a) * bprob.GetValueAt(b) - 
                      LOSCORES.PROBABility(a) * LOSCORES.PROBABility(a);
                S0->SetValueAt( a, b, val );
            }

        if(( err = S0->Scale( scale )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        *S0inv = *S0;

        if(( err = S0inv->CholeskyDecompose()) != 0 )
            throw myruntime_error( TranslatePSLError( err ));
        if(( err = S0inv->CDedLogDet( &ldet )) != 0 )
            throw myruntime_error( TranslatePSLError( err ));
        if(( err = S0inv->CDedInvert()) != 0 )
            throw myruntime_error( TranslatePSLError( err ));

        GetMenu()->SetS0ScaleFac( GetS0ScaleFac());
        GetMenu()->SetS0( S0 ); S0 = NULL;
        GetMenu()->SetInvS0( S0inv ); S0inv = NULL;
        GetMenu()->SetLDetS0( ldet );
        GetMenu()->SetMu0( mu0 ); mu0 = NULL;
        GetMenu()->SetKappa0( 1.0 / GetS0ScaleFac());
        GetMenu()->SetNu0( GetS0ScaleFac());
        GetMenu()->SetDim( dim );
        GetMenu()->SetCtx( ctx );

        //probability factors
        CalcPriorProbFact();

    } catch( myexception const& ex ) {
        errstr = preamb + ex.what();
    }

    if( mu0 ) { delete mu0; mu0 = NULL; }
    if( S0 ) { delete S0; S0 = NULL; }
    if( S0inv ) { delete S0inv; S0inv = NULL; }

    if( !errstr.empty())
        throw myruntime_error( errstr );
    return 0;
}

// =========================================================================
// Read: read set of vectors from file; vectors are supposed to be divided 
//  into groups
//
void HDPCalcPrior::Read()
{
    if( !GetFilename())
        throw myruntime_error( "No filename." );

    Ivector dids;//dish indices
    HDPbase::ReadGroups( GetFilename(), &dids );
}

// =========================================================================
// PrintParameters: print prior parameters
//
void HDPCalcPrior::PrintParameters( const char* filename )
{
    HDPbase::PrintParameters( filename );
}

// -------------------------------------------------------------------------
// PrintParameters: print prior parameters
//
void HDPCalcPrior::PrintParameters( FILE* fp )
{
    if( fp == NULL )
        return;
    HDPbase::PrintPriorParams( fp );
}
