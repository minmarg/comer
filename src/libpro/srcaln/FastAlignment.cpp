/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "data.h"

#include "mystring.h"
#include "myexcept.h"

#include "FastAlignment.h"


////////////////////////////////////////////////////////////////////////////
// CLASS FastAlignment
//
// Constructor:
// -------------------------------------------------------------------------

FastAlignment::FastAlignment(
    const AbstractScoreMatrix*  scsystem,
    const char* queryres,
    const char* sbjctres )
:
    F( NULL ),
    IND( NULL ),
    querySize( 0 ),
    subjectSize( 0 ),

    scoreSystem( scsystem ),

    queryresidues( queryres ),
    sbjctresidues( sbjctres ),

    PATH( NULL ),

    alnscore( 0.0 ),
    alnsteps( 0 )
{
    if( scsystem == NULL )
        throw myruntime_error( mystring( "FastAlignment: Uninitialized score system." ));

	querySize = scsystem->GetQuerySize();
	subjectSize = scsystem->GetSubjectSize();

    if( MAXCOLUMNS < querySize || querySize < 1 ||
        MAXCOLUMNS < subjectSize || subjectSize < 1 )
        throw myruntime_error( mystring( "FastAlignment: Number of columns exceeds the maximum allowed." ));

	//one extra is reserved for dynamic programming matrix...
	F 	    = ( double(**)[szACluster] )malloc( sizeof( double* ) * ( subjectSize + 1 ));
	IND     = (     int(**)[szIndices] )malloc( sizeof( int*    ) * ( subjectSize + 1 ));

    PATH    = (     int (*)[szDims] )   malloc( sizeof( int ) * ( subjectSize + querySize ) * szDims );


    if( !F || !IND || !PATH )
        throw myruntime_error( mystring( "FastAlignment: Not enough memory." ));

    memset( PATH, 0, sizeof( int ) * ( subjectSize + querySize ) * szDims );

	for( int m = 0; m < subjectSize + 1; m++ )
    {
        F[m]        = ( double(*)[szACluster] )malloc( sizeof( double ) * ( querySize + 1 ) * szACluster );
        IND[m] 	    = ( int    (*)[szIndices] )malloc( sizeof( int    ) * ( querySize + 1 ) * szIndices );

        if( !F[m] || !IND[m] )
            throw myruntime_error( mystring( "FastAlignment: Not enough memory." ));

        memset( F[m],   0, sizeof( double ) * ( querySize + 1 ) * szACluster );
        memset( IND[m], 0, sizeof( int )    * ( querySize + 1 ) * szIndices );
	}
}

// -------------------------------------------------------------------------
// Destructor: deallocate memory used by this class
// -------------------------------------------------------------------------

FastAlignment::~FastAlignment()
{
    if( PATH )
        free( PATH );

	for( int m = 0; m < subjectSize + 1; m++ ) {
		if( F ) free( F[m] );
		if( IND ) free( IND[m] );
	}

	if( F ) free( F );
	if( IND ) free( IND );
}

// -------------------------------------------------------------------------
// Run: computes dynamic programming matrix to make an alignment between two
//     series
// -------------------------------------------------------------------------

void FastAlignment::Run()
{
    Align();
    MakeAlignmentPath();
}

// -------------------------------------------------------------------------
// Align: computes dynamic programming matrix
// -------------------------------------------------------------------------

void FastAlignment::Align()
{
    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();

    if( !scmatrix )
        throw myruntime_error( "FastAlignment: Unable to align profiles." );

    double      score;                  //score to align the query and subject at a position
    double      newscore;               //new score to be acquired
    double      maxvalue;               //max value so far for a position
    int         rowindex;
    int         colindex;
    int         anlength;

    //operate the query positions
    for( int n = 1; n < GetQuerySize() + 1; n++ )
    {
        //operate the subject positions
        for( int m = 1; m < GetSubjectSize() + 1; m++ )
        {
            score  = scmatrix->GetImageScore( m-1, n-1 );
            newscore = F[m-1][n-1][aDValues] + score;

            //determine maximum value for this position by using information from the upper-left area
            if( F[m-1][n][aMaxvals] < newscore &&
                F[m][n-1][aMaxvals] < newscore ) {
                maxvalue = newscore;
                rowindex = m;
                colindex = n;
            } else
                if( F[m-1][n][aMaxvals] < F[m][n-1][aMaxvals] ) {
                    maxvalue = F[m][n-1][aMaxvals];
                    rowindex = IND[m][n-1][aRow];
                    colindex = IND[m][n-1][aColumn];
                } else {
                    maxvalue = F[m-1][n][aMaxvals];
                    rowindex = IND[m-1][n][aRow];
                    colindex = IND[m-1][n][aColumn];
                }

            //decide whether to proceed the ungapped alignment or to start a new one
            if( 0.0 < newscore ) {
                F[m][n][aDValues] = newscore;
                IND[m][n][aLength] = IND[m-1][n-1][aLength] + 1;
            } else {
                //length of the previous maximum-scored ungapped alignment
                anlength = IND[rowindex][colindex][aLength];
                //if alignment length is less than the sum of gaps to the start of this alignment
                if( anlength <= m + 1 - rowindex + n + 1 - colindex ) {
                    F[m][n][aDValues] = 0.0;
                    IND[m][n][aLength] = 0;
                } else {
                    //do not assign 0,
                    //instead proceed with the maximum-valued ungapped alignment met so far
                    F[m][n][aDValues] = maxvalue;
                    IND[m][n][aLength] = 0;//anlength;
                }
                //to indicate a break between individual ungapped alignments
                IND[m][n][aMarker] = 1;
            }

            //fill in info for maximum value found for this position
            F[m][n][aMaxvals] = maxvalue;
            IND[m][n][aRow] = rowindex;
            IND[m][n][aColumn] = colindex;
        }
    }
}

// -------------------------------------------------------------------------
// MakeAlignmentPath: makes alignment path using dynamic programming matrix
//     constructed
// -------------------------------------------------------------------------

void FastAlignment::MakeAlignmentPath()
{
    double      score = 0;              //maximum score of the dynamic programming matrix
    int         row = 0;                //row index of the maximum score
    int         col = 0;                //column index of the maximum score
    int         step = 0;               //step index

    if( !F || !IND )
        throw myruntime_error( mystring( "FastAlignment: Unable to make alignment path." ));

    if( !PATH )
        throw myruntime_error( mystring( "FastAlignment: Unable to make alignment path: Uninitialized." ));

    //indices of the max value are saved in variable
    row = IND[ GetSubjectSize() ][ GetQuerySize() ][aRow];
    col = IND[ GetSubjectSize() ][ GetQuerySize() ][aColumn];

    if( GetSubjectSize() < row || row < 0 ||
        GetQuerySize()   < col || col < 0 )
        throw myruntime_error( mystring( "FastAlignment: Wrong indices." ));

    //maximum value itself
    score = F[ row ][ col ][aDValues];

    if( score <= 0.0 )
        return;


    while( 0 < row && 0 < col ) {
        while( 0 < row && 0 < col && IND[row][col][aMarker] && 0.0 < F[row][col][aDValues] ) {
            row = IND[row][col][aRow];
            col = IND[row][col][aColumn];
        }

        if( row <= 0 || col <= 0 || F[row][col][aDValues] <= 0.0 )
            break;

        PATH[step][aFirst] = col--; //index for query
        PATH[step][aSecnd] = row--; //index for subject
        step++;
    }

    SetScore( score );
    SetAlnSteps( step );
}

// -------------------------------------------------------------------------
// Print: prints alignment to string stream, space of which must
//     have been allocated before
//

void FastAlignment::Print( char* sp )
{
    if( !sp )
        return;
    *sp = 0;  //to ensure printing to the end of the stream
    Print( &string_print, sp );
}

// -------------------------------------------------------------------------
// Print: prints alignment to file
//

void FastAlignment::Print( FILE* fp )
{
    Print( &file_print, fp );
}

// -------------------------------------------------------------------------
// Print: prints alignment to stream pointed by vpn which can be either
//     file or string
// -------------------------------------------------------------------------

void FastAlignment::Print( TPrintFunction print_func, void* vpn )
{
    if( vpn == NULL )
        return;

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();

    const char* queryres = GetQueryResidues();
    const char* sbjctres = GetSubjectResidues();

    if( queryres == NULL || sbjctres == NULL )
        return;

    if( GetQuerySize() <= 0 || GetSubjectSize() <= 0 )
        throw myruntime_error( mystring( "FastAlignment: Unable to print alignment." ));

    if( !PATH )
        throw myruntime_error( mystring( "FastAlignment: No alignment path." ));

    if( !scmatrix )
        throw myruntime_error( mystring( "FastAlignment: No score matrix." ));


    char    qaa;    //query amino acid
    char    saa;    //subject amino acid

    int     qin;    //index for query 
    int     sin;    //index for subject

    int     identities = 0, positives = 0;  //useful statistics


    print_func( vpn, "\n" );

    print_func( vpn, " Query length = %d, Sbjct length = %d\n", GetQuerySize(), GetSubjectSize());
    print_func( vpn, " Score = %.2f\n", GetScore());

    if( identities )
        print_func( vpn, " Identities = %d/%d (%d%%)", identities, GetAlnSteps(), int( identities * 100 / GetAlnSteps() ));
    if( positives ) {
        if( identities )
            print_func( vpn, "," );
        print_func( vpn, " Positives = %d/%d (%d%%)", positives, GetAlnSteps(), int( positives * 100 / GetAlnSteps() ));
    }


    print_func( vpn, "\n\n" );

    mystring    query;      //query sequence of alignment
    mystring    sbjct;      //subject sequence of alignment
    mystring    match;      //match between query and subject strings
    int         laststep;   //last step

    laststep = GetAlnSteps() - 1;

    for( int step = GetAlnSteps() - 1; step >= 0; step-- )
    {
        qaa = 0; saa = 0;
        //-1 -- since prog. matrix indices are one greater...
        query += ( qaa = DehashCode( queryres[ qin = PATH[step][aFirst]-1 ] ));
        sbjct += ( saa = DehashCode( sbjctres[ sin = PATH[step][aSecnd]-1 ] ));

        if( qaa && saa && qaa == saa ) {
            match += qaa;
            identities++;
        } else
            if( qaa && saa && scmatrix->GetImageScore( sin, qin ) > 0 ) {
                match += '+';
                positives++;
            } else
                match += ' ';

        if( 0 <= step && (
            PATH[ step-1 ][aFirst] != PATH[ step ][aFirst] + 1 ||
            PATH[ step-1 ][aSecnd] != PATH[ step ][aSecnd] + 1 ))
        {
            PrintIsolatedSegment( print_func, vpn, query, sbjct, match, step, laststep );
            query.erase();
            sbjct.erase();
            match.erase();
            laststep = step - 1;
        }
    }
}

// -------------------------------------------------------------------------
// PrintIsolatedSegment: prints ungapped alignment isolated in the dynamic
//     programing matrix
// -------------------------------------------------------------------------

void FastAlignment::PrintIsolatedSegment(
    TPrintFunction  print_func,
    void*           vpn,
    mystring&       query,
    mystring&       sbjct,
    mystring&       match,
    int             step_from,
    int             step_to )
{
    const int   width = OUTPUTWIDTH;

    for( int n = step_to; n >= step_from; n -= width )
    {
        print_func( vpn, "Query: %5d %s %-5d\n",
                PATH[n][aFirst],
                query.substr( step_to - n, width ).c_str(),
                PATH[ ( n - width + 1 > step_from )? n - width + 1: step_from ][ aFirst ] );

        print_func( vpn, "%12c %s\n", 32, match.substr( step_to - n, width ).c_str());

        print_func( vpn, "Sbjct: %5d %s %-5d\n\n",
                PATH[n][aSecnd],
                sbjct.substr( step_to - n, width ).c_str(),
                PATH[ ( n - width + 1 > step_from )? n - width + 1: step_from ][ aSecnd ] );
    }

    print_func( vpn, "\n" );
}

// -------------------------------------------------------------------------
// PrintScoringMatrix: output the computed scoring system
// -------------------------------------------------------------------------

void FastAlignment::PrintDPMatrix( FILE* fp ) const
{
    int     l = 0;

    if( fp == NULL )
        return;

    fprintf( fp,"%12c Dynamic programming matrix\n", 32 );

    fprintf( fp, "%9c", 32 );

    for( int m = 0; m < GetSubjectSize(); m++ )
        fprintf( fp, "%4c", DehashCode( GetSubjectResidues()[m] ));


    for( int n = 0; n < GetQuerySize(); n++ ) {

        fprintf( fp, "\n%5d %c   ", n + 1, DehashCode( GetQueryResidues()[n] ));

        for( int m = 0; m < GetSubjectSize(); m++ )
            fprintf( fp, "(%1d,%2d,%3d,%3d;%3d,%3d) ",
                IND[m][n][aMarker],
                IND[m][n][aLength],
                IND[m][n][aColumn],
                IND[m][n][aRow],
                ( int )rint( F[m][n][aDValues] ),
                ( int )rint( F[m][n][aMaxvals] )
            );
    }
    fprintf( fp, "\n" );
}
