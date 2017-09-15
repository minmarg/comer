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
#include "libpro/srcpro/datapro.h"

#include "mystring.h"
#include "myexcept.h"

#include "AlignmentProbe.h"


////////////////////////////////////////////////////////////////////////////
// CLASS AlignmentProbe
//
// Constructor:
//
AlignmentProbe::AlignmentProbe(
    const int querylen, const char* queryres, 
    const int sbjctlen, const char* sbjctres )
:
    F_( NULL ),
    IND_( NULL ),
    querySize_( querylen ),
    subjectSize_( sbjctlen ),

    queryresidues_( queryres ),
    sbjctresidues_( sbjctres ),

    PATH_( NULL ),

    maxsrowind_(0),
    maxscolind_(0),
    alnscore_(0.0),
    alnsteps_(0)
{
    if( MAXCOLUMNS < querySize_ || querySize_ < 1 ||
        MAXCOLUMNS < subjectSize_ || subjectSize_ < 1 )
        throw myruntime_error( mystring( "AlignmentProbe: Maximum sequence length exceeded." ));

    //one extra is reserved for dynamic programming matrix...
    F_ = ( double(**)[szACluster])malloc( sizeof(double*)*(subjectSize_+1));
    IND_ = ( int(**)[szIndices])malloc( sizeof(int*)*(subjectSize_+1));

    PATH_ = ( int(*)[szDims])malloc( sizeof(int)*(subjectSize_+querySize_)*szDims);

    if( !F_ || !IND_ || !PATH_ )
        throw myruntime_error("AlignmentProbe: Not enough memory.");

    memset( PATH_, 0, sizeof(int)*(subjectSize_+querySize_)*szDims );

    for( int m = 0; m < subjectSize_+1; m++ )
    {
        F_[m] = ( double(*)[szACluster])malloc( sizeof(double)*(querySize_+1)*szACluster );
        IND_[m] = ( int(*)[szIndices])malloc( sizeof(int)*(querySize_+1)*szIndices );

        if( !F_[m] || !IND_[m] )
            throw myruntime_error("AlignmentProbe: Not enough memory.");

        memset( F_[m], 0, sizeof(double)*(querySize_+1)*szACluster );
        memset( IND_[m], 0, sizeof(int)*(querySize_+1)*szIndices );
    }
}

// -------------------------------------------------------------------------
// Destructor: deallocate memory
//
AlignmentProbe::~AlignmentProbe()
{
    int m;

    if( PATH_ ) {
        free( PATH_ );
        PATH_ = NULL;
    }

    for( m = 0; m < subjectSize_+1; m++ ) {
        if( F_ ) free( F_[m]);
        if( IND_ ) free( IND_[m]);
    }

    if( F_ ) { free( F_ ); F_ = NULL; }
    if( IND_ ) { free( IND_ ); IND_ = NULL; }
}

// -------------------------------------------------------------------------
// Run: align and make alignment path between two sequencies
//
void AlignmentProbe::Run()
{
    Align();
    MakeAlignmentPath();
}

// -------------------------------------------------------------------------
// Align: compute DP matrix with backtracing information
//
void AlignmentProbe::Align()
{
    double  score;
    double  newscore;//new score
    double  maxvalue;//max value at position so far 
    int     rowindex = 0;
    int     colindex = 0;
    int     anlength = 0;
    int     n, m, dir;

    maxsrowind_ = maxscolind_ = 0;

    //query sequence
    for( n = 1; n < GetQuerySize()+1; n++ )
    {   //subject sequence
        for( m = 1; m < GetSubjectSize()+1; m++ )
        {
            score = LOSCORES.PrecomputedEntry( sbjctresidues_[m-1], queryresidues_[n-1]);
            newscore = F_[m-1][n-1][aDValues] + score;

            //determine maximum value for this position from the upper-left area
            if( F_[m-1][n][aDValues] < newscore &&
                F_[m][n-1][aDValues] < newscore ) {
                maxvalue = newscore;
                dir = 3;
                maxsrowind_ = m;
                maxscolind_ = n;
            } else
                if( F_[m-1][n][aDValues] < F_[m][n-1][aDValues]) {
                    maxvalue = F_[m][n-1][aDValues];
                    dir = 2;
                } else {
                    maxvalue = F_[m-1][n][aDValues];
                    dir = 1;
                }

            //decide whether to proceed the ungapped alignment or to start a new one
            if( 0.0 < maxvalue ) {
                F_[m][n][aDValues] = maxvalue;
            } else {
                F_[m][n][aDValues] = 0;//maxvalue;
            }
            IND_[m][n][aDir] = dir;
        }
    }
}

// -------------------------------------------------------------------------
// MakeAlignmentPath: infer alignment path from the DP matrix
//
void AlignmentProbe::MakeAlignmentPath()
{
    double  score = 0;//maximum score of the DP matrix
    int     row = 0;//row index of the maximum score
    int     col = 0;//column index of the maximum score
    int     step = 0;//step index
    int     state = 3;//matched

    if( !F_ || !IND_ )
        throw myruntime_error("AlignmentProbe: MakeAlignmentPath: Null matrices.");

    if( !PATH_ )
        throw myruntime_error("AlignmentProbe: MakeAlignmentPath: Null path vector." );

    //index of the maximum value
    row = maxsrowind_;
    col = maxscolind_;
    state = IND_[row][col][aDir];

    if( GetSubjectSize() < row || row < 0 ||
        GetQuerySize() < col || col < 0 )
        throw myruntime_error( "AlignmentProbe: MakeAlignmentPath: Invalid indices.");

    //maximum value
    score = F_[row][col][aDValues];

    if( score <= 0.0 )
        return;

    identities_ = 0;
    for( ;0 < row && 0 < col && 0.0 < F_[row][col][aDValues]; )
    {
        state = IND_[row][col][aDir];
        if( state == 3 )
            if( queryresidues_[col-1] == sbjctresidues_[row-1])
                identities_++;

        PATH_[step][aFirst] = col;//query index
        PATH_[step][aSecnd] = row;//subject index
        step++;

        if( state == 3 ) {
            row--; col--;
        }
        else if( state == 2 )
            col--;
        else if( state == 1 )
            row--;
        else
            throw myruntime_error("AlignmentProbe: MakeAlignmentPath: Invalid state.");
    }

    SetScore( score );
    SetAlnSteps( step );
}

// -------------------------------------------------------------------------
// Print: print alignment to string stream; NOTE: memory allocation assumed 
//
void AlignmentProbe::Print( char* sp )
{
    if( !sp )
        return;
    *sp = 0;//to append to the end of the stream
    Print( &string_print, sp );
}

// -------------------------------------------------------------------------
// Print: print alignment to file
//
void AlignmentProbe::Print( FILE* fp )
{
    Print( &file_print, fp );
}

// -------------------------------------------------------------------------
// Print: print alignment to stream vpn (either file or string)
//
void AlignmentProbe::Print( TPrintFunction print_func, void* vpn )
{
    if( vpn == NULL )
        return;

    const char* queryres = GetQueryResidues();
    const char* sbjctres = GetSubjectResidues();

    if( queryres == NULL || sbjctres == NULL )
        return;

    if( GetQuerySize() <= 0 || GetSubjectSize() <= 0 )
        throw myruntime_error("AlignmentProbe: Print: Sequences of zero length.");

    if( !PATH_ )
        throw myruntime_error("AlignmentProbe: Print: Null alignment path.");

    char    qaa;//query residue
    char    saa;//subject residue

    int     qin;//query index
    int     sin;//subject index

    int     step;
    int     identities = 0, positives = 0, gaps = 0;

    mystring    query; query.reserve( KBYTE );
    mystring    sbjct; sbjct.reserve( KBYTE );
    mystring    match; match.reserve( KBYTE );

    for( step = GetAlnSteps()-1; 0 <= step; step-- )
    {
        qaa = 0; saa = 0;
        if( step < GetAlnSteps()-1  &&  PATH_[step][aFirst] == PATH_[step+1][aFirst])
              query.append( '-' );
              //-1 for DPM indices are one greater...
        else  query.append( qaa = DehashCode( qin = queryres[ PATH_[step][aFirst]-1 ] ));

        if( step < GetAlnSteps()-1  &&  PATH_[step][aSecnd] == PATH_[step+1][aSecnd])
              sbjct.append( '-' );
        else  sbjct.append( saa = DehashCode( sin = sbjctres[ PATH_[step][aSecnd]-1 ] ));

        if( qaa && saa && qaa == saa ) {
            match.append( qaa );
            identities++;
        } else if( qaa && saa && 0 < LOSCORES.PrecomputedEntry( sin, qin )) {
                    match.append( '+' );
                    positives++;
                } else {
                    match.append( ' ' );
                    if( !qaa || !saa )
                        gaps++;
                }
    }

    print_func( vpn, "\n" );

    print_func( vpn, "  Length: Query = %d, Sbjct = %d\n", GetQuerySize(), GetSubjectSize());
    print_func( vpn, "\n" );
    print_func( vpn, " Score = %.2f\n", GetScore());

    if( identities )
        print_func( vpn," Identities = %d/%d (%d%%)", 
                    identities, GetAlnSteps(), int(identities*100/GetAlnSteps()));
    if( positives ) {
        if( identities )
            print_func( vpn, "," );
        print_func( vpn," Positives = %d/%d (%d%%)", 
                    positives, GetAlnSteps(), int(positives*100/GetAlnSteps()));
    }
    if( gaps ) {
        if( identities || positives )
            print_func( vpn, "," );
        print_func( vpn, " Gaps = %d/%d (%d%%)", 
                    gaps, GetAlnSteps(), int(gaps*100/GetAlnSteps()));
    }

    print_func( vpn, "\n\n" );

    const int   width = OUTPUTWIDTH;
    int         begp, endp, ind, n;

    for( n = GetAlnSteps()-1; 0 <= n; n -= width )
    {
        begp = PATH_[n][aFirst];
        endp = PATH_[ ind = ( n - width + 1 > 0 )? n - width + 1: 0 ][ aFirst ];
        if( endp == begp && ind != n )
            endp++;
        if( n < GetAlnSteps()-1 && begp == PATH_[n+1][aFirst])
            begp++;

        print_func( vpn, "Query: %5d %s %-5d\n",
                begp,
                ( query.substr( GetAlnSteps()-n-1, width )).c_str(),
                endp );

        print_func( vpn, "%12c %s\n", 32, ( match.substr( GetAlnSteps()-n-1, width )).c_str());

        begp = PATH_[n][aSecnd];
        endp = PATH_[ ind = ( n - width + 1 > 0 )? n - width + 1: 0 ][ aSecnd ];
        if( endp == begp && ind != n )
            endp++;
        if( n < GetAlnSteps()-1 && begp == PATH_[n+1][aSecnd])
            begp++;

        print_func( vpn, "Sbjct: %5d %s %-5d\n\n",
                begp,
                ( sbjct.substr( GetAlnSteps()-n-1, width )).c_str(),
                endp );
    }
}

// -------------------------------------------------------------------------
// PrintDPMatrix: output DP matrix and related info
//
void AlignmentProbe::PrintDPMatrix( FILE* fp ) const
{
    int l = 0;

    if( fp == NULL )
        return;

    fprintf( fp,"%12c Dynamic programming matrix\n", 32 );

    fprintf( fp, "%9c", 32 );

    for( int m = 0; m < GetSubjectSize(); m++ )
        fprintf( fp, "%4c", DehashCode( GetSubjectResidues()[m] ));


    for( int n = 0; n < GetQuerySize(); n++ ) {

        fprintf( fp, "\n%5d %c   ", n + 1, DehashCode( GetQueryResidues()[n] ));

        for( int m = 0; m < GetSubjectSize(); m++ )
            fprintf( fp, "(%1d;%3d) ",
                IND_[m][n][aDir],
                (int)rint(F_[m][n][aDValues])
            );
    }
    fprintf( fp, "\n" );
}
