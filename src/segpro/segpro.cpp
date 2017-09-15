/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/



#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <errno.h>
#include <getopt.h>
#include <stdlib.h>

#include "rc.h"
#include "segpro.h"

#include "libpro/srcpro/MOptions.h"
#include "libpro/srcpro/Serializer.h"
#include "libpro/srcpro/GapScheme.h"
#include "libpro/srcpro/DistributionMatrix.h"
#include "libpro/srcpro/SEGProfile.h"
#include "libpro/srcaln/InputMultipleAlignment.h"


// #include "segdata.h"

// Functional-interface-----------------------------------------------------

void PrintSequence( mystring&, SEGProfile& segpro, const char* description, size_t width );
void MaskMultipleAlignment( mystring& input, mystring& output, SEGProfile& segpro, size_t width );
void PutMASequence( FILE*, const char* title, const char* sequence, size_t length, char, size_t width );

// =========================================================================

int main( int argc, char *argv[] )
{
    //string values of the parameters
    mystring        input;
    mystring        output;
    mystring        ascii;
    mystring        seqfilename;
    mystring        fastawidth;
    mystring        identity;
    mystring        segwindow;
    mystring        seglowent;
    mystring        seghighent;
    mystring        segdistance;
    bool            suppress = true;    //suppress warnings
    int             c;

    SetGlobalProgName( argv[0], version );

    if( argc <= 1 ) {
        fprintf( stdout, "%s", usage( argv[0], instructions, version, verdate ).c_str());
        return EXIT_SUCCESS;
    }

    while( 1 ) {
#ifdef USE_GETOPT_LONG
        int option_index = 0;
        static struct option long_options[] = {
            {"i",       required_argument, 0, 'i'},
            {"o",       required_argument, 0, 'o'},

            {"p",       required_argument, 0, 'p'},
            {"s",       required_argument, 0, 's'},
            {"L",       required_argument, 0, 'L'},

            {"t",       required_argument, 0, 't'},

            {"v",       no_argument,       0, 'v'},

            {"w",       required_argument, 0, 'w'},
            {"f",       required_argument, 0, 'f'},
            {"F",       required_argument, 0, 'F'},
            {"D",       required_argument, 0, 'D'},
            { 0, 0, 0, 0 }
        };
        if(( c = getopt_long_only( argc, argv, "hvi:o:p:s:L:t:w:f:F:D:", long_options, &option_index )) == -1 )
            break;
#else
        if(( c = getopt( argc, argv, "hvi:o:p:s:L:t:w:f:F:D:" )) == -1 )
            break;
#endif

        switch( c ) {
            case ':':   fprintf( stdout, "Argument must be followed by the value.\n%s",
                                                usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_FAILURE;
            case '?':   fprintf( stdout, "%s",  usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_FAILURE;
            case 'h':   fprintf( stdout, "%s",  usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_SUCCESS;

            case 'i':   input       = optarg;       break;
            case 'o':   output      = optarg;       break;

            case 'p':   ascii       = optarg;       break;
            case 's':   seqfilename = optarg;       break;
            case 'L':   fastawidth  = optarg;       break;
            case 't':   identity    = optarg;       break;

            case 'w':   segwindow   = optarg;       break;
            case 'f':   seglowent   = optarg;       break;
            case 'F':   seghighent  = optarg;       break;
            case 'D':   segdistance = optarg;       break;

            case 'v':   suppress   = false;         break;
            default:    break;
        }
    }

    char*   p;
    double  tmpval;
    int     intvalue;
    char    strbuf[BUF_MAX];

    mystring valSUBMAT = defSUBMAT;
    double  sidlevel = IDENTITY_LEVEL;
    size_t  segwinlenval = DEFAULT_SEGPRO_WIN_LENGTH;
    double  seglowentval = DEFAULT_SEGPRO_LOW_ENTROPY;
    double  seghighentval = DEFAULT_SEGPRO_HIGH_ENTROPY;
    double  segdistanceval = DEFAULT_SEGPRO_EQUALITY_DISTANCE;
    size_t  outputwidth = DEFAULT_SEGPRO_PRINT_WIDTH;


    // -- sub. matrix --
    try {
        SetLOSCORES( valSUBMAT, DEFAULT_PSCORES_FILENAME );
    } catch( myexception const& ex ) {
        error( ex.what());
        return EXIT_FAILURE;
    }

    if( input.empty()) {
        error( "Input is not specified." );
        return EXIT_FAILURE;
    }
    if( output.empty()) {
        error( "Name of output file is not specified." );
        return EXIT_FAILURE;
    }

    if( !fastawidth.empty()) {
        intvalue = strtol( fastawidth.c_str(), &p, 10 );

        if( errno || *p ) {
            error( "Line width is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue < 1 && intvalue != -1 ) {
            error( "Line width must be strictly positive." );
            return EXIT_FAILURE;
        }
        outputwidth = ( size_t )intvalue;
    }

    if( !identity.empty()) {
        intvalue = strtol( identity.c_str(), &p, 10 );

        if( errno || *p ) {
            error( "Identity threshold specified is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue < 1 || intvalue > 100 ) {
            error( "Identity threshold is not within interval [1-100]." );
            return EXIT_FAILURE;
        }
        sidlevel = ( double )intvalue / 100.0;
    }

    // SEG options --

    if( !segwindow.empty()) {
        intvalue = strtol( segwindow.c_str(), &p, 10 );

        if( errno || *p ) {
            error( "Window length is invalid." );
            return EXIT_FAILURE;
        }
        if( intvalue < 2 ) {
            error( "Window length must be >1." );
            return EXIT_FAILURE;
        }
        segwinlenval = intvalue;
    }

    if( !seglowent.empty()) {
        tmpval = strtod( seglowent.c_str(), &p );
        if( errno || *p ) {
            error( "Low entropy threshold is invalid." );
            return EXIT_FAILURE;
        }
        if( tmpval < 0.0 ) {
            error( "Low entropy threshold is negative." );
            return EXIT_FAILURE;
        }
        seglowentval = tmpval;
    }

    if( !seghighent.empty()) {
        tmpval = strtod( seghighent.c_str(), &p );
        if( errno || *p ) {
            error( "High entropy threshold is invalid." );
            return EXIT_FAILURE;
        }
        if( tmpval < 0.0 ) {
            error( "High entropy threshold is negative." );
            return EXIT_FAILURE;
        }
        seghighentval = tmpval;
    }

    if( !segdistance.empty()) {
        tmpval = strtod( segdistance.c_str(), &p );
        if( errno || *p ) {
            error( "Distance of equivalence is invalid." );
            return EXIT_FAILURE;
        }
        if( tmpval < 0.0 ) {
            error( "Distance of equivalence is negative." );
            return EXIT_FAILURE;
        }
        segdistanceval = tmpval;
    }

    // **

    SetQuiet( suppress );
    message( NULL );

    sprintf( strbuf, "SEG:  %4d, %.2f, %.2f; Distance %.2f",
            segwinlenval, seglowentval, seghighentval, segdistanceval );
    message( strbuf );

    int                     ret = EXIT_SUCCESS;
    InputMultipleAlignment* inmaln = NULL;


    try {
        bool            fasta = false;
        Serializer      serializer;
        FrequencyMatrix frequencies;
        LogOddsMatrix   pssm;
        GapScheme       gaps;


            try {
                // first try reading input in the binary format
//                 serializer.DeserializeProfile( frequencies, pssm, gaps, input.c_str());//obsolete
                serializer.ReadProfile( input.c_str(), frequencies, pssm, gaps );

            } catch( myexception const& ex )
            {
                message( "Trying fasta..." );
                inmaln = new InputMultipleAlignment;

                if( !inmaln ) {
                    error( "Not enough memory." );
                    return EXIT_FAILURE;
                }

                inmaln->SetIdentityLevel( sidlevel );
                inmaln->ReadAlignment( input.c_str());
                inmaln->ConstructProfile();

                inmaln->ExportFrequenciesTo( frequencies );
                inmaln->ExportPSSMTo( pssm );
                inmaln->ExportGapWeightsTo( gaps );

                fasta = true;
            }

        //SEG logic
        SEGProfile  segpro(
                frequencies,
                pssm,
                segwinlenval,
                seglowentval,
                seghighentval
        );
        segpro.SetDistance( segdistanceval );
        segpro.Run();
        segpro.MaskSeggedPositions( frequencies, pssm, gaps );

        try {
            PrintSequence( seqfilename, segpro, pssm.GetDescription(), outputwidth );
        } catch( myexception const& ex ) {
            error( ex.what());
        }
        // --


        if( !ascii.empty())
            //inmaln->OutputProfile( ascii.c_str());
            OutputProfile( ascii.c_str(), frequencies, pssm, gaps );

        if( fasta ) {
            MaskMultipleAlignment( input, output, segpro, outputwidth );
        }
        else
//             serializer.SerializeProfile( frequencies, pssm, gaps, output.c_str());//obsolete
            serializer.WriteProfile( output.c_str(), frequencies, pssm, gaps );


    } catch( myexception const& ex )
    {
        error( ex.what());
        ret = EXIT_FAILURE;
    }

    if( inmaln )
        delete inmaln;

    return ret;
}

// -------------------------------------------------------------------------
// PrintSequence: prints segged sequence to file
// -------------------------------------------------------------------------

void PrintSequence( mystring& filename, SEGProfile& segpro, const char* description, size_t width )
{
    if( filename.empty())
        return;

    FILE*   fp = fopen( filename.c_str(), "w" );

    if( fp == NULL )
        throw myruntime_error( mystring( "Failed to open file for writing of profile sequence." ));

    if( description != NULL )
        fprintf( fp, ">%s\n", description );
    else
        fprintf( fp, ">\n");

    segpro.PrintSeggedSequence( fp, width );

    fclose( fp );
}

// -------------------------------------------------------------------------
// MaskMultipleAlignment: masks multiple alignment given as input and
//     outputs it
// -------------------------------------------------------------------------

void MaskMultipleAlignment( mystring& input, mystring& output, SEGProfile& segpro, size_t width )
{
    if( input.empty() || output.empty())
        return;

    InputMultipleAlignment  inmultaln;

    inmultaln.SetKeepTitles( true );
    inmultaln.SetIgnoreGapsInQuery( false );
    inmultaln.ReadAlignment( input.c_str());

    if( inmultaln.size() == 0 || inmultaln.SequenceAt( 0 ) == NULL )
        return;

    FILE*   fp = fopen( output.c_str(), "w" );

    if( fp == NULL )
        throw myruntime_error( mystring( "Failed to open file for output." ));

    const PosDescriptionVector*     first = inmultaln.SequenceAt( 0 );
    char*       omitmask = ( char* )malloc( first->size() * sizeof( char ));
    char*       seggedone = ( char* )malloc( first->size() * sizeof( char ));
    mystring    errstr;
    size_t      n;
    const char  msym = 'x';

    if( omitmask == NULL || seggedone == NULL ) {
        fclose( fp );
        throw myruntime_error( mystring( "Not enough memory." ));
    }

    for( n = 0; n < first->size(); n++ )
        if( first->ResidueAt( n ) == GAP )
            omitmask[n] = true;
        else
            omitmask[n] = false;


    for( n = 0; n < inmultaln.size(); n++ )
    {
        const PosDescriptionVector* sequence = inmultaln.SequenceAt( n );

        if( sequence == NULL ) {
            errstr = "Memory access error";
            break;
        }

        if( sequence->size() != first->size()) {
            errstr = "Wrong sequence length.";
            break;
        }

        memcpy( seggedone, sequence->GetResidues(), sequence->size() * sizeof( char ));

        try {
            segpro.MaskSequence(
                seggedone,
                sequence->size(),
                msym,
                omitmask,
                GAP,
                ASTERISK
            );
            PutMASequence(
                fp,
                sequence->GetDescription(),
                seggedone,
                sequence->size(),
                msym,
                width
            );
        } catch( myexception const& ex ) {
            errstr = ex.what();
            break;
        }
    }


    free( seggedone );
    free( omitmask );
    fclose( fp );

    if( !errstr.empty())
        throw myruntime_error( errstr );
}

// -------------------------------------------------------------------------
// PutMASequence: puts multiple alignment sequence into file given
// -------------------------------------------------------------------------

void PutMASequence(
    FILE* fp,
    const char* title,
    const char* sequence,
    size_t length,
    char msym,
    size_t width )
{
    if( sequence == NULL )
        return;

    size_t  n = 0;
    size_t  p = 0;
    char*   residues = ( char* )malloc( sizeof( char ) * ( length + 1 ));
    char*   buffer = ( char* )malloc( sizeof( char ) * ((( width < length )? width: length ) + 1 ));

    if( residues == NULL || buffer == NULL )
        throw myruntime_error( mystring( "Not enough memory." ));

    for( n = 0; n < length; n ++ )
        if( sequence[n] != msym )
            residues[n] = DehashCode( sequence[n] );
        else
            residues[n] = sequence[n];


    if( title != NULL )
        fprintf( fp, ">%s\n", title );
    else
        fprintf( fp, ">\n");

    for( n = 0; n < length; n += p ) {
        memcpy( buffer, residues + n, p = ( n + width < length )? width: length - n );
        buffer[p] = 0;
        fprintf( fp, "%s\n", buffer );
    }

    free( residues );
    free( buffer );
}

