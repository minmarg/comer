/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdlib.h>
#include <string.h>

#include "version.h"

#include "Database.h"
#include "GapScheme.h"
#include "DistributionMatrix.h"

#include "Serializer.h"



const char* Database::db_signature[] = {
    "PROFILE DATABASE VERSION ",
    dbversion,
    NULL
};

const char* Database::db_extensions[] = {//extensions of the database files
    ".prd",
    ".frq",
    ".prb",
    NULL
};

// -------------------------------------------------------------------------
// CLASS HitInformation
//
// Constructor
//
Database::Database( const char* name, int fstype )
:
    dbname( name ),
    data_dir( NULL ),
    profiles( NULL ),
    no_profs( 0 ),

    probnorm( 0.0 ),
    no_vectors( 0 ),
    no_sequences( 0 ),
    db_size( 0 ),
    nextsym_( 0 ),
    name_buffer( NULL ),

    usingseg( false ),
    segwinlen( 0 ),
    seglowentropy( -1.0 ),
    seghighentropy( -1.0 ),
    segdistance( 0.0 ),

    store( NULL ),
    fstrtype_( fstype )
{
    Init();
}

// Constructor
//
Database::Database( const char* output, const char* directory, int fstype )
:
    dbname( output ),
    data_dir( directory ),
    profiles( NULL ),
    no_profs( 0 ),

    probnorm( 0.0 ),
    no_vectors( 0 ),
    no_sequences( 0 ),
    db_size( 0 ),
    nextsym_( 0 ),
    name_buffer( NULL ),

    usingseg( false ),
    segwinlen( 0 ),
    seglowentropy( -1.0 ),
    seghighentropy( -1.0 ),
    segdistance( 0.0 ),

    store( NULL ),
    fstrtype_( fstype )
{
    Init();
}

// Constructor
//
Database::Database( const char* output, char* arguments[], int no_args, int fstype )
:
    dbname( output ),
    data_dir( NULL ),
    profiles( arguments ),
    no_profs( no_args ),

    probnorm( 0.0 ),
    no_vectors( 0 ),
    no_sequences( 0 ),
    db_size( 0 ),
    nextsym_( 0 ),
    name_buffer( NULL ),

    usingseg( false ),
    segwinlen( 0 ),
    seglowentropy( -1.0 ),
    seghighentropy( -1.0 ),
    segdistance( 0.0 ),

    store( NULL ),
    fstrtype_( fstype )
{
    Init();
}

// Default constructor
//
Database::Database( int fstype )
:
    dbname( NULL ),
    data_dir( NULL ),
    profiles( NULL ),
    no_profs( 0 ),
    probnorm( 0.0 ),
    no_vectors( 0 ),
    no_sequences( 0 ),
    db_size( 0 ),
    nextsym_( 0 ),
    name_buffer( NULL ),

    usingseg( false ),
    segwinlen( 0 ),
    seglowentropy( -1.0 ),
    seghighentropy( -1.0 ),
    segdistance( 0.0 ),

    store( NULL ),
    fstrtype_( fstype )
{
    throw( myruntime_error( mystring( "Default object initialization is not allowed." )));
}

// Destructor
//
Database::~Database()
{
    Close();
    if( name_buffer )
        free( name_buffer );
    if( store )
        delete store;
}

// -------------------------------------------------------------------------
// Init: necessary initialization
// -------------------------------------------------------------------------

void Database::Init()
{
    for( int n = 0; n < cntFiles; n++ )
        db_fp[n] = NULL;

    if( !GetDbName())
        return;

    name_buffer = ( char* )malloc( strlen( GetDbName()) + strlen( db_extensions[MAIN]) + 1 );

    if( !name_buffer )
        throw myruntime_error( mystring( "Database: Not enough memory." ));

    strcpy( name_buffer, GetDbName());

    store = new FrequencyStore( GetFreqStrType());

    if( !store )
        throw myruntime_error( mystring( "Database: Not enough memory." ));
}

// -------------------------------------------------------------------------
// GetMainDbName: obtains the main database name
// -------------------------------------------------------------------------

const char* Database::GetMainDbName()
{
    strcpy( name_buffer + strlen( GetDbName()), db_extensions[MAIN] );
    return  name_buffer;
}

// -------------------------------------------------------------------------
// GetFreqDbName: obtains the name of frequency file
// -------------------------------------------------------------------------

const char* Database::GetFreqDbName()
{
    strcpy( name_buffer + strlen( GetDbName()), db_extensions[FREQ] );
    return  name_buffer;
}

// -------------------------------------------------------------------------
// GetProbsDbName: get profile pair probabilities filename
//
const char* Database::GetProbsDbName()
{
    strcpy( name_buffer + strlen( GetDbName()), db_extensions[PPRO] );
    return  name_buffer;
}

// -------------------------------------------------------------------------
// Open: opens the database and initializes file descriptor related to it
// -------------------------------------------------------------------------

void Database::Open()
{
    if( !GetMainDbName())
        throw myruntime_error( mystring( "Unable to open database." ));

    if( db_fp[MAIN] != NULL )
        throw myruntime_error( mystring( "Unable to open database." ));

//     db_fp[MAIN] = fopen( GetMainDbName(), "rb" );//obsolete
    db_fp[MAIN] = fopen( GetMainDbName(), "r" );

    if( db_fp[MAIN] == NULL )
        throw myruntime_error( mystring( "Failed to open database." ));

    try {
        SetNextSym( 0 );
//         GetHeader( db_fp[MAIN] );//obsolete
        ReadTextHeader( db_fp[MAIN] );

        //read profile pair probabilities
        proprobs_.ReadProbs( GetProbsDbName());

    } catch( myexception const& ex )
    {
        Close();
        throw myruntime_error( ex.what(), ex.eclass());
    }
}

// -------------------------------------------------------------------------
// Close: closes the database
// -------------------------------------------------------------------------

void Database::Close( TFile which )
{
    if( cntFiles <= which || which < 0 ) {
        for( int n = 0; n < cntFiles; n++ )
            if( db_fp[n] != NULL ) {
                fclose( db_fp[n] );
                db_fp[n] = NULL;
            }
    } else
        if( db_fp[which] != NULL ) {
            fclose( db_fp[which] );
            db_fp[which] = NULL;
        }
}

// -------------------------------------------------------------------------
// Next: reads next bunch of profile matrices from the database; returns
//     flag indicating the end of the database
// -------------------------------------------------------------------------

bool Database::Next( FrequencyMatrix& freq, LogOddsMatrix& pssm, GapScheme& gaps,
                    int gapopencost, int gapextncost, bool fixedcosts )
{
    int nsym = 0;

    if( GetNextSym() == EOF )
        return false;

    if( GetDbDesc() == NULL )
        throw myruntime_error( mystring( "Unable to get information: Database is not opened." ));


    try {
//         serializer.DeserializeProfile( freq, pssm, gaps, GetDbDesc());//obsolete
        serializer.ReadProfile( GetDbDesc(), freq, pssm, gaps );

    } catch( myexception const& ex )
    {
        Close();
        throw myruntime_error( ex.what(), ex.eclass());
    }

    freq.CheckForAllZeros();    //necessary to ensure valid scores

    //set gap costs in order to properly compute position-specific costs
    gaps.Prepare( gapopencost, gapextncost, fixedcosts );


    if( GetNextSym() != EOF ) {
        nsym = fgetc( GetDbDesc());
        ungetc( nsym, GetDbDesc());
        SetNextSym( nsym );
    }
    return true;
}

// -------------------------------------------------------------------------
// Make: creates the database by processing profiles
// -------------------------------------------------------------------------

void Database::Make()
{
    if( !GetMainDbName())
        throw myruntime_error( mystring( "Unable to make database." ));

    if( db_fp[MAIN] != NULL || db_fp[FREQ] != NULL )
        throw myruntime_error( mystring( "Unable to make database." ));

    message( "Preprocessing..." );

    try {

        ResetProbNormTerm();

        //clear marginal counts before updating them
        proprobs_.ClearMargs();

        //1st pass: find out attributes of the database to be created
        ProcessInput( &Database::ProcessFile, NULL );

        //calculate probabilities after marginal counts have been updated
        proprobs_.CalcProbs();


    //     db_fp[MAIN] = fopen( GetMainDbName(), "wb" );//obsolete
        db_fp[MAIN] = fopen( GetMainDbName(), "w" );

        if( db_fp[MAIN] == NULL )
            throw myruntime_error("Unable to make database.");

    //     db_fp[FREQ] = fopen( GetFreqDbName(), "wb" );//obsolete
        db_fp[FREQ] = fopen( GetFreqDbName(), "w" );

        if( db_fp[FREQ] == NULL )
            throw myruntime_error("Unable to make database.");


        message( "Writing..." );
//         PutHeader( db_fp[MAIN] );//obsolete
        WriteTextHeader( db_fp[MAIN] );
        //2nd pass: write down the attributes and profiles themselves
        ProcessInput( &Database::WriteFile, db_fp );

        message( "Writing frequencies..." );
        WriteFrequencies( db_fp[FREQ] );

        message( "Writing probabilities..." );
        proprobs_.WriteProbs( GetProbsDbName());

    } catch( myexception const& ex )
    {
        Close();
        throw myruntime_error( ex.what(), ex.eclass());
    }

    message( "Finished." );
    Close();
}

// -------------------------------------------------------------------------
// ProcessInput: process profiles given explicitly or found in directory
// -------------------------------------------------------------------------

void Database::ProcessInput( PMETHOD processing_method, void* args )
{
    if( GetDataDir()) {
        struct stat     info;
        dirent*         entry = NULL;
        DIR*            direct = opendir( GetDataDir());
        mystring        dirname = GetDataDir() + mystring( DIRSEP );
        mystring        fullname;

        if( !direct )
            throw myruntime_error( mystring( "Unable to process the specified directory." ));

        while(( entry = readdir( direct )))
        {
            fullname = dirname + entry->d_name;

            if( stat( fullname.c_str(), &info ) == -1 ) {
                error( mystring( mystring( "Cannot stat '" ) + mystring( entry->d_name ) + "'; skipping." ).c_str());
                continue;
            }
            if(( S_IFMT & info.st_mode ) != S_IFREG )
//             if( !S_ISREG( info.st_mode ))
                continue;

            ( this->*processing_method )( fullname.c_str()/*entry->d_name*/, args );
        }
        closedir( direct );

    } else
        if( GetProfiles()) {
            bool    option = false;

            for( int n = 0; n < GetNoProfs(); n++ ) {
                if( option ) {
                    option = false;
                    continue;
                }
                if( GetProfile( n )[0] == '-' && GetProfile( n )[2] == 0 ) {
                    if( GetProfile( n )[1] != 'v' )
                        option = true;  //assume it is an option
                    continue;
                }
                ( this->*processing_method )( GetProfile( n ), args );
            }
        }
}

// -------------------------------------------------------------------------
// PreprocessProfile: Preprocesses profile: performs segment
//     masking if needed
// -------------------------------------------------------------------------

void Database::PreprocessProfile( FrequencyMatrix& freq, LogOddsMatrix& pssm, GapScheme& gaps )
{
    if( GetUsingSeg()) {
        //SEG logic
        SEGProfile  segpro(
                freq,
                pssm,
                GetSegWinLength(),
                GetSegLowEntropy(),
                GetSegHighEntropy()
        );
        segpro.SetDistance( GetSegDistance());
        segpro.Run();
        segpro.MaskSeggedPositions( freq, pssm, gaps );
    }
}

// -------------------------------------------------------------------------
// ProcessFile: read data from file and process it
// -------------------------------------------------------------------------

void Database::ProcessFile( const char* filename, void* )
{
    FrequencyMatrix     freq;
    LogOddsMatrix       pssm;
    GapScheme           gaps;
    double              prob = 0.0;

    try {
//         serializer.DeserializeProfile( freq, pssm, gaps, filename );//obsolete
        serializer.ReadProfile( filename, freq, pssm, gaps );
        PreprocessProfile( freq, pssm, gaps );
        IncNoSequences();
        IncDbSize( pssm.GetColumns());

        //update marginal counts
        proprobs_.UpdateMargs( pssm.GetEffNoSequences(), pssm.GetColumns());


        //store all different columns of frequency matrix in
        //  order to write them later to a file
        freq.CheckForAllZeros();    //necessary to ensure valid scores

        for( int n = 0; n < freq.GetColumns(); n++ )
        {
            //omit profile columns constructed for Xs
            if( freq.GetResidueAt( n ) == X )
                continue;

            FrequencyVector vect(   *freq.GetVectorAt( n ), pssm.GetFrequencyWeightAt( n ), 
                                     pssm.GetMIDExpNoObservationsAt( n, PS_M ),
                                    *pssm.GetVectorAt( n ), pssm.GetInformationAt( n ));
            IncNoVectors();
            FrequencyVector quest_vect = ( const char* )store->Store( vect );

            if( quest_vect.GetVector() == vect.GetVector()) {
                //vector has been stored
                //compute probability of this frequency vector to occur
                prob = vect.ComputeProbability();
                //adjust probability-normalizing term
                IncProbNormTerm( prob );
            }
            else {
                //if the vector has not been stored, update probabilities of the original...
                quest_vect.UpdateProbability();
                //and destroy the duplicate
                vect.Destroy();
            }
        }

    } catch( myexception const& ex )
    {
        error( mystring( ex.what() + mystring( " Skipping '" ) + filename + mystring( "'." )).c_str());
    }
}

// -------------------------------------------------------------------------
// WriteFile: write contents of file, i.e. profile matrices, to the
//     database
// -------------------------------------------------------------------------

void Database::WriteFile( const char* filename, void* args )
{
    FILE**  FP = ( FILE** )args;

    if( !FP || !FP[MAIN] || !FP[FREQ] )
        return;

    FrequencyMatrix     freq;
    LogOddsMatrix       pssm;
    GapScheme           gaps;

    try {
//         serializer.DeserializeProfile( freq, pssm, gaps, filename );//obsolete
        serializer.ReadProfile( filename, freq, pssm, gaps );
        PreprocessProfile( freq, pssm, gaps );

//         serializer.SerializeProfile( freq, pssm, gaps, FP[MAIN] );//obsolete
        serializer.WriteProfile( FP[MAIN], freq, pssm, gaps );

    } catch( myexception const& ex )
    {
        error( mystring( ex.what() + mystring( " Skipping '" ) + filename + mystring( "'." )).c_str());
    }
}

// -------------------------------------------------------------------------
// WriteFrequencies: write data vectors to the database file
// -------------------------------------------------------------------------

void Database::WriteFrequencies( FILE* fd )
{
    if( !fd || !store )
        return;

    const SimpleVector* frequencies = store->GetFrequencies();
    size_t              size = 0;
    size_t              n;

    if( frequencies )
        size = frequencies->GetSize();

    //write total number of frequency vectors and number of distinct
    //frequency vectors, respectively
//     Serializer::Write( fd, ( char* )&no_vectors, sizeof( no_vectors ), 1 );//obsolete
//     Serializer::Write( fd, ( char* )&size, sizeof( size ), 1 );//obsolete

    fprintf( fd, "%ld %ld\n", GetNoVectors(), size );

    for( n = 0; n < size; n++ ) {
        FrequencyVector vector(( const char* )frequencies->GetValueAt( n ));
        //normalize probability before writing frequency vector
        vector.NormalizeProbability( GetProbNormTerm());
//         serializer.SerializeFrequencies( vector, fd );//obsolete
        serializer.WriteVector( fd, vector );
    }

    char    strbuf[UCHAR_MAX];

    sprintf( strbuf, "distinct vectors,     %9d", size );                      message( strbuf, false );
    sprintf( strbuf, "1st-level collisions, %9d", store->GetCollisions1());    message( strbuf, false );
    sprintf( strbuf, "2nd-level collisions, %9d", store->GetCollisions2());    message( strbuf );
}

// -------------------------------------------------------------------------
// ReadInFrequencies: read frequencies in the internal storage
// -------------------------------------------------------------------------

void Database::ReadInFrequencies()
{
    if( !store )
        return;

    size_t          length, rbts;
    const size_t    locsize = KBYTE;
    char            locbuffer[locsize+1] = {0};
    char*           p;
    size_t          n;
    int             emsg;
    int             intval;
    long long int   llintval;
    size_t      size = 0;       //number of distinct vectors
    double      reps = 0.0;     //number of occurences of vectors in the database
    mystring    errstr;
    int         eclass;

    if( !GetFreqDbName())
        throw myruntime_error( mystring( "Unable to open one of database files." ));

    if( db_fp[FREQ] != NULL )
        throw myruntime_error( mystring( "Unable to open one of database files." ));

//     db_fp[FREQ] = fopen( GetFreqDbName(), "rb" );//obsolete
    db_fp[FREQ] = fopen( GetFreqDbName(), "r" );

    if( db_fp[FREQ] == NULL )
        throw myruntime_error( mystring( "Failed to open one of database files." ));

    try {
//         Serializer::Read( db_fp[FREQ], ( char* )&no_vectors, sizeof( no_vectors ), 1 );//obsolete
//         Serializer::Read( db_fp[FREQ], ( char* )&size, sizeof( size ), 1 );//obsolete

        if(( emsg = skip_comments( db_fp[FREQ], locbuffer, locsize, &length )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( feof( db_fp[FREQ] ) || !length )
            throw myruntime_error( "Wrong format of database file of vectors." );

        p = locbuffer;

        //number of vectors
        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong format of data vectors." );

        if(( emsg = read_llinteger( p, length - size_t( p - locbuffer ), &llintval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( llintval < 1 )
            throw myruntime_error( "Wrong format of data vectors" );

        SetNoVectors( llintval );
        p += rbts;

        //number of distinct vectors
        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong format of data vectors." );

        if(( emsg = read_llinteger( p, length - size_t( p - locbuffer ), &llintval, &rbts )) != 0 )
            throw myruntime_error( TranslateReadError( emsg ));

        if( llintval < 1 )
            throw myruntime_error( "Wrong format of data vectors" );

        size = llintval;
        p += rbts;


        //allocate needed memory at once to avoid excess and sequential reallocation
        store->SetNoFrequencyVectors( size );

        for( n = 0; n < size; n++ ) {
            FrequencyVector vector;
            //method's smart to allocate required memory
//             serializer.DeserializeFrequencies( vector, db_fp[FREQ] );//obsolete
            serializer.ReadVector( db_fp[FREQ], vector );

            if( store->Store( vector ) != vector.GetVector())
                throw myruntime_error( mystring( "Failed to store data vectors." ));

            reps += vector.GetProbability() + 1.0;
            vector.FinalProbability( GetNoVectors());
        }

        if( ! AreVectorsConsistent( reps ))
            throw myruntime_error( mystring( "Data vectors corrupted." ));

    } catch( myexception const& ex )
    {
        errstr = ex.what();
        eclass = ex.eclass();
    }

    Close( FREQ );

    if( !errstr.empty())
        throw myruntime_error( errstr.c_str(), eclass );
}



// -------------------------------------------------------------------------
// PutHeader: puts header information of the database
// -------------------------------------------------------------------------

void Database::PutHeader( FILE* fp )
{
    TFVectorProbabilities   distrib = GetDistributionType();

    PutSignature( fp );
    Serializer::Write( fp, ( char* )&no_sequences, sizeof( no_sequences ), 1 );
    Serializer::Write( fp, ( char* )&db_size, sizeof( db_size ), 1 );
    Serializer::Write( fp, ( char* )&distrib, sizeof( distrib ), 1 );
}

// -------------------------------------------------------------------------
// GetHeader: gets header information from the database
// -------------------------------------------------------------------------

void Database::GetHeader( FILE* fp )
{
    TFVectorProbabilities   distrib;

    GetSignature( fp );
    Serializer::Read( fp, ( char* )&no_sequences, sizeof( no_sequences ), 1 );
    Serializer::Read( fp, ( char* )&db_size, sizeof( db_size ), 1 );
    Serializer::Read( fp, ( char* )&distrib, sizeof( distrib ), 1 );

    switch( distrib ) {
        case DISCRETE:
        case PROVECTOR:
        case MULTINOMIAL:
            break;

        default:
            throw myruntime_error( mystring( "Database: Unknown vector distribution type." ));
    }

    SetDistributionType( distrib );
}



// =========================================================================

const char* patstrDBVER = "COMER profile database v";
const char* patstrNoPROFS = "PROFS:";
const char* patstrSIZE = "SIZE:";
const char* patstrDISTR = "DISTR:";

// =========================================================================
// WriteTextHeader: write database header to file descriptor
//
void Database::WriteTextHeader( FILE* fp )
{
    if( !fp )
        return;

    TFVectorProbabilities   distrib = GetDistributionType();
    mystring                distrstr = GetDistributionText( distrib );

    if( distrstr.empty())
        throw myruntime_error( mystring( "Database: Unknown vector distribution type." ));

    fprintf( fp, "%s%s\n", patstrDBVER, dataversion );
    fprintf( fp, "%-9s%ld\n", patstrNoPROFS, GetNoSequences());
    fprintf( fp, "%-9s%ld\n", patstrSIZE, GetDbSize());
    fprintf( fp, "%-9s%s\n", patstrDISTR, distrstr.c_str());
}

// -------------------------------------------------------------------------
// ReadTextHeader: read database header from file descriptor
//
void Database::ReadTextHeader( FILE* fp )
{
    if( !fp )
        return;

    size_t          length, rbts;
    const size_t    locsize = KBYTE;
    char            locbuffer[locsize+1] = {0};
    char*           p;
    int             emsg;
    int             intval;
    long long int   llintval;
    mystring        distrstr;
    TFVectorProbabilities   distrib;


    //read version number
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong database format." );

    if(( p = strstr( locbuffer, patstrDBVER )) == NULL )
        throw myruntime_error( "Wrong database format." );

    p += strlen( patstrDBVER );

    if( length <= size_t( p - locbuffer ))
        throw myruntime_error( mystring( "Wrong database format." ));

    if(( p = strstr( p, dataversion )) == NULL )
        throw myruntime_error( "Wrong database version number." );


    //read number of profiles
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong database format." );

    if(( p = strstr( locbuffer, patstrNoPROFS )) == NULL )
        throw myruntime_error( "Wrong database format." );

    p += strlen( patstrNoPROFS );

    if( length <= size_t( p - locbuffer ))
        throw myruntime_error( "Wrong database format." );

    if(( emsg = read_llinteger( p, length - size_t( p - locbuffer ), &llintval, &rbts )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( llintval < 1 )
        throw myruntime_error( "Wrong database format: Invalid number of profiles." );

    SetNoSequences( llintval );


    //read database size
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong database format." );

    if(( p = strstr( locbuffer, patstrSIZE )) == NULL )
        throw myruntime_error( "Wrong database format." );

    p += strlen( patstrSIZE );

    if( length <= size_t( p - locbuffer ))
        throw myruntime_error( "Wrong database format." );

    if(( emsg = read_llinteger( p, length - size_t( p - locbuffer ), &llintval, &rbts )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( llintval < 1 )
        throw myruntime_error( "Wrong database format: Invalid size." );

    SetDbSize( llintval );


    //read distribution description
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != 0 )
        throw myruntime_error( TranslateReadError( emsg ));

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong database format." );

    if(( p = strstr( locbuffer, patstrDISTR )) == NULL )
        throw myruntime_error( "Wrong database format." );

    p += strlen( patstrDISTR );
    for( ; *p == ' ' || *p == '\t' ; p++ );

    if( length <= size_t( p - locbuffer ))
        throw myruntime_error( mystring( "Wrong database format." ));

    distrstr = p;
    distrib = GetDistributionType( distrstr );

    if( DTypeUNKNOWN <= distrib )
        throw myruntime_error( "Wrong database format: Unrecognized distribution." );

    SetDistributionType( distrib );
}



// -------------------------------------------------------------------------
// PutSignature: puts signature of the database
// -------------------------------------------------------------------------

void Database::PutSignature( FILE* fp )
{
    for( int n = 0; db_signature[n]; n++ )
        Serializer::Write( fp, ( char* )db_signature[n], 1, strlen( db_signature[n] ));
}

// -------------------------------------------------------------------------
// GetSignature: gets signature from the file and checks whether it is valid
// -------------------------------------------------------------------------

void Database::GetSignature( FILE* fp )
{
    char    loc_signature[CHAR_MAX];
    //
    for( int n = 0; db_signature[n]; n++ ) {
        Serializer::Read( fp, loc_signature, 1, strlen( db_signature[n] ));

        if( strncmp( db_signature[n], loc_signature, strlen( db_signature[n] ))) {
            if( n == DATAVER )
                warning( "Database format version number is invalid." );
            else
                throw myruntime_error( mystring( "Database format is invalid." ));
        }
    }
}

