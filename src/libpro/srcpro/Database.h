/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __Database__
#define __Database__

#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "rc.h"
#include "mystring.h"
#include "myexcept.h"
#include "SEGProfile.h"
#include "Serializer.h"
#include "FrequencyStore.h"
#include "DBProfileProbs.h"


// _________________________________________________________________________
// Class Database
//

class Database {

    //enumerations for signatures...
    enum {  BEGTEXT,    //signature text
            DATAVER,    //version number
            END
    };
    //enumeration for database extensions...
    enum TFile{
            MAIN,    //main profile database file
            FREQ,    //file of frequency vectors
            PPRO,    //file of profile pair probabilities
            cntFiles
    };
public:
    //typedefs...
    typedef void ( Database::*PMETHOD )( const char*, void* );
    //
    Database( const char* name, int = FrequencyStore::TFreqSimpleVector );
    Database( const char* output, const char* directory, int = FrequencyStore::TFreqSimpleVector );
    Database( const char* output, char* arguments[], int no_args, int = FrequencyStore::TFreqSimpleVector );
    ~Database();

    void                    ReadInFrequencies();        //read frequencies in the internal storage

    void                    Open();                     //open database
    void                    Close( TFile = cntFiles );  //close database
                            //read one profile information...
    bool                    Next( FrequencyMatrix&, LogOddsMatrix&, GapScheme&, int goc, int gec, bool );
                            //make the database by processing and gluing profiles
    void                    Make();

    const char*             GetMainDbName();            //get main database name
    const char*             GetFreqDbName();            //get name of frequency file
    const char*             GetProbsDbName();           //get probabilities filename

    const char*             GetDbName() const           { return dbname; }
    size_t                  GetNoVectors() const        { return no_vectors; }
    size_t                  GetNoSequences() const      { return no_sequences; }
    Uint64                  GetDbSize() const           { return db_size; }

    bool            GetUsingSeg() const         { return usingseg; }
    size_t          GetSegWinLength() const     { return segwinlen; }
    double          GetSegLowEntropy() const    { return seglowentropy; }
    double          GetSegHighEntropy() const   { return seghighentropy; }
    double          GetSegDistance() const      { return segdistance; }

    void            SetSegParameters( size_t winlen, double lowent, double highent, double distance )  {
                        SetUsingSeg( true );
                        SetSegWinLength( winlen );
                        SetSegLowEntropy( lowent );
                        SetSegHighEntropy( highent );
                        SetSegDistance( distance );
                    }

    const FrequencyStore*   GetStore() const { return store; }

    const DBProfileProbs*   GetProProbs() const { return &proprobs_; }


    static mystring                 GetDistributionText( int type );
    static TFVectorProbabilities    GetDistributionType( const mystring& distrstr );
    static TFVectorProbabilities    GetDistributionType()           { return FrequencyStore::GetDistributionType(); }
    static void SetDistributionType( TFVectorProbabilities type )   { FrequencyStore::SetDistributionType( type ); }

protected:
    explicit Database( int = FrequencyStore::TFreqSimpleVector );

    void    Init();                                     //necessary initialization

    void    ProcessInput( PMETHOD, void* );             //process profiles given explicitly or found in directory
    void    PreprocessProfile( FrequencyMatrix&, LogOddsMatrix&, GapScheme& );
    void    ProcessFile( const char* filename, void* ); //process profile information
    void    WriteFile( const char* filename, void* );   //append contents of file, i.e. profile,  to the database
    void    WriteFrequencies( FILE* fd );               //write frequency vectors to file descriptor


    void        SetUsingSeg( bool value )           { usingseg = value; }
    void        SetSegWinLength( size_t value )     { segwinlen = value; }
    void        SetSegLowEntropy( double value )    { seglowentropy = value; }
    void        SetSegHighEntropy( double value )   { seghighentropy = value; }
    void        SetSegDistance( double value )      { segdistance = value; }


    FILE*       GetDbDesc() const       { return db_fp[MAIN]; }
    const char* GetDataDir() const      { return data_dir; }
    const char* const* GetProfiles() const  { return profiles; }
    const char* GetProfile( int ) const; 
    const int   GetNoProfs() const      { return no_profs; }

    double  GetProbNormTerm() const         { return probnorm; }
    void    IncProbNormTerm( double prob )  { probnorm += prob; }
    void    ResetProbNormTerm()             { probnorm = 0.0; }

    void    SetNoVectors( size_t value )    { no_vectors = value; }
    void    IncNoVectors()                  { no_vectors++; }
    void    SetNoSequences( size_t value )  { no_sequences = value; }
    void    IncNoSequences()                { no_sequences++; }
    void    SetDbSize( Uint64 value )       { db_size = value; }
    void    IncDbSize( size_t amount )      { db_size += amount; }

    int     GetFreqStrType() const { return fstrtype_; }

                                                        //verify the consistency of the vectors
    bool        AreVectorsConsistent( double reps ) const;

    void        WriteTextHeader( FILE* );
    void        ReadTextHeader( FILE* );

    void        PutHeader( FILE* );                     //put header of the database
    void        GetHeader( FILE* );                     //get header from the database

    void        PutSignature( FILE* );                  //put signature of the database
    void        GetSignature( FILE* );                  //get signature of the database

private:
    int     GetNextSym() const          { return nextsym_; }
    void    SetNextSym( int value )     { nextsym_ = value; }

private:
    const char*         dbname;             //name of database
    const char*         data_dir;           //directory that contains profiles to be appended to construct the database
    const char* const*  profiles;           //names of profiles used to construct the database
    const int           no_profs;           //number of profiles

    double              probnorm;           //probability-normalizing term computed as sum of all probabilities met so far
    size_t              no_vectors;         //overall number of frequency vectors within the database
    size_t              no_sequences;       //number of sequences within the database
    Uint64              db_size;            //size of the database

    int                 nextsym_;           //next symbol from file
    FILE*               db_fp[cntFiles];    //file descriptors of the database
    Serializer          serializer;         //object to serialize data to file

    char*               name_buffer;        //auxiliary buffer


    bool                usingseg;           //seg in use
    size_t              segwinlen;          //seg window length
    double              seglowentropy;      //seg low entropy threshold
    double              seghighentropy;     //seg  high entropy threshold
    double              segdistance;        //seg  distance of equivalence


    FrequencyStore*     store;              //store of frequency vectors
    const int           fstrtype_;          //frequency structure type

    DBProfileProbs      proprobs_;          //profile pair probabilities

    static const char*  db_signature[];     // database signature
    static const char*  db_extensions[];    // extensions of the database files
};


// /////////////////////////////////////////////////////////////////////////
// INLINES
//
inline
const char* Database::GetProfile( int n ) const
{
#ifdef __DEBUG__
    if( no_profs <= n )
        throw myruntime_error( mystring( "Database: Memory access error." ));
#endif
    return profiles[ n ];
}

// -------------------------------------------------------------------------
// AreVectorsConsistent: verifies if the number of vectors is valid
// -------------------------------------------------------------------------

inline
bool Database::AreVectorsConsistent( double reps ) const
{
    if( !GetStore())
        return false;

    return GetStore()->IsConsistent( reps, no_vectors );
}

// -------------------------------------------------------------------------
// GetDistributionText: return distribution description given type
//
inline
mystring Database::GetDistributionText( int distr )
{
    return( FrequencyStore::GetDistributionText( distr ));
}

// GetDistributionType: return distribution type given description
//
inline
TFVectorProbabilities Database::GetDistributionType( const mystring& distrstr )
{
    return( FrequencyStore::GetDistributionType( distrstr ));
}


#endif//__Database__
