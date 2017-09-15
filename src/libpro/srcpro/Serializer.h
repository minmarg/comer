/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __Serializer__
#define __Serializer__

#define SER_ANSI
// #define SERSCALEFACTOR ( 65536 )
#define SERSCALEFACTOR ( 10000 )

#include <stdio.h>

typedef class DistributionMatrix            FrequencyMatrix;
typedef class ExtendedDistributionMatrix    LogOddsMatrix;
class GapScheme;

class FrequencyVector;

class Serializer
{
public:
    Serializer();
    ~Serializer();

    void    SerializeProfile( const FrequencyMatrix&, const LogOddsMatrix&, const GapScheme& );
    void    SerializeProfile( const FrequencyMatrix&, const LogOddsMatrix&, const GapScheme&, FILE* );
    void    SerializeProfile( const FrequencyMatrix&, const LogOddsMatrix&, const GapScheme&, const char* fname );
    void    DeserializeProfile( FrequencyMatrix&, LogOddsMatrix&, GapScheme& );
    void    DeserializeProfile( FrequencyMatrix&, LogOddsMatrix&, GapScheme&, FILE* );
    void    DeserializeProfile( FrequencyMatrix&, LogOddsMatrix&, GapScheme&, const char* fname );

    void    WriteProfile( FILE*, const FrequencyMatrix&, const LogOddsMatrix&, const GapScheme&, int scale = SERSCALEFACTOR);
    void    WriteProfile( const char* fname, const FrequencyMatrix&, const LogOddsMatrix&, const GapScheme&, int scale  = SERSCALEFACTOR );
    void    ReadProfile( FILE*, FrequencyMatrix&, LogOddsMatrix&, GapScheme& );
    void    ReadProfile( const char* fname, FrequencyMatrix&, LogOddsMatrix&, GapScheme& );


    void    SerializeFrequencies( const FrequencyVector& );
    void    SerializeFrequencies( const FrequencyVector&, FILE* );
    void    DeserializeFrequencies( FrequencyVector& );
    void    DeserializeFrequencies( FrequencyVector&, FILE* );

    void    WriteVector( FILE*, const FrequencyVector& );
    void    ReadVector( FILE*, FrequencyVector& );


    static void     Write( FILE*, char*, size_t, size_t );      //static write
    static void     Read( FILE*, char*, size_t, size_t );       //static read

    void            Write( char*, size_t, size_t );             //write information
    void            Read( char*, size_t, size_t );              //read information

protected:
    void    PutSignature();
    void    GetSignature();

private:
    enum {
            BEGTEXT,
            DATAVER,
            BINSIGN,
            END
    };

#ifdef SER_ANSI
    FILE*   fp;                         // pointer to file structure
#else
    int     fd;                         // file descriptor
#endif
    static const char*  signature[];    // signature to put it or get it in/from the file as an indicator of validness
};


#endif//__Serializer__
