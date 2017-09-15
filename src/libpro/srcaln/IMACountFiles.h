/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __IMACountFiles__
#define __IMACountFiles__

#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "rc.h"
#include "mystring.h"
#include "myexcept.h"
#include "IMACounts.h"

// _________________________________________________________________________
// Class IMACountFiles
//

class IMACountFiles
{
public:
    //typedefs...
    typedef void ( IMACountFiles::*PMETHOD )( const char*, void* );
    //
    IMACountFiles( const char* name );
    IMACountFiles( const char* output, const char* directory );
    IMACountFiles( const char* output, char* arguments[], int no_args );
    ~IMACountFiles();

    void                Make();

    const char*         GetOutputName() const       { return outname_; }

protected:
    explicit IMACountFiles();

    void                ProcessInput( PMETHOD, void* );             //iterate over all files
    void                ProcessFile( const char* filename, void* ); //process one file

    void                OutputScores();

    IMACounts&          GetOvrCounts()              { return overallcounts_; }
    const IMACounts&    GetOvrCounts() const        { return overallcounts_; }

    const char*         GetDataDir() const          { return data_dir_; }
    const char* const*  GetFiles() const            { return files_; }
    const char*         GetFileAt( int ) const;
    const int           GetNoFiles() const          { return no_files_; }

    int                 GetNoProcessedFiles() const { return no_processedfiles_; }
    void                IncNoProcessedFiles()       { no_processedfiles_++; }
    void                ResetNoProcessedFiles()     { no_processedfiles_ = 0; }

private:
    IMACounts           overallcounts_;     //overall observed counts
    const char*         outname_;           //output name
    const char*         data_dir_;          //directory that contains observed frequency files to be processed
    const char* const*  files_;             //names of files to be processed
    const int           no_files_;          //number of files
    int                 no_processedfiles_; //number of processed files
};


// /////////////////////////////////////////////////////////////////////////
// INLINES
//
inline
const char* IMACountFiles::GetFileAt( int n ) const
{
#ifdef __DEBUG__
    if( no_files_ <= n )
        throw myruntime_error( mystring( "IMACountFiles: Memory access error." ));
#endif
    return files_[ n ];
}


#endif//__IMACountFiles__
