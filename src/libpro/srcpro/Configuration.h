/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __Configuration__
#define __Configuration__

#include "compdef.h"
#include "debug.h"

#include "mystring.h"
#include "myexcept.h"


extern int abs( int j );

class Configuration;


// Version history:
//

static const char*  configversion = "0.1";

void SetUngappedParams( Configuration& config );

// Known gap penalty schemes
enum TProcomSchemes {
    ProcomUngapped,
    ProcomGapped,
    NoSchemes
};

static const char*  s_section_ungapped  = "UNGAPPED";
static const char*  s_section_gapped    = "GAPPED";
static const char*  s_space             = "_";

static const char*  s_gap_open_cost_key = "Gap-open-cost";
static const char*  s_gap_extn_cost_key = "Gap-extension-cost";
static const char*  s_lambda_key        = "Lambda";
static const char*  s_K_key             = "K";
static const char*  s_H_key             = "H";
static const char*  s_alpha_key         = "Alpha";
static const char*  s_beta_key          = "Beta";
static const char*  s_scalef_key        = "Scale-factor";


// _________________________________________________________________________
// Class Configuration
//
// Class of manipulation with toolset configuration files
//
class Configuration
{
public:
    enum { NOGAPVAL = 0 };

public:
    Configuration();
    Configuration( const char* fullname );
    Configuration( const char* fullname, int, int );
    Configuration( const Configuration& );
    ~Configuration();
                                                        //assignment
    const Configuration&    operator=( const Configuration& );

    void            ReadUngapped();                     //read parameters for ungapped configuration from file
    void            WriteUngapped();                    //write parameters for ungapped configuration to file
    void            Read();                             //read parameter configuration from file
    void            Write() const;                      //write parameter configuration to file

    bool            GetAutoGapOpenCost() const          { return auto_gap_cost;     }
    int             GetGapOpenCost() const              { return gap_open_cost;     }
    int             GetGapExtendCost() const            { return gap_extend_cost;   }
    double          GetLambda() const                   { return lambda;            }
    double          GetK() const                        { return K;                 }
    double          GetH() const                        { return H;                 }
    double          GetAlpha() const                    { return alpha;             }
    double          GetBeta() const                     { return beta;              }
    double          GetScaleFactor() const              { return scale;             }

    void            SetAutoGapOpenCost( bool value )    { auto_gap_cost     = value;        }
    void            SetGapOpenCost( int value )         { gap_open_cost     = abs( value ); }
    void            SetGapExtendCost( int value )       { gap_extend_cost   = abs( value ); }
    void            SetLambda( double value )           { lambda            = value; }
    void            SetK( double value )                { K                 = value; }
    void            SetH( double value )                { H                 = value; }
    void            SetAlpha( double value )            { alpha             = value; }
    void            SetBeta( double value )             { beta              = value; }
    void            SetScaleFactor( double value )      { scale             = value; }

    const char*     GetFilename() const                 { return filename;          }
    void            SetFilename( const char* name )     { filename = name;          }

    bool            IsUngapped() const                  { return GetGapOpenCost() == NOGAPVAL && GetGapExtendCost() == NOGAPVAL; }

protected:
    void            ReadSection( const char* section, bool read_costs, bool read_scale );           //read data from section
    void            WriteSection( const char* section, bool write_costs, bool write_scale ) const;  //write data of section

    static void     CheckCost( int value );
    static void     CheckLambda( double value );
    static void     CheckK( double value );
    static void     CheckH( double value );
    static void     CheckAlpha( double value );
    static void     CheckBeta( double value );
    static void     CheckScaleFactor( double value );

    static int      GetDefault_gap_open_cost()          { return -1; }
    static int      GetDefault_gap_extend_cost()        { return -1; }
    static double   GetDefault_lambda()                 { return -1.0; }
    static double   GetDefault_K()                      { return -1.0; }
    static double   GetDefault_H()                      { return -1.0; }
    static double   GetDefault_alpha()                  { return 0.0; }
    static double   GetDefault_beta()                   { return 0.0; }
    static double   GetDefault_scale()                  { return 1.0; }

private:
    enum {
            BEGTEXT,
            CONFVER,
            BINSIGN,
            END
    };

    const char*         filename;

    bool    auto_gap_cost;      //automatically computed gap open costs
    int     gap_open_cost;      //gap open penalty
    int     gap_extend_cost;    //gap extension penalty
    double  lambda;             //statistical parameter Lambda
    double  K;                  //statistical parameter K
    double  H;                  //statistical parameter H (entropy)
    double  alpha;              //statistical parameter alpha (slope of edge-effect correction with linear regression)
    double  beta;               //statistical parameter beta (intercept of edge-effect correction with linear regression)
    double  scale;              //scale factor to scale score matrix

    static const char*  s_autocostsym;  //symbol for automatic gap open cost
    static const char*  signature[];    //signature to put/get it in/from the configuration file
};

// INLINES
// -------------------------------------------------------------------------
// CheckCost: checks for the gap cost value
//
inline
void Configuration::CheckCost( int value )
{
    if( value < 0 )
        throw myruntime_error( mystring( "Configuration: Gap cost specified is wrong." ));
}

// CheckLambda: checks for the lambda value
//
inline
void Configuration::CheckLambda( double value )
{
    if( value <= 0.0 /*|| 1.0 <= value */)
        throw myruntime_error( mystring( "Configuration: Parameter Lambda is wrong." ));
}

// CheckK: checks for the K value
//
inline
void Configuration::CheckK( double value )
{
    if( value <= 0.0 || 1.0 <= value )
        throw myruntime_error( mystring( "Configuration: Parameter K is wrong." ));
}

// CheckH: checks for the H value
//
inline
void Configuration::CheckH( double value )
{
    if( value <= 0.0 || 10.0 <= value )
        throw myruntime_error( mystring( "Configuration: Parameter H is wrong." ));
}

// CheckAlpha: checks for the alpha value
//
inline
void Configuration::CheckAlpha( double value )
{
    if( value <= 0.0 || 10.0 <= value )
        throw myruntime_error( mystring( "Configuration: Parameter Alpha is wrong." ));
}

// CheckBeta: checks for the beta value
//
inline
void Configuration::CheckBeta( double value )
{
    if( 0.0 < value && value < 0.0 )//second artificially included
        throw myruntime_error( mystring( "Configuration: Parameter Beta is wrong." ));
}

// CheckScaleFactor: checks for the scale factor value
//
inline
void Configuration::CheckScaleFactor( double value )
{
    if( value <= 0.0 || 10.0 <= value )
        throw myruntime_error( mystring( "Configuration: Scale-factor value is wrong." ));
}


#endif//__Configuration__
