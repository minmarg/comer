/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "rc.h"
#include "ConfigFile.h"
#include "datapro.h"
#include "Configuration.h"

//symbol for automatic gap open cost
const char* Configuration::s_autocostsym = "A";

//configformat1.0| -0x2a
const char* Configuration::signature[] = {
    "\nCONFIGURATION VERSION ",
    configversion,
    "\x39\x45\x44\x3c\x3f\x3d\x3c\x45\x48\x43\x37\x4a\x07\x04\x06\x00",
    NULL
};


// -------------------------------------------------------------------------
// SetUngappedParams: modifies Configuration to contain ungapped
//     parameters
//
void SetUngappedParams( Configuration& config )
{
    config.SetGapOpenCost( Configuration::NOGAPVAL );
    config.SetGapExtendCost( Configuration::NOGAPVAL );
    config.SetLambda(   LOSCORES.StatisParam( Ungapped, Lambda ));
    config.SetK(        LOSCORES.StatisParam( Ungapped, K ));
    config.SetH(        LOSCORES.StatisParam( Ungapped, H ));
    config.SetAlpha(    LOSCORES.StatisParam( Ungapped, alpha ));
    config.SetBeta(     LOSCORES.StatisParam( Ungapped, beta ));
    config.SetScaleFactor( 1.0 );
}

// -------------------------------------------------------------------------
// Constructor
//
Configuration::Configuration( const char* fullname )
:   filename( fullname )
{
    auto_gap_cost   = false;
    gap_open_cost   = NOGAPVAL;
    gap_extend_cost = NOGAPVAL;
    lambda          = GetDefault_lambda();
    K               = GetDefault_K();
    H               = GetDefault_H();
    alpha           = GetDefault_alpha();
    beta            = GetDefault_beta();
    scale           = GetDefault_scale();
}

// Constructor
//
Configuration::Configuration( const char* fullname, int cost_open, int cost_extend )
:   filename( fullname )
{
    auto_gap_cost   = false;
    gap_open_cost   = cost_open;
    gap_extend_cost = cost_extend;
    lambda          = GetDefault_lambda();
    K               = GetDefault_K();
    H               = GetDefault_H();
    alpha           = GetDefault_alpha();
    beta            = GetDefault_beta();
    scale           = GetDefault_scale();
}

// -------------------------------------------------------------------------
// Default construction
//
Configuration::Configuration()
:   filename( NULL )
{
    auto_gap_cost   = false;
    gap_open_cost   = GetDefault_gap_open_cost();
    gap_extend_cost = GetDefault_gap_extend_cost();
    lambda          = GetDefault_lambda();
    K               = GetDefault_K();
    H               = GetDefault_H();
    alpha           = GetDefault_alpha();
    beta            = GetDefault_beta();
    scale           = GetDefault_scale();
}

// Copy constructor
//
Configuration::Configuration( const Configuration& one )
{
    operator=( one );
}

// -------------------------------------------------------------------------
// Destructor
//
Configuration::~Configuration()
{
}

// -------------------------------------------------------------------------
// Assignment operator
//
const Configuration& Configuration::operator=( const Configuration& one )
{
    SetFilename( one.GetFilename());

    SetAutoGapOpenCost( one.GetAutoGapOpenCost());
    SetGapOpenCost( one.GetGapOpenCost());
    SetGapExtendCost( one.GetGapExtendCost());
    SetLambda( one.GetLambda());
    SetK( one.GetK());
    SetH( one.GetH());
    SetAlpha( one.GetAlpha());
    SetBeta( one.GetBeta());
    SetScaleFactor( one.GetScaleFactor());

    return *this;
}

// -------------------------------------------------------------------------
// ReadUngapped: reads parameters for ungapped configuration from file
//
void Configuration::ReadUngapped()
{
    SetGapOpenCost( NOGAPVAL );
    SetGapExtendCost( NOGAPVAL );
    Read();
}

// -------------------------------------------------------------------------
// WriteUngapped: write parameters of ungapped configuration to file
//
void Configuration::WriteUngapped()
{
    SetGapOpenCost( NOGAPVAL );
    SetGapExtendCost( NOGAPVAL );
    Write();
}

// -------------------------------------------------------------------------
// Read: reads parameter configuration from file
// NOTE: gap open cost and gap extension cost must be specified before
//     call of this method, because section name is determined by analysing
//     these value
//
void Configuration::Read()
{
    char    section[BUF_MAX];
    bool    auto_open = GetAutoGapOpenCost();
    int     cost_open = GetGapOpenCost();
    int     cost_extend = GetGapExtendCost();
    bool    read_costs = false;//true;
    bool    read_scale = false;

    if(( !auto_open &&
        cost_open == GetDefault_gap_open_cost()) ||
        cost_extend == GetDefault_gap_extend_cost())
        throw myruntime_error( mystring( "Configuration: Gap costs must be specified before reading of parameters." ));

    if( cost_open == NOGAPVAL && cost_extend == NOGAPVAL ) {
        strcpy( section, s_section_ungapped );
        read_costs = false;     //for ungapped configuration, do not read gap cost information
        read_scale = true;      //for ungapped configuration, scale factor must be read
    } else
        if( cost_open == NOGAPVAL || cost_extend == NOGAPVAL )
            throw myruntime_error( mystring( "Configuration: Gap costs are wrongly specified before reading of parameters." ));
        else {
            if( auto_open )
                sprintf( section, "%s%s%s%s%d", s_section_gapped, s_space, s_autocostsym, s_space, cost_extend );
            else
                sprintf( section, "%s%s%d%s%d", s_section_gapped, s_space, cost_open, s_space, cost_extend );
        }

    ReadSection( section, read_costs, read_scale );

    if(( !auto_open &&
        cost_open != GetGapOpenCost()) ||
        cost_extend != GetGapExtendCost())
        throw myruntime_error( mystring( "Configuration: Gap costs specified in configuration do not correspond section name." ));
}

// -------------------------------------------------------------------------
// Write: write parameter configuration to file; checks whether the costs
//     specified are valid at the same time
//
void Configuration::Write() const
{
    char    section[BUF_MAX];
    bool    auto_open = GetAutoGapOpenCost();
    int     cost_open = GetGapOpenCost();
    int     cost_extend = GetGapExtendCost();
    bool    write_costs = true;
    bool    write_scale = false;

    if(( !auto_open &&
        cost_open == GetDefault_gap_open_cost()) ||
        cost_extend == GetDefault_gap_extend_cost())
        throw myruntime_error( mystring( "Configuration: Write failed: Gap costs are unspecified." ));

    if( cost_open == NOGAPVAL && cost_extend == NOGAPVAL ) {
        strcpy( section, s_section_ungapped );
        write_costs = false;    //for ungapped configuration, do not write gap cost information
        write_scale = true;     //for ungapped configuration, write scale factor
    } else
        if( cost_open == NOGAPVAL || cost_extend == NOGAPVAL )
            throw myruntime_error( mystring( "Configuration: Write failed: Gap costs are wrongly specified." ));
        else {
            if( auto_open )
                sprintf( section, "%s%s%s%s%d", s_section_gapped, s_space, s_autocostsym, s_space, cost_extend );
            else
                sprintf( section, "%s%s%d%s%d", s_section_gapped, s_space, cost_open, s_space, cost_extend );
        }

    WriteSection( section, write_costs, write_scale );
}

// -------------------------------------------------------------------------
// ReadSection: reads data from section; 
//
void Configuration::ReadSection( const char* section, bool read_costs, bool read_scale )
{
    const int   tmpsz = 10;
    char    tmpstrval[tmpsz];
    int     tmpintval;
    double  tmpdblval;
    bool    auto_open = GetAutoGapOpenCost();
    bool    found = false;

    ConfigFile  config( section, GetFilename());



    if( read_costs ) {
        // Gap open cost
        if( auto_open )
            found = config.GetString( s_gap_open_cost_key, tmpstrval, tmpsz );
        else
            found = config.GetInt( s_gap_open_cost_key, &tmpintval );

        if( !found )
            throw myruntime_error( mystring( "Configuration: No key Gap-open-cost specified in configuration." ));

        if( auto_open ) {
            if( strlen( tmpstrval ) != strlen( s_autocostsym ) &&
                strcmp( tmpstrval, s_autocostsym ))
                throw myruntime_error(
                    mystring( "Configuration: Value of key Gap-open-cost does not correspond section name." ));
        } else
            SetGapOpenCost(  ( tmpintval < 0 )? -tmpintval: tmpintval );

        // Gap extension cost
        found = config.GetInt( s_gap_extn_cost_key, &tmpintval );

        if( !found )
            throw myruntime_error( mystring( "Configuration: No key Gap-extension-cost specified in configuration." ));

        SetGapExtendCost(( tmpintval < 0 )? -tmpintval: tmpintval );
    }



    // parameter Lambda
    found = config.GetDouble( s_lambda_key, &tmpdblval );

    if( !found )
        throw myruntime_error( mystring( "Configuration: No parameter Lambda found in configuration." ));

    CheckLambda( tmpdblval );
    SetLambda( tmpdblval );



    // parameter K
    found = config.GetDouble( s_K_key, &tmpdblval );

    if( !found )
        throw myruntime_error( mystring( "Configuration: No parameter K found in configuration." ));

    CheckK( tmpdblval );
    SetK( tmpdblval );


    // parameter H 
    found = config.GetDouble( s_H_key, &tmpdblval );

    if( !found )
        throw myruntime_error( mystring( "Configuration: No parameter H found in configuration." ));

    CheckH( tmpdblval );
    SetH( tmpdblval );



    // parameter alpha 
    found = config.GetDouble( s_alpha_key,  &tmpdblval );

    if( !found )
        throw myruntime_error( mystring( "Configuration: No parameter Alpha found in configuration." ));

    CheckAlpha( tmpdblval );
    SetAlpha( tmpdblval );



    // parameter beta
    found = config.GetDouble( s_beta_key,   &tmpdblval );

    if( !found )
        throw myruntime_error( mystring( "Configuration: No parameter Beta found in configuration." ));

    CheckBeta( tmpdblval );
    SetBeta( tmpdblval );


    if( ! read_scale )
        return;


    // scale factor
    found = config.GetDouble( s_scalef_key, &tmpdblval );

    if( !found )
        throw myruntime_error( mystring( "Configuration: No key Scale-factor found in configuration." ));

    CheckScaleFactor( tmpdblval );
    SetScaleFactor( tmpdblval );
}

// -------------------------------------------------------------------------
// WriteSection: write data of section
//
void Configuration::WriteSection( const char* section, bool write_costs, bool write_scale ) const
{
    ConfigFile  config( section, GetFilename());
    bool        auto_open = GetAutoGapOpenCost();

    if( write_costs ) {
        // Gap open cost
        if( auto_open )
            config.WriteString( s_gap_open_cost_key, s_autocostsym );
        else {
            CheckCost( GetGapOpenCost());
            config.WriteInt( s_gap_open_cost_key, GetGapOpenCost());
        }

        // Gap extension cost
        CheckCost( GetGapExtendCost());
        config.WriteInt( s_gap_extn_cost_key, GetGapExtendCost());
    }

    // parameter Lambda
    CheckLambda( GetLambda());
    config.WriteDouble( s_lambda_key, GetLambda());

    // parameter K
    CheckK( GetK());
    config.WriteDouble( s_K_key, GetK());

    // parameter H 
    CheckH( GetH());
    config.WriteDouble( s_H_key, GetH());

    // parameter alpha 
    CheckAlpha( GetAlpha());
    config.WriteDouble( s_alpha_key, GetAlpha());

    // parameter beta
    CheckBeta( GetBeta());
    config.WriteDouble( s_beta_key, GetBeta());

    if( ! write_scale )
        return;

    // scale factor
    CheckScaleFactor( GetScaleFactor());
    config.WriteDouble( s_scalef_key, GetScaleFactor());
}

