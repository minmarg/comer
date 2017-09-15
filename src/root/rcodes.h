/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __rcodes__
#define __rcodes__


#define PRT_OK          (   0 )
#define PRT_SUCCESS     (   0 )
#define PRT_ERR_DOMAIN  (   1 )
#define PRT_MAXITERATS  (   3 )
#define PRT_ERR_ADDRESS (   5 )

inline const char* TranslatePRTError( int code )
{
    switch( code ) {
        case PRT_SUCCESS:       return "Converged";
        case PRT_ERR_DOMAIN:    return "Domain error: Not bracketed root in interval";
        case PRT_MAXITERATS:    return "Maximum number of iterations reached";
        case PRT_ERR_ADDRESS:   return "Memory access error";
    }
    return "Unknown";
}

#endif//__rcodes__
