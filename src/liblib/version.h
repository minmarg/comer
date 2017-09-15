/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __version_h__
#define __version_h__

// Version history:
//
// 0.2, . information of statistical parameters included
// 0.3, . weights for observed frequencies appended
// 0.4, . alpha for frequency vector appended
// 0.5  . gap scheme format format changed; gap schemes now include delete-state information
// ----
//
// 1.0  . Text format introduced
//
// dbversion:
// 0.5, . type of frequency vector distribution added to the database
// 0.6, . information content is appended to frequency vectors of database
// 0.7, . effective thickness is appended to frequency vectors of database
// 0.8  . profile format changed: gap schemes now include delete-state information
// ----
//
// 1.0  . Text format introduced
// 1.1  . Background probabilities added


static const char*  dataversion = "1.1";
static const char*  dbversion = "1.1";


#endif//__version_h__
