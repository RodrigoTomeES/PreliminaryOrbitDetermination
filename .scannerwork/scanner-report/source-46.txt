//$Header$
//------------------------------------------------------------------------------
//                           timediff.h
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides the definitions of timediff function.
*
* @note
*/
//------------------------------------------------------------------------------


#ifndef TIMEDIFF_H
#define TIMEDIFF_H

#include "MatLabUtilites.h"

void timediff(double UT1_UTC,double TAI_UTC,double * UT1_TAI, double * UTC_GPS, double * UT1_GPS, double * TT_UTC, double * GPS_UTC);

#endif /* TIMEDIFF_H */