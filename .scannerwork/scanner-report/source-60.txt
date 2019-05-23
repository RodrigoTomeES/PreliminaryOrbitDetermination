//$Header$
//------------------------------------------------------------------------------
//                                   gmst.h
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides the definitions of gmst function.
*
* @note
*/
//------------------------------------------------------------------------------

#ifndef GMST_H
#define GMST_H

#include "MatLabUtilites.h"
#include "Frac.h"

//------------------------------------------------------------------------------
//  double gmst(double Mjd_UT1)
//------------------------------------------------------------------------------
/**
* Greenwich Mean Sidereal Time
*
* @param  -
* Input:
*  Mjd_UT1    Modified Julian Date UT1
* @return -
* Output:
*  gmstime     GMST in [rad]
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
double gmst(double Mjd_UT1);

#endif /* GMST_H */