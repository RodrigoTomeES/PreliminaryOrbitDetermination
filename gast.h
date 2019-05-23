//$Header$
//------------------------------------------------------------------------------
//                                   gast.h
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides the definitions of gast function.
*
* @note
*/
//------------------------------------------------------------------------------

#ifndef GAST_H
#define GAST_H

#include "MatLabUtilites.h"
#include "IERS.h"
#include "timediff.h"
#include "gmst.h"
#include "EqnEquinox.h"

//------------------------------------------------------------------------------
//  double gast(double Mjd_UT1, double ** eop, double filas,double columnas)
//------------------------------------------------------------------------------
/**
* Greenwich Apparent Sidereal Time
*
* @param  -
* Input:
*   Mjd_UT1   Modified Julian Date UT1
* @return -
* Output:
*   gstime    GAST in [rad]
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
double gast(double Mjd_UT1, double ** eop, double filas,double columnas);

#endif /* GAST_H */