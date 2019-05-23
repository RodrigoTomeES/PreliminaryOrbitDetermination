//$Header$
//------------------------------------------------------------------------------
//                                   IERS.h
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides the definitions of IERS function.
*
* @note
*/
//------------------------------------------------------------------------------

#ifndef IERS_H
#define IERS_H

#include <stdlib.h>
#include "MatLabUtilites.h"

//------------------------------------------------------------------------------
//  void IERS(double ** eop, double filas,double columnas,double Mjd_UTC,
//            char interp, double * UT1_UTC, double * TAI_UTC, double * x_pole,
//            double * y_pole, double * ddpsi, double * ddeps);
//------------------------------------------------------------------------------
/**
* Management of IERS time and polar motion data
*
* @param  - double ** eop, double filas,double columnas,double Mjd_UTC, char interp
* @return - double * UT1_UTC, double * TAI_UTC, double * x_pole, double * y_pole,
*           double * ddpsi, double * ddeps
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
void IERS(double ** eop, double filas,double columnas,double Mjd_UTC, char interp, double * UT1_UTC, double * TAI_UTC, double * x_pole, double * y_pole, double * ddpsi, double * ddeps);

#endif /* IERS_H */