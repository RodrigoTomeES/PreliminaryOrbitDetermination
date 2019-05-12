//$Header$
//------------------------------------------------------------------------------
//                           NutAngles.h
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides the definitions of NutAngles function.
*
* @note
*/
//------------------------------------------------------------------------------


#ifndef NUT_ANGLES_H
#define NUT_ANGLES_H

#include "MatLabUtilites.h"

void NutAngles(double Mjd_TT, double NutAngles_const[]);

#endif /* NUT_ANGLES_H */