//$Header$
//------------------------------------------------------------------------------
//                           Frac.h
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides the definitions of Frac function.
*
* @note
*/
//------------------------------------------------------------------------------

#ifndef FRAC_H
#define FRAC_H

#include "MatLabUtilites.h"

//------------------------------------------------------------------------------
//  void Frac(double x, double res[][])
//------------------------------------------------------------------------------
/**
* Fractional part of a number (y=x-[x])
*
* @param - x
* @return - res
* @exception - none
* @see - none
* @note - none
*/
//------------------------------------------------------------------------------
double Frac(double x);

#endif /* FRAC_H */