//$Header$
//------------------------------------------------------------------------------
//                           Mjday.h
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides the definitions of Mjday function.
*
* @note
*/
//------------------------------------------------------------------------------


#ifndef MJDAY_H
#define MJDAY_H

#include <stdio.h>
#include "MatLabUtilites.h"

double Mjday(double year, double month, double day, double hour, double min, double sec);

#endif /* MJDAY_H */