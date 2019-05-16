//$Header$
//------------------------------------------------------------------------------
//                           Newtonnu.h
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides the definitions of newtonnu function.
*
* @note
*/
//------------------------------------------------------------------------------

#ifndef NEWTONNU_H
#define NEWTONNU_H

#include <math.h>
#include "MatLabUtilites.h"

void newtonnu(double ecc, double nu, double * e0, double * m);

#endif /* NEWTONNU_H */