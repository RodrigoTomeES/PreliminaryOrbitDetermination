//$Header$
//------------------------------------------------------------------------------
//                           anglesg.h
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides the definitions of anglesg function.
*
* @note
*/
//------------------------------------------------------------------------------


#ifndef ANGLESG_H
#define ANGLESG_H

#include "MatLabUtilites.h"

void anglesg(double Alpha1, double Alpha2, double Alpha3, double Delta1, double Delta2, double Delta3, double JD1, double JD2, double JD3, double RS1[], double RS2[], double RS3[], double R2[], double V2[]);

#endif /* ANGLESG_H */