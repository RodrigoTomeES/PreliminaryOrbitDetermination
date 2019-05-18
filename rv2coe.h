//$Header$
//------------------------------------------------------------------------------
//                           rv2coe.h
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides the definitions of rv2coe function.
*
* @note
*/
//------------------------------------------------------------------------------


#ifndef RC2COE_H
#define RC2COE_H

#include <string.h>
#include "MatLabUtilites.h"
#include "angl.h"
#include "newtonnu.h"
#define ROWS 3
#define COLS 3

void rv2coe(double r[], double v[], double *p, double *a, double *ecc, double *incl, double *omega, double *argp, double *nu, double *m, double *arglat, double *truelon, double *lonper);

#endif /* RC2COE_H */