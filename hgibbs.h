//$Header$
//------------------------------------------------------------------------------
//                           hgibbs.h
//------------------------------------------------------------------------------
// POD: Preleminary Orbit Determination.
//
// Legal: MIT  License
//
// Author: David Lacalle & Rodrigo Tomé
// Created: 2019/04/27
//
/**
* Provides the definitions of hgibbs function.
*
* @note
*/
//------------------------------------------------------------------------------


#ifndef HGIBBS_H
#define HGIBBS_H

#include <string.h>
#include "MatLabUtilites.h"
#include "unit.h"
#include "angl.h"

void hgibbs(double r1[], double r2[], double r3[], double MJD1, double MJD2, double MJD3, double v2[], double *theta, double *theta1, double *copa, char error[]);

#endif /* HGIBBS_H */